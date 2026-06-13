//! The right-side **Rocket** workbench panel — the coupled
//! **design → simulate** loop over `valenx-rocket-demo`.
//!
//! Mirrors the other single-file workbenches (`frames_workbench`,
//! `fem_workbench`, …): a resizable [`egui::SidePanel`] gated on
//! `crate::ValenxApp::show_rocket_workbench`, toggled from the View menu.
//!
//! This is the front-end for the worked example that ties two valenx
//! engines into one pipeline:
//!
//! - **Trajectory** — [`valenx_astro`] flies the medium-lift two-stage
//!   preset to orbit and reports the orbit reached, the `Δv` budget,
//!   max-Q and the peak axial g-load.
//! - **Structure** — [`valenx_fem::beam`] sizes the interstage thrust
//!   struts against *that* flight: the peak g-load becomes an inertial
//!   load `F = m · a_max`, the per-strut stress is `σ = F / (N · A)`, and
//!   the safety factor is `σ_yield / σ`.
//!
//! The panel is **reactive**: editing any design knob re-runs the coupled
//! `design_and_simulate` pipeline (a bounded RK4 ascent that completes in
//! well under a frame) and refreshes the orbit + a live SAFE /
//! OVER-STRESSED verdict. A sizing aid reports the minimum per-strut area
//! needed to hit a chosen target safety factor.
//!
//! Honest scope: research / preliminary-design grade (see
//! [`valenx_rocket_demo`]). The trajectory is the fixed medium-lift preset
//! — only the interstage structure is parameterised here.

use std::path::PathBuf;

use eframe::egui;
use egui_plot::{Line, Plot, PlotPoints};

use crate::types::LoadedMesh;
use crate::ValenxApp;
use valenx_astro::{
    simulate_ascent, AscentConfig, DragModel, GuidanceMode, GuidanceProgram, Stage, Vehicle,
    WindModel,
};
use valenx_rocket_demo::{design_and_simulate, RocketDesign, RocketReport};

/// A cached Valenx LV-1 ascent: the altitude-vs-time series for the
/// in-panel plot, plus a one-glance summary line.
struct Lv1Flight {
    /// `[time_s, altitude_km]` samples for the ascent plot.
    alt_pts: Vec<[f64; 2]>,
    /// Multi-line summary (orbit / Δv / max-Q / peak g).
    summary: String,
}

/// The best feasible design an ascent-optimization run found.
#[derive(Clone, Copy)]
struct OptBest {
    payload_kg: f64,
    pitch_kick_deg: f64,
    vertical_rise_s: f64,
    periapsis_km: f64,
    apoapsis_km: f64,
    safety_factor: f64,
}

/// Result of an ascent-optimization run: the best design plus the
/// best-so-far convergence series for the plot.
struct OptResult {
    best: Option<OptBest>,
    reached_orbit: usize,
    n_evals: usize,
    /// `[eval_index, best_payload_kg_so_far]` for the convergence plot.
    convergence: Vec<[f64; 2]>,
}

/// Persistent form + result state for the Rocket workbench.
pub struct RocketWorkbenchState {
    /// The interstage design knobs fed to the coupled pipeline.
    design: RocketDesign,
    /// Target structural safety factor for the sizing-aid readout.
    target_sf: f64,
    /// The design the cached [`Self::report`] was computed for. Drives the
    /// reactive recompute: when it differs from `design` the pipeline runs.
    last_design: Option<RocketDesign>,
    /// The most recent coupled design→simulate result.
    report: Option<RocketReport>,
    /// The last pipeline error, shown in red. `None` on success.
    error: Option<String>,
    /// Cached Valenx LV-1 full-ascent flight, for the in-panel plot. `None`
    /// until first computed (lazily on first draw, or on the Fly button).
    lv1: Option<Lv1Flight>,
    /// Deferred request to load the 3-D rocket mesh into the central
    /// viewport (set by the button / first open; serviced after the panel).
    show_3d_request: bool,
    /// One-time guard so the 3-D rocket auto-loads on first open only.
    loaded_3d_once: bool,
    /// Last ascent-optimization run (best design + convergence), if any.
    opt: Option<OptResult>,
}

impl Default for RocketWorkbenchState {
    fn default() -> Self {
        Self {
            design: RocketDesign::default(),
            target_sf: 1.5,
            last_design: None,
            report: None,
            error: None,
            lv1: None,
            show_3d_request: false,
            loaded_3d_once: false,
            opt: None,
        }
    }
}

/// The Valenx LV-1 launch vehicle — a from-scratch two-stage kerolox
/// small-lift launcher (mirrors `valenx_rocket_demo::valenx_lv1`, kept
/// here so the panel flies it directly via `valenx-astro`).
fn lv1_vehicle() -> Vehicle {
    Vehicle {
        stages: vec![
            Stage {
                name: "LV-1 first stage (kerolox)".into(),
                dry_mass: 6_000.0,
                propellant_mass: 90_000.0,
                thrust_sl: 1_500_000.0,
                thrust_vac: 1_650_000.0,
                isp_sl: 283.0,
                isp_vac: 311.0,
            },
            Stage {
                name: "LV-1 second stage (kerolox, vacuum)".into(),
                dry_mass: 1_500.0,
                propellant_mass: 12_000.0,
                thrust_sl: 180_000.0,
                thrust_vac: 180_000.0,
                isp_sl: 345.0,
                isp_vac: 345.0,
            },
        ],
        payload_mass: 2_000.0,
        reference_area: 2.5,
        drag: DragModel::generic_launch_vehicle(),
    }
}

/// The LV-1 ascent profile — a gravity turn tuned (20 s vertical rise,
/// 12° pitch kick) to fly the vehicle into a bound orbit.
fn lv1_config() -> AscentConfig {
    AscentConfig {
        launch_altitude_m: 0.0,
        guidance: GuidanceProgram {
            vertical_rise_time: 20.0,
            pitch_kick_deg: 12.0,
            kick_duration: 5.0,
        },
        time_step: 0.1,
        max_time: 1_800.0,
        sample_interval: 2.0,
        mode: GuidanceMode::OpenLoopGravityTurn,
        wind: WindModel::None,
    }
}

/// Search the design space — payload × pitch-kick × vertical-rise — across
/// `n_evals` real `valenx-astro` ascent sims, keeping only designs that
/// reach a bound orbit with interstage safety factor ≥ `target_sf`, and
/// maximising payload to orbit. Deterministic (seeded) so it is testable.
fn optimize_ascent(target_sf: f64, n_evals: usize) -> OptResult {
    // Interstage capacity: 8 Al-2024-T3 struts of 15 cm² each carry the
    // wet upper stage + payload; F = m·a_max sets the load at peak g.
    let capacity_n = 324.0e6 * 8.0 * 1.5e-3;
    let upper_wet = {
        let v = lv1_vehicle();
        v.stages[1].dry_mass + v.stages[1].propellant_mass
    };

    // A tiny deterministic LCG so the search is reproducible (no rng dep).
    let mut seed: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut unit = || {
        seed = seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (seed >> 33) as f64 / ((1u64 << 31) as f64)
    };

    let mut best: Option<OptBest> = None;
    let mut best_payload = 0.0_f64;
    let mut reached_orbit = 0usize;
    let mut convergence = Vec::with_capacity(n_evals);

    for i in 0..n_evals {
        let payload = 500.0 + unit() * 7_500.0; // 0.5–8 t
        let pitch = 8.0 + unit() * 10.0; // 8–18°
        let rise = 14.0 + unit() * 31.0; // 14–45 s

        let mut v = lv1_vehicle();
        v.payload_mass = payload;
        let mut c = lv1_config();
        c.guidance.pitch_kick_deg = pitch;
        c.guidance.vertical_rise_time = rise;
        c.max_time = 900.0; // enough to insert; keeps each eval fast

        if let Ok(r) = simulate_ascent(&v, &c) {
            if r.reached_orbit && r.periapsis_km() > 100.0 {
                reached_orbit += 1;
                let load =
                    (upper_wet + payload) * r.max_acceleration_g * valenx_astro::constants::G0;
                let sf = if load > 0.0 {
                    capacity_n / load
                } else {
                    f64::INFINITY
                };
                if sf >= target_sf && payload > best_payload {
                    best_payload = payload;
                    best = Some(OptBest {
                        payload_kg: payload,
                        pitch_kick_deg: pitch,
                        vertical_rise_s: rise,
                        periapsis_km: r.periapsis_km(),
                        apoapsis_km: r.apoapsis_km(),
                        safety_factor: sf,
                    });
                }
            }
        }
        convergence.push([i as f64, best_payload]);
    }

    OptResult {
        best,
        reached_orbit,
        n_evals,
        convergence,
    }
}

/// Fly the Valenx LV-1 (via `valenx_astro::simulate_ascent`) and build the
/// cached altitude-vs-time plot series + the summary line.
fn fly_lv1() -> Lv1Flight {
    match simulate_ascent(&lv1_vehicle(), &lv1_config()) {
        Ok(r) => {
            let alt_pts = r
                .samples
                .iter()
                .map(|s| [s.time, s.altitude_m / 1000.0])
                .collect();
            let summary = format!(
                "{}  ·  {:.0} × {:.0} km orbit\nΔv {:.0} m/s · max-Q {:.0} kPa · peak {:.1} g",
                if r.reached_orbit {
                    "✔ reached orbit"
                } else {
                    "✖ suborbital"
                },
                r.periapsis_km(),
                r.apoapsis_km(),
                r.ideal_delta_v,
                r.max_dynamic_pressure / 1000.0,
                r.max_acceleration_g,
            );
            Lv1Flight { alt_pts, summary }
        }
        Err(e) => Lv1Flight {
            alt_pts: Vec::new(),
            summary: format!("ascent error: {e}"),
        },
    }
}

/// Minimum per-strut cross-sectional area (m²) needed to reach `target_sf`
/// against `load_n`, sharing the load across `strut_count` struts:
/// `A = SF · F / (σ_yield · N)` (the inverse of the demo's safety-factor
/// formula `SF = σ_yield · N · A / F`). `None` for non-positive inputs.
fn required_area_per_strut_m2(
    load_n: f64,
    yield_pa: f64,
    strut_count: usize,
    target_sf: f64,
) -> Option<f64> {
    let n = strut_count as f64;
    if load_n > 0.0 && yield_pa > 0.0 && n > 0.0 && target_sf > 0.0 {
        Some(target_sf * load_n / (yield_pa * n))
    } else {
        None
    }
}

/// Run the coupled design→simulate pipeline for the current design and
/// cache the result. Extracted from the draw closure so it is unit-testable.
fn recompute(s: &mut RocketWorkbenchState) {
    match design_and_simulate(&s.design) {
        Ok(r) => {
            s.report = Some(r);
            s.error = None;
        }
        Err(e) => {
            s.report = None;
            s.error = Some(format!("trajectory error: {e}"));
        }
    }
    s.last_design = Some(s.design);
}

/// Build the 3-D Valenx LV-1 rocket mesh and load it into the central
/// viewport (replacing any current STL / mesh) so it can be orbited.
fn load_lv1_rocket_3d(app: &mut ValenxApp) {
    let mesh = crate::rocket_mesh::lv1_rocket_mesh();
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<rocket>/valenx-lv1"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Draw the Rocket workbench right-side panel. A no-op when the
/// `show_rocket_workbench` toggle is off.
pub fn draw_rocket_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_rocket_workbench {
        return;
    }

    egui::SidePanel::right("valenx_rocket_workbench")
        .resizable(true)
        .default_width(380.0)
        .width_range(330.0..=620.0)
        .show(ctx, |ui| {
            ui.heading("Rocket — design → simulate");
            ui.label(
                egui::RichText::new("coupled ascent + structural check · valenx-rocket-demo")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.rocket;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    // ── Valenx LV-1 — watch it fly to orbit ───────────────
                    ui.label(egui::RichText::new("Valenx LV-1 — ascent to orbit").strong());
                    ui.label(
                        egui::RichText::new(
                            "a from-scratch two-stage launcher, flown live by valenx-astro",
                        )
                        .weak()
                        .small(),
                    );
                    let fly_clicked = ui
                        .button(egui::RichText::new("▶ Fly the Valenx LV-1").strong())
                        .clicked();
                    if s.lv1.is_none() || fly_clicked {
                        s.lv1 = Some(fly_lv1());
                    }
                    // 3-D rocket model → central viewport (auto-loads once
                    // on first open; the button reloads / re-frames it).
                    let show_3d = ui
                        .button(egui::RichText::new("Show the 3-D rocket model").strong())
                        .on_hover_text(
                            "Loads a 3-D model of the LV-1 into the centre viewport — orbit / zoom it.",
                        )
                        .clicked();
                    if show_3d || !s.loaded_3d_once {
                        s.loaded_3d_once = true;
                        s.show_3d_request = true;
                    }
                    if let Some(f) = &s.lv1 {
                        ui.label(egui::RichText::new(&f.summary).monospace().small());
                        if !f.alt_pts.is_empty() {
                            ui.label(
                                egui::RichText::new("altitude (km) vs time (s)")
                                    .weak()
                                    .small(),
                            );
                            Plot::new("lv1_ascent_plot").height(210.0).show(ui, |pui| {
                                pui.line(
                                    Line::new(PlotPoints::from(f.alt_pts.clone()))
                                        .name("altitude (km)"),
                                );
                            });
                        }
                    }
                    // ── AI optimizer — maximize payload to orbit ──────────
                    ui.add_space(6.0);
                    ui.label(
                        egui::RichText::new("AI optimizer — maximize payload to orbit").strong(),
                    );
                    ui.label(
                        egui::RichText::new(
                            "searches payload × pitch-kick × vertical-rise across many real \
                             valenx-astro sims, keeping interstage SF ≥ the target below.",
                        )
                        .weak()
                        .small(),
                    );
                    if ui
                        .button(egui::RichText::new("Run AI optimization (200 sims)").strong())
                        .on_hover_text(
                            "Flies 200 candidate designs through the real ascent engine and \
                             converges on the heaviest payload that still reaches orbit safely.",
                        )
                        .clicked()
                    {
                        s.opt = Some(optimize_ascent(s.target_sf, 200));
                    }
                    if let Some(o) = &s.opt {
                        match &o.best {
                            Some(b) => ui.label(
                                egui::RichText::new(format!(
                                    "ran {} sims · {} reached orbit\n\
                                     best payload {:.0} kg → {:.0} × {:.0} km  (SF {:.2})\n\
                                     @ pitch {:.1}° · rise {:.0} s",
                                    o.n_evals,
                                    o.reached_orbit,
                                    b.payload_kg,
                                    b.periapsis_km,
                                    b.apoapsis_km,
                                    b.safety_factor,
                                    b.pitch_kick_deg,
                                    b.vertical_rise_s,
                                ))
                                .monospace()
                                .small(),
                            ),
                            None => ui.colored_label(
                                egui::Color32::from_rgb(220, 160, 60),
                                format!(
                                    "ran {} sims · none met SF ≥ {:.2}",
                                    o.n_evals, s.target_sf
                                ),
                            ),
                        };
                        if o.convergence.len() > 1 {
                            ui.label(
                                egui::RichText::new("best payload (kg) vs sim #")
                                    .weak()
                                    .small(),
                            );
                            Plot::new("lv1_opt_plot").height(170.0).show(ui, |pui| {
                                pui.line(
                                    Line::new(PlotPoints::from(o.convergence.clone()))
                                        .name("best payload (kg)"),
                                );
                            });
                        }
                    }
                    ui.add_space(8.0);
                    ui.separator();
                    ui.label(
                        egui::RichText::new(
                            "Trajectory: fixed medium-lift two-stage preset. \
                             Edit the interstage struts below — the panel re-flies \
                             the ascent and re-checks the structure live.",
                        )
                        .weak()
                        .small(),
                    );

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Interstage structure").strong());
                    ui.horizontal(|ui| {
                        ui.label("supported mass");
                        ui.add(
                            egui::DragValue::new(&mut s.design.supported_mass_kg)
                                .speed(100.0)
                                .range(100.0..=200_000.0)
                                .suffix(" kg"),
                        )
                        .on_hover_text("Upper stage + payload carried by the interstage.");
                    });
                    ui.horizontal(|ui| {
                        ui.label("strut count N");
                        ui.add(
                            egui::DragValue::new(&mut s.design.strut_count)
                                .speed(0.2)
                                .range(1..=64),
                        );
                    });
                    ui.horizontal(|ui| {
                        ui.label("strut area A");
                        // Edit in cm² for usability; store in m². Only write back
                        // on an actual edit so the reactive dirty-check below stays
                        // stable (no per-frame float drift).
                        let mut area_cm2 = s.design.strut_area_m2 * 1.0e4;
                        if ui
                            .add(
                                egui::DragValue::new(&mut area_cm2)
                                    .speed(0.1)
                                    .range(0.01..=2000.0)
                                    .suffix(" cm²"),
                            )
                            .changed()
                        {
                            s.design.strut_area_m2 = area_cm2 * 1.0e-4;
                        }
                    });
                    ui.horizontal(|ui| {
                        ui.label("material yield σy");
                        let mut yield_mpa = s.design.material_yield_pa / 1.0e6;
                        if ui
                            .add(
                                egui::DragValue::new(&mut yield_mpa)
                                    .speed(5.0)
                                    .range(1.0..=2000.0)
                                    .suffix(" MPa"),
                            )
                            .on_hover_text("≈ 324 MPa for Al-2024-T3.")
                            .changed()
                        {
                            s.design.material_yield_pa = yield_mpa * 1.0e6;
                        }
                    });

                    // ── Reactive recompute ────────────────────────────────
                    // Re-fly + re-check whenever the design changed (and on the
                    // first draw, when `last_design` is None).
                    if s.last_design != Some(s.design) {
                        recompute(s);
                    }

                    if let Some(e) = &s.error {
                        ui.add_space(4.0);
                        ui.colored_label(egui::Color32::from_rgb(220, 90, 90), e);
                    }

                    if let Some(r) = s.report {
                        ui.separator();
                        ui.label(egui::RichText::new("Trajectory — valenx-astro").strong());
                        let orbit = if r.reached_orbit {
                            "reached orbit"
                        } else {
                            "no stable orbit"
                        };
                        ui.label(
                            egui::RichText::new(format!(
                                "  {orbit}\n  \
                                 apoapsis  : {:.0} km\n  \
                                 periapsis : {:.0} km\n  \
                                 Δv budget : {:.0} m/s\n  \
                                 max-Q     : {:.1} kPa\n  \
                                 peak g    : {:.1} g",
                                r.apoapsis_km,
                                r.periapsis_km,
                                r.delta_v_budget_ms,
                                r.max_q_pa / 1000.0,
                                r.max_acceleration_g,
                            ))
                            .monospace()
                            .small(),
                        );

                        ui.add_space(4.0);
                        ui.label(
                            egui::RichText::new("Structure — valenx-fem (loaded by peak g)")
                                .strong(),
                        );
                        ui.label(
                            egui::RichText::new(format!(
                                "  peak axial load : {:.0} kN\n  \
                                 strut stress    : {:.0} MPa\n  \
                                 safety factor   : {:.2}",
                                r.peak_axial_load_n / 1000.0,
                                r.strut_stress_pa / 1.0e6,
                                r.structural_safety_factor,
                            ))
                            .monospace()
                            .small(),
                        );

                        // Live verdict banner — margin of safety = SF − 1.
                        let ms_pct = (r.structural_safety_factor - 1.0) * 100.0;
                        let (txt, col) = if r.structurally_safe {
                            (
                                format!("✔ SAFE  ·  margin of safety {ms_pct:+.0}%"),
                                egui::Color32::from_rgb(80, 220, 120),
                            )
                        } else {
                            (
                                format!("✖ OVER-STRESSED  ·  margin {ms_pct:+.0}%"),
                                egui::Color32::from_rgb(220, 90, 90),
                            )
                        };
                        ui.add_space(2.0);
                        ui.colored_label(col, egui::RichText::new(txt).strong());

                        // ── Sizing aid ────────────────────────────────────
                        ui.add_space(6.0);
                        ui.label(egui::RichText::new("Sizing aid").strong());
                        ui.horizontal(|ui| {
                            ui.label("target SF");
                            ui.add(egui::Slider::new(&mut s.target_sf, 1.0..=3.0));
                        });
                        if let Some(a_req) = required_area_per_strut_m2(
                            r.peak_axial_load_n,
                            s.design.material_yield_pa,
                            s.design.strut_count,
                            s.target_sf,
                        ) {
                            let have = s.design.strut_area_m2;
                            let meets = have >= a_req;
                            ui.label(
                                egui::RichText::new(format!(
                                    "  need ≥ {:.2} cm²/strut for SF ≥ {:.2}  (have {:.2} cm²)",
                                    a_req * 1.0e4,
                                    s.target_sf,
                                    have * 1.0e4,
                                ))
                                .monospace()
                                .small(),
                            );
                            let (txt, col) = if meets {
                                (
                                    "✔ current struts meet the target".to_string(),
                                    egui::Color32::from_rgb(80, 220, 120),
                                )
                            } else {
                                (
                                    format!(
                                        "✖ enlarge struts ×{:.2} to meet the target",
                                        a_req / have.max(1e-12)
                                    ),
                                    egui::Color32::from_rgb(220, 160, 60),
                                )
                            };
                            ui.colored_label(col, egui::RichText::new(txt).small());
                        }
                    }
                });
        });

    // Deferred (outside the panel borrow): load the 3-D rocket mesh into
    // the central viewport when requested (the button, or first open).
    if app.rocket.show_3d_request {
        app.rocket.show_3d_request = false;
        load_lv1_rocket_3d(app);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = RocketWorkbenchState::default();
        assert!(s.report.is_none());
        assert!(s.error.is_none());
        assert!(s.last_design.is_none());
        assert!(s.lv1.is_none());
        assert!(s.opt.is_none());
        assert!((s.target_sf - 1.5).abs() < 1e-12);
    }

    #[test]
    fn optimizer_finds_a_feasible_payload_with_monotone_convergence() {
        let r = optimize_ascent(1.5, 120);
        assert_eq!(r.n_evals, 120);
        assert_eq!(r.convergence.len(), 120);
        assert!(r.reached_orbit > 0, "some candidates should reach orbit");
        let b = r.best.expect("a feasible design is found");
        assert!((500.0..=8000.0).contains(&b.payload_kg));
        assert!(b.periapsis_km > 100.0, "periapsis {} km", b.periapsis_km);
        assert!(b.safety_factor >= 1.5, "SF {}", b.safety_factor);
        // Best-so-far is monotonically non-decreasing across the search.
        for w in r.convergence.windows(2) {
            assert!(w[1][1] >= w[0][1] - 1e-9, "convergence is non-decreasing");
        }
        // Deterministic (seeded): a second run gives the same best payload.
        let r2 = optimize_ascent(1.5, 120);
        assert_eq!(
            r2.best.map(|b| b.payload_kg as u64),
            Some(b.payload_kg as u64)
        );
    }

    #[test]
    fn lv1_flight_reaches_orbit_with_a_plottable_ascent() {
        // The in-panel LV-1 flight produces a non-empty altitude series and
        // a summary confirming it reaches orbit (the astro engine is
        // validated separately in valenx-rocket-demo).
        let f = fly_lv1();
        assert!(!f.alt_pts.is_empty(), "ascent yields samples to plot");
        assert!(
            f.summary.contains("reached orbit"),
            "summary: {}",
            f.summary
        );
        // The series climbs (last sample altitude well above the first).
        let first = f.alt_pts.first().unwrap()[1];
        let last = f.alt_pts.last().unwrap()[1];
        assert!(last > first + 100.0, "climbs from {first} to {last} km");
    }

    #[test]
    fn recompute_default_reaches_orbit_and_is_safe() {
        let mut s = RocketWorkbenchState::default();
        recompute(&mut s);
        assert!(s.error.is_none(), "unexpected error: {:?}", s.error);
        let r = s.report.expect("a successful run yields a report");
        assert!(r.reached_orbit, "default preset reaches orbit");
        assert!(r.periapsis_km > 100.0, "periapsis {} km", r.periapsis_km);
        // Default 8×10 cm² Al struts survive the peak g-load.
        assert!(r.structurally_safe, "SF {}", r.structural_safety_factor);
        // The dirty-check anchor is updated so the panel won't recompute again
        // until the design changes.
        assert_eq!(s.last_design, Some(s.design));
    }

    #[test]
    fn shrinking_struts_flips_to_over_stressed() {
        let mut s = RocketWorkbenchState::default();
        recompute(&mut s);
        let safe_sf = s.report.unwrap().structural_safety_factor;
        // Shrink each strut 100× → stress 100× → safety factor far below 1.
        s.design.strut_area_m2 = 1.0e-5;
        recompute(&mut s);
        let r = s.report.unwrap();
        assert!(
            !r.structurally_safe,
            "tiny struts over-stressed: SF {}",
            r.structural_safety_factor
        );
        assert!(r.structural_safety_factor < safe_sf);
    }

    #[test]
    fn heavier_payload_raises_the_axial_load() {
        // The structural load is coupled to the mass: doubling the supported
        // mass doubles the peak axial load (same fixed trajectory ⇒ same g).
        let mut s = RocketWorkbenchState::default();
        recompute(&mut s);
        let load1 = s.report.unwrap().peak_axial_load_n;
        s.design.supported_mass_kg *= 2.0;
        recompute(&mut s);
        let load2 = s.report.unwrap().peak_axial_load_n;
        assert!(
            (load2 / load1 - 2.0).abs() < 1e-9,
            "load ∝ mass: {load1} → {load2}"
        );
    }

    #[test]
    fn required_area_matches_closed_form() {
        // A = SF·F/(σ·N): 2·100/(10·5) = 4.
        assert_eq!(required_area_per_strut_m2(100.0, 10.0, 5, 2.0), Some(4.0));
        // Non-positive inputs → None (no panic, no divide-by-zero).
        assert!(required_area_per_strut_m2(0.0, 10.0, 5, 2.0).is_none());
        assert!(required_area_per_strut_m2(100.0, 10.0, 0, 2.0).is_none());

        // Round-trip against the demo: sizing each strut to exactly the
        // computed required area yields a design whose safety factor equals
        // the target (within tolerance).
        let mut s = RocketWorkbenchState::default();
        recompute(&mut s);
        let r = s.report.unwrap();
        let target = 2.0;
        let a = required_area_per_strut_m2(
            r.peak_axial_load_n,
            s.design.material_yield_pa,
            s.design.strut_count,
            target,
        )
        .unwrap();
        s.design.strut_area_m2 = a;
        recompute(&mut s);
        assert!(
            (s.report.unwrap().structural_safety_factor - target).abs() < 1e-9,
            "sizing to the required area hits the target SF"
        );
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    /// Render the whole workbench panel once in a headless egui context.
    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_rocket_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_rocket_workbench);
        draw_workbench(&mut app);
        // Hidden ⇒ nothing computed.
        assert!(app.rocket.report.is_none());
    }

    #[test]
    fn workbench_computes_and_draws_on_first_open() {
        // Showing the panel auto-runs the coupled pipeline on the first draw
        // (last_design is None) and renders the result without panicking.
        let mut app = ValenxApp::default();
        app.show_rocket_workbench = true;
        draw_workbench(&mut app);
        assert!(app.rocket.report.is_some(), "first draw computes a report");
    }

    #[test]
    fn workbench_draws_over_stressed_state_without_panic() {
        let mut app = ValenxApp::default();
        app.show_rocket_workbench = true;
        app.rocket.design.strut_area_m2 = 1.0e-5; // tiny struts → over-stressed
        draw_workbench(&mut app);
        let r = app.rocket.report.expect("computed");
        assert!(!r.structurally_safe);
    }

    #[test]
    fn workbench_draws_an_error_state_without_panic() {
        // A surfaced error line must render without panicking.
        let mut app = ValenxApp::default();
        app.show_rocket_workbench = true;
        app.rocket.error = Some("trajectory error: bad case".to_string());
        draw_workbench(&mut app);
    }
}
