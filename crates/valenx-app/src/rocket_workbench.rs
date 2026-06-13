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

use eframe::egui;

use crate::ValenxApp;
use valenx_rocket_demo::{design_and_simulate, RocketDesign, RocketReport};

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
}

impl Default for RocketWorkbenchState {
    fn default() -> Self {
        Self {
            design: RocketDesign::default(),
            target_sf: 1.5,
            last_design: None,
            report: None,
            error: None,
        }
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
        assert!((s.target_sf - 1.5).abs() < 1e-12);
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
