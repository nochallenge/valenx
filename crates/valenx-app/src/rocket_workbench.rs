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
use std::sync::mpsc::{Receiver, TryRecvError};

use eframe::egui;
use egui_plot::{Line, Plot, PlotPoints, Points};

use crate::types::LoadedMesh;
use crate::ValenxApp;
use valenx_astro::{
    simulate_ascent, AscentConfig, DragModel, GuidanceMode, GuidanceProgram, Stage,
    TrajectorySample, Vehicle, WindModel,
};
use valenx_rocket_demo::{auto_design, design_and_simulate, RocketDesign, RocketReport};

/// A cached Valenx LV-1 ascent: the altitude-vs-time series for the
/// in-panel plot, plus a one-glance summary line.
struct Lv1Flight {
    /// `[time_s, altitude_km]` samples for the ascent plot.
    alt_pts: Vec<[f64; 2]>,
    /// Multi-line summary (orbit / Δv / max-Q / peak g).
    summary: String,
    /// `[time_s, altitude_km]` of the staging + MECO flight events — drawn
    /// as markers on the ascent plot.
    event_pts: Vec<[f64; 2]>,
    /// Full per-instant trajectory state, retained so the playback scrubber
    /// can report altitude / speed / Mach / dynamic pressure / g-load at any
    /// mission time without re-flying the ascent.
    samples: Vec<TrajectorySample>,
    /// Mission-elapsed time of the final sample (s) — the scrubber's upper
    /// bound. Zero for a failed ascent (empty `samples`).
    t_max: f64,
}

/// What the AI optimizer drives toward. Every objective keeps the same hard
/// constraint — reach a bound orbit (periapsis > 100 km) with interstage
/// safety factor ≥ the target — and differs only in what it then rewards.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
enum OptObjective {
    /// Heaviest payload delivered to orbit.
    #[default]
    MaxPayload,
    /// Highest apoapsis (orbital energy / reach).
    MaxApoapsis,
    /// Gentlest ride — lowest peak axial g-load.
    MinPeakG,
}

impl OptObjective {
    /// Short label for the radio + result readout.
    fn label(self) -> &'static str {
        match self {
            OptObjective::MaxPayload => "max payload",
            OptObjective::MaxApoapsis => "max apoapsis",
            OptObjective::MinPeakG => "min peak-g",
        }
    }
}

/// Which trajectory channel the LV-1 ascent plot draws against time. The
/// scrubber readout always lists every quantity; this only selects the curve.
#[derive(Clone, Copy, PartialEq, Eq, Default)]
enum PlotQuantity {
    /// Geometric altitude (km).
    #[default]
    Altitude,
    /// Speed relative to the rotating atmosphere (m/s).
    Speed,
    /// Mach number relative to the local air.
    Mach,
    /// Dynamic pressure (kPa).
    DynPressure,
    /// Sensed (non-gravitational) acceleration (g).
    AccelG,
}

impl PlotQuantity {
    /// This channel's value for one sample, in the plotted display units.
    fn value(self, s: &TrajectorySample) -> f64 {
        match self {
            PlotQuantity::Altitude => s.altitude_m / 1000.0,
            PlotQuantity::Speed => s.speed_relative,
            PlotQuantity::Mach => s.mach,
            PlotQuantity::DynPressure => s.dynamic_pressure / 1000.0,
            PlotQuantity::AccelG => s.acceleration_g,
        }
    }

    /// The plot series / y-axis label, including units.
    fn axis_label(self) -> &'static str {
        match self {
            PlotQuantity::Altitude => "altitude (km)",
            PlotQuantity::Speed => "airspeed (m/s)",
            PlotQuantity::Mach => "Mach",
            PlotQuantity::DynPressure => "dynamic pressure (kPa)",
            PlotQuantity::AccelG => "sensed accel (g)",
        }
    }
}

/// The best feasible design an ascent-optimization run found.
#[derive(Clone, Copy)]
struct OptBest {
    payload_kg: f64,
    pitch_kick_deg: f64,
    vertical_rise_s: f64,
    periapsis_km: f64,
    apoapsis_km: f64,
    /// Peak axial g-load over the ascent (the MinPeakG objective).
    peak_g: f64,
    safety_factor: f64,
    /// The optimized metric's value in its own display units (kg / km / g).
    objective_value: f64,
}

/// Result of an ascent-optimization run: the best design plus the
/// best-so-far convergence series for the plot.
struct OptResult {
    best: Option<OptBest>,
    /// The objective this run optimized — drives the readout + plot labels.
    objective: OptObjective,
    reached_orbit: usize,
    n_evals: usize,
    /// `[eval_index, best_objective_value_so_far]` for the convergence plot,
    /// in the objective's display units. `NaN` until the first feasible
    /// design is found (filtered out before plotting).
    convergence: Vec<[f64; 2]>,
}

/// A live progress snapshot from a running background optimization.
#[derive(Clone, Copy, Default)]
struct OptProgress {
    /// Evaluations completed so far.
    done: usize,
    /// How many of those reached a bound orbit.
    reached_orbit: usize,
    /// Best objective value so far, in display units (`NaN` until feasible).
    best_value: f64,
}

/// A message from the background optimizer thread to the UI.
enum OptMsg {
    /// A throttled progress tick.
    Progress(OptProgress),
    /// The finished result (sent once, last).
    Done(OptResult),
}

/// A background optimization in flight: the UI polls `rx` each frame for
/// progress + the final result, so thousands of sims never block the UI.
struct OptJob {
    rx: Receiver<OptMsg>,
    /// The objective being optimized (fixed at spawn; for the progress label).
    objective: OptObjective,
    /// Total evaluations the run will perform.
    n_evals: usize,
    /// The most recent progress tick received.
    last_progress: OptProgress,
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
    /// Which objective the AI optimizer drives toward (radio-selected).
    opt_objective: OptObjective,
    /// A background optimization in flight (None when idle). Polled each
    /// frame so the UI never blocks while thousands of sims run.
    opt_job: Option<OptJob>,
    /// Last automated end-to-end design (engine + ascent), if any — produced
    /// by the one-click "Auto-design", where the AI does the whole search.
    auto: Option<auto_design::BestDesign>,
    /// Playback-scrubber position (mission-elapsed seconds) for the LV-1
    /// ascent inspector. Clamped to the current flight's `t_max` each draw.
    scrub_t: f64,
    /// Which trajectory channel the ascent plot draws against time.
    plot_y: PlotQuantity,
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
            opt_objective: OptObjective::default(),
            opt_job: None,
            auto: None,
            scrub_t: 0.0,
            plot_y: PlotQuantity::default(),
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

/// The no-progress convenience wrapper used by tests (the live UI uses
/// [`spawn_opt_job`] → [`optimize_ascent_with`]).
#[cfg(test)]
fn optimize_ascent(objective: OptObjective, target_sf: f64, n_evals: usize) -> OptResult {
    optimize_ascent_with(objective, target_sf, n_evals, |_, _, _| {})
}

/// Search the design space — payload × pitch-kick × vertical-rise — across
/// `n_evals` real `valenx-astro` ascent sims, keeping only designs that
/// reach a bound orbit with interstage safety factor ≥ `target_sf`, and
/// rewarding whichever `objective` is selected (heaviest payload, highest
/// apoapsis, or gentlest peak-g ride). Deterministic (seeded) so it is
/// testable — and the candidate sequence is identical across objectives, so
/// each objective is a true arg-best over the same feasible set. `on_progress`
/// is called once per evaluation with `(evals_done, reached_orbit_count,
/// best_objective_value_so_far)` to drive a live counter on a background run.
fn optimize_ascent_with(
    objective: OptObjective,
    target_sf: f64,
    n_evals: usize,
    mut on_progress: impl FnMut(usize, usize, f64),
) -> OptResult {
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
    // The objective is always framed so HIGHER is better (MinPeakG negates g),
    // so one `>` comparison drives every objective.
    let mut best_score = f64::NEG_INFINITY;
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
                if sf >= target_sf {
                    // Reward — higher is better for every objective.
                    let score = match objective {
                        OptObjective::MaxPayload => payload,
                        OptObjective::MaxApoapsis => r.apoapsis_km(),
                        OptObjective::MinPeakG => -r.max_acceleration_g,
                    };
                    if score > best_score {
                        best_score = score;
                        best = Some(OptBest {
                            payload_kg: payload,
                            pitch_kick_deg: pitch,
                            vertical_rise_s: rise,
                            periapsis_km: r.periapsis_km(),
                            apoapsis_km: r.apoapsis_km(),
                            peak_g: r.max_acceleration_g,
                            safety_factor: sf,
                            objective_value: match objective {
                                OptObjective::MaxPayload => payload,
                                OptObjective::MaxApoapsis => r.apoapsis_km(),
                                OptObjective::MinPeakG => r.max_acceleration_g,
                            },
                        });
                    }
                }
            }
        }
        // Best-so-far in the objective's display units; NaN until the first
        // feasible design (the plot filters non-finite points).
        let display = best.map(|b| b.objective_value).unwrap_or(f64::NAN);
        convergence.push([i as f64, display]);
        on_progress(i + 1, reached_orbit, display);
    }

    OptResult {
        best,
        objective,
        reached_orbit,
        n_evals,
        convergence,
    }
}

/// Spawn the optimizer on a background thread, streaming throttled progress
/// and the final result back over a channel. Returns immediately; the UI
/// polls the returned [`OptJob`] via [`poll_opt_job`] and never blocks.
fn spawn_opt_job(objective: OptObjective, target_sf: f64, n_evals: usize) -> OptJob {
    let (tx, rx) = std::sync::mpsc::channel::<OptMsg>();
    std::thread::spawn(move || {
        let prog_tx = tx.clone();
        // Report ~1% steps (and the final tick) to keep the channel light.
        let step = (n_evals / 100).max(1);
        let result = optimize_ascent_with(objective, target_sf, n_evals, |done, reached, best| {
            if done % step == 0 || done == n_evals {
                let _ = prog_tx.send(OptMsg::Progress(OptProgress {
                    done,
                    reached_orbit: reached,
                    best_value: best,
                }));
            }
        });
        let _ = tx.send(OptMsg::Done(result));
    });
    OptJob {
        rx,
        objective,
        n_evals,
        last_progress: OptProgress::default(),
    }
}

/// Drain a running background optimization's channel: store the finished
/// [`OptResult`] on the state when it arrives, otherwise return a live
/// progress snapshot `(progress, n_evals, objective)` (requesting a repaint so
/// the counter keeps ticking). `None` when idle or just-finished.
fn poll_opt_job(
    s: &mut RocketWorkbenchState,
    ctx: &egui::Context,
) -> Option<(OptProgress, usize, OptObjective)> {
    let job = s.opt_job.as_mut()?;
    let mut finished: Option<OptResult> = None;
    let mut disconnected = false;
    loop {
        match job.rx.try_recv() {
            Ok(OptMsg::Progress(p)) => job.last_progress = p,
            Ok(OptMsg::Done(r)) => {
                finished = Some(r);
                break;
            }
            Err(TryRecvError::Empty) => break,
            Err(TryRecvError::Disconnected) => {
                disconnected = true;
                break;
            }
        }
    }
    if let Some(r) = finished {
        s.opt = Some(r);
        s.opt_job = None;
        None
    } else if disconnected {
        // The worker vanished without a result (e.g. it panicked) — clear the
        // job so the UI doesn't sit on a stuck progress bar forever.
        s.opt_job = None;
        None
    } else {
        let snapshot = (job.last_progress, job.n_evals, job.objective);
        ctx.request_repaint();
        Some(snapshot)
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
            let event_pts = r
                .events
                .iter()
                .filter(|e| {
                    let k = e.kind.to_ascii_lowercase();
                    k.contains("staging") || k.contains("meco")
                })
                .map(|e| [e.time, e.altitude_m / 1000.0])
                .collect();
            let t_max = r.samples.last().map(|s| s.time).unwrap_or(0.0);
            let samples = r.samples;
            Lv1Flight {
                alt_pts,
                summary,
                event_pts,
                samples,
                t_max,
            }
        }
        Err(e) => Lv1Flight {
            alt_pts: Vec::new(),
            summary: format!("ascent error: {e}"),
            event_pts: Vec::new(),
            samples: Vec::new(),
            t_max: 0.0,
        },
    }
}

/// Linearly interpolate the recorded trajectory state at mission time `t`
/// (seconds). `t` is clamped to the sampled window: before the first sample
/// it returns the first, after the last it returns the last. Returns `None`
/// only when there are no samples (a failed ascent). Each field is
/// interpolated linearly between the two bracketing samples — exact enough
/// for the densely down-sampled ascent series the inspector scrubs over.
fn sample_at(samples: &[TrajectorySample], t: f64) -> Option<TrajectorySample> {
    let first = samples.first()?;
    let last = samples.last()?;
    if t <= first.time {
        return Some(*first);
    }
    if t >= last.time {
        return Some(*last);
    }
    let hi = samples
        .iter()
        .position(|s| s.time >= t)
        .unwrap_or(samples.len() - 1)
        .max(1);
    let a = &samples[hi - 1];
    let b = &samples[hi];
    let span = b.time - a.time;
    let frac = if span > 0.0 { (t - a.time) / span } else { 0.0 };
    let lerp = |x: f64, y: f64| x + (y - x) * frac;
    Some(TrajectorySample {
        time: t,
        altitude_m: lerp(a.altitude_m, b.altitude_m),
        downrange_m: lerp(a.downrange_m, b.downrange_m),
        speed_inertial: lerp(a.speed_inertial, b.speed_inertial),
        speed_relative: lerp(a.speed_relative, b.speed_relative),
        mach: lerp(a.mach, b.mach),
        mass: lerp(a.mass, b.mass),
        dynamic_pressure: lerp(a.dynamic_pressure, b.dynamic_pressure),
        acceleration_g: lerp(a.acceleration_g, b.acceleration_g),
    })
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

/// Density of Al-2024-T3 (kg/m³) — the interstage strut material — for the
/// lightest-interstage mass estimate.
const AL2024_DENSITY_KG_M3: f64 = 2_780.0;
/// Representative interstage strut length (m) used to turn the minimum
/// load-bearing area into an estimated minimum strut mass.
const INTERSTAGE_STRUT_LENGTH_M: f64 = 2.0;

/// Minimum **total** load-bearing strut cross-section (m²) needed to reach
/// `target_sf` against `load_n`: `A_total = SF · F / σ_yield`. This is
/// independent of how many struts share it — the per-strut area scales as
/// `1/N`, so the total is fixed — which is why it is the honest "lightest
/// interstage" figure of merit. `None` for non-positive inputs.
fn min_total_strut_area_m2(load_n: f64, yield_pa: f64, target_sf: f64) -> Option<f64> {
    if load_n > 0.0 && yield_pa > 0.0 && target_sf > 0.0 {
        Some(target_sf * load_n / yield_pa)
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
            if crate::workbench_ui::header(
                ui,
                "Rocket — design → simulate",
                "coupled ascent + structural check · valenx-rocket-demo",
            ) {
                app.show_rocket_workbench = false;
            }

            let s = &mut app.rocket;
            // Poll any background optimization before drawing (non-blocking);
            // `opt_running` is a live progress snapshot while one is in flight.
            let opt_running = poll_opt_job(s, ui.ctx());
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    // ── Auto-design — the AI does the whole search ────────
                    ui.label(
                        egui::RichText::new("Auto-design — let the AI design it for you")
                            .strong(),
                    );
                    ui.label(
                        egui::RichText::new(
                            "one click: optimize the engine + the trajectory together and \
                             return the best rocket — no tuning. (You can still tune by hand \
                             below.)",
                        )
                        .weak()
                        .small(),
                    );
                    if ui
                        .button(egui::RichText::new("🤖 Auto-design the best rocket").strong())
                        .on_hover_text(
                            "Runs the full automated search — an optimized, cooled engine plus \
                             the best ascent — and returns the heaviest payload to orbit that \
                             stays structurally sound.",
                        )
                        .clicked()
                    {
                        s.auto = auto_design::auto_design(2_000);
                    }
                    if let Some(d) = &s.auto {
                        ui.label(
                            egui::RichText::new(format!(
                                "✦ AI design: {:.0} kg → {:.0} × {:.0} km\n\
                                 engine {:.0} bar · ε {:.1} · vac Isp {:.0} s · cooling {:.2}\n\
                                 pitch {:.1}° · rise {:.0} s · peak {:.1} g · interstage SF {:.2}",
                                d.payload_kg,
                                d.periapsis_km,
                                d.apoapsis_km,
                                d.engine.chamber_pressure / 1.0e5,
                                d.engine.expansion_ratio,
                                d.engine_vacuum.isp,
                                d.engine_cooling.cooling_margin,
                                d.pitch_kick_deg,
                                d.vertical_rise_s,
                                d.peak_g,
                                d.structural_sf,
                            ))
                            .monospace()
                            .small(),
                        );
                    }
                    ui.add_space(6.0);
                    ui.separator();

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
                    if s.lv1.is_some() {
                        // Summary + flight window (immutable view, released before
                        // the scrubber borrows `scrub_t` mutably below).
                        let (t_max, has_plot) = {
                            let f = s.lv1.as_ref().unwrap();
                            ui.label(egui::RichText::new(&f.summary).monospace().small());
                            (f.t_max.max(0.0), !f.alt_pts.is_empty())
                        };
                        s.scrub_t = s.scrub_t.clamp(0.0, t_max);

                        if has_plot {
                            // Pick which telemetry channel the plot draws; the
                            // numeric readout below always shows every quantity.
                            ui.horizontal_wrapped(|ui| {
                                ui.label(egui::RichText::new("plot:").small());
                                ui.radio_value(&mut s.plot_y, PlotQuantity::Altitude, "alt");
                                ui.radio_value(&mut s.plot_y, PlotQuantity::Speed, "speed");
                                ui.radio_value(&mut s.plot_y, PlotQuantity::Mach, "Mach");
                                ui.radio_value(&mut s.plot_y, PlotQuantity::DynPressure, "q");
                                ui.radio_value(&mut s.plot_y, PlotQuantity::AccelG, "g");
                            });
                            let q = s.plot_y;
                            ui.label(
                                egui::RichText::new(format!(
                                    "{} vs time (s) · dots = staging + MECO · drag to scrub",
                                    q.axis_label()
                                ))
                                .weak()
                                .small(),
                            );

                            // ▶ Playback scrubber — drag through the flight and
                            // read the full vehicle state at that instant.
                            if t_max > 0.0 {
                                ui.add(
                                    egui::Slider::new(&mut s.scrub_t, 0.0..=t_max)
                                        .text("▶ playback · t (s)"),
                                );
                            }

                            let f = s.lv1.as_ref().unwrap();
                            let scrub = sample_at(&f.samples, s.scrub_t);
                            if let Some(p) = scrub {
                                ui.label(
                                    egui::RichText::new(format!(
                                        "t {:>6.1} s   alt {:>7.2} km   downrange {:>7.2} km\n\
                                         v_rel {:>6.0} m/s   v_inert {:>6.0} m/s   Mach {:>5.2}\n\
                                         q {:>6.1} kPa   accel {:>5.2} g   mass {:>7.2} t",
                                        p.time,
                                        p.altitude_m / 1000.0,
                                        p.downrange_m / 1000.0,
                                        p.speed_relative,
                                        p.speed_inertial,
                                        p.mach,
                                        p.dynamic_pressure / 1000.0,
                                        p.acceleration_g,
                                        p.mass / 1000.0,
                                    ))
                                    .monospace()
                                    .small(),
                                );
                            }

                            // Line + event + now markers all follow the selected
                            // channel, derived from the one retained sample series.
                            let line_pts: Vec<[f64; 2]> = if q == PlotQuantity::Altitude {
                                f.alt_pts.clone()
                            } else {
                                f.samples
                                    .iter()
                                    .map(|smp| [smp.time, q.value(smp)])
                                    .collect()
                            };
                            let event_pts: Vec<[f64; 2]> = f
                                .event_pts
                                .iter()
                                .filter_map(|ep| {
                                    sample_at(&f.samples, ep[0]).map(|smp| [ep[0], q.value(&smp)])
                                })
                                .collect();
                            let now_pt = scrub.map(|p| [p.time, q.value(&p)]);
                            Plot::new("lv1_ascent_plot").height(210.0).show(ui, |pui| {
                                pui.line(Line::new(PlotPoints::from(line_pts)).name(q.axis_label()));
                                if !event_pts.is_empty() {
                                    pui.points(
                                        Points::new(PlotPoints::from(event_pts))
                                            .radius(5.0)
                                            .name("staging / MECO"),
                                    );
                                }
                                if let Some(pt) = now_pt {
                                    pui.points(
                                        Points::new(PlotPoints::from(vec![pt]))
                                            .radius(7.0)
                                            .color(egui::Color32::from_rgb(255, 196, 0))
                                            .name("now"),
                                    );
                                }
                            });
                        }
                    }
                    // ── AI optimizer — multi-objective design search ──────
                    ui.add_space(6.0);
                    ui.label(
                        egui::RichText::new("AI optimizer — multi-objective design search")
                            .strong(),
                    );
                    ui.label(
                        egui::RichText::new(
                            "searches payload × pitch-kick × vertical-rise across many real \
                             valenx-astro sims, keeping interstage SF ≥ the target below.",
                        )
                        .weak()
                        .small(),
                    );
                    ui.horizontal(|ui| {
                        ui.label("objective:");
                        ui.radio_value(
                            &mut s.opt_objective,
                            OptObjective::MaxPayload,
                            "max payload",
                        );
                        ui.radio_value(
                            &mut s.opt_objective,
                            OptObjective::MaxApoapsis,
                            "max apoapsis",
                        );
                        ui.radio_value(&mut s.opt_objective, OptObjective::MinPeakG, "min peak-g");
                    });
                    let busy = opt_running.is_some();
                    let run = ui
                        .add_enabled(
                            !busy,
                            egui::Button::new(
                                egui::RichText::new("Run AI optimization (2000 sims)").strong(),
                            ),
                        )
                        .on_hover_text(
                            "Flies 2000 candidate designs through the real ascent engine on a \
                             background thread — the UI stays responsive — and converges on the \
                             best design for the selected objective.",
                        );
                    if run.clicked() && !busy {
                        s.opt_job = Some(spawn_opt_job(s.opt_objective, s.target_sf, 2_000));
                    }
                    if let Some((p, n, obj)) = opt_running {
                        let frac = (p.done as f32 / n.max(1) as f32).clamp(0.0, 1.0);
                        // Live best-so-far in the objective's units (blank until
                        // the first feasible design is found).
                        let best = if p.best_value.is_finite() {
                            match obj {
                                OptObjective::MaxPayload => {
                                    format!(" · best {:.0} kg", p.best_value)
                                }
                                OptObjective::MaxApoapsis => {
                                    format!(" · best {:.0} km", p.best_value)
                                }
                                OptObjective::MinPeakG => format!(" · best {:.1} g", p.best_value),
                            }
                        } else {
                            String::new()
                        };
                        ui.add(egui::ProgressBar::new(frac).text(format!(
                            "{} · {}/{} sims · {} orbit{}",
                            obj.label(),
                            p.done,
                            n,
                            p.reached_orbit,
                            best,
                        )));
                    }
                    if let Some(o) = &s.opt {
                        match &o.best {
                            Some(b) => ui.label(
                                egui::RichText::new(format!(
                                    "{}: ran {} sims · {} reached orbit\n\
                                     payload {:.0} kg → {:.0} × {:.0} km\n\
                                     peak {:.1} g · SF {:.2} · pitch {:.1}° · rise {:.0} s",
                                    o.objective.label(),
                                    o.n_evals,
                                    o.reached_orbit,
                                    b.payload_kg,
                                    b.periapsis_km,
                                    b.apoapsis_km,
                                    b.peak_g,
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
                        // Plot best-so-far in the objective's own units —
                        // skip the NaN evals before the first feasible design.
                        let conv: Vec<[f64; 2]> = o
                            .convergence
                            .iter()
                            .copied()
                            .filter(|p| p[1].is_finite())
                            .collect();
                        if conv.len() > 1 {
                            let (axis, series) = match o.objective {
                                OptObjective::MaxPayload => {
                                    ("best payload (kg) vs sim #", "best payload (kg)")
                                }
                                OptObjective::MaxApoapsis => {
                                    ("best apoapsis (km) vs sim #", "best apoapsis (km)")
                                }
                                OptObjective::MinPeakG => ("best peak-g vs sim #", "best peak g"),
                            };
                            ui.label(egui::RichText::new(axis).weak().small());
                            Plot::new("lv1_opt_plot").height(170.0).show(ui, |pui| {
                                pui.line(Line::new(PlotPoints::from(conv.clone())).name(series));
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

                        // ── Lightest interstage ───────────────────────────
                        // Minimum total load-bearing area to hit the target
                        // SF (A_total = SF·F/σy — strut-count independent),
                        // plus an estimated minimum strut mass.
                        if let Some(a_total) = min_total_strut_area_m2(
                            r.peak_axial_load_n,
                            s.design.material_yield_pa,
                            s.target_sf,
                        ) {
                            let mass_min =
                                a_total * INTERSTAGE_STRUT_LENGTH_M * AL2024_DENSITY_KG_M3;
                            ui.add_space(4.0);
                            ui.label(
                                egui::RichText::new(format!(
                                    "lightest interstage @ SF {:.2}: ≥ {:.1} cm² total \
                                     (≈ {:.0} kg · {:.1} m Al-2024 struts)",
                                    s.target_sf,
                                    a_total * 1.0e4,
                                    mass_min,
                                    INTERSTAGE_STRUT_LENGTH_M,
                                ))
                                .monospace()
                                .small(),
                            )
                            .on_hover_text(
                                "Minimum total load-bearing cross-section to meet the target \
                                 safety factor — A_total = SF·F/σy, independent of strut count. \
                                 Mass assumes solid Al-2024-T3 struts of the length shown.",
                            );
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
        assert_eq!(s.opt_objective, OptObjective::MaxPayload);
        assert!((s.target_sf - 1.5).abs() < 1e-12);
    }

    #[test]
    fn optimizer_finds_a_feasible_payload_with_monotone_convergence() {
        let r = optimize_ascent(OptObjective::MaxPayload, 1.5, 120);
        assert_eq!(r.n_evals, 120);
        assert_eq!(r.convergence.len(), 120);
        assert_eq!(r.objective, OptObjective::MaxPayload);
        assert!(r.reached_orbit > 0, "some candidates should reach orbit");
        let b = r.best.expect("a feasible design is found");
        assert!((500.0..=8000.0).contains(&b.payload_kg));
        assert!(b.periapsis_km > 100.0, "periapsis {} km", b.periapsis_km);
        assert!(b.safety_factor >= 1.5, "SF {}", b.safety_factor);
        // Best-so-far payload is monotonically non-decreasing (ignoring the
        // NaN entries before the first feasible design).
        let finite: Vec<f64> = r
            .convergence
            .iter()
            .map(|p| p[1])
            .filter(|y| y.is_finite())
            .collect();
        for w in finite.windows(2) {
            assert!(w[1] >= w[0] - 1e-9, "convergence is non-decreasing");
        }
        // Deterministic (seeded): a second run gives the same best payload.
        let r2 = optimize_ascent(OptObjective::MaxPayload, 1.5, 120);
        assert_eq!(
            r2.best.map(|b| b.payload_kg as u64),
            Some(b.payload_kg as u64)
        );
    }

    #[test]
    fn each_objective_optimizes_its_own_metric() {
        // Identical seeded candidate sequence across objectives ⇒ each best is
        // a true arg-best over the SAME feasible set, so these comparisons are
        // exact (modulo float epsilon), not statistical.
        let n = 160;
        let pay = optimize_ascent(OptObjective::MaxPayload, 1.5, n)
            .best
            .expect("payload feasible");
        let apo = optimize_ascent(OptObjective::MaxApoapsis, 1.5, n)
            .best
            .expect("apoapsis feasible");
        let g = optimize_ascent(OptObjective::MinPeakG, 1.5, n)
            .best
            .expect("min-g feasible");

        // The max-payload design leads on payload; the max-apoapsis design
        // leads on apoapsis; the min-peak-g design is the gentlest ride.
        assert!(
            pay.payload_kg >= apo.payload_kg - 1e-6,
            "max-payload leads on payload: {} vs {}",
            pay.payload_kg,
            apo.payload_kg
        );
        assert!(
            apo.apoapsis_km >= pay.apoapsis_km - 1e-6,
            "max-apoapsis leads on apoapsis: {} vs {}",
            apo.apoapsis_km,
            pay.apoapsis_km
        );
        assert!(
            g.peak_g <= pay.peak_g + 1e-6,
            "min-peak-g is gentlest: {} vs {}",
            g.peak_g,
            pay.peak_g
        );

        // MinPeakG convergence is monotonically non-increasing (best = lowest).
        let r = optimize_ascent(OptObjective::MinPeakG, 1.5, n);
        assert_eq!(r.objective, OptObjective::MinPeakG);
        let finite: Vec<f64> = r
            .convergence
            .iter()
            .map(|p| p[1])
            .filter(|y| y.is_finite())
            .collect();
        assert!(finite.len() > 1, "min-g run finds feasible designs");
        for w in finite.windows(2) {
            assert!(w[1] <= w[0] + 1e-9, "min-g convergence is non-increasing");
        }
    }

    #[test]
    fn progress_callback_fires_each_eval_monotonically() {
        let mut ticks: Vec<(usize, usize)> = Vec::new();
        let r = optimize_ascent_with(OptObjective::MaxPayload, 1.5, 80, |done, reached, _best| {
            ticks.push((done, reached));
        });
        assert_eq!(ticks.len(), 80, "one progress tick per evaluation");
        assert_eq!(ticks.first().unwrap().0, 1);
        assert_eq!(ticks.last().unwrap().0, 80);
        // `done` increments by one; reached_orbit only grows; the final tick's
        // orbit count matches the result.
        for w in ticks.windows(2) {
            assert_eq!(w[1].0, w[0].0 + 1, "done increments by one");
            assert!(w[1].1 >= w[0].1, "reached_orbit only grows");
        }
        assert_eq!(ticks.last().unwrap().1, r.reached_orbit);
    }

    #[test]
    fn background_job_matches_synchronous_run() {
        // The threaded optimizer is deterministic and identical to the direct
        // call (same seed, same candidate sequence) — proves the non-blocking
        // path changes only *when* the work runs, not *what* it computes.
        let job = spawn_opt_job(OptObjective::MaxApoapsis, 1.5, 120);
        let mut result = None;
        while let Ok(msg) = job.rx.recv() {
            if let OptMsg::Done(r) = msg {
                result = Some(r);
                break;
            }
        }
        let bg = result.expect("background job sends a Done result");
        let sync = optimize_ascent(OptObjective::MaxApoapsis, 1.5, 120);
        assert_eq!(bg.objective, OptObjective::MaxApoapsis);
        assert_eq!(
            bg.best.map(|b| b.apoapsis_km as u64),
            sync.best.map(|b| b.apoapsis_km as u64),
            "threaded result equals the synchronous one"
        );
        assert_eq!(bg.reached_orbit, sync.reached_orbit);
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
    fn sample_at_clamps_and_interpolates() {
        let s = |time: f64, altitude_m: f64, dynamic_pressure: f64| TrajectorySample {
            time,
            altitude_m,
            downrange_m: 0.0,
            speed_inertial: 0.0,
            speed_relative: 0.0,
            mach: 0.0,
            mass: 0.0,
            dynamic_pressure,
            acceleration_g: 0.0,
        };
        let samples = vec![
            s(0.0, 0.0, 10.0),
            s(10.0, 1000.0, 50.0),
            s(20.0, 4000.0, 30.0),
        ];
        // Empty input → None.
        assert!(sample_at(&[], 1.0).is_none());
        // Before the window clamps to the first sample, after it to the last.
        assert!(sample_at(&samples, -5.0).unwrap().altitude_m.abs() < 1e-12);
        assert!((sample_at(&samples, 99.0).unwrap().altitude_m - 4000.0).abs() < 1e-9);
        // Midway between the first two samples is the linear midpoint.
        let mid = sample_at(&samples, 5.0).unwrap();
        assert!((mid.time - 5.0).abs() < 1e-9);
        assert!(
            (mid.altitude_m - 500.0).abs() < 1e-9,
            "alt {}",
            mid.altitude_m
        );
        assert!(
            (mid.dynamic_pressure - 30.0).abs() < 1e-9,
            "q {}",
            mid.dynamic_pressure
        );
        // Landing exactly on a recorded sample returns that sample's state.
        let on = sample_at(&samples, 20.0).unwrap();
        assert!((on.altitude_m - 4000.0).abs() < 1e-9);
    }

    #[test]
    fn fly_lv1_keeps_full_samples_for_the_scrubber() {
        let f = fly_lv1();
        assert!(
            !f.samples.is_empty(),
            "the scrubber needs the per-instant series"
        );
        assert!(f.t_max > 0.0, "t_max {}", f.t_max);
        // The scrubbed state at the end matches the last recorded sample.
        let last = *f.samples.last().unwrap();
        let at_end = sample_at(&f.samples, f.t_max).unwrap();
        assert!((at_end.altitude_m - last.altitude_m).abs() < 1e-6);
        // Mid-flight the vehicle is actually moving.
        let mid = sample_at(&f.samples, f.t_max * 0.5).unwrap();
        assert!(
            mid.speed_relative > 0.0 || mid.speed_inertial > 0.0,
            "should be moving mid-flight"
        );
    }

    #[test]
    fn plot_quantity_value_picks_the_right_channel_in_display_units() {
        let s = TrajectorySample {
            time: 1.0,
            altitude_m: 2000.0,
            downrange_m: 0.0,
            speed_inertial: 0.0,
            speed_relative: 300.0,
            mach: 0.9,
            mass: 0.0,
            dynamic_pressure: 25_000.0,
            acceleration_g: 3.0,
        };
        // Altitude and dynamic pressure are converted to km / kPa; the rest
        // are passed through in their native units.
        assert!((PlotQuantity::Altitude.value(&s) - 2.0).abs() < 1e-9);
        assert!((PlotQuantity::Speed.value(&s) - 300.0).abs() < 1e-9);
        assert!((PlotQuantity::Mach.value(&s) - 0.9).abs() < 1e-9);
        assert!((PlotQuantity::DynPressure.value(&s) - 25.0).abs() < 1e-9);
        assert!((PlotQuantity::AccelG.value(&s) - 3.0).abs() < 1e-9);
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

    #[test]
    fn min_total_area_is_strut_count_independent() {
        // A_total = SF·F/σy: 2·100/10 = 20, with no N term.
        assert_eq!(min_total_strut_area_m2(100.0, 10.0, 2.0), Some(20.0));
        // Non-positive inputs → None (no panic, no divide-by-zero).
        assert!(min_total_strut_area_m2(0.0, 10.0, 2.0).is_none());
        assert!(min_total_strut_area_m2(100.0, 0.0, 2.0).is_none());

        // The total equals N × the per-strut required area for every N — the
        // per-strut formula divides by N, the total multiplies it back, so
        // the lightest interstage is genuinely strut-count independent.
        let (load, yield_pa, target) = (123_456.0, 324.0e6, 1.8);
        let total = min_total_strut_area_m2(load, yield_pa, target).unwrap();
        for n in [1usize, 4, 8, 16, 32] {
            let per = required_area_per_strut_m2(load, yield_pa, n, target).unwrap();
            assert!(
                (per * n as f64 - total).abs() < 1e-12,
                "N={n}: {per}×{n} should equal total {total}"
            );
        }
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

    #[test]
    fn workbench_draws_with_scrubber_set_without_panic() {
        let mut app = ValenxApp::default();
        app.show_rocket_workbench = true;
        // First draw computes + caches the LV-1 flight.
        draw_workbench(&mut app);
        // Scrub to mid-flight and redraw — the readout + "now" marker render.
        let t_max = app.rocket.lv1.as_ref().map(|f| f.t_max).unwrap_or(0.0);
        app.rocket.scrub_t = t_max * 0.5;
        draw_workbench(&mut app);
        assert!(app.rocket.lv1.is_some());
    }

    #[test]
    fn workbench_draws_each_plot_channel_without_panic() {
        let mut app = ValenxApp::default();
        app.show_rocket_workbench = true;
        draw_workbench(&mut app); // first draw computes the flight
        for q in [
            PlotQuantity::Altitude,
            PlotQuantity::Speed,
            PlotQuantity::Mach,
            PlotQuantity::DynPressure,
            PlotQuantity::AccelG,
        ] {
            app.rocket.plot_y = q;
            draw_workbench(&mut app);
        }
        assert!(app.rocket.lv1.is_some());
    }
}
