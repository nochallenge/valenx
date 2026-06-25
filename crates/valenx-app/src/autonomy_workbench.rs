//! The right-side **Autonomy V&V Workbench** panel — a native front-end over
//! the in-house `valenx-autonomy-vnv` crate (scenario-based verification &
//! validation of an autonomous vehicle carrying simulated sensors, layered on
//! `valenx-sensors`).
//!
//! Mirrors the other workbenches (`ocean_workbench`, `fluids_workbench`,
//! `sensors_workbench`): a [`crate::workbench_chrome::workbench_shell`] panel
//! gated on [`crate::ValenxApp::show_autonomy_workbench`], toggled from the View
//! menu and openable by the agent bridge under the workbench id `"autonomy"`
//! (see [`crate::project_tabs::TabKind`] / [`crate::agent_commands`]).
//!
//! The user edits an **initial vehicle state** (planar position, heading,
//! speed), the **scenario duration** (time step + step count), a single
//! spherical **obstacle** (position + radius), and a few **requirement
//! thresholds** (a `MinClearance` distance, a `NoCollision` vehicle radius, and
//! a `StayInBounds` box of half-extents about the origin). Clicking **Run**
//! builds a [`Scenario`] — the vehicle coasts in a straight line at its initial
//! velocity — drives the `valenx-sensors` [`Harness`] step-by-step into a
//! [`Trace`], then scores a [`RequirementSet`] against it into a [`VnvReport`].
//!
//! The 2-D top-down view (ENU `+x` east → right, `+y` north → up) draws the
//! bounds box, the obstacle circle, and the ground-truth trajectory polyline
//! from the trace; a table lists each requirement PASS/FAIL with its signed
//! margin (green for slack, red for a violation), and the readouts show the
//! overall pass/fail, the minimum clearance achieved, and the worst (smallest)
//! margin across the requirement set.
//!
//! Honesty: this is the **V&V *methodology* / harness** over the analytic,
//! model-grade `valenx-sensors` models (a kinematic vehicle, an analytic scene,
//! graphics-grade sensors) — it validates that the *framework logic* (does a
//! violating run fail with the right margin, a safe one pass) is correct, NOT a
//! certified safety case or real-world autonomy assurance. The crates are
//! explicit about this. Every error from `valenx-autonomy-vnv` surfaces verbatim
//! in-panel — the workbench never invents a number, and degenerate parameters
//! (e.g. `steps == 0`, `dt ≤ 0`, a non-positive obstacle radius, a non-finite
//! threshold) show an in-panel error, NOT a panic.

use eframe::egui;
use nalgebra::{UnitQuaternion, Vector3};
use valenx_autonomy_vnv::{
    evaluate, run_scenario, Aabb, CommandSeq, Requirement, RequirementSet, Scenario, VnvReport,
};
// The harness/scene primitives the scenario is built from. `valenx-autonomy-vnv`
// re-exports `Command`/`VehicleState`/`Scene`, but `Sphere` lives only in
// `valenx-sensors`, so source the scene-building types from there directly
// (valenx-sensors is already a valenx-app dependency).
use valenx_sensors::{Command, Scene, Sphere, VehicleState};

use crate::agent_commands::AgentValue;
use crate::ValenxApp;

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

/// Editable autonomy-scenario + requirement parameters shown in the workbench.
#[derive(Clone, Debug)]
pub struct AutonomyParams {
    // --- Initial vehicle state (planar, in the ENU world: +x east, +y north) ---
    /// Initial east position `x` (m).
    pub start_x: f64,
    /// Initial north position `y` (m).
    pub start_y: f64,
    /// Initial heading (deg), measured from +x (east) toward +y (north). Sets
    /// both the orientation about +z and the initial velocity direction.
    pub heading_deg: f64,
    /// Initial ground speed (m/s, ≥ 0) along the heading. The vehicle coasts
    /// (zero commanded acceleration), so this is held for the whole run.
    pub speed: f64,

    // --- Scenario duration ---
    /// Per-tick time step `dt` (s) — must be finite and > 0.
    pub dt: f64,
    /// Number of ticks — must be ≥ 1. Total duration is `dt · steps`.
    pub steps: usize,

    // --- Obstacle (one sphere in the scene) ---
    /// Obstacle centre east `x` (m).
    pub obstacle_x: f64,
    /// Obstacle centre north `y` (m).
    pub obstacle_y: f64,
    /// Obstacle radius (m) — must be finite and > 0.
    pub obstacle_radius: f64,

    // --- Requirement thresholds ---
    /// `MinClearance` required distance `d` (m, ≥ 0): the vehicle must stay at
    /// least `d` clear of the obstacle surface at every tick.
    pub min_clearance: f64,
    /// `NoCollision` vehicle collision radius (m, ≥ 0): the vehicle's hull
    /// (a sphere of this radius) must never touch the obstacle.
    pub collision_radius: f64,
    /// `StayInBounds` box half-extent in `x` (m, > 0), centred on the origin.
    pub bounds_half_x: f64,
    /// `StayInBounds` box half-extent in `y` (m, > 0), centred on the origin.
    pub bounds_half_y: f64,
}

impl Default for AutonomyParams {
    fn default() -> Self {
        // A vehicle starting at the origin heading east at 3 m/s for ~2.5 s, with
        // an obstacle 12 m ahead. Default thresholds make this a *passing*
        // scenario (it travels ~7.5 m, stopping well clear of the obstacle at
        // x=12 r=1, staying inside a 20 m box).
        Self {
            start_x: 0.0,
            start_y: 0.0,
            heading_deg: 0.0,
            speed: 3.0,
            dt: 0.1,
            steps: 25,
            obstacle_x: 12.0,
            obstacle_y: 0.0,
            obstacle_radius: 1.0,
            min_clearance: 2.0,
            collision_radius: 0.5,
            bounds_half_x: 20.0,
            bounds_half_y: 20.0,
        }
    }
}

// ---------------------------------------------------------------------------
// Simulation result
// ---------------------------------------------------------------------------

/// One requirement's verdict, flattened for the painter/table (so the view
/// borrows neither the report nor the requirement set).
#[derive(Clone, Debug)]
pub struct ReqRow {
    /// Short label of the requirement (e.g. `MinClearance(d=2)`).
    pub label: String,
    /// Whether it passed (`margin >= 0`).
    pub pass: bool,
    /// Signed margin: positive slack if satisfied, negative violation if not.
    pub margin: f64,
}

/// Cached V&V output for the painter + readouts.
#[derive(Default, Clone)]
pub struct AutonomyResult {
    /// The full ground-truth trajectory (initial state + every tick), as planar
    /// `(x, y)` points (m) for the top-down polyline.
    pub path: Vec<(f64, f64)>,
    /// Per-requirement verdicts in requirement-set order.
    pub rows: Vec<ReqRow>,
    /// Overall verdict: every requirement passed.
    pub overall_pass: bool,
    /// The worst (minimum) margin across all requirements, or `None` for an
    /// empty requirement set.
    pub worst_margin: Option<f64>,
    /// The minimum clearance (m) the vehicle achieved to the obstacle surface
    /// across the run (distance to the sphere surface, clamped ≥ 0 inside).
    pub min_clearance_achieved: f64,
    /// Obstacle centre `(x, y)` (m) and radius (m), echoed for the painter.
    pub obstacle: (f64, f64, f64),
    /// The `StayInBounds` box as `(min_x, min_y, max_x, max_y)` (m).
    pub bounds: (f64, f64, f64, f64),
}

// ---------------------------------------------------------------------------
// Workbench state
// ---------------------------------------------------------------------------

/// Persistent state for the Autonomy V&V workbench.
#[derive(Default)]
pub struct AutonomyWorkbenchState {
    /// User-editable parameters.
    pub params: AutonomyParams,
    /// Last successful result (populated after a successful run).
    pub result: Option<AutonomyResult>,
    /// Status / error line shown below the controls.
    pub status: String,
}

impl AutonomyWorkbenchState {
    /// Build the [`Scenario`] from the current parameters as a fail-loud
    /// `Result`.
    ///
    /// Returns `Err` (shown in-panel) rather than panicking when the user has
    /// entered degenerate values (non-finite start/heading/speed, a non-positive
    /// obstacle radius, `dt ≤ 0`, or `steps == 0`). The vehicle is given an
    /// orientation + initial velocity from its heading and coasts (zero
    /// commanded acceleration) so the path is a clean straight line; the
    /// obstacle is a single [`Sphere`] in the [`Scene`].
    pub fn build_scenario(&self) -> Result<Scenario, String> {
        let p = &self.params;

        // Validate the planar/heading/speed inputs up front with clear messages
        // (the scenario also validates non-finite state, but we *construct* the
        // velocity/orientation here so guard the inputs we consume first).
        for (name, v) in [
            ("start x", p.start_x),
            ("start y", p.start_y),
            ("heading", p.heading_deg),
            ("speed", p.speed),
        ] {
            if !v.is_finite() {
                return Err(format!("{name} must be finite, got {v}"));
            }
        }
        if p.speed < 0.0 {
            return Err(format!("speed must be ≥ 0, got {}", p.speed));
        }

        // Heading about +z (up): orientation rotates the body +x onto the world
        // heading, and the initial velocity points along it at `speed`.
        let heading = p.heading_deg.to_radians();
        let orientation = UnitQuaternion::from_axis_angle(&Vector3::z_axis(), heading);
        let dir = Vector3::new(heading.cos(), heading.sin(), 0.0);
        let initial_state = VehicleState {
            position: Vector3::new(p.start_x, p.start_y, 0.0),
            velocity: dir * p.speed,
            orientation,
            angular_rate: Vector3::zeros(),
        };

        // One spherical obstacle. `Sphere::new` fails loud on a non-positive /
        // non-finite radius — surface that error verbatim.
        let sphere = Sphere::new(
            Vector3::new(p.obstacle_x, p.obstacle_y, 0.0),
            p.obstacle_radius,
        )
        .map_err(|e| e.to_string())?;
        let mut scene = Scene::new();
        scene.push_sphere(sphere);

        // Coast for `steps` ticks of `dt`. `CommandSeq` (via `Scenario::validate`
        // inside `run_scenario`) rejects `dt ≤ 0` / `steps == 0`, so we don't
        // pre-empt those checks — they surface with the engine's own message.
        let commands = CommandSeq::Constant {
            command: Command::coast(),
            dt: p.dt,
            steps: p.steps,
        };

        Ok(
            Scenario::new("autonomy-vnv", initial_state, scene, commands)
                .with_param("start_x", p.start_x)
                .with_param("speed", p.speed)
                .with_param("obstacle_x", p.obstacle_x),
        )
    }

    /// Build the [`RequirementSet`] (MinClearance + NoCollision + StayInBounds)
    /// from the current thresholds, as a fail-loud `Result`.
    ///
    /// The `StayInBounds` box is centred on the origin with the given
    /// half-extents (a tall `z` span so the planar run is always within it in
    /// the vertical axis). `Aabb::new` fails loud on non-finite corners, and each
    /// requirement's own `validate` (run inside [`evaluate`]) rejects a negative
    /// clearance/radius — those errors surface verbatim.
    pub fn build_requirements(&self) -> Result<RequirementSet, String> {
        let p = &self.params;
        if !(p.bounds_half_x.is_finite() && p.bounds_half_x > 0.0) {
            return Err(format!(
                "bounds half-extent x must be finite and > 0, got {}",
                p.bounds_half_x
            ));
        }
        if !(p.bounds_half_y.is_finite() && p.bounds_half_y > 0.0) {
            return Err(format!(
                "bounds half-extent y must be finite and > 0, got {}",
                p.bounds_half_y
            ));
        }
        let bounds = Aabb::new(
            Vector3::new(-p.bounds_half_x, -p.bounds_half_y, -1000.0),
            Vector3::new(p.bounds_half_x, p.bounds_half_y, 1000.0),
        )
        .map_err(|e| e.to_string())?;

        Ok(RequirementSet::new(vec![
            Requirement::MinClearance { d: p.min_clearance },
            Requirement::NoCollision {
                radius: p.collision_radius,
            },
            Requirement::StayInBounds { bounds },
        ]))
    }

    /// Run the full V&V pipeline: build the scenario, drive the harness into a
    /// [`Trace`] via [`run_scenario`], build the requirement set, and score it
    /// with [`evaluate`] into a [`VnvReport`] — then flatten everything the
    /// painter/readouts need.
    ///
    /// Every failure is returned as an `Err(String)` — no panics, no invented
    /// numbers.
    pub fn run(&self) -> Result<AutonomyResult, String> {
        let p = &self.params;

        let scenario = self.build_scenario()?;
        let requirements = self.build_requirements()?;

        // Drive the harness → trace (fails loud on a bad command sequence /
        // non-finite state).
        let trace = run_scenario(&scenario).map_err(|e| e.to_string())?;

        // Score the requirement set → report (fails loud on a setup mismatch).
        let report: VnvReport = evaluate(&requirements, &trace).map_err(|e| e.to_string())?;

        // Planar trajectory (initial + every tick).
        let path: Vec<(f64, f64)> = trace
            .states_with_initial()
            .iter()
            .map(|s| (s.position.x, s.position.y))
            .collect();

        // Minimum clearance achieved to the obstacle surface across the run
        // (distance to the sphere surface, clamped ≥ 0 inside). Computed here for
        // the readout independently of the requirement margin.
        let oc = Vector3::new(p.obstacle_x, p.obstacle_y, 0.0);
        let min_clearance_achieved = trace
            .states_with_initial()
            .iter()
            .map(|s| ((s.position - oc).norm() - p.obstacle_radius).max(0.0))
            .fold(f64::INFINITY, f64::min);
        let min_clearance_achieved = if min_clearance_achieved.is_finite() {
            min_clearance_achieved
        } else {
            0.0
        };

        let rows: Vec<ReqRow> = report
            .outcomes
            .iter()
            .map(|o| ReqRow {
                label: o.label.clone(),
                pass: o.pass,
                margin: o.margin,
            })
            .collect();

        Ok(AutonomyResult {
            path,
            rows,
            overall_pass: report.overall_pass,
            worst_margin: report.worst_margin(),
            min_clearance_achieved,
            obstacle: (p.obstacle_x, p.obstacle_y, p.obstacle_radius),
            bounds: (
                -p.bounds_half_x,
                -p.bounds_half_y,
                p.bounds_half_x,
                p.bounds_half_y,
            ),
        })
    }

    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (see [`crate::agent_commands`]). Returned by `ListControls`
    /// so an agent can discover the name space. The captions match exactly what
    /// the workbench form draws (and what each control is `labelled_by`).
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            "start x (m)",
            "start y (m)",
            "heading (deg)",
            "speed (m/s)",
            "time step dt (s)",
            "steps",
            "obstacle x (m)",
            "obstacle y (m)",
            "obstacle radius (m)",
            "min clearance d (m)",
            "collision radius (m)",
            "bounds half-x (m)",
            "bounds half-y (m)",
        ]
    }

    /// Set one labelled control by its user-visible caption, for the agent
    /// `SetControl` bridge. Fail-loud: an unknown caption or a value of the wrong
    /// type returns `Err(String)` (the bridge turns it into a `warn` feed note) —
    /// never a panic, and no field is written on error. Every numeric caption
    /// reads [`AgentValue::as_f64`]; the integer `steps` count reads
    /// [`AgentValue::as_i64`]. Range/finiteness validation stays in
    /// [`run`](Self::run) (so e.g. a negative `steps` surfaces there, not here).
    pub fn agent_set(&mut self, name: &str, value: &AgentValue) -> Result<(), String> {
        let p = &mut self.params;
        match name {
            "start x (m)" => p.start_x = value.as_f64()?,
            "start y (m)" => p.start_y = value.as_f64()?,
            "heading (deg)" => p.heading_deg = value.as_f64()?,
            "speed (m/s)" => p.speed = value.as_f64()?,
            "time step dt (s)" => p.dt = value.as_f64()?,
            "steps" => {
                let n = value.as_i64()?;
                if n < 0 {
                    return Err(format!("steps must be >= 0, got {n}"));
                }
                p.steps = n as usize;
            }
            "obstacle x (m)" => p.obstacle_x = value.as_f64()?,
            "obstacle y (m)" => p.obstacle_y = value.as_f64()?,
            "obstacle radius (m)" => p.obstacle_radius = value.as_f64()?,
            "min clearance d (m)" => p.min_clearance = value.as_f64()?,
            "collision radius (m)" => p.collision_radius = value.as_f64()?,
            "bounds half-x (m)" => p.bounds_half_x = value.as_f64()?,
            "bounds half-y (m)" => p.bounds_half_y = value.as_f64()?,
            other => return Err(format!("unknown autonomy control: {other:?}")),
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Draw the Autonomy V&V workbench. A no-op unless toggled on via
/// View → Autonomy V&V.
///
/// Mirrors [`crate::ocean_workbench::draw_ocean_workbench`].
pub fn draw_autonomy_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_autonomy_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_autonomy_workbench",
        "Autonomy V&V (scenario verification)",
        autonomy_workbench_body,
    );
    if close {
        app.show_autonomy_workbench = false;
    }
}

// ---------------------------------------------------------------------------
// Workbench body
// ---------------------------------------------------------------------------

fn autonomy_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(
            "Scenario-based V&V of a kinematic autonomous vehicle + simulated sensors · \
             valenx-autonomy-vnv  [V&V methodology over model-grade sensors — NOT a certified \
             safety case]",
        )
        .weak()
        .small(),
    );
    ui.separator();

    let mut do_run = false;

    {
        let s = &mut app.autonomy;
        let p = &mut s.params;

        ui.label(egui::RichText::new("Initial vehicle state").strong());
        egui::Grid::new("autonomy_vehicle_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let lbl = ui.label("start x (m)");
                ui.add(egui::DragValue::new(&mut p.start_x).speed(0.1))
                    .labelled_by(lbl.id)
                    .on_hover_text("Initial east position (m) in the ENU world frame.");
                ui.end_row();

                let lbl = ui.label("start y (m)");
                ui.add(egui::DragValue::new(&mut p.start_y).speed(0.1))
                    .labelled_by(lbl.id)
                    .on_hover_text("Initial north position (m) in the ENU world frame.");
                ui.end_row();

                let lbl = ui.label("heading (deg)");
                ui.add(
                    egui::DragValue::new(&mut p.heading_deg)
                        .speed(1.0)
                        .range(-360.0..=360.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Heading (deg) from +x (east) toward +y (north). Sets the initial \
                     velocity direction; the vehicle coasts straight along it.",
                );
                ui.end_row();

                let lbl = ui.label("speed (m/s)");
                ui.add(
                    egui::DragValue::new(&mut p.speed)
                        .speed(0.1)
                        .range(0.0..=200.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text("Initial ground speed (m/s) along the heading; held (coasting).");
                ui.end_row();
            });

        ui.add_space(4.0);
        ui.label(egui::RichText::new("Scenario duration").strong());
        egui::Grid::new("autonomy_time_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let lbl = ui.label("time step dt (s)");
                ui.add(
                    egui::DragValue::new(&mut p.dt)
                        .speed(0.005)
                        .range(0.001..=10.0)
                        .suffix(" s"),
                )
                .labelled_by(lbl.id)
                .on_hover_text("Per-tick integration step (s). Must be > 0.");
                ui.end_row();

                let lbl = ui.label("steps");
                ui.add(
                    egui::DragValue::new(&mut p.steps)
                        .speed(1)
                        .range(1..=100_000),
                )
                .labelled_by(lbl.id)
                .on_hover_text("Number of ticks (≥ 1). Total run duration is dt · steps.");
                ui.end_row();
            });

        ui.add_space(4.0);
        ui.label(egui::RichText::new("Obstacle (sphere)").strong());
        egui::Grid::new("autonomy_obstacle_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let lbl = ui.label("obstacle x (m)");
                ui.add(egui::DragValue::new(&mut p.obstacle_x).speed(0.1))
                    .labelled_by(lbl.id)
                    .on_hover_text("Obstacle centre east position (m).");
                ui.end_row();

                let lbl = ui.label("obstacle y (m)");
                ui.add(egui::DragValue::new(&mut p.obstacle_y).speed(0.1))
                    .labelled_by(lbl.id)
                    .on_hover_text("Obstacle centre north position (m).");
                ui.end_row();

                let lbl = ui.label("obstacle radius (m)");
                ui.add(
                    egui::DragValue::new(&mut p.obstacle_radius)
                        .speed(0.05)
                        .range(0.001..=1000.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text("Obstacle sphere radius (m). Must be > 0.");
                ui.end_row();
            });

        ui.add_space(4.0);
        ui.label(egui::RichText::new("Requirements (acceptance thresholds)").strong());
        egui::Grid::new("autonomy_requirement_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let lbl = ui.label("min clearance d (m)");
                ui.add(
                    egui::DragValue::new(&mut p.min_clearance)
                        .speed(0.05)
                        .range(0.0..=1000.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "MinClearance: the vehicle must stay at least this far (m) clear of the \
                     obstacle surface at every tick.",
                );
                ui.end_row();

                let lbl = ui.label("collision radius (m)");
                ui.add(
                    egui::DragValue::new(&mut p.collision_radius)
                        .speed(0.05)
                        .range(0.0..=1000.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "NoCollision: the vehicle hull (a sphere of this radius, m) must never touch \
                     the obstacle.",
                );
                ui.end_row();

                let lbl = ui.label("bounds half-x (m)");
                ui.add(
                    egui::DragValue::new(&mut p.bounds_half_x)
                        .speed(0.5)
                        .range(0.01..=100_000.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "StayInBounds: half-width (m) of the allowed box about the origin in x.",
                );
                ui.end_row();

                let lbl = ui.label("bounds half-y (m)");
                ui.add(
                    egui::DragValue::new(&mut p.bounds_half_y)
                        .speed(0.5)
                        .range(0.01..=100_000.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "StayInBounds: half-height (m) of the allowed box about the origin in y.",
                );
                ui.end_row();
            });

        ui.add_space(6.0);
        ui.horizontal(|ui| {
            if ui
                .button(egui::RichText::new("Run").strong())
                .on_hover_text("Drive the scenario into a trace and score the requirements.")
                .clicked()
            {
                do_run = true;
            }
        });
    }

    // --- Execute (outside borrow) -------------------------------------------
    if do_run {
        run_and_store(app);
    }

    // --- Status line ---------------------------------------------------------
    let s = &app.autonomy;
    if !s.status.is_empty() {
        ui.add_space(6.0);
        let color = if s.status.starts_with('\u{26A0}') {
            egui::Color32::from_rgb(220, 120, 60)
        } else {
            egui::Color32::from_rgb(90, 180, 110)
        };
        ui.label(egui::RichText::new(&s.status).color(color).strong());
    }

    // --- Visualisation -------------------------------------------------------
    ui.add_space(6.0);
    ui.separator();
    draw_autonomy_viz(s, ui);
}

/// Run the pipeline and fold the result (or error) into the workbench status.
/// Factored out so the Run button (and tests) can use it.
pub(crate) fn run_and_store(app: &mut ValenxApp) {
    let s = &mut app.autonomy;
    match s.run() {
        Ok(res) => {
            let verdict = if res.overall_pass { "PASS" } else { "FAIL" };
            s.status = format!(
                "\u{2714} {verdict} · {} / {} requirements met · min clearance {:.2} m · worst margin {}",
                res.rows.iter().filter(|r| r.pass).count(),
                res.rows.len(),
                res.min_clearance_achieved,
                res.worst_margin
                    .map(|m| format!("{m:+.2}"))
                    .unwrap_or_else(|| "—".to_string()),
            );
            s.result = Some(res);
        }
        Err(e) => {
            s.status = format!("\u{26A0} {e}");
            s.result = None;
        }
    }
}

// ---------------------------------------------------------------------------
// 2-D top-down visualisation (scene + trajectory) + requirements table
// ---------------------------------------------------------------------------

fn draw_autonomy_viz(s: &AutonomyWorkbenchState, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Top-down view — scene + ground-truth trajectory").strong());
    ui.label(
        egui::RichText::new(
            "grey box = StayInBounds · orange circle = obstacle · cyan path = vehicle trajectory \
             (+x east → right, +y north → up)",
        )
        .weak()
        .small(),
    );

    let available = ui.available_size();
    let (rect, _) = ui.allocate_exact_size(
        egui::vec2(available.x.min(460.0), 300.0),
        egui::Sense::hover(),
    );

    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(14, 22, 34));

    let Some(res) = &s.result else {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "press \"Run\" to verify the scenario",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    };

    // World extent to show: cover the bounds box, the obstacle (with its
    // radius), and the whole trajectory, with a small margin so nothing clips.
    let (bx0, by0, bx1, by1) = res.bounds;
    let (ox, oy, orad) = res.obstacle;
    let mut min_x = bx0.min(ox - orad);
    let mut max_x = bx1.max(ox + orad);
    let mut min_y = by0.min(oy - orad);
    let mut max_y = by1.max(oy + orad);
    for &(x, y) in &res.path {
        min_x = min_x.min(x);
        max_x = max_x.max(x);
        min_y = min_y.min(y);
        max_y = max_y.max(y);
    }
    // Pad by 5% of the larger span (and guard a degenerate zero span).
    let span = ((max_x - min_x).max(max_y - min_y)).max(1.0);
    let pad = span * 0.05;
    min_x -= pad;
    max_x += pad;
    min_y -= pad;
    max_y += pad;
    let span_x = (max_x - min_x).max(f64::EPSILON);
    let span_y = (max_y - min_y).max(f64::EPSILON);

    let inner = rect.shrink(14.0);
    // Map world (x, y) → painter pixel. World +y (north) is up, so it maps to a
    // *smaller* screen-y. A single uniform scale keeps the aspect ratio square.
    let scale = (inner.width() as f64 / span_x).min(inner.height() as f64 / span_y);
    let to_px = |x: f64, y: f64| -> egui::Pos2 {
        let px = inner.left() as f64 + (x - min_x) * scale;
        let py = inner.bottom() as f64 - (y - min_y) * scale;
        egui::pos2(px as f32, py as f32)
    };

    // StayInBounds box.
    let box_rect = egui::Rect::from_two_pos(to_px(bx0, by0), to_px(bx1, by1));
    painter.rect_stroke(
        box_rect,
        0.0,
        egui::Stroke::new(1.5, egui::Color32::from_gray(110)),
    );

    // Obstacle circle (radius scaled to pixels via the uniform scale).
    let obstacle_center = to_px(ox, oy);
    let obstacle_r_px = (orad * scale) as f32;
    painter.circle_filled(
        obstacle_center,
        obstacle_r_px.max(1.0),
        egui::Color32::from_rgba_unmultiplied(220, 140, 60, 90),
    );
    painter.circle_stroke(
        obstacle_center,
        obstacle_r_px.max(1.0),
        egui::Stroke::new(1.5, egui::Color32::from_rgb(220, 140, 60)),
    );

    // Trajectory polyline.
    if res.path.len() >= 2 {
        let pts: Vec<egui::Pos2> = res.path.iter().map(|&(x, y)| to_px(x, y)).collect();
        // Colour the path by the overall verdict: green pass, red fail.
        let col = if res.overall_pass {
            egui::Color32::from_rgb(110, 200, 130)
        } else {
            egui::Color32::from_rgb(230, 90, 90)
        };
        painter.add(egui::Shape::line(pts.clone(), egui::Stroke::new(2.0, col)));
        // Start (filled) + end (ring) markers.
        if let Some(first) = pts.first() {
            painter.circle_filled(*first, 3.5, egui::Color32::from_rgb(120, 190, 240));
        }
        if let Some(last) = pts.last() {
            painter.circle_stroke(
                *last,
                3.5,
                egui::Stroke::new(2.0, egui::Color32::from_rgb(120, 190, 240)),
            );
        }
    } else if let Some(&(x, y)) = res.path.first() {
        painter.circle_filled(to_px(x, y), 3.5, egui::Color32::from_rgb(120, 190, 240));
    }

    // --- Requirements table --------------------------------------------------
    ui.add_space(6.0);
    ui.label(egui::RichText::new("Requirements").strong());
    let pass_col = egui::Color32::from_rgb(90, 190, 110);
    let fail_col = egui::Color32::from_rgb(230, 90, 90);
    egui::Grid::new("autonomy_req_table")
        .num_columns(3)
        .striped(true)
        .show(ui, |ui| {
            ui.label(egui::RichText::new("requirement").weak());
            ui.label(egui::RichText::new("verdict").weak());
            ui.label(egui::RichText::new("margin").weak());
            ui.end_row();
            for row in &res.rows {
                let col = if row.pass { pass_col } else { fail_col };
                ui.label(&row.label);
                ui.label(
                    egui::RichText::new(if row.pass { "PASS" } else { "FAIL" })
                        .color(col)
                        .strong(),
                );
                ui.label(egui::RichText::new(format!("{:+.3}", row.margin)).color(col));
                ui.end_row();
            }
        });

    // --- Readouts ------------------------------------------------------------
    ui.add_space(4.0);
    egui::Grid::new("autonomy_stats")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            let row = |ui: &mut egui::Ui, k: &str, v: String| {
                ui.label(k);
                ui.label(v);
                ui.end_row();
            };
            ui.label("overall verdict");
            ui.label(
                egui::RichText::new(if res.overall_pass { "PASS" } else { "FAIL" })
                    .color(if res.overall_pass { pass_col } else { fail_col })
                    .strong(),
            );
            ui.end_row();
            row(
                ui,
                "min clearance achieved (m)",
                format!("{:.3}", res.min_clearance_achieved),
            );
            row(
                ui,
                "worst margin",
                res.worst_margin
                    .map(|m| format!("{m:+.3}"))
                    .unwrap_or_else(|| "— (no requirements)".to_string()),
            );
            row(
                ui,
                "requirements met",
                format!(
                    "{} / {}",
                    res.rows.iter().filter(|r| r.pass).count(),
                    res.rows.len()
                ),
            );
        });
}

// ---------------------------------------------------------------------------
// Tests (unit + headless_ui_tests, mirroring ocean_workbench)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_run_passes_and_path_is_populated() {
        let s = AutonomyWorkbenchState::default();
        let res = s.run().expect("default autonomy run should succeed");
        // initial state + `steps` ticks.
        assert_eq!(
            res.path.len(),
            s.params.steps + 1,
            "path = initial state + one point per tick"
        );
        assert_eq!(
            res.rows.len(),
            3,
            "MinClearance + NoCollision + StayInBounds"
        );
        assert!(
            res.overall_pass,
            "default scenario stops well clear of the obstacle and stays in bounds"
        );
        // Every requirement passes with non-negative margin.
        assert!(res.rows.iter().all(|r| r.pass && r.margin >= 0.0));
        assert!(res.worst_margin.unwrap() >= 0.0);
    }

    #[test]
    fn collision_scenario_fails() {
        // Aim the vehicle straight at the obstacle and run long enough / fast
        // enough to drive into it: the NoCollision + MinClearance margins go
        // negative and the overall verdict fails.
        let mut s = AutonomyWorkbenchState::default();
        s.params.speed = 5.0;
        s.params.steps = 60; // travels 30 m past the obstacle at x=12
        let res = s
            .run()
            .expect("run should still succeed (it's a *failing* verdict, not an error)");
        assert!(!res.overall_pass, "driving into the obstacle must FAIL");
        // The worst margin is negative (a real violation).
        assert!(
            res.worst_margin.unwrap() < 0.0,
            "worst margin should be negative on collision, got {:?}",
            res.worst_margin
        );
        // The min clearance achieved is ~0 (it passed through the sphere).
        assert!(res.min_clearance_achieved < 1e-6);
    }

    #[test]
    fn out_of_bounds_scenario_fails_stay_in_bounds() {
        // Shrink the bounds box so the (passing-clearance) run leaves it.
        let mut s = AutonomyWorkbenchState::default();
        s.params.bounds_half_x = 2.0; // vehicle travels to ~7.5 m east → out
        let res = s.run().expect("run should succeed with a failing verdict");
        assert!(
            !res.overall_pass,
            "leaving the bounds box must FAIL overall"
        );
        let sib = res
            .rows
            .iter()
            .find(|r| r.label.starts_with("StayInBounds"))
            .expect("a StayInBounds row");
        assert!(!sib.pass, "StayInBounds must fail");
        assert!(sib.margin < 0.0, "its margin must be negative");
    }

    #[test]
    fn margin_sign_convention_holds() {
        // For every row, pass iff margin >= 0 (the crate's uniform convention,
        // surfaced faithfully).
        let s = AutonomyWorkbenchState::default();
        let res = s.run().expect("run");
        for r in &res.rows {
            assert_eq!(r.pass, r.margin >= 0.0, "pass must equal margin >= 0");
        }
    }

    // ---- degenerate-param tests — must return Err, NOT panic ----

    #[test]
    fn zero_steps_returns_err() {
        let mut s = AutonomyWorkbenchState::default();
        s.params.steps = 0;
        assert!(s.run().is_err(), "steps = 0 must return Err, not panic");
    }

    #[test]
    fn non_positive_dt_returns_err() {
        let mut s = AutonomyWorkbenchState::default();
        s.params.dt = 0.0;
        assert!(s.run().is_err(), "dt = 0 must return Err, not panic");
        s.params.dt = -0.1;
        assert!(s.run().is_err(), "dt < 0 must return Err, not panic");
    }

    #[test]
    fn non_positive_obstacle_radius_returns_err() {
        let mut s = AutonomyWorkbenchState::default();
        s.params.obstacle_radius = 0.0;
        assert!(
            s.run().is_err(),
            "obstacle radius = 0 must return Err, not panic"
        );
        s.params.obstacle_radius = -1.0;
        assert!(
            s.run().is_err(),
            "obstacle radius < 0 must return Err, not panic"
        );
    }

    #[test]
    fn non_finite_start_returns_err() {
        let mut s = AutonomyWorkbenchState::default();
        s.params.start_x = f64::NAN;
        assert!(s.run().is_err(), "start x = NaN must return Err, not panic");
    }

    #[test]
    fn negative_speed_returns_err() {
        let mut s = AutonomyWorkbenchState::default();
        s.params.speed = -1.0;
        assert!(s.run().is_err(), "speed < 0 must return Err, not panic");
    }

    #[test]
    fn non_positive_bounds_returns_err() {
        let mut s = AutonomyWorkbenchState::default();
        s.params.bounds_half_x = 0.0;
        assert!(
            s.run().is_err(),
            "bounds half-x = 0 must return Err, not panic"
        );
    }

    #[test]
    fn negative_min_clearance_returns_err() {
        // A negative required clearance is invalid config — the requirement's
        // own validate (inside evaluate) rejects it. Must be Err, not panic.
        let mut s = AutonomyWorkbenchState::default();
        s.params.min_clearance = -1.0;
        assert!(
            s.run().is_err(),
            "negative min clearance must return Err, not panic"
        );
    }

    // ---- agent_set / agent_control_names (the SetControl bridge) ----

    #[test]
    fn agent_set_sets_params_and_rejects_unknown_and_typemismatch() {
        let mut s = AutonomyWorkbenchState::default();

        // A representative float param, verified via state.
        s.agent_set("speed (m/s)", &AgentValue::Float(7.5))
            .expect("set speed");
        assert!((s.params.speed - 7.5).abs() < 1e-12);
        // The integer step-count (an Int and a whole-valued Float both work).
        s.agent_set("steps", &AgentValue::Int(40))
            .expect("set steps");
        assert_eq!(s.params.steps, 40);
        s.agent_set("steps", &AgentValue::Float(50.0))
            .expect("whole float -> usize");
        assert_eq!(s.params.steps, 50);
        // An Int widens to an f64 field.
        s.agent_set("obstacle radius (m)", &AgentValue::Int(2))
            .expect("set obstacle radius");
        assert!((s.params.obstacle_radius - 2.0).abs() < 1e-12);

        // Unknown caption -> Err (not a panic).
        assert!(s.agent_set("nope", &AgentValue::Int(1)).is_err());
        // Type mismatch: a numeric caption fed a string -> Err.
        assert!(s
            .agent_set("speed (m/s)", &AgentValue::Str("fast".into()))
            .is_err());
        // Type mismatch: the integer caption fed a fractional float -> Err.
        assert!(s.agent_set("steps", &AgentValue::Float(3.5)).is_err());
        // A negative steps value -> Err (not a silently-wrapped usize).
        assert!(s.agent_set("steps", &AgentValue::Int(-1)).is_err());

        // Every advertised control name is settable with a value of its type.
        for name in AutonomyWorkbenchState::agent_control_names() {
            let v = if *name == "steps" {
                AgentValue::Int(5)
            } else {
                AgentValue::Float(1.0)
            };
            assert!(
                s.agent_set(name, &v).is_ok(),
                "advertised control '{name}' must be settable"
            );
        }
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_autonomy_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    fn has_named_node(nodes: &[(NodeId, Node)], name: &str) -> bool {
        nodes.iter().any(|(_, n)| n.name() == Some(name))
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_autonomy_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_autonomy_workbench(&mut app, ctx);
        });
        // No panic = pass.
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_autonomy_workbench = true;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_passing_result_without_panic() {
        let mut app = ValenxApp::default();
        app.show_autonomy_workbench = true;
        let res = app.autonomy.run().expect("default run should pass");
        assert!(res.overall_pass, "default scenario should pass");
        app.autonomy.result = Some(res);
        app.autonomy.status = "\u{2714} PASS test result".to_string();
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_failing_collision_result_without_panic() {
        let mut app = ValenxApp::default();
        app.show_autonomy_workbench = true;
        // A collision scenario → a failing (but valid) report.
        app.autonomy.params.speed = 5.0;
        app.autonomy.params.steps = 60;
        let res = app.autonomy.run().expect("failing run is still Ok");
        assert!(!res.overall_pass, "collision scenario should fail");
        app.autonomy.result = Some(res);
        app.autonomy.status = "\u{2714} FAIL test result".to_string();
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_error_status_without_panic() {
        let mut app = ValenxApp::default();
        app.show_autonomy_workbench = true;
        // Trigger an error state (steps = 0 → Err).
        app.autonomy.params.steps = 0;
        let result = app.autonomy.run();
        app.autonomy.status = match result {
            Err(e) => format!("\u{26A0} {e}"),
            Ok(_) => "\u{26A0} simulated error for testing".to_string(),
        };
        app.autonomy.result = None;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn numeric_controls_are_labelled_by_named() {
        let mut app = ValenxApp::default();
        app.show_autonomy_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        // 4 vehicle + 2 duration + 3 obstacle + 4 requirement = 13 DragValues,
        // all exposed as SpinButton nodes that MUST carry an accessible name.
        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 13,
            "expected at least 13 numeric controls (DragValues), got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by a caption (AI-drivable name)"
        );

        // Check specific captions are present as named accessibility nodes.
        for caption in [
            "start x (m)",
            "start y (m)",
            "heading (deg)",
            "speed (m/s)",
            "time step dt (s)",
            "steps",
            "obstacle x (m)",
            "obstacle y (m)",
            "obstacle radius (m)",
            "min clearance d (m)",
            "collision radius (m)",
            "bounds half-x (m)",
            "bounds half-y (m)",
        ] {
            assert!(
                has_named_node(&nodes, caption),
                "caption '{caption}' must be a named node in the a11y tree"
            );
        }

        // The Run button must be a named, invokable node.
        assert!(
            nodes.iter().any(|(_, n)| {
                n.role() == Role::Button && n.name().is_some_and(|s| s.contains("Run"))
            }),
            "the Run button must be a named, invokable node"
        );
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Every numeric DragValue is a SpinButton and must be `labelled_by` its
        // caption (egui clears a DragValue's own Name), so an AI / screen reader
        // can find the control by its caption text. Beyond merely being
        // non-empty, each `labelled_by` target must RESOLVE to the caption node
        // — i.e. the spin button is correctly associated with a real named
        // label, not a dangling id.
        let mut app = ValenxApp::default();
        app.show_autonomy_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let by_id: std::collections::HashMap<NodeId, &Node> =
            nodes.iter().map(|(id, n)| (*id, n)).collect();

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // 4 vehicle + 2 duration + 3 obstacle + 4 requirement = 13 DragValues,
        // all unconditionally rendered (no mode gating).
        assert!(
            spin_buttons.len() >= 13,
            "expected the numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by its caption (AI-drivable name)"
        );
        // Each spin button's labelled_by target resolves to a named caption node.
        assert!(
            spin_buttons.iter().all(|n| {
                n.labelled_by()
                    .iter()
                    .any(|id| by_id.get(id).is_some_and(|t| t.name().is_some()))
            }),
            "every DragValue's labelled_by must point at a named caption node"
        );
        // A couple of captions must exist as named nodes in the a11y tree.
        for caption in ["start x (m)", "speed (m/s)", "steps"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }

    #[test]
    fn degenerate_steps_shows_error_not_panic() {
        // When steps == 0 the workbench must surface the error in-panel, not
        // panic.
        let mut state = AutonomyWorkbenchState::default();
        state.params.steps = 0;
        assert!(
            state.run().is_err(),
            "steps = 0 must produce Err, not panic"
        );
    }

    #[test]
    fn agent_bridge_autonomy_id_resolves_and_sets_flag() {
        // Verify the two mechanisms the agent bridge uses for
        //   `OpenWorkbench { id: "autonomy" }`:
        //   1. TabKind::from_id("autonomy") → Some(TabKind::Autonomy)
        //   2. set_workbench_flag(app, "autonomy", true) → show_autonomy_workbench
        use crate::project_tabs::{set_workbench_flag, TabKind};

        // 1. Lookup.
        assert_eq!(
            TabKind::from_id("autonomy"),
            Some(TabKind::Autonomy),
            "\"autonomy\" must resolve to TabKind::Autonomy"
        );
        // Case-insensitive + whitespace-tolerant.
        assert_eq!(TabKind::from_id("AUTONOMY"), Some(TabKind::Autonomy));
        assert_eq!(TabKind::from_id("  autonomy  "), Some(TabKind::Autonomy));

        // 2. Flag toggle.
        let mut app = ValenxApp::default();
        assert!(!app.show_autonomy_workbench);
        set_workbench_flag(&mut app, "autonomy", true);
        assert!(
            app.show_autonomy_workbench,
            "set_workbench_flag(\"autonomy\", true) must set show_autonomy_workbench"
        );
        set_workbench_flag(&mut app, "autonomy", false);
        assert!(!app.show_autonomy_workbench);
    }
}
