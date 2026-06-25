//! The right-side **Aerodynamics / Wind Tunnel** workbench panel.
//!
//! [`valenx_aero`] is a native 3-D external-
//! aerodynamics CFD engine (a virtual wind tunnel: immersed-boundary
//! incompressible RANS, k-ε / k-ω SST turbulence, drag / lift / moment
//! coefficients, Cp / velocity / vorticity fields, an angle-of-attack
//! polar sweep). It shipped as a library + agent API with **no UI**.
//!
//! This module is that UI — a polished egui side panel that feels like
//! commercial wind-tunnel software. It mirrors the CAD-side
//! [`crate::mesh_toolbox`] and the [`crate::genetics_workbench`]
//! idioms: a resizable right-hand [`egui::SidePanel`], toggled from the
//! View menu, off by default. The panel walks a clear wind-tunnel
//! workflow in eight collapsible sections — Body, Wind conditions,
//! Ground & wheels, Tunnel & mesh, Solver, Run, Results, Flow
//! visualization — and the long solve runs on a background thread so
//! the window never freezes.
//!
//! The workbench owns one [`AeroWorkbenchState`] (a field on
//! [`crate::ValenxApp`]); the real work is split across the
//! [`crate::aero`] sub-modules.

use eframe::egui;

use crate::aero::compute::{AeroRunHandle, PolarSweepResult};
use crate::aero::model::{
    BodySource, CutAxis, FlowField, GridResolution, RunMode, TurbulenceChoice, WindTunnelForm,
};
use crate::aero::panels;
use crate::ValenxApp;
use valenx_aero::{AeroReport, AeroResult};

/// All Wind-Tunnel workbench form + result state.
///
/// One instance lives on [`crate::ValenxApp`] (the `aero` field),
/// exactly as the CAD-side `MeshToolboxState` and the
/// `GeneticsWorkbenchState` do.
#[derive(Default)]
pub struct AeroWorkbenchState {
    /// Every wind-tunnel form input — see [`WindTunnelForm`].
    pub form: WindTunnelForm,

    /// The live background run, if one is in flight. `None` when idle.
    pub run: Option<AeroRunHandle>,
    /// A coarse status line for the Run section (the current stage, or
    /// the last completion message).
    pub status: String,
    /// The convergence history of the current / last run — `(iteration,
    /// residual)` pairs feeding the live residual plot.
    pub residual_history: Vec<(f64, f64)>,
    /// Sweep points accumulated live during an angle-of-attack run —
    /// `(angle_deg, cd, cl)`.
    pub sweep_progress: Vec<(f64, f64, f64)>,

    /// The last completed steady solve, if any.
    pub last_result: Option<Box<AeroResult>>,
    /// The human-readable report for [`Self::last_result`].
    pub last_report: Option<Box<AeroReport>>,
    /// The last completed angle-of-attack polar, if any.
    pub last_polar: Option<PolarSweepResult>,
    /// A short label for the body the last run tested.
    pub last_body_label: String,
    /// The last error message, shown in red. Cleared on a new run.
    pub error: Option<String>,
    /// Status line from the last flow-visualization push (which field
    /// was sent to the viewport, or why it could not be).
    pub viz_status: Option<String>,
    /// Which flow field is currently pushed into the 3-D viewport, if
    /// any — drives the "active overlay" read-out + the clear button.
    pub last_field_overlay: Option<crate::aero::model::FlowField>,
    /// Undo / redo over the wind-tunnel form. A snapshot lands on the
    /// stack when the user presses Run, so `Ctrl+Z` reverses the
    /// settings of the last completed solve.
    pub history: crate::undo::History<WindTunnelForm>,
}

/// Parse a test-body name (agent `SetControl` bridge) into a [`BodySource`].
/// Case-insensitive; accepts the combo words + short forms. Fail-loud.
fn parse_body_source(s: &str) -> Result<BodySource, String> {
    match s.trim().to_ascii_lowercase().as_str() {
        "current cad model" | "cad" | "current" => Ok(BodySource::CurrentCadModel),
        "imported stl file" | "stl" | "imported stl" => Ok(BodySource::ImportedStl),
        "built-in demo box" | "demo box" | "demo" | "box" => Ok(BodySource::DemoBox),
        other => Err(format!(
            "unknown test body '{other}' (expected 'current cad model', 'imported stl file', \
             or 'demo box')"
        )),
    }
}

/// Parse a grid-resolution name into a [`GridResolution`]. Case-insensitive.
fn parse_grid_resolution(s: &str) -> Result<GridResolution, String> {
    match s.trim().to_ascii_lowercase().as_str() {
        "coarse" | "coarse (fast)" | "fast" => Ok(GridResolution::Coarse),
        "medium" | "medium (balanced)" | "balanced" => Ok(GridResolution::Medium),
        "fine" | "fine (accurate)" | "accurate" => Ok(GridResolution::Fine),
        other => Err(format!(
            "unknown grid resolution '{other}' (expected 'coarse', 'medium', or 'fine')"
        )),
    }
}

/// Parse a turbulence-model name into a [`TurbulenceChoice`]. Case-insensitive;
/// accepts both the hyphen + greek-letter spellings.
fn parse_turbulence(s: &str) -> Result<TurbulenceChoice, String> {
    match s.trim().to_ascii_lowercase().as_str() {
        "k-epsilon" | "kepsilon" | "k-e" | "ke" => Ok(TurbulenceChoice::KEpsilon),
        "k-omega sst" | "k-omega" | "komega sst" | "sst" | "kw" => Ok(TurbulenceChoice::KOmegaSST),
        other => Err(format!(
            "unknown turbulence model '{other}' (expected 'k-epsilon' or 'k-omega SST')"
        )),
    }
}

/// Parse a run-mode name into a [`RunMode`]. Case-insensitive.
fn parse_run_mode(s: &str) -> Result<RunMode, String> {
    match s.trim().to_ascii_lowercase().as_str() {
        "steady point" | "steady" | "point" => Ok(RunMode::SteadyPoint),
        "angle-of-attack sweep" | "angle sweep" | "sweep" | "aoa sweep" => Ok(RunMode::AngleSweep),
        other => Err(format!(
            "unknown run mode '{other}' (expected 'steady point' or 'angle-of-attack sweep')"
        )),
    }
}

/// Parse a flow-field name into a [`FlowField`]. Case-insensitive; accepts the
/// button labels + short forms.
fn parse_flow_field(s: &str) -> Result<FlowField, String> {
    match s.trim().to_ascii_lowercase().as_str() {
        "surface cp" | "cp" => Ok(FlowField::SurfaceCp),
        "skin friction cf" | "skin friction" | "cf" => Ok(FlowField::SkinFriction),
        "velocity magnitude" | "velocity" | "|u|" => Ok(FlowField::VelocityMagnitude),
        "pressure cut plane" | "pressure slice" | "pressure" => Ok(FlowField::PressureSlice),
        "wake / q-criterion" | "q-criterion" | "vortex q" | "wake" | "q" => Ok(FlowField::VortexQ),
        other => Err(format!(
            "unknown flow field '{other}' (expected 'surface cp', 'skin friction cf', \
             'velocity magnitude', 'pressure cut plane', or 'wake / q-criterion')"
        )),
    }
}

/// Parse a cut-plane-axis name into a [`CutAxis`]. Case-insensitive; accepts the
/// bare axis letter or the combo words.
fn parse_cut_axis(s: &str) -> Result<CutAxis, String> {
    match s.trim().to_ascii_lowercase().as_str() {
        "x" | "x (frontal)" | "frontal" => Ok(CutAxis::X),
        "y" | "y (centreline)" | "centreline" => Ok(CutAxis::Y),
        "z" | "z (plan)" | "plan" => Ok(CutAxis::Z),
        other => Err(format!(
            "unknown cut-plane axis '{other}' (expected 'x', 'y', or 'z')"
        )),
    }
}

impl AeroWorkbenchState {
    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (see [`crate::agent_commands`]), matching exactly what the
    /// wind-tunnel form draws (and what each control is `labelled_by`). The six
    /// enum selectors (`Test body`, `Grid resolution`, `Turbulence model`,
    /// `Run mode`, `Field`, `Cut-plane axis`) are set by option name; the three
    /// checkbox/radio toggles by boolean; the rest are numeric. (The demo-box
    /// size has three per-axis captions; they only affect the demo-box body.)
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            // -- Body --
            "Test body",
            "length X",
            "width Y",
            "height Z",
            // -- Wind conditions --
            "Free-stream speed",
            "Yaw (sideslip)",
            "Pitch (angle of attack)",
            "Air density",
            "Air viscosity",
            "Turbulence intensity",
            "Air temperature",
            "Aspect ratio  AR",
            "Span efficiency  e",
            "Apply Prandtl-Glauert compressibility correction",
            // -- Ground & wheels --
            "moving ground",
            "Show rolling-wheel speed (reference)",
            "Wheel radius",
            // -- Tunnel & mesh --
            "Grid resolution",
            "Override the auto-sized domain",
            "Upstream clearance",
            "Downstream (wake)",
            "Lateral clearance",
            // -- Solver --
            "Turbulence model",
            "Max iterations",
            "Convergence tolerance",
            "Run mode",
            "Start angle",
            "End angle",
            "Sweep points",
            // -- Flow visualization --
            "Field",
            "Cut-plane axis",
            "Cut position",
        ]
    }

    /// Set one labelled control by its user-visible caption, for the agent
    /// `SetControl` bridge. The six enum captions read [`AgentValue::as_str`] and
    /// match a small option-name set; the toggles read [`AgentValue::as_bool`];
    /// numeric fields read [`AgentValue::as_f64`] / [`AgentValue::as_i64`].
    ///
    /// Fail-loud: an unknown caption, a wrong type, an unrecognised enum name, or
    /// a negative count returns `Err` (the bridge turns it into a `warn` note) —
    /// never a panic, no field written on error. Out-of-range numbers are caught
    /// at solve time, which surfaces an in-panel error rather than a panic.
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        let f = &mut self.form;
        match name {
            // -- Body --
            "Test body" => f.body_source = parse_body_source(value.as_str()?)?,
            "length X" => f.demo_box_size[0] = value.as_f64()?,
            "width Y" => f.demo_box_size[1] = value.as_f64()?,
            "height Z" => f.demo_box_size[2] = value.as_f64()?,
            // -- Wind conditions --
            "Free-stream speed" => f.speed_ms = value.as_f64()?,
            "Yaw (sideslip)" => f.yaw_deg = value.as_f64()?,
            "Pitch (angle of attack)" => f.pitch_deg = value.as_f64()?,
            "Air density" => f.air_density = value.as_f64()?,
            "Air viscosity" => f.air_viscosity = value.as_f64()?,
            "Turbulence intensity" => f.turbulence_intensity = value.as_f64()?,
            "Air temperature" => f.temperature_k = value.as_f64()?,
            "Aspect ratio  AR" => f.aspect_ratio = value.as_f64()?,
            "Span efficiency  e" => f.span_efficiency = value.as_f64()?,
            "Apply Prandtl-Glauert compressibility correction" => {
                f.apply_compressibility = value.as_bool()?
            }
            // -- Ground & wheels --
            "moving ground" => f.moving_ground = value.as_bool()?,
            "Show rolling-wheel speed (reference)" => f.rotating_wheels = value.as_bool()?,
            "Wheel radius" => f.wheel_radius = value.as_f64()?,
            // -- Tunnel & mesh --
            "Grid resolution" => f.resolution = parse_grid_resolution(value.as_str()?)?,
            "Override the auto-sized domain" => f.override_domain = value.as_bool()?,
            "Upstream clearance" => f.upstream_lengths = value.as_f64()?,
            "Downstream (wake)" => f.downstream_lengths = value.as_f64()?,
            "Lateral clearance" => f.lateral_lengths = value.as_f64()?,
            // -- Solver --
            "Turbulence model" => f.turbulence = parse_turbulence(value.as_str()?)?,
            "Max iterations" => {
                let n = value.as_i64()?;
                if n < 0 {
                    return Err(format!("Max iterations must be >= 0, got {n}"));
                }
                f.max_iterations = n as usize;
            }
            "Convergence tolerance" => f.tolerance = value.as_f64()?,
            "Run mode" => f.run_mode = parse_run_mode(value.as_str()?)?,
            "Start angle" => f.sweep_start_deg = value.as_f64()?,
            "End angle" => f.sweep_end_deg = value.as_f64()?,
            "Sweep points" => {
                let n = value.as_i64()?;
                if n < 0 {
                    return Err(format!("Sweep points must be >= 0, got {n}"));
                }
                f.sweep_points = n as usize;
            }
            // -- Flow visualization --
            "Field" => f.flow_field = parse_flow_field(value.as_str()?)?,
            "Cut-plane axis" => f.cut_axis = parse_cut_axis(value.as_str()?)?,
            "Cut position" => f.cut_fraction = value.as_f64()?,
            other => return Err(format!("unknown Wind Tunnel control: {other:?}")),
        }
        Ok(())
    }

    /// `true` while a background solve is running.
    pub fn is_running(&self) -> bool {
        self.run.is_some()
    }

    /// Record the current form state on the undo stack. The Run
    /// section calls this when the user spawns a solve so a later
    /// `Ctrl+Z` rewinds them back to the prior settings.
    pub fn record_form(&mut self) {
        self.history.record(self.form.clone());
    }

    /// Undo the last form-state snapshot.
    pub fn undo_edit(&mut self) -> bool {
        let current = self.form.clone();
        if let Some(prev) = self.history.undo(current) {
            self.form = prev;
            self.error = None;
            true
        } else {
            false
        }
    }
    /// Redo the most recently undone form-state snapshot.
    pub fn redo_edit(&mut self) -> bool {
        let current = self.form.clone();
        if let Some(next) = self.history.redo(current) {
            self.form = next;
            self.error = None;
            true
        } else {
            false
        }
    }
    /// `true` if Ctrl+Z would change the form state.
    pub fn can_undo(&self) -> bool {
        self.history.can_undo()
    }
    /// `true` if Ctrl+Y would change the form state.
    pub fn can_redo(&self) -> bool {
        self.history.can_redo()
    }
}

/// Draw the Aerodynamics / Wind Tunnel right-side panel.
///
/// Mirrors [`crate::genetics_workbench::draw_genetics_workbench`]: a
/// no-op when the `show_aero_workbench` toggle is off, otherwise a
/// resizable [`egui::SidePanel`] mounted before the central viewport so
/// egui docks it to the right (alongside the Mesh Toolbox / Genetics
/// workbench when several are open).
pub fn draw_aero_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_aero_workbench {
        return;
    }

    // Drain the background run's progress + completion before drawing,
    // so the panel always shows the freshest residuals / results.
    pump_aero_run(app, ctx);

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_aero_workbench",
        "Wind Tunnel",
        aero_workbench_body,
    );
    if close {
        app.show_aero_workbench = false;
    }
}

/// The Wind Tunnel (Aero) workbench body — body / wind / ground / tunnel /
/// solver / run / results / visualization sections. Extracted from
/// [`draw_aero_workbench`] so it can be hosted by the classic
/// [`crate::workbench_chrome::workbench_shell`] *or* the opt-in dockable
/// tile layout ([`crate::dock_layout`]). Drains the background-run poll up
/// front (cheap `Context` clone) so the dock path shows fresh residuals.
pub(crate) fn aero_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let ctx = ui.ctx().clone();
    pump_aero_run(app, &ctx);
    ui.label(
        egui::RichText::new("Virtual wind tunnel — 3-D external-aerodynamics CFD")
            .weak()
            .small(),
    );
    ui.label(
        egui::RichText::new("backed by `valenx-aero`")
            .weak()
            .small(),
    );
    ui.separator();

    // Fade-in animation on workbench open — when the user
    // toggles the workbench on via Ctrl+3 / View → Wind Tunnel
    // the panel body fades in over 0.18 s rather than popping
    // in instantly. The animation auto-resets when the panel
    // closes.
    let anim_id = egui::Id::new("valenx_aero_workbench_open");
    let t = ui.ctx().animate_bool_with_time(anim_id, true, 0.18);
    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            ui.scope(|ui| {
                ui.set_opacity(t.clamp(0.0, 1.0));
                panels::draw_body_section(app, ui);
                panels::draw_wind_section(app, ui);
                panels::draw_ground_section(app, ui);
                panels::draw_tunnel_section(app, ui);
                panels::draw_solver_section(app, ui);
                panels::draw_run_section(app, ui);
                panels::draw_results_section(app, ui);
                panels::draw_visualization_section(app, ui);
            });
        });
}

/// Poll the background wind-tunnel run: drain progress messages into
/// the live plots, and on completion move the result into the
/// workbench state.
fn pump_aero_run(app: &mut ValenxApp, ctx: &egui::Context) {
    use crate::aero::compute::{AeroOutcome, AeroProgress};

    let Some(handle) = app.aero.run.as_mut() else {
        return;
    };

    // Repaint promptly while a run is live so the residual plot and the
    // status line keep up with the channel.
    ctx.request_repaint_after(std::time::Duration::from_millis(80));

    // Drain progress.
    for msg in handle.poll() {
        match msg {
            AeroProgress::Stage(s) => app.aero.status = s,
            AeroProgress::Residual(it, r) => {
                app.aero.residual_history.push((it as f64, r));
            }
            AeroProgress::SweepPoint(angle, cd, cl) => {
                app.aero.sweep_progress.push((angle, cd, cl));
            }
        }
    }

    // Check for completion.
    if let Some(outcome) = handle.take_outcome() {
        app.aero.run = None;
        match outcome {
            AeroOutcome::Steady(result, report) => {
                app.aero.status = if result.converged {
                    format!(
                        "Converged — {} iterations, residual {:.2e}",
                        result.flow.iterations, result.flow.residual
                    )
                } else {
                    format!(
                        "Stopped at the iteration cap ({} iters, residual {:.2e}) \
                         — coefficients are provisional",
                        result.flow.iterations, result.flow.residual
                    )
                };
                // Mirror the residual history off the result in case
                // any messages were missed between frames.
                if app.aero.residual_history.is_empty() {
                    app.aero.residual_history = result
                        .flow
                        .residual_history
                        .iter()
                        .enumerate()
                        .map(|(i, &r)| ((i + 1) as f64, r))
                        .collect();
                }
                app.aero.last_result = Some(result);
                app.aero.last_report = Some(report);
                app.aero.last_polar = None;
            }
            AeroOutcome::Sweep(curve) => {
                app.aero.status = format!("Sweep complete — {} angles", curve.points.len());
                app.aero.last_polar = Some(PolarSweepResult { curve: *curve });
                app.aero.last_result = None;
                app.aero.last_report = None;
            }
            AeroOutcome::Failed(e) => {
                app.aero.status = "Run failed".to_string();
                app.aero.error = Some(e);
            }
        }
    }
}

/// The agent-bridge product for the aerodynamics workbench
/// (`show_3d{kind="aero"}`).
///
/// Aero computes a *field on a body*, not geometry — so this product uses the
/// workbench's built-in **demo box** (`valenx_aero::geometry::box_body`, the
/// always-available canonical bluff body that needs no user-loaded mesh), runs a
/// small bounded steady RANS solve over it (`valenx_aero::run_windtunnel` with a
/// coarse grid and capped iterations so the builder stays cheap and
/// deterministic), then paints the **surface pressure coefficient `Cp`** onto
/// the voxelized body shell. [`crate::aero::viz::build_flow_viz`] yields a
/// `(valenx_mesh::Mesh, valenx_fields::Field)` pair (one flat-shaded quad per
/// surface face, the `Cp` value on its nodes); the field is mapped to
/// triangle-major per-vertex `vertex_colors` through the shared cool-to-warm
/// ramp ([`crate::products_registry::node_field_to_vertex_colors`]) so the body
/// renders as a `Cp` map (blue low → red high). Pure and app-state-free. The
/// readout reports the demo case and the drag / lift coefficients.
pub(crate) fn aero_product() -> crate::WorkspaceProduct {
    /// The geometry + colours + headline coefficients a successful demo solve
    /// yields — a named struct so the fallible builder's return type stays
    /// simple (no complex tuple).
    struct Built {
        mesh: valenx_mesh::Mesh,
        colors: Vec<[f32; 3]>,
        cd: f64,
        cl: f64,
    }

    // A small canonical demo body (car-ish proportions, ~1 m scale) + a bounded
    // coarse solve so the synchronous build is fast and deterministic.
    let speed = 20.0_f64;
    let built = (|| -> Result<Built, String> {
        use crate::aero::model::{CutAxis, FlowField};
        use valenx_aero::{
            geometry::box_body, run_windtunnel, AeroRequest, TunnelSizing, TurbulenceModel,
        };
        let body = box_body(
            nalgebra::Vector3::new(-0.6, -0.3, -0.2),
            nalgebra::Vector3::new(0.6, 0.3, 0.2),
        );
        let req = AeroRequest::new(speed)
            .with_turbulence(TurbulenceModel::KEpsilon)
            .with_sizing(TunnelSizing {
                cells_across_body: 10,
                max_cells: 400_000,
                ..TunnelSizing::default()
            })
            .with_max_iterations(20);
        let result = run_windtunnel(&body, &req).map_err(|e| e.to_string())?;
        let cd = result.coefficients.cd;
        let cl = result.coefficients.cl;
        let viz = crate::aero::viz::build_flow_viz(&result, FlowField::SurfaceCp, CutAxis::Y, 0.5)?;
        // The field range drives the colour ramp endpoints.
        let (min, max) = viz.field.range.unwrap_or_else(|| {
            let mut lo = f64::INFINITY;
            let mut hi = f64::NEG_INFINITY;
            for &v in &viz.field.data {
                if v.is_finite() {
                    lo = lo.min(v);
                    hi = hi.max(v);
                }
            }
            if lo <= hi {
                (lo, hi)
            } else {
                (0.0, 1.0)
            }
        });
        let colors = crate::products_registry::node_field_to_vertex_colors(
            &viz.mesh,
            &viz.field.data,
            min,
            max,
        );
        Ok(Built {
            mesh: viz.mesh,
            colors,
            cd,
            cl,
        })
    })();

    match built {
        Ok(b) => {
            let loaded = crate::products_registry::loaded_mesh_from(b.mesh, "<aero>/surface-cp");
            let camera = crate::products_registry::camera_for(&loaded.mesh);
            let lines = vec![
                format!("wind tunnel: demo box @ {speed:.0} m/s (k-\u{03B5})"),
                format!("Cd {:+.4} · Cl {:+.4}", b.cd, b.cl),
                "surface coloured by pressure coefficient Cp".to_string(),
            ];
            crate::WorkspaceProduct {
                title: "Aero (Cp on demo body)".into(),
                lines,
                mesh: Some(loaded),
                vertex_colors: Some(b.colors),
                camera,
                kind2d: None,
                last_export: None,
                image: None,
                image_texture: None,
                animation: None,
            }
        }
        Err(e) => {
            // Theoretically unreachable for the bounded demo solve; degrade to a
            // tiny placeholder triangle + a note rather than panicking.
            let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
            block.connectivity = vec![0, 1, 2];
            let mut placeholder = valenx_mesh::Mesh::new("valenx-aero-surface");
            placeholder.nodes = vec![
                nalgebra::Vector3::new(0.0, 0.0, 0.0),
                nalgebra::Vector3::new(1.0, 0.0, 0.0),
                nalgebra::Vector3::new(0.0, 1.0, 0.0),
            ];
            placeholder.element_blocks.push(block);
            placeholder.recompute_stats();
            let loaded =
                crate::products_registry::loaded_mesh_from(placeholder, "<aero>/surface-cp");
            let camera = crate::products_registry::camera_for(&loaded.mesh);
            crate::WorkspaceProduct {
                title: "Aero (Cp on demo body)".into(),
                lines: vec![
                    "aero surface-pressure field".to_string(),
                    format!("solve unavailable — showing placeholder ({e})"),
                ],
                mesh: Some(loaded),
                vertex_colors: None,
                camera,
                kind2d: None,
                last_export: None,
                image: None,
                image_texture: None,
                animation: None,
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fresh_state_is_idle_with_a_default_form() {
        let s = AeroWorkbenchState::default();
        assert!(!s.is_running());
        assert!(s.last_result.is_none());
        assert!(s.last_polar.is_none());
        assert!(s.error.is_none());
        assert!(s.residual_history.is_empty());
        // The default form is the road-car case.
        assert!((s.form.speed_ms - 30.0).abs() < 1e-12);
    }

    #[test]
    fn workbench_is_off_by_default_on_a_fresh_app() {
        // The Wind Tunnel workbench, like the Genetics workbench, is
        // hidden until the user turns it on from the View menu.
        let app = ValenxApp::default();
        assert!(!app.show_aero_workbench);
    }

    #[test]
    fn agent_set_sets_param_unknown_and_type_mismatch_err() {
        use crate::agent_commands::AgentValue;
        let mut s = AeroWorkbenchState::default();
        // A representative float set lands in the form.
        s.agent_set("Free-stream speed", &AgentValue::Float(45.0))
            .unwrap();
        assert!((s.form.speed_ms - 45.0).abs() < 1e-12);
        // An integer count field.
        s.agent_set("Max iterations", &AgentValue::Int(150))
            .unwrap();
        assert_eq!(s.form.max_iterations, 150);
        // An enum-by-name selector.
        s.agent_set("Turbulence model", &AgentValue::Str("k-epsilon".into()))
            .unwrap();
        assert_eq!(s.form.turbulence, TurbulenceChoice::KEpsilon);
        // A boolean toggle (radio-backed bool).
        s.agent_set("moving ground", &AgentValue::Bool(false))
            .unwrap();
        assert!(!s.form.moving_ground);
        // The double-space caption resolves exactly.
        s.agent_set("Aspect ratio  AR", &AgentValue::Float(8.0))
            .unwrap();
        assert!((s.form.aspect_ratio - 8.0).abs() < 1e-12);
        // Unknown caption -> Err.
        assert!(s.agent_set("no such control", &AgentValue::Int(1)).is_err());
        // Unknown enum name -> Err.
        assert!(s
            .agent_set("Turbulence model", &AgentValue::Str("les".into()))
            .is_err());
        // Type mismatch (string into the float speed) -> Err, field untouched.
        assert!(s
            .agent_set("Free-stream speed", &AgentValue::Str("fast".into()))
            .is_err());
        assert!(
            (s.form.speed_ms - 45.0).abs() < 1e-12,
            "rejected set leaves field untouched"
        );
    }
}

/// Headless egui UI-logic tests for the Wind Tunnel workbench host
/// panel.
#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use crate::aero::model::{BodySource, RunMode};

    /// Run the whole workbench panel once in a headless context.
    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_aero_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        // With the toggle off the workbench draws nothing and never
        // panics — the default state.
        let mut app = ValenxApp::default();
        assert!(!app.show_aero_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        // With the workbench shown, the whole eight-section side panel
        // mounts headlessly — fresh state, the demo box selected.
        let mut app = ValenxApp::default();
        app.show_aero_workbench = true;
        app.aero.form.body_source = BodySource::DemoBox;
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_in_sweep_mode_without_panic() {
        let mut app = ValenxApp::default();
        app.show_aero_workbench = true;
        app.aero.form.body_source = BodySource::DemoBox;
        app.aero.form.run_mode = RunMode::AngleSweep;
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_an_error_state_without_panic() {
        let mut app = ValenxApp::default();
        app.show_aero_workbench = true;
        app.aero.error = Some("invalid case: bad air".to_string());
        app.aero.status = "Run failed".to_string();
        draw_workbench(&mut app);
    }
}
