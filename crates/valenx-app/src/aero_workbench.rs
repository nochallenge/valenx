//! The right-side **Aerodynamics / Wind Tunnel** workbench panel.
//!
//! Round 6 Block 15 shipped [`valenx_aero`] — a native 3-D external-
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
use crate::aero::model::WindTunnelForm;
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

impl AeroWorkbenchState {
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
    pub fn can_undo(&self) -> bool { self.history.can_undo() }
    /// `true` if Ctrl+Y would change the form state.
    pub fn can_redo(&self) -> bool { self.history.can_redo() }
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

    egui::SidePanel::right("valenx_aero_workbench")
        .resizable(true)
        .default_width(400.0)
        .width_range(330.0..=680.0)
        .show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.heading("Wind Tunnel");
            });
            ui.label(
                egui::RichText::new(
                    "Virtual wind tunnel — 3-D external-aerodynamics CFD",
                )
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
            let t = ctx.animate_bool_with_time(anim_id, true, 0.18);
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
                app.aero.status =
                    format!("Sweep complete — {} angles", curve.points.len());
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
