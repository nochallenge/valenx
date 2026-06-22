//! The right-side **Astro / Launch** workbench panel.
//!
//! The `valenx-astro` crate is a native launch-vehicle ascent +
//! trajectory simulator — point a rocket at the sky, fly it to orbit,
//! and read back the engineering answer (the orbit reached, the `Δv`
//! budget, peak dynamic pressure, the staging timeline) — plus a family
//! of closed-form mission planners (Hohmann transfers, hoverslam
//! ignition altitude, two-impulse rendezvous, launch azimuth). It
//! shipped as a library + agent API with **no UI**.
//!
//! This module is that UI — a polished egui side panel that mirrors the
//! CFD-side [`crate::aero_workbench`]: a resizable right-hand
//! [`egui::SidePanel`], toggled from the View menu (Ctrl+4), off by
//! default, with a fade-in on open. The panel has two tabs — an
//! **Ascent to orbit** simulator (vehicle setup → Run → result summary +
//! flight-profile chart) and a set of **mission planners** — and the
//! real work is split across the [`crate::astro`] sub-modules.
//!
//! Unlike the wind tunnel, whose steady RANS solve runs on a background
//! thread, the launch ascent is a bounded fixed-step RK4 integration
//! that completes in well under a frame, so this workbench runs it
//! **synchronously on click** (see [`crate::astro::run`]).

use eframe::egui;

use crate::astro::model::{AscentForm, AstroTab, PlannerForm};
use crate::astro::panels;
use crate::ValenxApp;
use valenx_astro::AscentResult;

/// All Astro / Launch workbench form + result state.
///
/// One instance lives on [`crate::ValenxApp`] (the `astro` field),
/// exactly as the CFD-side `AeroWorkbenchState` does. Off by default —
/// no simulation runs until the user clicks Run.
#[derive(Default)]
pub struct AstroWorkbenchState {
    /// Which sub-view (Ascent / Planners) is selected.
    pub tab: AstroTab,

    /// The ascent-simulation form inputs (stages, payload, guidance, …).
    pub ascent: AscentForm,
    /// The closed-form mission-planner inputs.
    pub planner: PlannerForm,

    /// The last completed ascent run, if any. `None` until the user
    /// clicks Run (so the panel is off / idle by default). Boxed because
    /// [`AscentResult`] carries the full trajectory sample series.
    pub last_result: Option<Box<AscentResult>>,
    /// A coarse status line for the Run section (the last completion or
    /// failure message).
    pub status: String,
    /// The last error message, shown in red. Cleared on a new run.
    pub error: Option<String>,

    /// Undo / redo over the ascent form. A snapshot lands on the stack
    /// when the user presses Run, so `Ctrl+Z` reverses the settings of
    /// the last run (mirrors the aero workbench).
    pub history: crate::undo::History<AscentForm>,
}

impl AstroWorkbenchState {
    /// Record the current ascent form on the undo stack. The Run action
    /// calls this when the user runs a sim so a later `Ctrl+Z` rewinds
    /// them back to the prior settings.
    pub fn record_form(&mut self) {
        self.history.record(self.ascent.clone());
    }

    /// Undo the last form-state snapshot.
    pub fn undo_edit(&mut self) -> bool {
        let current = self.ascent.clone();
        if let Some(prev) = self.history.undo(current) {
            self.ascent = prev;
            self.error = None;
            true
        } else {
            false
        }
    }
    /// Redo the most recently undone form-state snapshot.
    pub fn redo_edit(&mut self) -> bool {
        let current = self.ascent.clone();
        if let Some(next) = self.history.redo(current) {
            self.ascent = next;
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

/// Draw the Astro / Launch right-side panel.
///
/// Mirrors [`crate::aero_workbench::draw_aero_workbench`]: a no-op when
/// the `show_astro_workbench` toggle is off, otherwise a resizable
/// [`egui::SidePanel`] mounted before the central viewport so egui docks
/// it to the right (alongside the Mesh Toolbox / Genetics / Wind Tunnel
/// workbenches when several are open).
pub fn draw_astro_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_astro_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_astro_workbench",
        "Astro / Launch",
        astro_workbench_body,
    );
    if close {
        app.show_astro_workbench = false;
    }
}

/// The Astro / Launch workbench body — the ascent-simulator + mission-
/// planner tabs. Extracted from [`draw_astro_workbench`] so it can be
/// hosted by the classic [`crate::workbench_chrome::workbench_shell`] *or*
/// the opt-in dockable tile layout ([`crate::dock_layout`]) without
/// duplicating logic.
pub(crate) fn astro_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new("Launch-vehicle ascent + trajectory simulator")
            .weak()
            .small(),
    );
    ui.label(
        egui::RichText::new("backed by `valenx-astro`")
            .weak()
            .small(),
    );
    ui.separator();

    // Fade-in animation on workbench open — when the user toggles
    // the workbench on via Ctrl+4 / View → Astro / Launch the
    // panel body fades in over 0.18 s rather than popping in
    // instantly. The animation auto-resets when the panel closes
    // (matches the wind tunnel).
    let anim_id = egui::Id::new("valenx_astro_workbench_open");
    let t = ui.ctx().animate_bool_with_time(anim_id, true, 0.18);
    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            ui.scope(|ui| {
                ui.set_opacity(t.clamp(0.0, 1.0));
                panels::draw_tab_selector(app, ui);
                match app.astro.tab {
                    AstroTab::Ascent => {
                        panels::draw_vehicle_section(app, ui);
                        panels::draw_guidance_section(app, ui);
                        panels::draw_run_section(app, ui);
                        panels::draw_results_section(app, ui);
                    }
                    AstroTab::Planners => {
                        panels::draw_planners_section(app, ui);
                    }
                }
            });
        });
}

/// Build the **Astro / Launch** result card for the Workbench+Agent bridge — a
/// DATA-ONLY [`crate::WorkspaceProduct`] (`mesh: None`) whose `lines` are genuine
/// computed results for the canonical default forms: (1) the ascent-to-orbit
/// simulation ([`valenx_astro::simulate_ascent`]) of the default launch vehicle —
/// outcome, apoapsis, periapsis, ideal Δv budget — and (2) the default Hohmann
/// transfer ([`valenx_astro::hohmann_transfer`], a LEO→GEO 300 → 35 786 km
/// transfer) — the two burns and the total Δv. Registered as the `"astro"`
/// producer in [`crate::products_registry::lookup`]; the tile renders it as a
/// text card, not a 3-D view.
///
/// Both computations are state-free (built from `AscentForm::default()` /
/// `PlannerForm::default()`) and cheap — the ascent is a bounded fixed-step RK4
/// integration that completes well under a frame, and the Hohmann transfer is
/// closed-form. On the (canonically-unreachable) solver error each section
/// carries that message instead of panicking.
pub(crate) fn astro_product() -> crate::WorkspaceProduct {
    let mut lines = Vec::new();

    // (1) Ascent to orbit — the genuine simulate_ascent result for the default
    // launch vehicle + config (same path the Run button drives).
    let ascent = AscentForm::default();
    let vehicle = ascent.build_vehicle();
    let config = ascent.build_config();
    match valenx_astro::simulate_ascent(&vehicle, &config) {
        Ok(r) => {
            lines.push(format!("ascent outcome   : {:?}", r.outcome));
            lines.push(format!("apoapsis         : {:.0} km", r.apoapsis_km()));
            lines.push(format!("periapsis        : {:.0} km", r.periapsis_km()));
            lines.push(format!("ideal Δv budget  : {:.0} m/s", r.ideal_delta_v));
        }
        Err(e) => lines.push(format!(
            "ascent failed: {}",
            crate::astro::model::friendly_error(&e)
        )),
    }

    // (2) Hohmann transfer — the genuine closed-form Δv between the default
    // departure / arrival circular altitudes (LEO 300 km → GEO 35 786 km).
    let planner = PlannerForm::default();
    let r1 = crate::astro::model::altitude_km_to_radius_m(planner.hohmann_from_km);
    let r2 = crate::astro::model::altitude_km_to_radius_m(planner.hohmann_to_km);
    lines.push(String::new());
    lines.push(format!(
        "Hohmann {:.0} → {:.0} km:",
        planner.hohmann_from_km, planner.hohmann_to_km
    ));
    match valenx_astro::hohmann_transfer(r1, r2) {
        Ok(t) => {
            lines.push(format!("  burn 1 Δv      : {:.0} m/s", t.delta_v1));
            lines.push(format!("  burn 2 Δv      : {:.0} m/s", t.delta_v2));
            lines.push(format!("  total Δv       : {:.0} m/s", t.total_delta_v));
        }
        Err(e) => lines.push(format!(
            "  cannot plan: {}",
            crate::astro::model::friendly_error(&e)
        )),
    }

    crate::WorkspaceProduct {
        title: "Astro / Launch".into(),
        lines,
        mesh: None,
        vertex_colors: None,
        camera: Default::default(),
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
        animation: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fresh_state_is_idle_with_a_default_form() {
        let s = AstroWorkbenchState::default();
        assert!(s.last_result.is_none());
        assert!(s.error.is_none());
        assert!(!s.can_undo());
        assert!(!s.can_redo());
        // The default tab is the ascent simulator.
        assert_eq!(s.tab, AstroTab::Ascent);
        // The default form is the medium-lift preset (two stages).
        assert_eq!(s.ascent.stages.len(), 2);
    }

    #[test]
    fn workbench_is_off_by_default_on_a_fresh_app() {
        // The Astro workbench, like the Wind Tunnel workbench, is hidden
        // until the user turns it on from the View menu.
        let app = ValenxApp::default();
        assert!(!app.show_astro_workbench);
    }

    #[test]
    fn undo_redo_round_trips_the_form() {
        let mut s = AstroWorkbenchState::default();
        s.record_form();
        s.ascent.payload_mass = 999.0;
        assert!(s.undo_edit());
        // Undo restored the recorded payload (the preset 10 000 kg).
        assert!((s.ascent.payload_mass - 10_000.0).abs() < 1e-9);
        assert!(s.redo_edit());
        assert!((s.ascent.payload_mass - 999.0).abs() < 1e-9);
    }
}

/// Headless egui UI-logic tests for the Astro / Launch workbench host
/// panel.
///
/// The whole panel is rendered into a windowless [`egui::Context`] in
/// each tab + state; nothing opens an OS window and nothing reaches
/// `rfd::FileDialog` (the panel has no file IO at all).
#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use crate::astro::model::GuidanceChoice;

    /// Run the whole workbench panel once in a headless context.
    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_astro_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        // With the toggle off the workbench draws nothing and never
        // panics — the default state.
        let mut app = ValenxApp::default();
        assert!(!app.show_astro_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_the_ascent_tab_without_panic() {
        // With the workbench shown on the Ascent tab, the whole vehicle
        // setup + guidance + run + results column mounts headlessly with
        // the fresh (preset) state.
        let mut app = ValenxApp::default();
        app.show_astro_workbench = true;
        app.astro.tab = AstroTab::Ascent;
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_the_planners_tab_without_panic() {
        // The Planners tab renders all four closed-form planner cards,
        // each computing its output live from the default inputs.
        let mut app = ValenxApp::default();
        app.show_astro_workbench = true;
        app.astro.tab = AstroTab::Planners;
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_post_run_results_without_panic() {
        // A completed ascent → the Results section's summary cards +
        // flight-profile chart must render. The run is synchronous, so
        // the result is present the same call.
        let mut app = ValenxApp::default();
        app.show_astro_workbench = true;
        crate::astro::run::run_ascent(&mut app);
        assert!(
            app.astro.last_result.is_some() || app.astro.error.is_some(),
            "a run should produce a result or an error"
        );
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_closed_loop_results_without_panic() {
        // Closed-loop insertion exercises the target-altitude path + the
        // flight-events list with the circularisation event.
        let mut app = ValenxApp::default();
        app.show_astro_workbench = true;
        app.astro.ascent.guidance = GuidanceChoice::ClosedLoopInsertion;
        app.astro.ascent.pitch_kick_deg = 12.9;
        crate::astro::run::run_ascent(&mut app);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_an_error_state_without_panic() {
        // An error line + a failure status must render without panicking.
        let mut app = ValenxApp::default();
        app.show_astro_workbench = true;
        app.astro.error = Some("invalid case: bad stage".to_string());
        app.astro.status = "Ascent run failed".to_string();
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_singular_rendezvous_without_panic() {
        // A half-period transfer fraction makes the rendezvous BVP
        // singular — the planner must show the error text, not panic.
        let mut app = ValenxApp::default();
        app.show_astro_workbench = true;
        app.astro.tab = AstroTab::Planners;
        app.astro.planner.rdv_transfer_fraction = 0.5; // nT = π -> singular
        draw_workbench(&mut app);
    }
}
