//! Headless render-smoke tests for the UI-coverage workbenches.
//!
//! Each workbench added in the coverage sprint (and the neural-interface
//! readout) is drawn once in a **headless egui frame**. If any panel's draw
//! path panics — bad layout, a missing-id clash, an out-of-bounds access —
//! the test fails. This is the automated stand-in for "open every panel and
//! confirm it renders", which the dev binary can't be driven for via
//! computer-use (it isn't an installed app).

#![allow(clippy::field_reassign_with_default)]

use eframe::egui;

use crate::ValenxApp;

/// Draw every coverage workbench once in a headless context.
fn run_frame(app: &mut ValenxApp) {
    let ctx = egui::Context::default();
    let _ = ctx.run(egui::RawInput::default(), |ctx| {
        crate::cad_workbench::draw_cad_workbench(app, ctx);
        crate::draft2d_workbench::draw_draft2d_workbench(app, ctx);
        crate::reinforcement_workbench::draw_reinforcement_workbench(app, ctx);
        crate::render_workbench::draw_render_workbench(app, ctx);
        crate::hvac_workbench::draw_hvac_workbench(app, ctx);
        crate::reverse_workbench::draw_reverse_workbench(app, ctx);
        crate::interior_workbench::draw_interior_workbench(app, ctx);
        crate::animate_workbench::draw_animate_workbench(app, ctx);
        crate::variant_effect_workbench::draw_variant_effect_workbench(app, ctx);
        crate::neuro_workbench::draw_neuro_workbench(app, ctx);
        crate::fem_workbench::draw_fem_workbench(app, ctx);
        crate::cfd_workbench::draw_cfd_workbench(app, ctx);
        crate::param_sketch_panel::draw_param_sketch_workbench(app, ctx);
    });
}

#[test]
fn all_coverage_workbenches_draw_when_shown_without_panic() {
    let mut app = ValenxApp::default();
    app.show_cad_workbench = true;
    app.show_draft2d_workbench = true;
    app.show_reinforcement_workbench = true;
    app.show_render_workbench = true;
    app.show_hvac_workbench = true;
    app.show_reverse_workbench = true;
    app.show_interior_workbench = true;
    app.show_animate_workbench = true;
    app.show_variant_effect_workbench = true;
    app.show_neuro_workbench = true;
    app.show_fem_workbench = true;
    app.show_cfd_workbench = true;
    app.show_param_sketch = true;
    // Two frames: render, then re-render against retained state.
    run_frame(&mut app);
    run_frame(&mut app);
}

#[test]
fn coverage_workbenches_are_noops_when_hidden() {
    // Default state — all toggles off — draws nothing and never panics.
    let mut app = ValenxApp::default();
    run_frame(&mut app);
}
