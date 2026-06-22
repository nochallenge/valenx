//! Part Design → **Parametric Sketch (constraints)** — a first-class,
//! discoverable host for the in-house `valenx-sketch` constraint
//! sketcher.
//!
//! ## Why this module exists
//!
//! valenx already ships a complete FreeCAD-Phase-1 parametric sketcher
//! (the `valenx-sketch` crate: a Newton-Raphson / Levenberg-Marquardt
//! constraint solver over 2-D `Point`/`Line`/`Circle`/`Arc`/… entities)
//! and a fully-wired UI for it in
//! [`crate::mesh_toolbox::draw_sketcher_panel`]. That panel does the
//! whole Fusion-grade flow — draw entities → apply dimensional +
//! geometric constraints → **Solve** (with under/over-constrained DOF
//! feedback) → **Extrude** into a `valenx_cad::Solid` loaded back into
//! the viewport.
//!
//! The problem was *discoverability*: that panel lived buried inside the
//! Mesh Toolbox as a `collapsing("Sketcher", …)` sub-section, several
//! clicks deep, so users couldn't find it. This module re-surfaces the
//! *exact same* pipeline as a standalone right-side workbench reachable
//! from a dedicated **Part Design → "Parametric Sketch (constraints)"**
//! top-bar menu item, sitting next to the lightweight "Sketch (draw on
//! canvas)" polygon canvas. The constraint sketcher is the more powerful
//! mode; the polygon canvas stays as-is.
//!
//! ## What is reused vs. new
//!
//! - **Reused, unchanged:** the entire panel body
//!   [`crate::mesh_toolbox::draw_sketcher_panel`] (entity creation,
//!   constraint palette, solve + DOF diagnostics, Pad/extrude), the
//!   constraint solver `valenx_sketch::solver::solve`, the live viewport
//!   overlay [`crate::sketch_overlay`], and the shared sketch state
//!   `app.mesh_toolbox.sketcher` ([`crate::mesh_toolbox::SketcherPanelState`]).
//!   Because both the menu's polygon-canvas mode and this constraint
//!   mode are *additive*, opening this panel does not disturb the
//!   feature-tree CAD workbench.
//! - **New:** only the visibility flag
//!   [`crate::ValenxApp::show_param_sketch`], this thin
//!   [`workbench_shell`]-hosted wrapper, a one-line discoverable intro,
//!   and the Part Design menu entry.
//!
//! [`workbench_shell`]: crate::workbench_chrome::workbench_shell

use crate::ValenxApp;
use eframe::egui;

/// egui id + persisted-chrome key for the Parametric Sketch workbench.
/// Used by [`crate::workbench_chrome::workbench_shell`] to key its
/// docked / floating / torn-off window state, and by the dock/tab
/// machinery if this panel is later added there.
pub const PARAM_SKETCH_ID: &str = "valenx_param_sketch";

/// Human title shown in the panel header and tab strip.
pub const PARAM_SKETCH_TITLE: &str = "Parametric Sketch";

/// Draw the Parametric Sketch (constraints) workbench as a first-class
/// right-side panel. A no-op when [`ValenxApp::show_param_sketch`] is
/// off.
///
/// The body delegates verbatim to
/// [`crate::mesh_toolbox::draw_sketcher_panel`] — the same constraint
/// sketcher UI the Mesh Toolbox hosts — preceded by a short
/// discoverable intro that names the entity → constraint → solve →
/// extrude flow. Hosting it through
/// [`crate::workbench_chrome::workbench_shell`] gives it the standard
/// collapse / float / tear-off / ✕-close chrome every other workbench
/// has, and keeps its sketch state shared with the Mesh Toolbox copy
/// (both read & write `app.mesh_toolbox.sketcher`), so a sketch started
/// in one is visible in the other.
pub fn draw_param_sketch_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_param_sketch {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        PARAM_SKETCH_ID,
        PARAM_SKETCH_TITLE,
        |app, ui| {
            // Discoverable intro: spell out the flow so a first-time
            // user knows this is the constraint-driven sketch → solid
            // path (Fusion-grade: under/over-constrained feedback).
            ui.label(
                egui::RichText::new(
                    "In-house parametric constraint sketcher · valenx-sketch. \
                     Flow: draw entities → apply geometric + dimensional \
                     constraints → Solve (with under/over-constrained DOF \
                     feedback) → Pad/extrude into a 3-D solid.",
                )
                .weak()
                .small(),
            )
            .on_hover_text(
                "Constraint-driven 2-D sketcher: a Newton-Raphson / \
                 Levenberg-Marquardt solver drives every variable to satisfy \
                 every constraint, then the closed profile extrudes into a \
                 valenx_cad solid. The live sketch renders on the XY plane in \
                 the viewport (toggle below).",
            );
            ui.separator();
            // Reuse the full, already-wired sketcher pipeline verbatim.
            crate::mesh_toolbox::draw_sketcher_panel(app, ui);
        },
    );

    if close {
        app.show_param_sketch = false;
    }
}

#[cfg(test)]
// `ValenxApp` is a large struct; flipping a single visibility flag after
// `::default()` is the codebase-wide idiom for these workbench tests
// (see e.g. `bolt_workbench`), clearer here than spelling out a struct
// literal with `..Default::default()`.
#[allow(clippy::field_reassign_with_default)]
mod tests {
    use super::*;

    /// Run a single headless egui frame and call `f` inside it. Mirrors
    /// the pattern in [`crate::coverage_ui_tests`] — a panic anywhere in
    /// the workbench draw path fails the test.
    fn run_frame(app: &mut ValenxApp, f: impl FnOnce(&mut ValenxApp, &egui::Context)) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| f(app, ctx));
    }

    /// The workbench is a clean no-op when its visibility flag is off —
    /// the full entry point must draw nothing, mutate nothing, and not
    /// panic.
    #[test]
    fn hidden_workbench_is_a_no_op() {
        let mut app = ValenxApp::default();
        assert!(!app.show_param_sketch, "defaults to hidden");
        let entities_before = app.mesh_toolbox.sketcher.sketch.entities.len();
        run_frame(&mut app, draw_param_sketch_workbench);
        assert_eq!(
            app.mesh_toolbox.sketcher.sketch.entities.len(),
            entities_before,
            "hidden panel must not mutate the sketch"
        );
    }

    /// When shown, the full workbench (shell + body) renders without
    /// panicking and operates on the shared `app.mesh_toolbox.sketcher`
    /// state — a sketch seeded there survives the render, proving the
    /// panel reuses the same sketch the Mesh Toolbox sees.
    #[test]
    fn shown_workbench_renders_shared_sketcher() {
        let mut app = ValenxApp::default();
        app.show_param_sketch = true;
        // Seed a tiny sketch into the SHARED state.
        let p = app.mesh_toolbox.sketcher.sketch.add_point(0.0, 0.0);
        let q = app.mesh_toolbox.sketcher.sketch.add_point(1.0, 0.0);
        app.mesh_toolbox
            .sketcher
            .sketch
            .add_line(p, q)
            .expect("seed line");
        run_frame(&mut app, draw_param_sketch_workbench);
        assert!(
            app.mesh_toolbox.sketcher.sketch.entities.len() >= 3,
            "shared sketch (2 points + 1 line) survives a render"
        );
    }

    /// End-to-end of the headline flow this panel surfaces, exercised
    /// through the real `valenx-sketch` API (the same calls the panel's
    /// buttons make): build a closed unit square, constrain two edges
    /// Horizontal / Vertical and pin a side length with Distance, Solve
    /// → expect Converged, then Pad/extrude → expect a non-empty solid.
    ///
    /// This is the "Fusion value" path: constraints → solve →
    /// well-/under-/over-constrained DOF feedback → extrude.
    #[test]
    fn constrained_square_solves_converged_then_extrudes_nonempty() {
        use valenx_sketch::constraint::Constraint;

        let mut sketch = valenx_sketch::Sketch::new();
        // Four corners of a (roughly) unit square.
        let a = sketch.add_point(0.0, 0.0);
        let b = sketch.add_point(1.0, 0.0);
        let c = sketch.add_point(1.0, 1.0);
        let d = sketch.add_point(0.0, 1.0);
        // Closed loop of four lines a→b→c→d→a.
        let bottom = sketch.add_line(a, b).expect("bottom edge");
        let right = sketch.add_line(b, c).expect("right edge");
        let _top = sketch.add_line(c, d).expect("top edge");
        let left = sketch.add_line(d, a).expect("left edge");

        // Geometric constraints: bottom horizontal, left/right vertical.
        sketch.add_constraint(Constraint::Horizontal(bottom));
        sketch.add_constraint(Constraint::Vertical(left));
        sketch.add_constraint(Constraint::Vertical(right));
        // Dimensional constraint: pin the bottom edge to length 2.0
        // (Distance between its two endpoints). The solver should move
        // b from x=1 to x=2 to satisfy it.
        sketch.add_constraint(Constraint::Distance {
            a,
            b,
            target: 2.0,
        });

        // Solve — expect convergence.
        let report = valenx_sketch::solver::solve(&mut sketch, valenx_sketch::SolverConfig::default())
            .expect("solver runs on the square");
        assert_eq!(
            report.status,
            valenx_sketch::SolverStatus::Converged,
            "constrained square should converge (report: {report:?})"
        );
        // The Distance constraint should have been satisfied: |a→b| ≈ 2.
        let (ax, ay) = sketch
            .point_at(a)
            .expect("point a")
            .read(&sketch.vars);
        let (bx, by) = sketch
            .point_at(b)
            .expect("point b")
            .read(&sketch.vars);
        let len = ((bx - ax).powi(2) + (by - ay).powi(2)).sqrt();
        assert!(
            (len - 2.0).abs() < 1e-3,
            "Distance(2.0) constraint satisfied after solve, got {len}"
        );

        // Pad/extrude the solved profile → a non-empty solid, then
        // tessellate it the way the panel's "Pad" button does.
        let solid = sketch.extrude(3.0).expect("extrude the closed square");
        let mesh = valenx_cad::solid_to_mesh(&solid, valenx_cad::DEFAULT_TESS_TOLERANCE)
            .expect("tessellate the padded solid");
        assert!(
            !mesh.nodes.is_empty(),
            "the extruded square should tessellate into a non-empty mesh"
        );
    }
}
