//! The right-side **Pulley System Workbench** panel — native
//! rope-and-pulley mechanical-advantage analysis over `valenx-pulley`.
//!
//! Mirrors the Belt Drive / Heat Transfer workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_pulley_workbench`,
//! toggled from the View menu. The form picks the arrangement (a single
//! fixed pulley, a single movable pulley, or a block-and-tackle of `n`
//! supporting rope segments), the load weight, a lumped efficiency `eta`
//! and the lift distance `s`; "Analyze" reports the ideal mechanical
//! advantage `MA = n`, the equal velocity ratio, the ideal and
//! friction-aware real effort, the actual mechanical advantage `MA * eta`
//! and the load that real effort raises, plus the energy story for one
//! lift through `s` — the rope length pulled `VR * s`, the useful output
//! work `W * s`, the operator's input work `W * s / eta` and the work lost
//! to friction — and "Show 3-D" loads a representative pulley assembly (a fixed
//! block over a movable block, drawn as sheave cylinders) into the central
//! viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_pulley::{
    actual_mechanical_advantage, effort_distance, ideal_effort, input_work, load_from_effort,
    output_work, real_effort, work_lost, PulleySystem,
};

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// The pulley arrangement the form is configured for. The numeric
/// mechanical advantage always comes from the resulting
/// [`PulleySystem`]'s supporting-rope count, never from this tag alone.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum Arrangement {
    /// A single fixed pulley (`MA = 1`, direction change only).
    Fixed,
    /// A single movable pulley (`MA = 2`).
    Movable,
    /// A block-and-tackle with `supporting_ropes` segments (`MA = n`).
    BlockAndTackle,
}

/// Persistent form + result state for the Pulley System Workbench.
pub struct PulleyWorkbenchState {
    /// The selected pulley arrangement.
    arrangement: Arrangement,
    /// Number of rope segments supporting the load — used only when the
    /// arrangement is a block-and-tackle.
    supporting_ropes: u32,
    /// Load (weight) to raise (N).
    load_n: f64,
    /// Lumped mechanical efficiency `eta` in `(0, 1]`.
    efficiency: f64,
    /// Distance to raise the load (m) — drives the rope-pulled length and
    /// the work / energy readout.
    lift_distance_m: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D pulley assembly (serviced after
    /// the panel draws).
    show_3d_request: bool,
}

impl Default for PulleyWorkbenchState {
    fn default() -> Self {
        // A 4-rope block-and-tackle lifting a 1200 N load at 90%
        // efficiency: MA = VR = 4, ideal effort 300 N, real effort
        // 1200 / 3.6 = 333.33 N, actual MA 3.60. Raising it 2 m pulls
        // 8 m of rope, doing 2400 J on the load for 2666.67 J of input
        // (266.67 J lost to friction).
        Self {
            arrangement: Arrangement::BlockAndTackle,
            supporting_ropes: 4,
            load_n: 1200.0,
            efficiency: 0.9,
            lift_distance_m: 2.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Pulley System Workbench right-side panel. A no-op when the
/// `show_pulley_workbench` toggle is off.
pub fn draw_pulley_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_pulley_workbench {
        return;
    }

    egui::SidePanel::right("valenx_pulley_workbench")
        .resizable(true)
        .default_width(360.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            if crate::workbench_ui::header(
                ui,
                "Pulley System",
                "native rope-and-pulley mechanical advantage · valenx-pulley",
            ) {
                app.show_pulley_workbench = false;
            }

            let s = &mut app.pulley;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Arrangement").strong());
                    ui.radio_value(&mut s.arrangement, Arrangement::Fixed, "fixed pulley (MA = 1)");
                    ui.radio_value(
                        &mut s.arrangement,
                        Arrangement::Movable,
                        "movable pulley (MA = 2)",
                    );
                    ui.radio_value(
                        &mut s.arrangement,
                        Arrangement::BlockAndTackle,
                        "block and tackle (MA = n)",
                    );
                    if s.arrangement == Arrangement::BlockAndTackle {
                        ui.horizontal(|ui| {
                            ui.label("supporting ropes n");
                            ui.add(
                                egui::DragValue::new(&mut s.supporting_ropes)
                                    .speed(1.0)
                                    .range(1..=24),
                            );
                        });
                    }

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Load & efficiency").strong());
                    ui.horizontal(|ui| {
                        ui.label("load W (N)");
                        ui.add(egui::DragValue::new(&mut s.load_n).speed(10.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("efficiency η (0–1]");
                        ui.add(egui::DragValue::new(&mut s.efficiency).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("lift distance s (m)");
                        ui.add(
                            egui::DragValue::new(&mut s.lift_distance_m)
                                .speed(0.1)
                                .range(0.0..=f64::INFINITY),
                        );
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_pulley(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative pulley assembly (a fixed block over a movable block, drawn as sheave cylinders) as a 3-D solid and load it into the central viewport to orbit",
                        )
                        .clicked()
                    {
                        s.show_3d_request = true;
                    }

                    if let Some(e) = &s.error {
                        ui.add_space(4.0);
                        ui.colored_label(egui::Color32::from_rgb(220, 90, 90), e);
                    }

                    if !s.result.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new("Mechanical advantage").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        });

    // Serviced after the panel draws (the `&mut app.pulley` borrow is
    // released here): build the assembly's 3-D solid and load it.
    if app.pulley.show_3d_request {
        app.pulley.show_3d_request = false;
        load_pulley_3d(app);
    }
}

/// Validate the form, analyse the system and format the readout.
fn run_pulley(s: &mut PulleyWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the [`PulleySystem`] the form describes, mapping any domain error
/// to a display string. Extracted so it is shared by the readout and the
/// 3-D gate.
fn system_of(s: &PulleyWorkbenchState) -> Result<PulleySystem, String> {
    match s.arrangement {
        Arrangement::Fixed => Ok(PulleySystem::fixed()),
        Arrangement::Movable => Ok(PulleySystem::movable()),
        Arrangement::BlockAndTackle => {
            PulleySystem::block_and_tackle(s.supporting_ropes).map_err(|e| e.to_string())
        }
    }
}

/// Analyse the configured system and format the full readout, mapping any
/// domain error to a display string. Extracted so it is unit-testable.
fn compute(s: &PulleyWorkbenchState) -> Result<String, String> {
    let sys = system_of(s)?;
    let ma = sys.ideal_mechanical_advantage();
    let vr = sys.velocity_ratio();
    let f_ideal = ideal_effort(sys, s.load_n).map_err(|e| e.to_string())?;
    let f_real = real_effort(sys, s.load_n, s.efficiency).map_err(|e| e.to_string())?;
    let ama = actual_mechanical_advantage(sys, s.efficiency).map_err(|e| e.to_string())?;
    // The load that real effort raises is the round-trip of `real_effort`.
    let raised = load_from_effort(sys, f_real, s.efficiency).map_err(|e| e.to_string())?;

    // Energy story for raising the load through `lift_distance_m`: the
    // operator pulls `VR * s` of rope, doing `W_in` of work for `W_out`
    // useful work on the load, with the difference lost to friction.
    let rope_pulled = effort_distance(sys, s.lift_distance_m).map_err(|e| e.to_string())?;
    let w_out = output_work(s.load_n, s.lift_distance_m).map_err(|e| e.to_string())?;
    let w_in = input_work(s.load_n, s.lift_distance_m, s.efficiency).map_err(|e| e.to_string())?;
    let w_loss = work_lost(s.load_n, s.lift_distance_m, s.efficiency).map_err(|e| e.to_string())?;

    Ok(format!(
        "arrangement     : {arrangement}\n\
         supporting ropes: {n}\n\
         load W          : {load:.1} N\n\
         efficiency η    : {eta:.3}\n\n\
         ideal MA (= n)  : {ma:.2}\n\
         velocity ratio  : {vr:.2}\n\
         ideal effort    : {f_ideal:.2} N\n\
         real effort     : {f_real:.2} N\n\
         actual MA (η·MA): {ama:.2}\n\
         raises load     : {raised:.1} N\n\n\
         lift distance s : {lift:.2} m\n\
         rope pulled     : {rope_pulled:.2} m\n\
         output work W_out: {w_out:.1} J\n\
         input work W_in : {w_in:.1} J\n\
         work lost (frict): {w_loss:.1} J",
        arrangement = arrangement_label(s.arrangement),
        n = sys.supporting_ropes(),
        load = s.load_n,
        eta = s.efficiency,
        lift = s.lift_distance_m,
    ))
}

/// A short human label for an [`Arrangement`].
fn arrangement_label(a: Arrangement) -> &'static str {
    match a {
        Arrangement::Fixed => "fixed pulley",
        Arrangement::Movable => "movable pulley",
        Arrangement::BlockAndTackle => "block and tackle",
    }
}

/// Append an `x`-axis sheave cylinder (a pulley wheel lying with its axle
/// along `x`) to the buffers: a tube of radius `r` and axial `length`
/// starting at `base`, with `seg` angular facets.
fn push_cyl_x(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    base: Vector3<f64>,
    length: f64,
    r: f64,
    seg: usize,
) {
    let (x0, x1) = (base.x, base.x + length);
    let lo = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(x0, base.y + r * a.cos(), base.z + r * a.sin()));
    }
    let hi = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(x1, base.y + r * a.cos(), base.z + r * a.sin()));
    }
    for j in 0..seg {
        let jn = (j + 1) % seg;
        tris.extend_from_slice(&[
            lo + j,
            hi + j,
            hi + jn,
            lo + j,
            hi + jn,
            lo + jn,
            lo + j,
            hi + jn,
            hi + j,
            lo + j,
            lo + jn,
            hi + jn,
        ]);
    }
}

/// Build the pulley assembly as a triangle [`Mesh`] — a fixed block (a row
/// of sheave cylinders high up) over a movable block (a row lower down),
/// representing the `n` supporting rope segments by splitting the sheaves
/// between the two blocks. Representative geometry (not to scale; the
/// mechanical-advantage numbers are the `valenx-pulley` result). `None`
/// for an invalid configuration.
fn pulley_solid_mesh(s: &PulleyWorkbenchState) -> Option<Mesh> {
    let sys = system_of(s).ok()?;
    let n = sys.supporting_ropes();

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Split the supporting segments between a fixed (upper) and movable
    // (lower) block: ceil(n/2) sheaves up, floor(n/2) (at least one) down.
    let fixed_sheaves = n.div_ceil(2).max(1);
    let movable_sheaves = (n / 2).max(1);
    let r = 0.18;
    let width = 0.12;
    let pitch = 0.5;

    // Fixed block: a row of sheaves along y at the top (+z).
    let fixed_span = (fixed_sheaves.saturating_sub(1)) as f64 * pitch;
    for i in 0..fixed_sheaves {
        let y = i as f64 * pitch - 0.5 * fixed_span;
        push_cyl_x(
            &mut nodes,
            &mut tris,
            Vector3::new(-0.5 * width, y, 0.9),
            width,
            r,
            28,
        );
    }

    // Movable block: a row of sheaves lower down (the load hangs here).
    let movable_span = (movable_sheaves.saturating_sub(1)) as f64 * pitch;
    for i in 0..movable_sheaves {
        let y = i as f64 * pitch - 0.5 * movable_span;
        push_cyl_x(
            &mut nodes,
            &mut tris,
            Vector3::new(-0.5 * width, y, 0.1),
            width,
            r,
            28,
        );
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-pulley");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D pulley assembly and load it into the central viewport.
fn load_pulley_3d(app: &mut ValenxApp) {
    let Some(mesh) = pulley_solid_mesh(&app.pulley) else {
        app.pulley.error =
            Some("pulley parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<pulley>/valenx-pulley"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = PulleyWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_advantage_and_effort() {
        let mut s = PulleyWorkbenchState::default();
        run_pulley(&mut s);
        assert!(
            s.error.is_none(),
            "default tackle should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("ideal MA"));
        assert!(s.result.contains("actual MA"));
        // n = 4 -> MA 4.00, ideal effort 1200 / 4 = 300.00 N, real effort
        // 1200 / (4 * 0.9) = 333.33 N, actual MA 4 * 0.9 = 3.60.
        assert!(s.result.contains("4.00"));
        assert!(s.result.contains("300.00"));
        assert!(s.result.contains("333.33"));
        assert!(s.result.contains("3.60"));
    }

    #[test]
    fn fixed_pulley_has_unit_advantage() {
        let mut s = PulleyWorkbenchState {
            arrangement: Arrangement::Fixed,
            load_n: 500.0,
            efficiency: 1.0,
            ..Default::default()
        };
        run_pulley(&mut s);
        assert!(s.error.is_none());
        // MA = 1, ideal effort == load == 500.0 N.
        assert!(s.result.contains("ideal MA (= n)  : 1.00"));
        assert!(s.result.contains("ideal effort    : 500.00 N"));
    }

    #[test]
    fn analyze_rejects_efficiency_above_one() {
        let mut s = PulleyWorkbenchState {
            efficiency: 1.5,
            ..Default::default()
        };
        run_pulley(&mut s);
        assert!(s.error.is_some());
    }

    /// Ground truth: the ideal mechanical advantage of a block-and-tackle
    /// equals its number of supporting rope segments, and the ideal effort
    /// is the load divided by that count — checked against hand values.
    #[test]
    fn ideal_advantage_equals_supporting_ropes() {
        let load = 600.0;
        for n in 1..=8u32 {
            let sys = PulleySystem::block_and_tackle(n).unwrap();
            assert!((sys.ideal_mechanical_advantage() - f64::from(n)).abs() < 1e-12);
            let f = ideal_effort(sys, load).unwrap();
            let expected = load / f64::from(n);
            assert!((f - expected).abs() < 1e-9, "n = {n}");
        }
    }

    /// Ground truth for the energy block surfaced from `valenx-pulley`'s
    /// `effort_distance` / `output_work` / `input_work` / `work_lost`.
    /// Defaults: n = 4 (VR = 4), W = 1200 N, η = 0.9, lifting s = 2 m.
    /// Hand values: rope pulled = VR·s = 4·2 = 8.00 m; output work
    /// W_out = W·s = 1200·2 = 2400 J; input work W_in = W_out/η =
    /// 2400/0.9 = 2666.67 J; work lost = W_in − W_out = 266.67 J.
    #[test]
    fn analyze_default_reports_rope_pulled_and_work() {
        let mut s = PulleyWorkbenchState::default();
        run_pulley(&mut s);
        assert!(
            s.error.is_none(),
            "default tackle should analyze: {:?}",
            s.error
        );
        assert!(
            s.result.contains("rope pulled     : 8.00 m"),
            "{}",
            s.result
        );
        assert!(
            s.result.contains("output work W_out: 2400.0 J"),
            "{}",
            s.result
        );
        assert!(
            s.result.contains("input work W_in : 2666.7 J"),
            "{}",
            s.result
        );
        assert!(
            s.result.contains("work lost (frict): 266.7 J"),
            "{}",
            s.result
        );
    }

    /// Cross-check the energy identity directly against `valenx-pulley`:
    /// the surfaced rope length is VR·s and W_in − W_out equals the work
    /// lost, with no loss at η = 1.
    #[test]
    fn energy_block_matches_crate_ground_truth() {
        let load = 1200.0_f64;
        let s_load = 2.0_f64;
        let eta = 0.9_f64;
        let sys = PulleySystem::block_and_tackle(4).unwrap();

        let rope = effort_distance(sys, s_load).unwrap();
        assert!((rope - 8.0).abs() < 1e-12, "rope = {rope}");

        let w_out = output_work(load, s_load).unwrap();
        let w_in = input_work(load, s_load, eta).unwrap();
        let w_loss = work_lost(load, s_load, eta).unwrap();
        assert!((w_out - 2400.0).abs() < 1e-9, "w_out = {w_out}");
        assert!((w_in - 2400.0 / 0.9).abs() < 1e-9, "w_in = {w_in}");
        assert!((w_loss - (w_in - w_out)).abs() < 1e-9, "w_loss = {w_loss}");

        // No friction loss for the perfect machine.
        let lossless = work_lost(load, s_load, 1.0).unwrap();
        assert!(lossless.abs() < 1e-9, "lossless = {lossless}");
    }

    #[test]
    fn pulley_mesh_for_default_is_nonempty_and_in_range() {
        let s = PulleyWorkbenchState::default();
        let mesh = pulley_solid_mesh(&s).expect("default tackle yields a solid");
        assert!(mesh.nodes.len() > 8, "expected several sheave cylinders");
        let count = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < count));
        }
    }

    #[test]
    fn pulley_mesh_none_for_invalid() {
        let s = PulleyWorkbenchState {
            arrangement: Arrangement::BlockAndTackle,
            supporting_ropes: 0,
            ..Default::default()
        };
        assert!(pulley_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_pulley_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_pulley_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_pulley_workbench = true;
        run_pulley(&mut app.pulley);
        draw_workbench(&mut app);
    }
}
