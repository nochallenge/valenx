//! The right-side **Statics Workbench** panel — native 2-D rigid-body
//! statics over `valenx-statics`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_statics_workbench`,
//! toggled from the View menu. The form lays out a simply-supported beam
//! (a pin at the left support, a roller at the right) carrying up to three
//! downward point loads; "Analyze" solves the support reactions in closed
//! form — `R_A` at the pin and `R_B` at the roller from `sum Fy = 0` and
//! `sum M = 0` — and reports the equivalent single load resultant and an
//! independent `sum Fx = sum Fy = sum M = 0` equilibrium check, and
//! "Show 3-D" loads a representative horizontal beam resting on two
//! supports with downward load arrows into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_statics::{PointLoad, SimpleBeam};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// How many of the three point-load slots are active.
#[derive(Debug, Clone, Copy, PartialEq)]
enum LoadCount {
    /// A single point load.
    One,
    /// Two point loads.
    Two,
    /// Three point loads.
    Three,
}

/// Persistent form + result state for the Statics Workbench.
pub struct StaticsWorkbenchState {
    /// Span between the supports (m); the pin sits at `0`, the roller at
    /// this position.
    span_m: f64,
    /// How many of the three load slots are active.
    load_count: LoadCount,
    /// Position of load 1 along the beam axis (m from the pin).
    load1_pos_m: f64,
    /// Downward magnitude of load 1 (N, positive = downward).
    load1_down_n: f64,
    /// Position of load 2 along the beam axis (m from the pin).
    load2_pos_m: f64,
    /// Downward magnitude of load 2 (N, positive = downward).
    load2_down_n: f64,
    /// Position of load 3 along the beam axis (m from the pin).
    load3_pos_m: f64,
    /// Downward magnitude of load 3 (N, positive = downward).
    load3_down_n: f64,
    /// Formatted reaction readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D beam solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for StaticsWorkbenchState {
    fn default() -> Self {
        // A 10 m simply-supported beam with a single 1000 N load at
        // mid-span: a textbook case that splits evenly, R_A = R_B = 500 N.
        Self {
            span_m: 10.0,
            load_count: LoadCount::One,
            load1_pos_m: 5.0,
            load1_down_n: 1000.0,
            load2_pos_m: 7.0,
            load2_down_n: 600.0,
            load3_pos_m: 9.0,
            load3_down_n: 400.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Statics Workbench right-side panel. A no-op when the
/// `show_statics_workbench` toggle is off.
pub fn draw_statics_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_statics_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_statics_workbench",
        "Statics",
        |app, ui| {
            ui.label(
                egui::RichText::new("native simply-supported beam reactions · valenx-statics")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.statics;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Beam").strong());
                    ui.horizontal(|ui| {
                        ui.label("span — pin→roller (m)");
                        ui.add(egui::DragValue::new(&mut s.span_m).speed(0.25));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Loads").strong());
                    ui.horizontal(|ui| {
                        ui.label("count");
                        ui.radio_value(&mut s.load_count, LoadCount::One, "1");
                        ui.radio_value(&mut s.load_count, LoadCount::Two, "2");
                        ui.radio_value(&mut s.load_count, LoadCount::Three, "3");
                    });

                    ui.horizontal(|ui| {
                        ui.label("load 1  x (m)");
                        ui.add(egui::DragValue::new(&mut s.load1_pos_m).speed(0.1));
                        ui.label("P (N)");
                        ui.add(egui::DragValue::new(&mut s.load1_down_n).speed(10.0));
                    });
                    if matches!(s.load_count, LoadCount::Two | LoadCount::Three) {
                        ui.horizontal(|ui| {
                            ui.label("load 2  x (m)");
                            ui.add(egui::DragValue::new(&mut s.load2_pos_m).speed(0.1));
                            ui.label("P (N)");
                            ui.add(egui::DragValue::new(&mut s.load2_down_n).speed(10.0));
                        });
                    }
                    if matches!(s.load_count, LoadCount::Three) {
                        ui.horizontal(|ui| {
                            ui.label("load 3  x (m)");
                            ui.add(egui::DragValue::new(&mut s.load3_pos_m).speed(0.1));
                            ui.label("P (N)");
                            ui.add(egui::DragValue::new(&mut s.load3_down_n).speed(10.0));
                        });
                    }

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_statics(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative horizontal beam resting on a pin and a roller support, with downward load arrows, as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Reactions").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_statics_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.statics` borrow is
    // released here): build the beam's 3-D solid and load it.
    if app.statics.show_3d_request {
        app.statics.show_3d_request = false;
        load_beam_3d(app);
    }
}

/// Validate the form, solve the beam and format the readout.
fn run_statics(s: &mut StaticsWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// The active point loads implied by the form's load count.
fn active_loads(s: &StaticsWorkbenchState) -> Vec<PointLoad> {
    let mut loads = vec![PointLoad::vertical(s.load1_pos_m, s.load1_down_n)];
    if matches!(s.load_count, LoadCount::Two | LoadCount::Three) {
        loads.push(PointLoad::vertical(s.load2_pos_m, s.load2_down_n));
    }
    if matches!(s.load_count, LoadCount::Three) {
        loads.push(PointLoad::vertical(s.load3_pos_m, s.load3_down_n));
    }
    loads
}

/// Assemble the [`SimpleBeam`] (pin at `0`, roller at `span`) with the
/// active loads. The quantities both the readout and the 3-D gate need.
/// Extracted so it is unit-testable and shared.
fn build_beam(s: &StaticsWorkbenchState) -> Result<SimpleBeam, String> {
    let mut beam = SimpleBeam::new(0.0, s.span_m).map_err(|e| e.to_string())?;
    for load in active_loads(s) {
        beam.add_load(load).map_err(|e| e.to_string())?;
    }
    Ok(beam)
}

/// Solve the beam and format the full readout, mapping any domain error to
/// a display string. Extracted so it is unit-testable.
fn compute(s: &StaticsWorkbenchState) -> Result<String, String> {
    let beam = build_beam(s)?;
    let r = beam.reactions();

    // Independent equilibrium check on the assembled load + reaction
    // system: sum Fx = sum Fy = sum M = 0.
    let sys = beam.equilibrium_system();
    let balanced = sys.is_in_equilibrium(1e-6);
    let check = if balanced { "OK" } else { "FAIL" };
    // The actual residual net moment of the solved load+reaction system
    // about the origin (the number behind the OK/FAIL verdict). At the
    // closed-form solution this is ~0 N·m; a non-zero value would flag a
    // numerical or modelling error in the assembled system.
    let residual_moment = sys.sum_moment_origin();

    // The single equivalent load resultant (magnitude at the centroid).
    let resultant = beam
        .vertical_load_resultant()
        .map(|res| format!("{:.1} N at x = {:.3} m", res.down, res.position))
        .unwrap_or_else(|| "none (net zero load)".to_string());

    Ok(format!(
        "span (pin→roller): {:.3} m\n\
         total load   : {:.1} N\n\
         load resultant: {resultant}\n\n\
         R_A pin    : {:.2} N up\n\
         R_B roller : {:.2} N up\n\
         sum R       : {:.2} N up\n\
         pin horiz   : {:.2} N\n\n\
         equilibrium (sum Fx,Fy,M = 0): {check}\n\
         residual sum M @ origin: {residual_moment:.3} N·m",
        beam.span(),
        beam.total_down(),
        r.pin_vertical,
        r.roller_vertical,
        r.total_vertical(),
        r.pin_horizontal,
    ))
}

/// Append an outward-facing box (centre `c`, half-extents `h`) to the
/// buffers.
fn push_box(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    c: Vector3<f64>,
    h: Vector3<f64>,
) {
    let base = nodes.len();
    let signs = [
        (-1.0, -1.0, -1.0),
        (1.0, -1.0, -1.0),
        (1.0, 1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
        (1.0, -1.0, 1.0),
        (1.0, 1.0, 1.0),
        (-1.0, 1.0, 1.0),
    ];
    for (sx, sy, sz) in signs {
        nodes.push(c + Vector3::new(sx * h.x, sy * h.y, sz * h.z));
    }
    let faces = [
        [1, 2, 6, 5],
        [0, 4, 7, 3],
        [3, 7, 6, 2],
        [0, 1, 5, 4],
        [4, 5, 6, 7],
        [0, 3, 2, 1],
    ];
    for f in faces {
        tris.extend_from_slice(&[
            base + f[0],
            base + f[1],
            base + f[2],
            base + f[0],
            base + f[2],
            base + f[3],
        ]);
    }
}

/// Append an upward-pointing triangular-prism support (a wedge whose apex
/// touches the beam at `apex`), extruded along `y` by half-depth `half_y`.
/// Models the pin / roller as the textbook knife-edge support triangle.
#[allow(clippy::too_many_arguments)]
fn push_support(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    apex_x: f64,
    apex_z: f64,
    half_width: f64,
    height: f64,
    half_y: f64,
) {
    let base = nodes.len();
    let base_z = apex_z - height;
    // Two triangular end faces (front at +y, back at -y), apex on top.
    // Front face: 0 apex, 1 base-left, 2 base-right.
    nodes.push(Vector3::new(apex_x, half_y, apex_z));
    nodes.push(Vector3::new(apex_x - half_width, half_y, base_z));
    nodes.push(Vector3::new(apex_x + half_width, half_y, base_z));
    // Back face: 3 apex, 4 base-left, 5 base-right.
    nodes.push(Vector3::new(apex_x, -half_y, apex_z));
    nodes.push(Vector3::new(apex_x - half_width, -half_y, base_z));
    nodes.push(Vector3::new(apex_x + half_width, -half_y, base_z));
    // Front triangle, back triangle, and the three quad side faces.
    let quads = [(1usize, 2usize, 5usize, 4usize), (0, 1, 4, 3), (2, 0, 3, 5)];
    tris.extend_from_slice(&[base, base + 1, base + 2]);
    tris.extend_from_slice(&[base + 3, base + 5, base + 4]);
    for (a, b, c, d) in quads {
        tris.extend_from_slice(&[base + a, base + b, base + c, base + a, base + c, base + d]);
    }
}

/// Append a downward load arrow at beam position `x_frac` (a fraction of
/// the rendered beam length): a thin vertical shaft above the beam topped
/// by a wider arrowhead box pointing down at the beam.
fn push_load_arrow(nodes: &mut Vec<Vector3<f64>>, tris: &mut Vec<usize>, x: f64, beam_top_z: f64) {
    // Shaft.
    push_box(
        nodes,
        tris,
        Vector3::new(x, 0.0, beam_top_z + 0.55),
        Vector3::new(0.03, 0.03, 0.35),
    );
    // Arrowhead just above the beam.
    push_box(
        nodes,
        tris,
        Vector3::new(x, 0.0, beam_top_z + 0.13),
        Vector3::new(0.09, 0.09, 0.13),
    );
}

/// Build the beam-on-supports scene as a triangle [`Mesh`] — a horizontal
/// beam bar resting on a pin and a roller support, with a downward arrow
/// over each active load. The geometry is representative (the beam is
/// drawn to a fixed rendered length with loads placed by span fraction;
/// the reaction numbers are the `valenx-statics` result). `None` for an
/// invalid configuration.
fn beam_solid_mesh(s: &StaticsWorkbenchState) -> Option<Mesh> {
    let beam = build_beam(s).ok()?;
    let span = beam.span();

    // Fixed rendered beam length; supports sit a little inboard of the
    // ends. Real load positions map onto [x_min, x_max] by span fraction.
    let length = 4.0;
    let x_min = -length / 2.0;
    let x_max = length / 2.0;
    let beam_center_z = 1.0;
    let beam_half_z = 0.08;
    let beam_top_z = beam_center_z + beam_half_z;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // The beam bar.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, beam_center_z),
        Vector3::new(length / 2.0, 0.25, beam_half_z),
    );

    // Pin support under the left end (apex touches the beam underside).
    push_support(
        &mut nodes,
        &mut tris,
        x_min + 0.2,
        beam_center_z - beam_half_z,
        0.35,
        0.9,
        0.22,
    );
    // Roller support under the right end.
    push_support(
        &mut nodes,
        &mut tris,
        x_max - 0.2,
        beam_center_z - beam_half_z,
        0.35,
        0.9,
        0.22,
    );

    // A downward arrow over each active load, placed by span fraction.
    for load in active_loads(s) {
        let frac = (load.position / span).clamp(0.0, 1.0);
        let x = x_min + frac * length;
        push_load_arrow(&mut nodes, &mut tris, x, beam_top_z);
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-statics");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D beam solid and load it into the central viewport.
fn load_beam_3d(app: &mut ValenxApp) {
    let Some(mesh) = beam_solid_mesh(&app.statics) else {
        app.statics.error = Some("beam parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<beam>/valenx-statics"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical statics workbench as a 3-D solid plus its
/// `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn statics_product() -> crate::WorkspaceProduct {
    let s = StaticsWorkbenchState::default();
    let mesh = beam_solid_mesh(&s).expect("canonical statics ⇒ beam solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<statics>/valenx-beam");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical statics ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Statics beam (reactions/shear/moment)".into(),
        lines,
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = StaticsWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_reactions_and_equilibrium() {
        let mut s = StaticsWorkbenchState::default();
        run_statics(&mut s);
        assert!(
            s.error.is_none(),
            "default beam should solve: {:?}",
            s.error
        );
        assert!(s.result.contains("R_A pin"));
        assert!(s.result.contains("R_B roller"));
        assert!(s.result.contains("equilibrium"));
        assert!(s.result.contains("OK"));
        // 1000 N at mid-span of a 10 m beam splits evenly: 500 N each.
        assert!(s.result.contains("500.00 N up"));
    }

    #[test]
    fn analyze_rejects_zero_span() {
        let mut s = StaticsWorkbenchState {
            span_m: 0.0,
            ..Default::default()
        };
        run_statics(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn ground_truth_lever_rule_off_center_load() {
        // Ground truth, hand-computed: a single 100 N load 8 m from the
        // pin on a 10 m beam. Moment about A gives R_B = P*a/L =
        // 100*8/10 = 80 N; vertical balance gives R_A = P - R_B = 20 N.
        let mut s = StaticsWorkbenchState {
            span_m: 10.0,
            load_count: LoadCount::One,
            load1_pos_m: 8.0,
            load1_down_n: 100.0,
            ..Default::default()
        };
        let beam = build_beam(&s).expect("valid beam");
        let r = beam.reactions();
        let r_b: f64 = r.roller_vertical;
        let r_a: f64 = r.pin_vertical;
        assert!((r_b - 80.0).abs() < 1e-9, "R_B {r_b}");
        assert!((r_a - 20.0).abs() < 1e-9, "R_A {r_a}");
        // And the formatted readout carries the hand-computed values.
        run_statics(&mut s);
        assert!(s.result.contains("80.00 N up"));
        assert!(s.result.contains("20.00 N up"));
    }

    #[test]
    fn ground_truth_residual_moment_is_zero_at_solution() {
        // Ground truth, hand-computed. Off-centre 100 N load 8 m from
        // the pin on a 10 m beam: R_by = 100*8/10 = 80, R_ay = 20.
        // The assembled load+reaction system (all forces at y = 0) has
        // residual moment about the origin sum(x_i * Fy_i):
        //   load   8 * (-100) = -800
        //   pin    0 *  (+20) =    0
        //   roller 10 * (+80) = +800
        //   sum                =    0.000 N·m  (the solution closes).
        let mut s = StaticsWorkbenchState {
            span_m: 10.0,
            load_count: LoadCount::One,
            load1_pos_m: 8.0,
            load1_down_n: 100.0,
            ..Default::default()
        };
        let beam = build_beam(&s).expect("valid beam");
        let resid = beam.equilibrium_system().sum_moment_origin();
        assert!(resid.abs() < 1e-9, "residual moment {resid}");
        // And the formatted readout carries the hand-computed residual.
        run_statics(&mut s);
        assert!(
            s.result.contains("residual sum M @ origin: 0.000 N·m"),
            "missing residual readout: {}",
            s.result
        );
    }

    #[test]
    fn three_loads_sum_into_reactions() {
        // With three loads active the total upward reaction equals the
        // total downward load (sum Fy = 0).
        let s = StaticsWorkbenchState {
            span_m: 10.0,
            load_count: LoadCount::Three,
            load1_pos_m: 2.0,
            load1_down_n: 300.0,
            load2_pos_m: 5.0,
            load2_down_n: 500.0,
            load3_pos_m: 8.0,
            load3_down_n: 200.0,
            ..Default::default()
        };
        let beam = build_beam(&s).expect("valid beam");
        assert_eq!(beam.loads.len(), 3);
        let total: f64 = beam.total_down();
        let r_sum: f64 = beam.reactions().total_vertical();
        assert!((total - 1000.0).abs() < 1e-9, "total {total}");
        assert!((r_sum - total).abs() < 1e-9, "sum R {r_sum}");
    }

    #[test]
    fn beam_mesh_for_default_is_nonempty_and_in_range() {
        let s = StaticsWorkbenchState::default();
        let mesh = beam_solid_mesh(&s).expect("default beam yields a solid");
        assert!(mesh.nodes.len() > 8, "expected beam + supports + arrow");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn beam_mesh_none_for_invalid() {
        let s = StaticsWorkbenchState {
            span_m: 0.0,
            ..Default::default()
        };
        assert!(beam_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_statics_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_statics_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_statics_workbench = true;
        run_statics(&mut app.statics);
        draw_workbench(&mut app);
    }
}
