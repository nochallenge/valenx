//! The right-side **Lever Workbench** panel — ideal rigid-lever mechanics
//! over `valenx-leverage`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_leverage_workbench`,
//! toggled from the View menu. The form sets the two arm lengths, an
//! applied effort force and an effort travel, plus the declared lever
//! class. "Analyze" applies the static moment-balance law
//! (`effort * effort_arm = load * load_arm`) and reports the mechanical
//! advantage, the inferred class, the balanced load and shared moment, the
//! load travel and the balance check; "Show 3-D" loads a representative
//! beam-on-fulcrum solid into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_leverage::{Lever, LeverClass};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Lever Workbench.
pub struct LeverageWorkbenchState {
    /// Distance from the fulcrum to the applied effort force (m, > 0).
    effort_arm_m: f64,
    /// Distance from the fulcrum to the resisting load force (m, > 0).
    load_arm_m: f64,
    /// Applied effort force (N).
    effort_n: f64,
    /// A target load to size the lever against, in the inverse direction:
    /// the readout reports the effort required to balance it (N, ≥ 0).
    target_load_n: f64,
    /// Travel of the effort point as the beam swings (m).
    effort_travel_m: f64,
    /// The declared lever class (a labelling aid; the analysis also infers
    /// the class from the arm ratio).
    declared_class: LeverClass,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D lever solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for LeverageWorkbenchState {
    fn default() -> Self {
        // A crowbar-style force multiplier: a 1.2 m effort arm over a
        // 0.3 m load arm gives MA = 4, so a 150 N effort balances a 600 N
        // load, and a 0.40 m effort swing moves the load 0.10 m. Sizing the
        // inverse direction, holding a 1000 N target load needs 250 N of
        // effort (1000 / MA).
        Self {
            effort_arm_m: 1.2,
            load_arm_m: 0.3,
            effort_n: 150.0,
            target_load_n: 1000.0,
            effort_travel_m: 0.4,
            declared_class: LeverClass::Second,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Lever Workbench right-side panel. A no-op when the
/// `show_leverage_workbench` toggle is off.
pub fn draw_leverage_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_leverage_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_leverage_workbench",
        "Lever",
        |app, ui| {
            ui.label(
                egui::RichText::new("ideal rigid-lever mechanics · valenx-leverage")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.leverage;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Geometry").strong());
                    ui.horizontal(|ui| {
                        ui.label("effort arm (m)");
                        ui.add(egui::DragValue::new(&mut s.effort_arm_m).speed(0.05));
                    });
                    ui.horizontal(|ui| {
                        ui.label("load arm (m)");
                        ui.add(egui::DragValue::new(&mut s.load_arm_m).speed(0.05));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Forces & motion").strong());
                    ui.horizontal(|ui| {
                        ui.label("effort (N)");
                        ui.add(egui::DragValue::new(&mut s.effort_n).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("target load (N)");
                        ui.add(egui::DragValue::new(&mut s.target_load_n).speed(1.0))
                            .on_hover_text(
                                "A load to size the lever against: the readout reports the effort required to balance it",
                            );
                    });
                    ui.horizontal(|ui| {
                        ui.label("effort travel (m)");
                        ui.add(egui::DragValue::new(&mut s.effort_travel_m).speed(0.01));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Declared class").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.declared_class, LeverClass::First, "1st");
                        ui.radio_value(&mut s.declared_class, LeverClass::Second, "2nd");
                        ui.radio_value(&mut s.declared_class, LeverClass::Third, "3rd");
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_leverage(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative lever — a beam balanced on a wedge fulcrum, with the effort and load arms to scale — as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Mechanics").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_leverage_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.leverage` borrow is
    // released here): build the lever's 3-D solid and load it.
    if app.leverage.show_3d_request {
        app.leverage.show_3d_request = false;
        load_lever_3d(app);
    }
}

/// Validate the form, evaluate the lever and format the readout.
fn run_leverage(s: &mut LeverageWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated [`Lever`] from the form geometry, the quantity both
/// the readout and the 3-D gate need. Extracted so it is unit-testable and
/// shared.
fn lever(s: &LeverageWorkbenchState) -> Result<Lever, String> {
    Lever::new(s.effort_arm_m, s.load_arm_m).map_err(|e| e.to_string())
}

/// Evaluate the lever and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &LeverageWorkbenchState) -> Result<String, String> {
    let lever = lever(s)?;
    let ma = lever.mechanical_advantage();
    let inferred = lever.class();
    let bal = lever.balance_load(s.effort_n).map_err(|e| e.to_string())?;
    // Inverse sizing: the effort required to balance the target load
    // (`required = target_load / MA`). The kinematic dual of the forward
    // `balance_load` above, answering "what effort do I need to hold this?"
    let req = lever
        .balance_effort(s.target_load_n)
        .map_err(|e| e.to_string())?;
    let load_travel = lever
        .load_displacement(s.effort_travel_m)
        .map_err(|e| e.to_string())?;
    // At the balanced pair the net moment is zero to floating-point
    // precision; report it and the balance flag as a built-in cross-check.
    let net = lever
        .net_moment(bal.effort, bal.load)
        .map_err(|e| e.to_string())?;
    let balanced = lever
        .is_balanced(bal.effort, bal.load, 1e-9)
        .map_err(|e| e.to_string())?;

    Ok(format!(
        "effort / load arm: {:.3} / {:.3} m\n\
         declared class   : {}\n\
         inferred class   : {}\n\n\
         mech advantage   : {:.3}\n\
         effort applied   : {:.2} N\n\
         balanced load    : {:.2} N\n\
         moment           : {:.3} N·m\n\n\
         target load      : {:.2} N\n\
         required effort  : {:.2} N\n\n\
         effort travel    : {:.3} m\n\
         load travel      : {:.3} m\n\n\
         net moment       : {:.3e} N·m\n\
         balanced         : {}",
        s.effort_arm_m,
        s.load_arm_m,
        s.declared_class.label(),
        inferred.label(),
        ma,
        bal.effort,
        bal.load,
        bal.moment,
        req.load,
        req.effort,
        s.effort_travel_m,
        load_travel,
        net,
        balanced,
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

/// Append a triangular-prism wedge fulcrum (apex up, length along `y`),
/// centred at `c` with half-length `half_len`, half-width `half_w` at the
/// base and height `height`, to the buffers.
fn push_wedge(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    c: Vector3<f64>,
    half_len: f64,
    half_w: f64,
    height: f64,
) {
    let base = nodes.len();
    // Two triangular end caps in the x-z plane (base corners + apex),
    // duplicated at -y and +y, apex riding the beam at z = c.z + height.
    let verts = [
        Vector3::new(c.x - half_w, c.y - half_len, c.z),
        Vector3::new(c.x + half_w, c.y - half_len, c.z),
        Vector3::new(c.x, c.y - half_len, c.z + height),
        Vector3::new(c.x - half_w, c.y + half_len, c.z),
        Vector3::new(c.x + half_w, c.y + half_len, c.z),
        Vector3::new(c.x, c.y + half_len, c.z + height),
    ];
    for v in verts {
        nodes.push(v);
    }
    // End caps (0,1,2) and (3,5,4), then the three quad side faces split
    // into triangles.
    let caps = [[0, 1, 2], [3, 5, 4]];
    for c in caps {
        tris.extend_from_slice(&[base + c[0], base + c[1], base + c[2]]);
    }
    let quads = [[0, 1, 4, 3], [1, 2, 5, 4], [2, 0, 3, 5]];
    for q in quads {
        tris.extend_from_slice(&[
            base + q[0],
            base + q[1],
            base + q[2],
            base + q[0],
            base + q[2],
            base + q[3],
        ]);
    }
}

/// Build the lever as a triangle [`Mesh`] — a slender beam balanced on a
/// wedge fulcrum, the effort arm and load arm drawn to their relative
/// lengths about the pivot. Representative geometry (the mechanics numbers
/// are the `valenx-leverage` result). `None` for an invalid configuration.
fn lever_solid_mesh(s: &LeverageWorkbenchState) -> Option<Mesh> {
    lever(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Place the fulcrum at the origin; the beam runs along x with the load
    // arm to -x and the effort arm to +x, scaled so the longer arm reads
    // about 1.0 in viewport units.
    let span = (s.effort_arm_m + s.load_arm_m).max(1e-6);
    let scale = 1.6 / span;
    let effort = s.effort_arm_m * scale;
    let load = s.load_arm_m * scale;
    let beam_centre_x = 0.5 * (effort - load);
    let beam_half_len = 0.5 * (effort + load);

    // Beam: long in x, slender in y/z, sitting on the fulcrum apex.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(beam_centre_x, 0.0, 0.32),
        Vector3::new(beam_half_len, 0.07, 0.03),
    );
    // Wedge fulcrum under the pivot, apex meeting the beam underside.
    push_wedge(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.0),
        0.1,
        0.18,
        0.29,
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-leverage");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D lever solid and load it into the central viewport.
fn load_lever_3d(app: &mut ValenxApp) {
    let Some(mesh) = lever_solid_mesh(&app.leverage) else {
        app.leverage.error =
            Some("lever parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<lever>/valenx-leverage"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical leverage workbench as a 3-D solid plus its
/// `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn leverage_product() -> crate::WorkspaceProduct {
    let s = LeverageWorkbenchState::default();
    let mesh = lever_solid_mesh(&s).expect("canonical leverage ⇒ lever solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<leverage>/valenx-lever");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical leverage ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Lever (mechanical advantage)".into(),
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
        let s = LeverageWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_advantage_load_and_balance() {
        let mut s = LeverageWorkbenchState::default();
        run_leverage(&mut s);
        assert!(
            s.error.is_none(),
            "default lever should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("mech advantage"));
        assert!(s.result.contains("balanced load"));
        assert!(s.result.contains("load travel"));
        // 1.2 m / 0.3 m = MA 4.000; 150 N effort -> 600.00 N load.
        assert!(s.result.contains("4.000"));
        assert!(s.result.contains("600.00"));
        // The reported pair is in static balance.
        assert!(s.result.contains("balanced"));
        assert!(s.result.contains(": true"));
    }

    #[test]
    fn analyze_reports_required_effort_for_target_load() {
        // Ground truth: required effort = target_load / MA. With the default
        // arms MA = 1.2 / 0.3 = 4, so holding a 1000 N target load needs
        // 1000 / 4 = 250 N of effort — the inverse-direction sizing the
        // forward `balanced load` line does not give.
        let s = LeverageWorkbenchState::default();
        let lever = lever(&s).expect("default arms are valid");
        let ma = lever.mechanical_advantage();
        let req = lever.balance_effort(s.target_load_n).unwrap();
        assert!(
            (req.effort - s.target_load_n / ma).abs() < 1e-9,
            "required effort {req_effort} must equal target_load / MA",
            req_effort = req.effort
        );
        let expected = s.target_load_n / ma;
        assert!(
            (expected - 250.0).abs() < 1e-9,
            "hand-computed required effort is 250 N, got {expected}"
        );
        // The readout surfaces the target load and the required effort, both
        // at 2-dp: 1000.00 N target -> 250.00 N effort.
        let out = compute(&s).expect("default lever computes");
        assert!(out.contains("required effort  : 250.00 N"), "got:\n{out}");
        assert!(out.contains("target load      : 1000.00 N"), "got:\n{out}");
    }

    #[test]
    fn analyze_rejects_zero_arm() {
        let mut s = LeverageWorkbenchState {
            load_arm_m: 0.0,
            ..Default::default()
        };
        run_leverage(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn mechanical_advantage_is_effort_arm_over_load_arm() {
        // Ground truth: MA is the arm ratio, hand-computed 1.2 / 0.3 = 4.0,
        // and at balance the load is the effort scaled by MA.
        let s = LeverageWorkbenchState::default();
        let lever = lever(&s).expect("default arms are valid");
        let ma = lever.mechanical_advantage();
        assert!((ma - 4.0).abs() < 1e-12);
        let bal = lever.balance_load(s.effort_n).unwrap();
        assert!((bal.load - s.effort_n * ma).abs() < 1e-9);
        // The balance law itself: effort * effort_arm == load * load_arm.
        assert!((bal.effort * lever.effort_arm - bal.load * lever.load_arm).abs() < 1e-9);
    }

    #[test]
    fn lever_mesh_for_default_is_nonempty_and_in_range() {
        let s = LeverageWorkbenchState::default();
        let mesh = lever_solid_mesh(&s).expect("default lever yields a solid");
        assert!(mesh.nodes.len() > 8, "expected beam + wedge fulcrum");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn lever_mesh_none_for_invalid() {
        let s = LeverageWorkbenchState {
            effort_arm_m: 0.0,
            ..Default::default()
        };
        assert!(lever_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_leverage_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_leverage_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_leverage_workbench = true;
        run_leverage(&mut app.leverage);
        draw_workbench(&mut app);
    }
}
