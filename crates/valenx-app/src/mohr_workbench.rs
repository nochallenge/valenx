//! The right-side **Mohr's Circle Workbench** panel — native 2-D
//! plane-stress transformation over `valenx-mohr`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_mohr_workbench`,
//! toggled from the View menu. The form sets a plane-stress state (the two
//! normal stresses `sx`, `sy` and the in-plane shear `txy`) plus a query
//! plane angle; "Analyze" computes the closed-form Mohr's-circle results —
//! mean normal stress (circle centre), radius, principal stresses
//! `s1 >= s2`, maximum in-plane and absolute (Tresca) shear, the principal
//! orientation, and the normal / shear stress transformed onto the query
//! plane — and "Show 3-D" loads a unit stress-element cube into the central
//! viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_mohr::StressState;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Mohr's Circle Workbench.
pub struct MohrWorkbenchState {
    /// Normal stress on the `x` face `sx` (MPa, positive in tension).
    sx: f64,
    /// Normal stress on the `y` face `sy` (MPa, positive in tension).
    sy: f64,
    /// In-plane shear stress `txy` (MPa).
    txy: f64,
    /// Query plane angle `theta` (degrees, counter-clockwise from the `x`
    /// axis to the outward normal of the cut plane).
    theta_deg: f64,
    /// Formatted results readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D stress element (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for MohrWorkbenchState {
    fn default() -> Self {
        // Classic worked example (Hibbeler / Gere): sx = -20, sy = 90,
        // txy = 60 MPa. Centre 35, radius ~81.39, so s1 ~ 116.39 and
        // s2 ~ -46.39 MPa, with the principal plane (carrying s1) at the
        // crate's quadrant-correct angle ~66.26 deg.
        Self {
            sx: -20.0,
            sy: 90.0,
            txy: 60.0,
            theta_deg: 30.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Mohr's Circle Workbench right-side panel. A no-op when the
/// `show_mohr_workbench` toggle is off.
pub fn draw_mohr_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_mohr_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_mohr_workbench",
        "Mohr's Circle",
        |app, ui| {
            ui.label(
                egui::RichText::new("native 2-D plane-stress transformation · valenx-mohr")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.mohr;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Plane-stress state").strong());
                    // Associate each numeric `DragValue` with its caption via
                    // `labelled_by`, so the spin button carries the caption as
                    // its accessibility / UI-Automation Name (egui clears a
                    // DragValue's own Name otherwise, leaving it anonymous to a
                    // screen reader / AI driver).
                    ui.horizontal(|ui| {
                        let l = ui.label("σx (MPa)");
                        ui.add(egui::DragValue::new(&mut s.sx).speed(1.0))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("σy (MPa)");
                        ui.add(egui::DragValue::new(&mut s.sy).speed(1.0))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("τxy (MPa)");
                        ui.add(egui::DragValue::new(&mut s.txy).speed(1.0))
                            .labelled_by(l.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Query plane").strong());
                    ui.horizontal(|ui| {
                        let l = ui.label("θ (deg, ccw)");
                        ui.add(egui::DragValue::new(&mut s.theta_deg).speed(1.0))
                            .labelled_by(l.id);
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_mohr(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a unit stress-element cube and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Results").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_mohr_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.mohr` borrow is
    // released here): build the stress element's 3-D solid and load it.
    if app.mohr.show_3d_request {
        app.mohr.show_3d_request = false;
        load_element_3d(app);
    }
}

/// Validate the form, evaluate the transformation and format the readout.
fn run_mohr(s: &mut MohrWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Evaluate the Mohr's-circle results and format the full readout, mapping
/// any domain error to a display string. Extracted so it is unit-testable.
fn compute(s: &MohrWorkbenchState) -> Result<String, String> {
    let state = StressState::new(s.sx, s.sy, s.txy).map_err(|e| e.to_string())?;
    let theta = s.theta_deg.to_radians();
    let on = state.stress_on_plane(theta).map_err(|e| e.to_string())?;

    let centre = state.mean_normal();
    let radius = state.radius();
    let p = state.principal_stresses();
    let tau_max = state.max_shear();
    let tau_abs = state.absolute_max_shear();
    let theta_p_deg = state.principal_angle().to_degrees();

    Ok(format!(
        "σx / σy        : {:.2} / {:.2} MPa\n\
         τxy            : {:.2} MPa\n\n\
         centre (σavg)  : {:.2} MPa\n\
         radius R       : {:.2} MPa\n\
         σ1 (max)       : {:.2} MPa\n\
         σ2 (min)       : {:.2} MPa\n\
         τ max in-plane : {:.2} MPa\n\
         τ abs max      : {:.2} MPa\n\
         principal angle: {:.2} deg\n\n\
         on plane θ = {:.2} deg\n\
         σ normal       : {:.2} MPa\n\
         τ shear        : {:.2} MPa",
        s.sx,
        s.sy,
        s.txy,
        centre,
        radius,
        p.s1,
        p.s2,
        tau_max,
        tau_abs,
        theta_p_deg,
        s.theta_deg,
        on.normal,
        on.shear,
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

/// Build the stress element as a triangle [`Mesh`] — a unit cube centred at
/// the origin representing the infinitesimal material element whose faces
/// carry the plane-stress state. Representative geometry (the stress numbers
/// are the `valenx-mohr` result). `None` for an invalid configuration.
fn element_solid_mesh(s: &MohrWorkbenchState) -> Option<Mesh> {
    StressState::new(s.sx, s.sy, s.txy).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // The stress-element cube.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(0.5, 0.5, 0.5),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-mohr");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D stress element and load it into the central viewport.
fn load_element_3d(app: &mut ValenxApp) {
    let Some(mesh) = element_solid_mesh(&app.mohr) else {
        app.mohr.error = Some("stress state is invalid — cannot build the 3-D element".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<element>/valenx-mohr"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical mohr workbench as a 3-D solid plus its
/// `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn mohr_product() -> crate::WorkspaceProduct {
    let s = MohrWorkbenchState::default();
    let mesh = element_solid_mesh(&s).expect("canonical mohr ⇒ stress-element solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<mohr>/valenx-element");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical mohr ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Mohr's circle (principal stresses)".into(),
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

    /// Tolerance for analytic float comparisons in the readout substrings.
    const EPS: f64 = 1e-9;

    #[test]
    fn default_state_is_idle() {
        let s = MohrWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_principals_and_shear() {
        let mut s = MohrWorkbenchState::default();
        run_mohr(&mut s);
        assert!(
            s.error.is_none(),
            "default state should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("centre (σavg)"));
        assert!(s.result.contains("σ1 (max)"));
        assert!(s.result.contains("σ2 (min)"));
        assert!(s.result.contains("τ abs max"));
        assert!(s.result.contains("on plane θ"));
        // Hibbeler case: centre 35.00, s1 ~ 116.39, s2 ~ -46.39 MPa.
        assert!(s.result.contains("35.00"));
        assert!(s.result.contains("116.39"));
        assert!(s.result.contains("-46.39"));
    }

    #[test]
    fn analyze_rejects_non_finite_stress() {
        let mut s = MohrWorkbenchState {
            txy: f64::NAN,
            ..Default::default()
        };
        run_mohr(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn principals_are_mean_plus_minus_radius_ground_truth() {
        // Ground truth: for sx = -20, sy = 90, txy = 60, the principal
        // stresses are the circle centre plus/minus the radius, with
        // R = sqrt(((sx - sy)/2)^2 + txy^2). Hand-computed:
        //   centre = (-20 + 90)/2 = 35
        //   R      = sqrt((-55)^2 + 60^2) = sqrt(3025 + 3600) = sqrt(6625)
        let state = StressState::new(-20.0, 90.0, 60.0).unwrap();
        let half_diff = (-20.0_f64 - 90.0) / 2.0;
        let r_hand = (half_diff * half_diff + 60.0_f64 * 60.0).sqrt();
        let centre_hand = 35.0;
        let p = state.principal_stresses();
        assert!((p.s1 - (centre_hand + r_hand)).abs() < EPS, "s1 = {}", p.s1);
        assert!((p.s2 - (centre_hand - r_hand)).abs() < EPS, "s2 = {}", p.s2);
        assert!((state.radius() - r_hand).abs() < EPS);
        assert!((state.mean_normal() - centre_hand).abs() < EPS);
    }

    #[test]
    fn element_mesh_for_default_is_nonempty_and_in_range() {
        let s = MohrWorkbenchState::default();
        let mesh = element_solid_mesh(&s).expect("default state yields a solid");
        assert_eq!(mesh.nodes.len(), 8, "a single stress-element cube");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn element_mesh_none_for_invalid() {
        let s = MohrWorkbenchState {
            sx: f64::INFINITY,
            ..Default::default()
        };
        assert!(element_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_mohr_workbench(app, ctx);
        });
    }

    /// As [`draw_workbench`], but with accesskit enabled, returning the emitted
    /// accessibility tree nodes — the same tree a screen reader / AI driver
    /// consumes. `accesskit` is re-exported by egui, so no extra dependency.
    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_mohr_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_mohr_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_mohr_workbench = true;
        run_mohr(&mut app.mohr);
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Each plane-stress / query-plane DragValue is a SpinButton; each must
        // be `labelled_by` its caption (egui clears a DragValue's own Name), so
        // an AI / screen reader can find the control by the caption text.
        let mut app = ValenxApp::default();
        app.show_mohr_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // σx, σy, τxy, θ.
        assert!(
            spin_buttons.len() >= 4,
            "expected the Mohr numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every Mohr DragValue must be labelled_by its caption (AI-drivable name)"
        );

        for caption in ["σx (MPa)", "τxy (MPa)", "θ (deg, ccw)"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }
}
