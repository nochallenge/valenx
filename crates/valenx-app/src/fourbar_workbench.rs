//! The right-side **Four-Bar Linkage Workbench** panel — native planar
//! four-bar mechanism kinematics over `valenx-kinematics`.
//!
//! Mirrors the Wind Turbine / Rail workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_fourbar_workbench`,
//! toggled from the View menu. The form drives a
//! [`valenx_kinematics::FourBar`]; "Analyze" reports the Grashof class and,
//! at a chosen crank angle, the solved pin positions, the coupler / rocker
//! angles and the transmission angle, and "Show 3-D linkage" loads the four
//! links at that pose as a solid into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_kinematics::{Assembly, FourBar, GrashofClass};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Four-Bar Linkage Workbench.
pub struct FourBarWorkbenchState {
    /// Ground / frame link length `r1` (O2-O4).
    ground: f64,
    /// Crank (input) link length `r2`.
    crank: f64,
    /// Coupler link length `r3`.
    coupler: f64,
    /// Rocker (output) link length `r4`.
    rocker: f64,
    /// Input crank angle `theta2` (degrees, CCW from +x).
    crank_angle_deg: f64,
    /// Take the crossed assembly branch (else the open branch).
    crossed: bool,
    /// Formatted readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D linkage solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for FourBarWorkbenchState {
    fn default() -> Self {
        // A Grashof crank-rocker: the shortest link (the crank) fully
        // rotates. s + l = 1 + 4 = 5 < p + q = 3 + 3.5 = 6.5.
        Self {
            ground: 4.0,
            crank: 1.0,
            coupler: 3.5,
            rocker: 3.0,
            crank_angle_deg: 60.0,
            crossed: false,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Four-Bar Linkage Workbench right-side panel. A no-op when the
/// `show_fourbar_workbench` toggle is off.
pub fn draw_fourbar_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_fourbar_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_fourbar_workbench",
        "Four-Bar Linkage",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native planar four-bar mechanism kinematics · valenx-kinematics",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.fourbar;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Link lengths").strong());
                    ui.horizontal(|ui| {
                        let g = ui.label("ground r1");
                        ui.add(egui::DragValue::new(&mut s.ground).speed(0.05))
                            .labelled_by(g.id);
                    });
                    ui.horizontal(|ui| {
                        let cr = ui.label("crank r2");
                        ui.add(egui::DragValue::new(&mut s.crank).speed(0.05))
                            .labelled_by(cr.id);
                    });
                    ui.horizontal(|ui| {
                        let co = ui.label("coupler r3");
                        ui.add(egui::DragValue::new(&mut s.coupler).speed(0.05))
                            .labelled_by(co.id);
                    });
                    ui.horizontal(|ui| {
                        let ro = ui.label("rocker r4");
                        ui.add(egui::DragValue::new(&mut s.rocker).speed(0.05))
                            .labelled_by(ro.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Pose").strong());
                    ui.horizontal(|ui| {
                        let th = ui.label("crank angle θ2 (°)");
                        ui.add(egui::DragValue::new(&mut s.crank_angle_deg).speed(1.0))
                            .labelled_by(th.id);
                    });
                    ui.checkbox(&mut s.crossed, "crossed assembly branch");

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_fourbar(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D linkage").strong())
                        .on_hover_text(
                            "Build the four links at the solved pose as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Kinematics").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_fourbar_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.fourbar` borrow is
    // released here): build the linkage's 3-D solid and load it.
    if app.fourbar.show_3d_request {
        app.fourbar.show_3d_request = false;
        load_linkage_3d(app);
    }
}

/// Build a validated [`FourBar`] from the form, mapping the domain error to a
/// display string.
fn build_fourbar(s: &FourBarWorkbenchState) -> Result<FourBar, String> {
    FourBar::new(s.ground, s.crank, s.coupler, s.rocker).map_err(|e| e.to_string())
}

/// The assembly branch the form currently selects.
fn assembly(s: &FourBarWorkbenchState) -> Assembly {
    if s.crossed {
        Assembly::Crossed
    } else {
        Assembly::Open
    }
}

/// Validate the form, solve the linkage and format the readout.
fn run_fourbar(s: &mut FourBarWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Compute the full readout, mapping any domain error to a display string.
/// Extracted so it is unit-testable.
fn compute(s: &FourBarWorkbenchState) -> Result<String, String> {
    let fb = build_fourbar(s)?;
    let class = match fb.grashof_class() {
        GrashofClass::Crank => "Grashof crank (a link fully rotates)",
        GrashofClass::Change => "change-point (collinear / folding)",
        GrashofClass::DoubleRocker => "non-Grashof double-rocker",
    };
    let pose = fb
        .solve(s.crank_angle_deg.to_radians(), assembly(s))
        .map_err(|e| e.to_string())?;
    let a = pose.pin_a;
    let b = pose.pin_b;
    // Transmission angle: the angle at the coupler-rocker pin B between the
    // coupler (B->A) and the rocker (B->O4), clamped for round-off.
    let u = [a[0] - b[0], a[1] - b[1]];
    let v = [s.ground - b[0], -b[1]];
    let un = (u[0] * u[0] + u[1] * u[1]).sqrt();
    let vn = (v[0] * v[0] + v[1] * v[1]).sqrt();
    let mu_deg = if un > 0.0 && vn > 0.0 {
        ((u[0] * v[0] + u[1] * v[1]) / (un * vn))
            .clamp(-1.0, 1.0)
            .acos()
            .to_degrees()
    } else {
        0.0
    };
    Ok(format!(
        "ground r1 : {:.3}\n\
         crank r2  : {:.3}\n\
         coupler r3: {:.3}\n\
         rocker r4 : {:.3}\n\n\
         Grashof   : {}\n\
         assembly  : {}\n\n\
         at crank θ2 = {:.1}°:\n\
         crank pin A  : ({:.3}, {:.3})\n\
         coupler pin B: ({:.3}, {:.3})\n\
         coupler θ3 : {:.1}°\n\
         rocker θ4  : {:.1}°\n\
         transmission μ: {:.1}°",
        s.ground,
        s.crank,
        s.coupler,
        s.rocker,
        class,
        if s.crossed { "crossed" } else { "open" },
        s.crank_angle_deg,
        a[0],
        a[1],
        b[0],
        b[1],
        pose.coupler_angle.to_degrees(),
        pose.rocker_angle.to_degrees(),
        mu_deg,
    ))
}

/// Append an outward-facing box (centre `c`, half-extents `h`) to the
/// buffers (used for the pivot / pin markers).
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

/// Append a rectangular bar (a box oriented along `p -> q` in the z = 0
/// plane) of half-width `half_w` (in-plane) and half-thickness `half_t`
/// (out-of-plane). Faces are emitted double-sided. Used for the links.
fn push_bar(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    p: Vector3<f64>,
    q: Vector3<f64>,
    half_w: f64,
    half_t: f64,
) {
    let dir = q - p;
    let len = dir.norm();
    if len < 1e-9 {
        return;
    }
    let d = dir / len;
    let perp = Vector3::new(-d.y, d.x, 0.0) * half_w;
    let zt = Vector3::new(0.0, 0.0, half_t);
    let base = nodes.len();
    nodes.push(p - perp - zt); // 0
    nodes.push(q - perp - zt); // 1
    nodes.push(q + perp - zt); // 2
    nodes.push(p + perp - zt); // 3
    nodes.push(p - perp + zt); // 4
    nodes.push(q - perp + zt); // 5
    nodes.push(q + perp + zt); // 6
    nodes.push(p + perp + zt); // 7
    let faces = [
        [0, 1, 2, 3],
        [4, 5, 6, 7],
        [0, 3, 7, 4],
        [1, 2, 6, 5],
        [0, 1, 5, 4],
        [3, 2, 6, 7],
    ];
    for f in faces {
        let (a, b, c, e) = (base + f[0], base + f[1], base + f[2], base + f[3]);
        tris.extend_from_slice(&[a, b, c, a, c, e, a, c, b, a, e, c]);
    }
}

/// Build the four-bar linkage at its solved pose as a triangle [`Mesh`] — the
/// ground, crank, coupler and rocker links as bars in the z = 0 plane plus
/// pivot / pin markers. `None` for an invalid linkage, or one that cannot
/// assemble at the chosen crank angle.
fn linkage_solid_mesh(s: &FourBarWorkbenchState) -> Option<Mesh> {
    let fb = build_fourbar(s).ok()?;
    let pose = fb.solve(s.crank_angle_deg.to_radians(), assembly(s)).ok()?;
    let o2 = Vector3::zeros();
    let o4 = Vector3::new(s.ground, 0.0, 0.0);
    let a = Vector3::new(pose.pin_a[0], pose.pin_a[1], 0.0);
    let b = Vector3::new(pose.pin_b[0], pose.pin_b[1], 0.0);

    let avg = (s.ground + s.crank + s.coupler + s.rocker) / 4.0;
    let hw = (avg * 0.035).max(1e-3);
    let ht = hw;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    // Links.
    push_bar(&mut nodes, &mut tris, o2, o4, hw * 0.7, ht * 0.7); // ground frame
    push_bar(&mut nodes, &mut tris, o2, a, hw, ht); // crank
    push_bar(&mut nodes, &mut tris, a, b, hw, ht); // coupler
    push_bar(&mut nodes, &mut tris, o4, b, hw, ht); // rocker
                                                    // Pivot / pin markers.
    let m = Vector3::new(hw * 1.6, hw * 1.6, ht * 1.6);
    for p in [o2, o4, a, b] {
        push_box(&mut nodes, &mut tris, p, m);
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-fourbar");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D linkage solid and load it into the central viewport.
fn load_linkage_3d(app: &mut ValenxApp) {
    let Some(mesh) = linkage_solid_mesh(&app.fourbar) else {
        app.fourbar.error =
            Some("linkage is invalid or cannot assemble at this crank angle — no 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<linkage>/valenx-kinematics"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical four-bar-linkage workbench as a 3-D solid
/// plus its `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn fourbar_product() -> crate::WorkspaceProduct {
    let s = FourBarWorkbenchState::default();
    let mesh = linkage_solid_mesh(&s).expect("canonical four-bar ⇒ linkage solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<fourbar>/valenx-linkage");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical four-bar ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Four-bar linkage (kinematics)".into(),
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
        let s = FourBarWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_grashof_and_pose() {
        let mut s = FourBarWorkbenchState::default();
        run_fourbar(&mut s);
        assert!(
            s.error.is_none(),
            "default crank-rocker solves: {:?}",
            s.error
        );
        assert!(s.result.contains("Grashof"));
        assert!(s.result.contains("transmission"));
        assert!(s.result.contains("rocker"));
    }

    #[test]
    fn analyze_rejects_nonpositive_length() {
        let mut s = FourBarWorkbenchState {
            ground: 0.0,
            ..Default::default()
        };
        run_fourbar(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn linkage_mesh_for_default_is_nonempty_and_in_range() {
        let s = FourBarWorkbenchState::default();
        let mesh = linkage_solid_mesh(&s).expect("default linkage yields a solid");
        assert!(mesh.nodes.len() > 8, "expected four links + pin markers");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn linkage_mesh_none_for_invalid() {
        let s = FourBarWorkbenchState {
            crank: 0.0,
            ..Default::default()
        };
        assert!(linkage_solid_mesh(&s).is_none());
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
            draw_fourbar_workbench(app, ctx);
        });
    }

    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_fourbar_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_fourbar_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_fourbar_workbench = true;
        run_fourbar(&mut app.fourbar);
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        let mut app = ValenxApp::default();
        app.show_fourbar_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);
        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 5,
            "expected numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by its caption (AI-drivable name)"
        );
        for caption in ["ground r1", "crank r2", "crank angle θ2 (°)"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }
}
