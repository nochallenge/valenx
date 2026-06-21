//! The right-side **Truss Workbench** panel — native planar pin-jointed
//! truss analysis over `valenx-truss`.
//!
//! Mirrors the Four-Bar / Wind Turbine workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_truss_workbench`,
//! toggled from the View menu. The form builds a parametric parallel-chord
//! **Warren truss** (a pin support at the left, a roller at the right, and a
//! load split across the interior bottom joints) and solves it with the
//! method of joints; "Analyze" reports the determinacy, the support
//! reactions and the peak tension / compression members, and "Show 3-D
//! truss" loads the bar frame into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_truss::{solve, Determinacy, Load, Member, Node, Support, Truss};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Truss Workbench.
pub struct TrussWorkbenchState {
    /// Number of panels (bays) along the span; at least two.
    panels: u32,
    /// Overall span `L` (length units).
    span: f64,
    /// Truss height `H` (the top chord rises this far above the bottom).
    height: f64,
    /// Total downward load, split evenly across the interior bottom joints.
    total_load: f64,
    /// Formatted readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D truss solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for TrussWorkbenchState {
    fn default() -> Self {
        // A 4-panel Warren truss, 12 long x 2.5 high, carrying 40 of load.
        Self {
            panels: 4,
            span: 12.0,
            height: 2.5,
            total_load: 40.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Truss Workbench right-side panel. A no-op when the
/// `show_truss_workbench` toggle is off.
pub fn draw_truss_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_truss_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_truss_workbench",
        "Truss",
        |app, ui| {
            ui.label(
                egui::RichText::new("native planar pin-jointed Warren truss · valenx-truss")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.truss;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Warren truss").strong());
                    ui.horizontal(|ui| {
                        ui.label("panels (bays)");
                        ui.add(egui::DragValue::new(&mut s.panels).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("span L");
                        ui.add(egui::DragValue::new(&mut s.span).speed(0.2));
                    });
                    ui.horizontal(|ui| {
                        ui.label("height H");
                        ui.add(egui::DragValue::new(&mut s.height).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("total load (down)");
                        ui.add(egui::DragValue::new(&mut s.total_load).speed(1.0));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_truss(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D truss").strong())
                        .on_hover_text(
                            "Build the truss members as a 3-D bar frame and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Statics").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_truss_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.truss` borrow is
    // released here): build the truss's 3-D solid and load it.
    if app.truss.show_3d_request {
        app.truss.show_3d_request = false;
        load_truss_3d(app);
    }
}

/// Build a parametric parallel-chord Warren truss: `n + 1` bottom-chord
/// joints (a pin at the left, a roller at the right, the total load split
/// across the interior joints) and `n` top-chord joints at the panel
/// mid-points, triangulated by the zig-zag diagonals. The result is
/// statically determinate (`m + r = 2N`).
fn build_truss(s: &TrussWorkbenchState) -> Result<Truss, String> {
    if s.panels < 2 {
        return Err("need at least 2 panels".into());
    }
    if !(s.span.is_finite() && s.span > 0.0) {
        return Err("span must be finite and > 0".into());
    }
    if !(s.height.is_finite() && s.height > 0.0) {
        return Err("height must be finite and > 0".into());
    }
    let n = s.panels as usize;
    let dx = s.span / n as f64;
    let per = s.total_load / (n - 1) as f64; // interior bottom joints

    let mut t = Truss::new();
    let mut bottom = Vec::with_capacity(n + 1);
    for i in 0..=n {
        let mut node = Node::new(i as f64 * dx, 0.0).map_err(|e| e.to_string())?;
        if i == 0 {
            node = node.with_support(Support::Pin);
        } else if i == n {
            node = node.with_support(Support::horizontal_roller());
        } else {
            node = node.with_load(Load::down(per).map_err(|e| e.to_string())?);
        }
        bottom.push(t.add_node(node));
    }
    let mut top = Vec::with_capacity(n);
    for i in 0..n {
        let x = (i as f64 + 0.5) * dx;
        top.push(t.add_node(Node::new(x, s.height).map_err(|e| e.to_string())?));
    }
    // Bottom chord.
    for i in 0..n {
        t.add_member(Member::new(bottom[i], bottom[i + 1]))
            .map_err(|e| e.to_string())?;
    }
    // Top chord.
    for i in 0..n.saturating_sub(1) {
        t.add_member(Member::new(top[i], top[i + 1]))
            .map_err(|e| e.to_string())?;
    }
    // Zig-zag diagonals: each top joint to the two bottom joints below it.
    for i in 0..n {
        t.add_member(Member::new(bottom[i], top[i]))
            .map_err(|e| e.to_string())?;
        t.add_member(Member::new(top[i], bottom[i + 1]))
            .map_err(|e| e.to_string())?;
    }
    Ok(t)
}

/// Validate the form, solve the truss and format the readout.
fn run_truss(s: &mut TrussWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build + solve the truss and format the readout. Extracted so it is
/// unit-testable.
fn compute(s: &TrussWorkbenchState) -> Result<String, String> {
    let t = build_truss(s)?;
    let det = match t.determinacy() {
        Determinacy::Determinate => "statically determinate",
        Determinacy::Mechanism => "mechanism (under-constrained)",
        Determinacy::Indeterminate => "statically indeterminate",
    };
    let sol = solve(&t).map_err(|e| e.to_string())?;
    let mut max_t = 0.0_f64;
    let mut max_t_m = 0usize;
    let mut max_c = 0.0_f64;
    let mut max_c_m = 0usize;
    for m in 0..t.members.len() {
        let f = sol.force(m); // tension positive, compression negative
        if f > max_t {
            max_t = f;
            max_t_m = m;
        }
        if f < max_c {
            max_c = f;
            max_c_m = m;
        }
    }
    Ok(format!(
        "panels     : {}\n\
         span L     : {:.2}\n\
         height H   : {:.2}\n\
         total load : {:.1} (down)\n\n\
         joints / members: {} / {}\n\
         {}\n\n\
         reactions Σ : Rx {:.2}, Ry {:.2}\n\n\
         max tension     : {:.2}  (member {})\n\
         max compression : {:.2}  (member {})\n\
         (sign: tension +, compression −)",
        s.panels,
        s.span,
        s.height,
        s.total_load,
        t.nodes.len(),
        t.members.len(),
        det,
        sol.total_reaction_fx(),
        sol.total_reaction_fy(),
        max_t,
        max_t_m,
        max_c,
        max_c_m,
    ))
}

/// Append an outward-facing box (centre `c`, half-extents `h`) to the
/// buffers (used for the joint markers).
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
/// plane) of half-width `half_w` (in-plane) and half-thickness `half_t`.
/// Faces are emitted double-sided. Used for the truss members.
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
    nodes.push(p - perp - zt);
    nodes.push(q - perp - zt);
    nodes.push(q + perp - zt);
    nodes.push(p + perp - zt);
    nodes.push(p - perp + zt);
    nodes.push(q - perp + zt);
    nodes.push(q + perp + zt);
    nodes.push(p + perp + zt);
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

/// Build the truss as a triangle [`Mesh`] — every member as a bar in the
/// z = 0 plane plus a marker box at each joint. `None` for an invalid truss.
fn truss_solid_mesh(s: &TrussWorkbenchState) -> Option<Mesh> {
    let t = build_truss(s).ok()?;
    let hw = (s.span * 0.008).max(1e-3);
    let ht = hw;
    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    for m in &t.members {
        let a = &t.nodes[m.a];
        let b = &t.nodes[m.b];
        push_bar(
            &mut nodes,
            &mut tris,
            Vector3::new(a.x, a.y, 0.0),
            Vector3::new(b.x, b.y, 0.0),
            hw,
            ht,
        );
    }
    let mk = Vector3::new(hw * 1.8, hw * 1.8, ht * 1.8);
    for nd in &t.nodes {
        push_box(&mut nodes, &mut tris, Vector3::new(nd.x, nd.y, 0.0), mk);
    }
    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-truss");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D truss solid and load it into the central viewport.
fn load_truss_3d(app: &mut ValenxApp) {
    let Some(mesh) = truss_solid_mesh(&app.truss) else {
        app.truss.error = Some("truss parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<truss>/valenx-truss"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"truss"}`** product: the parametric Warren
/// truss bar frame (members + joint markers) built from the canonical 4-panel
/// truss (12 long × 2.5 high, 40 of load), paired with the statics readout rows
/// (determinacy, reactions, peak tension / compression members), at a fixed 3/4
/// camera. Registered in [`crate::products_registry`]; the per-tool builder the
/// registry dispatches to. Pure — driven off [`TrussWorkbenchState::default`].
pub(crate) fn truss_product() -> crate::WorkspaceProduct {
    let s = TrussWorkbenchState::default();
    let mesh = truss_solid_mesh(&s).expect("canonical Warren truss ⇒ bar-frame solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<truss>/valenx-truss");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical Warren truss ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Warren truss (method of joints)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = TrussWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_forces_and_determinacy() {
        let mut s = TrussWorkbenchState::default();
        run_truss(&mut s);
        assert!(s.error.is_none(), "default truss solves: {:?}", s.error);
        assert!(s.result.contains("determinate"));
        assert!(s.result.contains("max tension"));
        assert!(s.result.contains("max compression"));
    }

    #[test]
    fn analyze_rejects_too_few_panels() {
        let mut s = TrussWorkbenchState {
            panels: 1,
            ..Default::default()
        };
        run_truss(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn truss_mesh_for_default_is_nonempty_and_in_range() {
        let s = TrussWorkbenchState::default();
        let mesh = truss_solid_mesh(&s).expect("default truss yields a solid");
        assert!(mesh.nodes.len() > 8, "expected members + joint markers");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn truss_mesh_none_for_invalid() {
        let s = TrussWorkbenchState {
            panels: 1,
            ..Default::default()
        };
        assert!(truss_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_truss_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_truss_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_truss_workbench = true;
        run_truss(&mut app.truss);
        draw_workbench(&mut app);
    }
}
