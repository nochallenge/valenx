//! The right-side **Springs Workbench** panel — native helical-spring
//! design over `valenx-springs`.
//!
//! Mirrors the CFD / FEM workbenches: a resizable [`egui::SidePanel`]
//! gated on `crate::ValenxApp::show_springs_workbench`, toggled from the
//! View menu. The form drives a [`valenx_springs::SpringSpec`]; the
//! "Analyze" button computes the design scalars — spring index, axial
//! stiffness, Wahl stress-correction factor, the outer/inner diameters,
//! the coil pitch, and the developed wire length — and renders them as a
//! monospace readout.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;

use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_springs::{
    compression_centerline, extension_centerline, spring_index, stiffness_n_per_mm,
    torsion_centerline, wahl_factor, SpringKind, SpringSpec,
};
use valenx_viz::{project_point, OrbitCamera, ViewDirection};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Springs Workbench.
pub struct SpringsWorkbenchState {
    /// Spring kind (compression / extension / torsion).
    kind: SpringKind,
    /// Wire diameter `d` (mm).
    wire_diameter_mm: f64,
    /// Mean coil diameter `D` (mm), measured to the wire centreline.
    mean_coil_diameter_mm: f64,
    /// Free (unloaded) length (mm).
    free_length_mm: f64,
    /// Number of active coils.
    n_active_coils: f64,
    /// Shear modulus `G` (MPa) — steel ≈ 79 300 MPa.
    shear_modulus_mpa: f64,
    /// Formatted design readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the spring's 3-D solid and load it into
    /// the central viewport (set by the "Show 3-D spring" button; serviced
    /// after the panel draws, like the Rocket workbench's 3-D model).
    show_3d_request: bool,
}

impl Default for SpringsWorkbenchState {
    fn default() -> Self {
        // Mirror valenx-springs' default compression spring.
        Self {
            kind: SpringKind::Compression,
            wire_diameter_mm: 1.0,
            mean_coil_diameter_mm: 10.0,
            free_length_mm: 30.0,
            n_active_coils: 8.0,
            shear_modulus_mpa: 79_300.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Springs Workbench right-side panel. A no-op when the
/// `show_springs_workbench` toggle is off.
pub fn draw_springs_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_springs_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_springs_workbench",
        "Springs",
        |app, ui| {
            ui.label(
                egui::RichText::new("native helical-spring design · valenx-springs")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.springs;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Kind").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.kind, SpringKind::Compression, "Compression");
                        ui.radio_value(&mut s.kind, SpringKind::Extension, "Extension");
                        ui.radio_value(&mut s.kind, SpringKind::Torsion, "Torsion");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Geometry (mm)").strong());
                    ui.horizontal(|ui| {
                        ui.label("wire d");
                        ui.add(egui::DragValue::new(&mut s.wire_diameter_mm).speed(0.05));
                    });
                    ui.horizontal(|ui| {
                        ui.label("mean coil D");
                        ui.add(egui::DragValue::new(&mut s.mean_coil_diameter_mm).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("free length");
                        ui.add(egui::DragValue::new(&mut s.free_length_mm).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("active coils n");
                        ui.add(egui::DragValue::new(&mut s.n_active_coils).speed(0.1));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Material").strong());
                    ui.horizontal(|ui| {
                        ui.label("shear modulus G (MPa)");
                        ui.add(egui::DragValue::new(&mut s.shear_modulus_mpa).speed(100.0));
                    });

                    // Live hint: the spring index C = D/d (4–12 is the
                    // manufacturability sweet spot).
                    if s.wire_diameter_mm > 0.0 {
                        let c = s.mean_coil_diameter_mm / s.wire_diameter_mm;
                        ui.label(
                            egui::RichText::new(format!("spring index C ≈ {c:.2}  (4–12 typical)"))
                                .weak()
                                .small(),
                        );
                    }

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_springs(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D spring").strong())
                        .on_hover_text(
                            "Build the wire as a solid 3-D tube and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Design").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }

                    // Live 3-D wireframe of the spring's centreline.
                    if let Some(pts) = preview_centerline(s) {
                        ui.add_space(6.0);
                        ui.label(egui::RichText::new("Centerline preview (isometric)").strong());
                        draw_centerline_preview(ui, &pts);
                    }
                });
        },
    );
    if close {
        app.show_springs_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.springs` borrow taken in
    // the closure above is released here): build the spring's 3-D solid and
    // load it into the central viewport, like the Rocket workbench's model.
    if app.springs.show_3d_request {
        app.springs.show_3d_request = false;
        load_spring_3d(app);
    }
}

/// Build a [`SpringSpec`] from the current form and return the selected
/// spring's 3-D centreline polyline — best-effort, `None` for an invalid
/// spec. Drives the live isometric preview.
fn preview_centerline(s: &SpringsWorkbenchState) -> Option<Vec<Vector3<f64>>> {
    if !(s.wire_diameter_mm.is_finite()
        && s.wire_diameter_mm > 0.0
        && s.mean_coil_diameter_mm.is_finite()
        && s.mean_coil_diameter_mm > s.wire_diameter_mm
        && s.free_length_mm.is_finite()
        && s.free_length_mm > 0.0
        && s.n_active_coils.is_finite()
        && s.n_active_coils > 0.0)
    {
        return None;
    }
    let spec = SpringSpec {
        kind: s.kind,
        wire_diameter_mm: s.wire_diameter_mm,
        mean_coil_diameter_mm: s.mean_coil_diameter_mm,
        free_length_mm: s.free_length_mm,
        n_active_coils: s.n_active_coils,
        shear_modulus_mpa: s.shear_modulus_mpa,
        ..SpringSpec::default_compression()
    };
    let pts = match s.kind {
        SpringKind::Compression => compression_centerline(&spec),
        SpringKind::Extension => extension_centerline(&spec),
        SpringKind::Torsion => torsion_centerline(&spec),
    };
    pts.ok().filter(|v| v.len() >= 2)
}

/// Sweep a circular wire cross-section (radius = wire diameter / 2) along
/// the spring's centreline to build a watertight 3-D tube and return it as
/// a triangle [`Mesh`]. `None` for an invalid spec (mirrors
/// `preview_centerline`). At each centreline point the cross-section is laid
/// out in the tangent frame `n1 = tangent x reference`, `n2 = tangent x n1`,
/// with `reference` switched from +Z to +X near the poles so the basis never
/// degenerates. Consecutive rings are joined by quad side-walls (two
/// triangles each) and the two open ends are closed with triangle-fan caps.
fn spring_solid_mesh(s: &SpringsWorkbenchState) -> Option<Mesh> {
    let path = preview_centerline(s)?;
    let r = s.wire_diameter_mm / 2.0;
    let seg = 12usize;

    let up = Vector3::new(0.0, 0.0, 1.0);
    let alt = Vector3::new(1.0, 0.0, 0.0);

    let mut nodes: Vec<Vector3<f64>> = Vec::with_capacity(path.len() * seg + 2);
    let mut tris: Vec<usize> = Vec::new();

    // One ring of `seg` vertices around each centreline point, in the plane
    // perpendicular to the local tangent.
    for i in 0..path.len() {
        let raw = if i + 1 < path.len() {
            path[i + 1] - path[i]
        } else {
            path[i] - path[i - 1]
        };
        let tan = if raw.norm() > 1e-12 {
            raw.normalize()
        } else {
            up
        };
        // Reference axis not (nearly) parallel to the tangent, so the
        // cross-section basis stays well-conditioned.
        let reference = if tan.dot(&up).abs() < 0.9 { up } else { alt };
        let n1 = tan.cross(&reference).normalize();
        let n2 = tan.cross(&n1).normalize();
        for j in 0..seg {
            let theta = j as f64 / seg as f64 * TAU;
            nodes.push(path[i] + (n1 * theta.cos() + n2 * theta.sin()) * r);
        }
    }

    // Side walls: connect ring i to ring i+1 with two triangles per segment.
    for i in 0..path.len() - 1 {
        let a = i * seg;
        let b = (i + 1) * seg;
        for j in 0..seg {
            let jn = (j + 1) % seg;
            tris.extend_from_slice(&[a + j, b + j, b + jn, a + j, b + jn, a + jn]);
        }
    }

    // End caps: a triangle fan from each end ring's centre point.
    let first_center = nodes.len();
    nodes.push(path[0]);
    for j in 0..seg {
        let jn = (j + 1) % seg;
        tris.extend_from_slice(&[first_center, jn, j]);
    }
    let last_ring = (path.len() - 1) * seg;
    let last_center = nodes.len();
    nodes.push(path[path.len() - 1]);
    for j in 0..seg {
        let jn = (j + 1) % seg;
        tris.extend_from_slice(&[last_center, last_ring + j, last_ring + jn]);
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-spring");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the spring's 3-D wire solid and load it into the central viewport
/// (replacing any current STL / mesh) so it can be orbited — mirrors the
/// Rocket workbench's `load_lv1_rocket_3d`. On an invalid spec, surfaces a
/// form error instead of loading anything.
fn load_spring_3d(app: &mut ValenxApp) {
    let Some(mesh) = spring_solid_mesh(&app.springs) else {
        app.springs.error =
            Some("spring parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<spring>/valenx-springs"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Draw a centreline polyline as an isometric wireframe in a fixed-height
/// canvas, with the camera auto-framed to the geometry's bounds.
fn draw_centerline_preview(ui: &mut egui::Ui, pts: &[Vector3<f64>]) {
    let (response, painter) = ui.allocate_painter(
        egui::vec2(ui.available_width(), 200.0),
        egui::Sense::hover(),
    );
    let rect = response.rect;

    let mut min = [f32::INFINITY; 3];
    let mut max = [f32::NEG_INFINITY; 3];
    for p in pts {
        for k in 0..3 {
            let v = p[k] as f32;
            min[k] = min[k].min(v);
            max[k] = max[k].max(v);
        }
    }

    let mut cam = OrbitCamera::default();
    cam.set_view(ViewDirection::Iso);
    cam.frame_bounds(min, max);

    let (w, h) = (rect.width(), rect.height());
    let stroke = egui::Stroke::new(1.5, egui::Color32::from_rgb(120, 200, 255));
    for pair in pts.windows(2) {
        let a = project_point(
            &cam,
            w,
            h,
            [pair[0].x as f32, pair[0].y as f32, pair[0].z as f32],
        );
        let b = project_point(
            &cam,
            w,
            h,
            [pair[1].x as f32, pair[1].y as f32, pair[1].z as f32],
        );
        if let (Some(a), Some(b)) = (a, b) {
            painter.line_segment(
                [
                    egui::pos2(rect.min.x + a.x, rect.min.y + a.y),
                    egui::pos2(rect.min.x + b.x, rect.min.y + b.y),
                ],
                stroke,
            );
        }
    }
}

/// Build a [`SpringSpec`] from the form, validate it, and format the
/// design-scalar readout. Extracted from the draw closure so it is
/// unit-testable.
fn run_springs(s: &mut SpringsWorkbenchState) {
    s.error = None;

    // Validate against the backend's preconditions (mean coil > wire,
    // free length > 0, active coils > 0) so the scalars stay finite + sane.
    if !(s.wire_diameter_mm.is_finite() && s.wire_diameter_mm > 0.0) {
        s.error = Some("wire diameter must be positive".into());
        return;
    }
    if !(s.mean_coil_diameter_mm.is_finite() && s.mean_coil_diameter_mm > s.wire_diameter_mm) {
        s.error = Some("mean coil diameter must exceed the wire diameter".into());
        return;
    }
    if !(s.free_length_mm.is_finite() && s.free_length_mm > 0.0) {
        s.error = Some("free length must be positive".into());
        return;
    }
    if !(s.n_active_coils.is_finite() && s.n_active_coils > 0.0) {
        s.error = Some("active-coil count must be positive".into());
        return;
    }
    if !(s.shear_modulus_mpa.is_finite() && s.shear_modulus_mpa > 0.0) {
        s.error = Some("shear modulus must be positive".into());
        return;
    }

    // Struct-update fills `end_treatment` (and any future field) from the
    // canonical default so we only surface the design-relevant inputs.
    let spec = SpringSpec {
        kind: s.kind,
        wire_diameter_mm: s.wire_diameter_mm,
        mean_coil_diameter_mm: s.mean_coil_diameter_mm,
        free_length_mm: s.free_length_mm,
        n_active_coils: s.n_active_coils,
        shear_modulus_mpa: s.shear_modulus_mpa,
        ..SpringSpec::default_compression()
    };

    let c = spring_index(&spec);
    let k = stiffness_n_per_mm(&spec);
    let kw = wahl_factor(c);

    s.result = format!(
        "kind          : {}\n\
         wire d        : {:.3} mm\n\
         mean coil D   : {:.3} mm\n\
         free length   : {:.3} mm\n\
         active coils  : {:.2}\n\
         shear mod G   : {:.0} MPa\n\n\
         spring index C: {:.3}  (D/d)\n\
         stiffness k   : {:.4} N/mm  (G·d⁴/(8·D³·n))\n\
         Wahl factor   : {:.4}  (stress correction)\n\
         outer diameter: {:.3} mm\n\
         inner diameter: {:.3} mm\n\
         pitch         : {:.3} mm/coil\n\
         wire length   : {:.2} mm  (developed helix)",
        s.kind.label(),
        s.wire_diameter_mm,
        s.mean_coil_diameter_mm,
        s.free_length_mm,
        s.n_active_coils,
        s.shear_modulus_mpa,
        c,
        k,
        kw,
        spec.outer_diameter_mm(),
        spec.inner_diameter_mm(),
        spec.pitch_mm(),
        spec.helix_length_mm(),
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = SpringsWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_compression_spring() {
        let mut s = SpringsWorkbenchState::default();
        run_springs(&mut s);
        assert!(s.error.is_none());
        assert!(!s.result.is_empty());
        // The readout names the core design scalars.
        assert!(s.result.contains("spring index C"));
        assert!(s.result.contains("stiffness k"));
        assert!(s.result.contains("Wahl factor"));
        assert!(s.result.contains("wire length"));
        // Default compression: C = D/d = 10/1 = 10.
        assert!(s.result.contains("10.000"));
    }

    #[test]
    fn analyze_rejects_coil_not_exceeding_wire() {
        let mut s = SpringsWorkbenchState {
            wire_diameter_mm: 10.0,
            mean_coil_diameter_mm: 8.0,
            ..Default::default()
        };
        run_springs(&mut s);
        assert!(s.error.is_some());
        assert!(s.result.is_empty());
    }

    #[test]
    fn analyze_rejects_nonpositive_inputs() {
        for bad in [
            SpringsWorkbenchState {
                wire_diameter_mm: 0.0,
                ..Default::default()
            },
            SpringsWorkbenchState {
                free_length_mm: 0.0,
                ..Default::default()
            },
            SpringsWorkbenchState {
                n_active_coils: 0.0,
                ..Default::default()
            },
            SpringsWorkbenchState {
                shear_modulus_mpa: 0.0,
                ..Default::default()
            },
        ] {
            let mut s = bad;
            run_springs(&mut s);
            assert!(s.error.is_some());
        }
    }

    #[test]
    fn preview_centerline_for_default_is_a_polyline() {
        let s = SpringsWorkbenchState::default();
        let pts = preview_centerline(&s).expect("default spec yields a centerline");
        assert!(pts.len() >= 2);
        // The compression helix has a non-zero axial (z) span.
        let zmax = pts.iter().map(|p| p.z).fold(f64::NEG_INFINITY, f64::max);
        let zmin = pts.iter().map(|p| p.z).fold(f64::INFINITY, f64::min);
        assert!(zmax - zmin > 0.0);
    }

    #[test]
    fn preview_centerline_none_for_invalid_spec() {
        // mean coil ≤ wire → the spec is invalid → no preview.
        let s = SpringsWorkbenchState {
            mean_coil_diameter_mm: 0.5,
            ..Default::default()
        };
        assert!(preview_centerline(&s).is_none());
    }

    #[test]
    fn spring_solid_mesh_for_default_is_a_nonempty_solid() {
        let s = SpringsWorkbenchState::default();
        let mesh = spring_solid_mesh(&s).expect("default spec yields a 3-D solid");
        // A swept tube: one ring of vertices per centreline point plus two
        // cap centres, and a populated triangle block.
        assert!(mesh.nodes.len() > 12, "expected a multi-ring tube");
        assert!(
            mesh.element_blocks
                .iter()
                .any(|b| !b.connectivity.is_empty()),
            "expected at least one triangle"
        );
        // The connectivity is whole triangles (3 indices each), all in range.
        let n = mesh.nodes.len() as u32;
        for b in &mesh.element_blocks {
            assert_eq!(b.connectivity.len() % 3, 0);
            assert!(b.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn spring_solid_mesh_none_for_invalid_spec() {
        let s = SpringsWorkbenchState {
            mean_coil_diameter_mm: 0.5,
            ..Default::default()
        };
        assert!(spring_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    /// Render the whole workbench panel once in a headless egui context.
    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_springs_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_springs_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_springs_workbench = true;
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_a_result_and_error_without_panic() {
        let mut app = ValenxApp::default();
        app.show_springs_workbench = true;
        run_springs(&mut app.springs);
        app.springs.error = Some("invalid spring parameters".to_string());
        draw_workbench(&mut app);
    }
}
