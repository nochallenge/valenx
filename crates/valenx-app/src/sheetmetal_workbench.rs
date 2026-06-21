//! The right-side **Sheet Metal Workbench** panel — bend allowance and
//! bend deduction over `valenx-sheet-metal`.
//!
//! Mirrors the springs / gears / geomatics / piping / collision
//! workbenches: a resizable [`egui::SidePanel`] gated on
//! `crate::ValenxApp::show_sheetmetal_workbench`, toggled from the View
//! menu. The form takes a bend's thickness, inside radius, angle, and
//! k-factor; the "Compute" button reports the neutral-axis radius, the
//! bend allowance (the neutral-axis arc length), and the bend deduction
//! (the flat-blank correction), as a monospace readout.

use eframe::egui;
use nalgebra::Vector3;

use valenx_sheet_metal::{Bend, Sheet};
use valenx_viz::{project_point, OrbitCamera, ViewDirection};

use crate::ValenxApp;

/// Persistent form + result state for the Sheet Metal Workbench.
pub struct SheetmetalWorkbenchState {
    /// Sheet thickness (mm).
    thickness_mm: f64,
    /// Inside (concave-side) bend radius (mm).
    inside_radius_mm: f64,
    /// Bend angle (degrees) — the swept angle of the bend.
    bend_angle_deg: f64,
    /// K-factor — the neutral-axis fraction (0.33–0.5 typical).
    k_factor: f64,
    /// Formatted bend readout (empty until the first compute).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
}

impl Default for SheetmetalWorkbenchState {
    fn default() -> Self {
        // A 90° bend in 1 mm stock, 1 mm inside radius, k = 0.44 — the
        // canonical worked example (bend allowance 1.44·π/2 ≈ 2.2619 mm).
        Self {
            thickness_mm: 1.0,
            inside_radius_mm: 1.0,
            bend_angle_deg: 90.0,
            k_factor: 0.44,
            result: String::new(),
            error: None,
        }
    }
}

/// Draw the Sheet Metal Workbench right-side panel. A no-op when the
/// `show_sheetmetal_workbench` toggle is off.
pub fn draw_sheetmetal_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_sheetmetal_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_sheetmetal_workbench",
        "Sheet Metal",
        |app, ui| {
            ui.label(
                egui::RichText::new("native bend allowance / deduction · valenx-sheet-metal")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.sheetmetal;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Bend parameters").strong());
                    ui.horizontal(|ui| {
                        ui.label("thickness (mm)");
                        ui.add(egui::DragValue::new(&mut s.thickness_mm).speed(0.05));
                    });
                    ui.horizontal(|ui| {
                        ui.label("inside radius (mm)");
                        ui.add(egui::DragValue::new(&mut s.inside_radius_mm).speed(0.05));
                    });
                    ui.horizontal(|ui| {
                        ui.label("bend angle (\u{00B0})");
                        ui.add(egui::DragValue::new(&mut s.bend_angle_deg).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("k-factor");
                        ui.add(egui::DragValue::new(&mut s.k_factor).speed(0.01));
                    });
                    ui.label(
                        egui::RichText::new("k-factor = neutral-axis fraction (0.33–0.5 typical)")
                            .weak()
                            .small(),
                    );

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("\u{25B6} Compute").strong())
                        .clicked()
                    {
                        run_sheetmetal(s);
                    }

                    if let Some(e) = &s.error {
                        ui.add_space(4.0);
                        ui.colored_label(egui::Color32::from_rgb(220, 90, 90), e);
                    }

                    if !s.result.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new("Bend").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }

                    // Live bent-sheet side profile (face-on).
                    if let Some(pts) = preview_bend(s) {
                        ui.add_space(6.0);
                        ui.label(egui::RichText::new("Bend profile preview").strong());
                        draw_bend_preview(ui, &pts);
                    }
                });
        },
    );
    if close {
        app.show_sheetmetal_workbench = false;
    }
}

/// The closed cross-section outline of a sheet of `thickness` bent through
/// `bend_angle_deg` at the given `inside_radius`, with straight flanges of
/// length `flange`. The bend centre is the origin: the inner surface is the
/// arc at `inside_radius`, the outer surface the concentric arc at
/// `inside_radius + thickness`, joined by straight tangent flanges and end
/// caps. Returned closed (last == first) in the XY plane (z = 0).
fn bend_profile(
    thickness: f64,
    inside_radius: f64,
    bend_angle_deg: f64,
    flange: f64,
) -> Vec<Vector3<f64>> {
    let theta = bend_angle_deg.to_radians();
    let half = theta / 2.0;
    let r_i = inside_radius;
    let r_o = inside_radius + thickness;
    const ARC_STEPS: usize = 24;

    // Radial unit vector at parameter ψ, measured from the +Y bisector.
    let radial = |psi: f64| Vector3::new(psi.sin(), psi.cos(), 0.0);
    // Outward flange directions at the ±half tangent points (⟂ to the radius).
    let dir_left = Vector3::new(-half.cos(), -half.sin(), 0.0);
    let dir_right = Vector3::new(half.cos(), -half.sin(), 0.0);

    let inner_l = radial(-half) * r_i;
    let inner_r = radial(half) * r_i;
    let outer_l = radial(-half) * r_o;
    let outer_r = radial(half) * r_o;

    let mut pts: Vec<Vector3<f64>> = Vec::new();
    // Left flange inner edge → inner arc (−half → +half) → right flange inner edge.
    pts.push(inner_l + dir_left * flange);
    pts.push(inner_l);
    for i in 1..ARC_STEPS {
        let psi = -half + (i as f64 / ARC_STEPS as f64) * theta;
        pts.push(radial(psi) * r_i);
    }
    pts.push(inner_r);
    pts.push(inner_r + dir_right * flange);
    // End cap across the thickness at the right flange tip.
    pts.push(outer_r + dir_right * flange);
    // Right flange outer edge → outer arc (+half → −half) → left flange outer edge.
    pts.push(outer_r);
    for i in 1..ARC_STEPS {
        let psi = half - (i as f64 / ARC_STEPS as f64) * theta;
        pts.push(radial(psi) * r_o);
    }
    pts.push(outer_l);
    pts.push(outer_l + dir_left * flange);
    // Close: end cap across the thickness at the left flange tip.
    let first = pts[0];
    pts.push(first);
    pts
}

/// Build the bent-sheet side profile for the live preview, best-effort `None`
/// when the form dimensions are invalid (the same preconditions as
/// [`run_sheetmetal`], minus the k-factor which only affects the scalars).
fn preview_bend(s: &SheetmetalWorkbenchState) -> Option<Vec<Vector3<f64>>> {
    if !(s.thickness_mm.is_finite() && s.thickness_mm > 0.0) {
        return None;
    }
    if !(s.inside_radius_mm.is_finite() && s.inside_radius_mm >= 0.0) {
        return None;
    }
    if !(s.bend_angle_deg.is_finite() && s.bend_angle_deg > 0.0 && s.bend_angle_deg <= 180.0) {
        return None;
    }
    // A flange length proportional to the bend so the profile is well framed.
    let flange = 3.0 * (s.inside_radius_mm + s.thickness_mm);
    Some(bend_profile(
        s.thickness_mm,
        s.inside_radius_mm,
        s.bend_angle_deg,
        flange,
    ))
}

/// Draw a closed XY polyline as a face-on (front-view, eye on +Z looking
/// along −Z) wireframe in a fixed-height canvas, the camera auto-framed to
/// the profile. A segment is painted only when both endpoints project in
/// front of the camera, so the render path never panics.
fn draw_bend_preview(ui: &mut egui::Ui, pts: &[Vector3<f64>]) {
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
    cam.set_view(ViewDirection::Front);
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

/// Build a [`Bend`] from the form, compute the bend allowance / deduction,
/// and format the readout. Extracted from the draw closure so it is
/// unit-testable. The bend-line endpoints do not affect either scalar, so
/// a unit line is used.
fn run_sheetmetal(s: &mut SheetmetalWorkbenchState) {
    s.error = None;

    if !(s.thickness_mm.is_finite() && s.thickness_mm > 0.0) {
        s.error = Some("thickness must be positive".into());
        return;
    }
    if !(s.inside_radius_mm.is_finite() && s.inside_radius_mm >= 0.0) {
        s.error = Some("inside radius must be \u{2265} 0".into());
        return;
    }
    if !(s.bend_angle_deg.is_finite() && s.bend_angle_deg > 0.0 && s.bend_angle_deg <= 180.0) {
        s.error = Some("bend angle must be in (0, 180]".into());
        return;
    }
    if !(s.k_factor.is_finite() && (0.0..=1.0).contains(&s.k_factor)) {
        s.error = Some("k-factor must be in [0, 1]".into());
        return;
    }

    let bend = Bend::new(
        [0.0, 0.0],
        [1.0, 0.0],
        s.bend_angle_deg.to_radians(),
        s.inside_radius_mm,
    );
    let r_neutral = s.inside_radius_mm + s.k_factor * s.thickness_mm;
    let ba = bend.bend_allowance(s.thickness_mm, s.k_factor);
    let bd = bend.bend_deduction(s.thickness_mm, s.k_factor);

    s.result = format!(
        "thickness      : {:.3} mm\n\
         inside radius  : {:.3} mm\n\
         bend angle     : {:.1}\u{00B0}\n\
         k-factor       : {:.3}\n\n\
         neutral radius : {:.4} mm   (r_i + k\u{00B7}t)\n\
         bend allowance : {:.4} mm   (r_n\u{00B7}\u{03B8}; neutral-axis arc)\n\
         bend deduction : {:.4} mm   (2\u{00B7}OSSB \u{2212} BA; flat-blank)",
        s.thickness_mm, s.inside_radius_mm, s.bend_angle_deg, s.k_factor, r_neutral, ba, bd,
    );
}

/// The agent-bridge product for the sheet-metal workbench
/// (`show_3d{kind="sheetmetal"}`).
///
/// Folds a **canonical demo part** — a 100 × 50 mm rectangular blank in 1 mm
/// stock with a single 90° bend down the middle — into a 3-D solid via
/// [`valenx_sheet_metal::Sheet::to_solid`], then extracts its `Tri3`
/// [`valenx_mesh::Mesh`] through the public `valenx_cad::solid_to_mesh` (a
/// mesh-backed [`valenx_cad::Solid`] short-circuits to its cached mesh — no
/// re-tessellation). The bend angle / inside radius mirror the workbench's
/// canonical worked example. Pure and app-state-free. The readout reports the
/// blank size, bend, and triangle count.
///
/// (The workbench panel itself is a 2-D bend-allowance calculator; the folded
/// solid is the natural 3-D companion the backing crate already produces.)
pub(crate) fn sheetmetal_product() -> crate::WorkspaceProduct {
    // Canonical demo blank (mm): 100 × 50 × 1 mm with one 90° bend at x = 50.
    const W: f64 = 100.0;
    const H: f64 = 50.0;
    const T: f64 = 1.0;
    const R: f64 = 1.0;
    let built = (|| -> Result<valenx_mesh::Mesh, String> {
        let sheet = Sheet::rectangle(W, H, T)
            .map_err(|e| e.to_string())?
            .add_bend(Bend::new(
                [W / 2.0, 0.0],
                [W / 2.0, H],
                std::f64::consts::FRAC_PI_2,
                R,
            ));
        let solid = sheet.to_solid().map_err(|e| e.to_string())?;
        valenx_cad::solid_to_mesh(&solid, valenx_cad::DEFAULT_TESS_TOLERANCE)
            .map_err(|e| e.to_string())
    })();
    let (mesh, lines) = match built {
        Ok(mesh) => {
            let tris = mesh.total_elements();
            let lines = vec![
                format!("sheet-metal blank {W:.0} \u{00D7} {H:.0} \u{00D7} {T:.0} mm"),
                "single 90\u{00B0} bend · 1 mm inside radius".to_string(),
                format!("folded solid: {tris} triangles"),
            ];
            (mesh, lines)
        }
        Err(e) => {
            // Theoretically unreachable for the canonical blank; degrade to a
            // tiny placeholder triangle + a note rather than panicking.
            let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
            block.connectivity = vec![0, 1, 2];
            let mut placeholder = valenx_mesh::Mesh::new("valenx-sheet-metal");
            placeholder.nodes = vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(0.0, 1.0, 0.0),
            ];
            placeholder.element_blocks.push(block);
            placeholder.recompute_stats();
            (
                placeholder,
                vec![
                    "folded sheet-metal part".to_string(),
                    format!("fold unavailable — showing placeholder ({e})"),
                ],
            )
        }
    };
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<sheetmetal>/folded");
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Sheet-metal part (folded)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
        kind2d: None,
        last_export: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = SheetmetalWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn compute_default_90deg_bend() {
        let mut s = SheetmetalWorkbenchState::default();
        run_sheetmetal(&mut s);
        assert!(s.error.is_none());
        assert!(!s.result.is_empty());
        assert!(s.result.contains("neutral radius"));
        assert!(s.result.contains("bend allowance"));
        assert!(s.result.contains("bend deduction"));
        // 90° bend, r_i = 1, t = 1, k = 0.44: r_n = 1.44, BA = 1.44·π/2 ≈
        // 2.2619 mm. Recompute via the backend to confirm the default.
        let bend = Bend::new([0.0, 0.0], [1.0, 0.0], 90_f64.to_radians(), 1.0);
        assert!((bend.bend_allowance(1.0, 0.44) - 2.261_95).abs() < 1e-3);
        assert!(s.result.contains("2.2619"));
    }

    #[test]
    fn compute_rejects_bad_inputs() {
        for bad in [
            SheetmetalWorkbenchState {
                thickness_mm: 0.0,
                ..Default::default()
            },
            SheetmetalWorkbenchState {
                bend_angle_deg: 0.0,
                ..Default::default()
            },
            SheetmetalWorkbenchState {
                bend_angle_deg: 200.0,
                ..Default::default()
            },
            SheetmetalWorkbenchState {
                k_factor: 1.5,
                ..Default::default()
            },
        ] {
            let mut s = bad;
            run_sheetmetal(&mut s);
            assert!(s.error.is_some());
            assert!(s.result.is_empty());
        }
    }

    #[test]
    fn bend_profile_is_closed_with_correct_radii_and_angle() {
        let (t, ri, ang) = (1.0, 1.0, 90.0);
        let pts = bend_profile(t, ri, ang, 6.0);
        assert!(pts.len() > 8);
        assert_eq!(pts.first(), pts.last()); // closed outline
        assert!(pts.iter().all(|p| p.z == 0.0)); // planar (z = 0)
        let dist = |p: &Vector3<f64>| (p.x * p.x + p.y * p.y).sqrt();
        // The inner surface hugs the inside radius: the closest points to the
        // bend centre (origin) sit at r_i.
        let dmin = pts.iter().map(dist).fold(f64::INFINITY, f64::min);
        assert!((dmin - ri).abs() < 1e-9);
        // The outer surface sits at r_i + t — at least one point is there.
        let r_o = ri + t;
        assert!(pts.iter().any(|p| (dist(p) - r_o).abs() < 1e-9));
        // The inner arc subtends the bend angle about the centre.
        let arc_psi: Vec<f64> = pts
            .iter()
            .filter(|&p| (dist(p) - ri).abs() < 1e-9)
            .map(|p| p.x.atan2(p.y))
            .collect();
        let amin = arc_psi.iter().cloned().fold(f64::INFINITY, f64::min);
        let amax = arc_psi.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!((amax - amin - ang.to_radians()).abs() < 1e-9);
    }

    #[test]
    fn preview_bend_some_for_default_none_for_invalid() {
        let s = SheetmetalWorkbenchState::default();
        let pts = preview_bend(&s).expect("default 90° bend previews");
        assert_eq!(pts.first(), pts.last());
        assert!(pts.len() > 8);
        assert!(pts.iter().all(|p| p.z == 0.0));
        for bad in [
            SheetmetalWorkbenchState {
                thickness_mm: 0.0,
                ..Default::default()
            },
            SheetmetalWorkbenchState {
                bend_angle_deg: 0.0,
                ..Default::default()
            },
            SheetmetalWorkbenchState {
                bend_angle_deg: 200.0,
                ..Default::default()
            },
        ] {
            assert!(preview_bend(&bad).is_none());
        }
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
            draw_sheetmetal_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_sheetmetal_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_sheetmetal_workbench = true;
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_a_result_and_error_without_panic() {
        let mut app = ValenxApp::default();
        app.show_sheetmetal_workbench = true;
        run_sheetmetal(&mut app.sheetmetal);
        app.sheetmetal.error = Some("invalid bend parameters".to_string());
        draw_workbench(&mut app);
    }
}
