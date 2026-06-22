//! The right-side **Fasteners Workbench** panel — ISO 4017 hex-bolt
//! dimensions over `valenx-fasteners`.
//!
//! Mirrors the springs / gears / … / fields workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_fasteners_workbench`,
//! toggled from the View menu. The form picks a standard metric bolt size
//! from the ISO 4017 hex table; the "Compute" button reports the width
//! across flats, head height, pitch diameter, and tensile stress area, as
//! a monospace readout.

use eframe::egui;
use nalgebra::Vector3;

use valenx_fasteners::bolt::{iso4017_hex_table, BoltSpec};
use valenx_viz::{project_point, OrbitCamera, ViewDirection};

use crate::mesh_prims::MeshBuilder;
use crate::ValenxApp;

/// Zinc-plated steel grey for the bolt head + shank.
const STEEL: [f32; 3] = [0.62, 0.64, 0.68];
/// Angular segments around the cylindrical shank.
const SHANK_SEGMENTS: usize = 32;

/// Persistent form + result state for the Fasteners Workbench.
pub struct FastenersWorkbenchState {
    /// Selected ISO metric nominal designation (e.g. `"M6"`).
    nominal: String,
    /// Formatted dimension readout (empty until the first compute).
    result: String,
    /// Validation / lookup error, if any.
    error: Option<String>,
}

impl Default for FastenersWorkbenchState {
    fn default() -> Self {
        Self {
            nominal: "M6".to_string(),
            result: String::new(),
            error: None,
        }
    }
}

/// Draw the Fasteners Workbench right-side panel. A no-op when the
/// `show_fasteners_workbench` toggle is off.
pub fn draw_fasteners_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_fasteners_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_fasteners_workbench",
        "Fasteners",
        |app, ui| {
            ui.label(
                egui::RichText::new("ISO 4017 hex-bolt dimensions · valenx-fasteners")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.fasteners;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Bolt size (ISO 4017 hex)").strong());
                    let table = iso4017_hex_table();
                    egui::ComboBox::from_id_source("valenx_fasteners_nominal")
                        .selected_text(s.nominal.clone())
                        .show_ui(ui, |ui| {
                            for spec in &table {
                                ui.selectable_value(
                                    &mut s.nominal,
                                    spec.nominal.clone(),
                                    spec.nominal.as_str(),
                                );
                            }
                        });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("\u{25B6} Compute").strong())
                        .clicked()
                    {
                        run_fasteners(s);
                    }

                    if let Some(e) = &s.error {
                        ui.add_space(4.0);
                        ui.colored_label(egui::Color32::from_rgb(220, 90, 90), e);
                    }

                    if !s.result.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new("Dimensions").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }

                    // Live hex-head outline (across-flats, face-on).
                    if let Some(pts) = preview_hex(s) {
                        ui.add_space(6.0);
                        ui.label(egui::RichText::new("Hex-head preview").strong());
                        draw_hex_preview(ui, &pts);
                    }
                });
        },
    );
    if close {
        app.show_fasteners_workbench = false;
    }
}

/// A regular hexagon with the given width **across flats** (the distance
/// between opposite parallel flats), centred at the origin in the XY plane.
/// The circumradius is `af / √3`; six vertices at 60° steps are returned
/// closed (7 points, last == first) so a `windows(2)` walk draws the outline.
fn hexagon(across_flats: f64) -> Vec<Vector3<f64>> {
    let r = across_flats / 3.0_f64.sqrt();
    let mut pts: Vec<Vector3<f64>> = (0..6)
        .map(|i| {
            let th = (i as f64) * std::f64::consts::FRAC_PI_3; // 0°, 60°, … 300°
            Vector3::new(r * th.cos(), r * th.sin(), 0.0)
        })
        .collect();
    let first = pts[0];
    pts.push(first);
    pts
}

/// The six **counter-clockwise** `[x, y]` corners of a regular hexagon of the
/// given width **across flats**, centred at the origin — the closed profile
/// [`MeshBuilder::extrude`] sweeps into the hex-head prism. Same circumradius
/// `R = af / √3` as [`hexagon`], but returned open (6 distinct vertices, no
/// closing duplicate) and as 2-D pairs, which is what `extrude` wants.
fn hex_profile(across_flats: f64) -> Vec<[f64; 2]> {
    let r = across_flats / 3.0_f64.sqrt();
    (0..6)
        .map(|i| {
            let th = (i as f64) * std::f64::consts::FRAC_PI_3; // 0°, 60°, … 300°
            [r * th.cos(), r * th.sin()]
        })
        .collect()
}

/// Build the selected bolt as a real 3-D solid **with per-vertex colours**: an
/// **ISO 4017 hex head** (a regular hexagon prism of width-across-flats `s` and
/// head height `k`, extruded along +z) sitting atop a **cylindrical shank**
/// (diameter = nominal `d`, length = the spec's `length_mm`), unioned into one
/// [`valenx_mesh::Mesh`] via [`MeshBuilder`]. Steel grey throughout.
///
/// The shank runs from `z = 0` (free end) up to `z = L` (under-head bearing
/// face); the head sits on top from `z = L` to `z = L + k`. Geometry is in
/// millimetres, matching the dimension readout.
///
/// Threads are **not** modelled — the shank is a clean smooth cylinder (the
/// thread is reported numerically as pitch diameter / tensile stress area in
/// the readout, the same simplification `valenx_fasteners::bolt::to_solid`
/// makes). Returns `(mesh, colors)` with `colors.len() == 3 × triangle_count`,
/// ready for [`crate::WorkspaceProduct::vertex_colors`], or `None` when the
/// spec's nominal diameter / length is non-positive.
fn bolt_mesh_colored(spec: &BoltSpec) -> Option<(valenx_mesh::Mesh, Vec<[f32; 3]>)> {
    let d = spec.thread.nominal_diameter;
    let length = spec.length_mm;
    let s = spec.width_across_flats_mm();
    let k = spec.head_height_mm();
    if !(d.is_finite() && d > 0.0 && length.is_finite() && length > 0.0 && s > 0.0 && k > 0.0) {
        return None;
    }

    let mut b = MeshBuilder::new();
    // Cylindrical shank, axis along +z, spanning z ∈ [0, L] (centre at L/2).
    b.cylinder(
        [0.0, 0.0, length * 0.5],
        [0.0, 0.0, 1.0],
        d * 0.5,
        length,
        SHANK_SEGMENTS,
        STEEL,
    );
    // Hex head: a regular-hexagon prism of across-flats `s`, extruded from the
    // bearing face (z = L) up by the head height `k`.
    b.extrude(&hex_profile(s), length, length + k, STEEL);

    let (mut mesh, colors) = b.into_mesh_and_colors();
    mesh.id = "valenx-fasteners-bolt".to_string();
    Some((mesh, colors))
}

/// Build the hex-head outline for the live preview from the selected bolt's
/// width across flats, best-effort `None` when the nominal size is not in the
/// ISO 4017 table (or has a non-positive across-flats).
fn preview_hex(s: &FastenersWorkbenchState) -> Option<Vec<Vector3<f64>>> {
    let table = iso4017_hex_table();
    let spec = table.iter().find(|b| b.nominal == s.nominal)?;
    let af = spec.width_across_flats_mm();
    if !(af.is_finite() && af > 0.0) {
        return None;
    }
    Some(hexagon(af))
}

/// Draw a closed XY polyline as a face-on (front-view, looking down −Z)
/// wireframe in a fixed-height canvas, the camera auto-framed to the
/// outline. A segment is painted only when both endpoints project in front
/// of the camera, so the render path never panics.
fn draw_hex_preview(ui: &mut egui::Ui, pts: &[Vector3<f64>]) {
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

/// Look up the selected bolt in the ISO 4017 table and format its
/// dimensions. Extracted from the draw closure so it is unit-testable.
fn run_fasteners(s: &mut FastenersWorkbenchState) {
    s.error = None;

    let table = iso4017_hex_table();
    let spec = match table.iter().find(|b| b.nominal == s.nominal) {
        Some(b) => b,
        None => {
            s.error = Some(format!("unknown bolt size '{}'", s.nominal));
            return;
        }
    };

    s.result = format!(
        "designation         : {}  (ISO 4017 hex)\n\
         width across flats  : {:.3} mm\n\
         head height         : {:.3} mm\n\
         pitch diameter      : {:.4} mm\n\
         tensile stress area : {:.3} mm\u{00B2}",
        s.nominal,
        spec.width_across_flats_mm(),
        spec.head_height_mm(),
        spec.pitch_diameter_mm(),
        spec.tensile_stress_area_mm2(),
    );
}

/// Build the **Fasteners** product for the Workbench+Agent bridge — a real
/// 3-D [`crate::WorkspaceProduct`]: the canonical default bolt (`M6`) built as
/// an ISO 4017 hex-head + cylindrical-shank solid ([`bolt_mesh_colored`]),
/// paired with the genuine ISO dimension readout rows ([`run_fasteners`]) and
/// framed at the shared 3/4 hero camera. Registered as the `"fasteners"`
/// producer in [`crate::products_registry::lookup`]; the tile renders the
/// orbitable bolt with its dimensions alongside. Threads are not modelled (the
/// shank is a smooth cylinder); the thread is reported numerically in `lines`.
pub(crate) fn fasteners_product() -> crate::WorkspaceProduct {
    let mut s = FastenersWorkbenchState::default();
    run_fasteners(&mut s);
    let spec = iso4017_hex_table()
        .into_iter()
        .find(|b| b.nominal == s.nominal)
        .expect("default fasteners designation (M6) is in the ISO 4017 table");
    let (mesh, colors) =
        bolt_mesh_colored(&spec).expect("canonical M6 bolt ⇒ a non-empty 3-D solid");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<fasteners>/valenx-bolt");
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Fasteners (ISO 4017 hex bolt)".into(),
        lines: crate::products_registry::lines_from_readout(&s.result),
        mesh: Some(loaded),
        vertex_colors: Some(colors),
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
        let s = FastenersWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn compute_default_m6() {
        let mut s = FastenersWorkbenchState::default();
        run_fasteners(&mut s);
        assert!(s.error.is_none());
        assert!(!s.result.is_empty());
        assert!(s.result.contains("width across flats"));
        assert!(s.result.contains("pitch diameter"));
        assert!(s.result.contains("tensile stress area"));
        // M6: pitch diameter 6 − 0.6495·1 = 5.3505 mm; tensile stress area ≈
        // 20.12 mm². Recompute via the backend table to confirm the default.
        let m6 = iso4017_hex_table()
            .into_iter()
            .find(|b| b.nominal == "M6")
            .expect("M6 is in the ISO 4017 table");
        assert!((m6.pitch_diameter_mm() - 5.3505).abs() < 1e-4);
        assert!((m6.tensile_stress_area_mm2() - 20.12).abs() < 0.05);
        assert!(s.result.contains("5.3505"));
    }

    #[test]
    fn compute_rejects_unknown_size() {
        let mut s = FastenersWorkbenchState {
            nominal: "M999".to_string(),
            ..Default::default()
        };
        run_fasteners(&mut s);
        assert!(s.error.is_some());
        assert!(s.result.is_empty());
    }

    #[test]
    fn hexagon_has_six_vertices_and_correct_across_flats() {
        let pts = hexagon(10.0);
        assert_eq!(pts.len(), 7); // 6 vertices + closing duplicate
        assert_eq!(pts.first(), pts.last());
        assert!(pts.iter().all(|p| p.z == 0.0));
        // All six vertices share the circumradius R = AF/√3.
        let r = 10.0 / 3.0_f64.sqrt();
        assert!(pts
            .iter()
            .all(|p| ((p.x * p.x + p.y * p.y).sqrt() - r).abs() < 1e-9));
        // The six unique vertices are pairwise distinct.
        let verts = &pts[..6];
        for (i, a) in verts.iter().enumerate() {
            for b in &verts[i + 1..] {
                assert!((a - b).norm() > 1e-6);
            }
        }
        // Width across flats = 2 × apothem (centre → edge-midpoint distance).
        let mid = (pts[0] + pts[1]) * 0.5;
        let af = 2.0 * (mid.x * mid.x + mid.y * mid.y).sqrt();
        assert!((af - 10.0).abs() < 1e-9);
    }

    #[test]
    fn preview_hex_default_m6_matches_table_across_flats() {
        let s = FastenersWorkbenchState::default(); // M6
        let pts = preview_hex(&s).expect("M6 is in the ISO 4017 table");
        assert_eq!(pts.len(), 7);
        let af = iso4017_hex_table()
            .into_iter()
            .find(|b| b.nominal == "M6")
            .expect("M6")
            .width_across_flats_mm();
        let mid = (pts[0] + pts[1]) * 0.5;
        let measured_af = 2.0 * (mid.x * mid.x + mid.y * mid.y).sqrt();
        assert!((measured_af - af).abs() < 1e-9);
    }

    #[test]
    fn preview_hex_none_for_unknown_size() {
        let s = FastenersWorkbenchState {
            nominal: "M999".to_string(),
            ..Default::default()
        };
        assert!(preview_hex(&s).is_none());
    }

    #[test]
    fn hex_profile_has_six_distinct_ccw_corners() {
        let p = hex_profile(12.0);
        assert_eq!(p.len(), 6, "open hex profile: 6 distinct corners, no dup");
        // Across-flats = 2 × apothem (centre → edge-midpoint distance).
        let mid = [(p[0][0] + p[1][0]) * 0.5, (p[0][1] + p[1][1]) * 0.5];
        let af = 2.0 * (mid[0] * mid[0] + mid[1] * mid[1]).sqrt();
        assert!((af - 12.0).abs() < 1e-9, "profile across-flats matches");
    }

    #[test]
    fn bolt_mesh_is_a_nonempty_steel_solid_with_aligned_colours() {
        // The default M6 bolt builds a real 3-D solid (head prism + shank
        // cylinder) with a per-vertex colour buffer aligned to the renderer's
        // coloured path (3 entries / triangle), all steel grey.
        let spec = iso4017_hex_table()
            .into_iter()
            .find(|b| b.nominal == "M6")
            .expect("M6 is in the ISO 4017 table");
        let (mesh, colors) = bolt_mesh_colored(&spec).expect("M6 ⇒ a solid");
        assert!(!mesh.nodes.is_empty(), "bolt mesh has vertices");
        assert!(mesh.total_elements() > 0, "bolt mesh has triangles");
        assert_eq!(
            colors.len(),
            mesh.total_elements() * 3,
            "vertex_colors must equal 3 × triangle count (coloured path)"
        );
        assert!(colors.iter().all(|&c| c == STEEL), "all steel grey");
        // The solid spans the shank length + head height in z, and the hex head
        // is wider (across-corners) than the shank radius.
        let max_z = mesh.nodes.iter().map(|p| p.z).fold(f64::MIN, f64::max);
        let min_z = mesh.nodes.iter().map(|p| p.z).fold(f64::MAX, f64::min);
        assert!(min_z >= -1e-9, "shank free end at z=0");
        assert!(
            (max_z - (spec.length_mm + spec.head_height_mm())).abs() < 1e-6,
            "top of head at z = L + k"
        );
    }

    #[test]
    fn bolt_mesh_none_for_degenerate_spec() {
        let mut spec = iso4017_hex_table()
            .into_iter()
            .find(|b| b.nominal == "M6")
            .expect("M6");
        spec.length_mm = 0.0;
        assert!(bolt_mesh_colored(&spec).is_none(), "zero length ⇒ None");
    }

    #[test]
    fn fasteners_product_is_a_real_3d_bolt() {
        // The agent-bridge product now carries a real mesh + aligned colours and
        // keeps its ISO dimension readout rows.
        let product = fasteners_product();
        let loaded = product.mesh.as_ref().expect("fasteners product has a mesh");
        assert!(loaded.mesh.total_elements() > 0, "non-empty bolt mesh");
        let colors = product
            .vertex_colors
            .as_ref()
            .expect("fasteners product carries vertex_colors");
        assert_eq!(
            colors.len(),
            loaded.mesh.total_elements() * 3,
            "product colours aligned to the coloured path"
        );
        assert!(
            product.lines.iter().any(|l| l.contains("width across flats")),
            "keeps the ISO dimension readout"
        );
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
            draw_fasteners_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_fasteners_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_fasteners_workbench = true;
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_a_result_and_error_without_panic() {
        let mut app = ValenxApp::default();
        app.show_fasteners_workbench = true;
        run_fasteners(&mut app.fasteners);
        app.fasteners.error = Some("invalid bolt".to_string());
        draw_workbench(&mut app);
    }
}
