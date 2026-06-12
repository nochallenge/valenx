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

use valenx_fasteners::bolt::iso4017_hex_table;
use valenx_viz::{project_point, OrbitCamera, ViewDirection};

use crate::ValenxApp;

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

    egui::SidePanel::right("valenx_fasteners_workbench")
        .resizable(true)
        .default_width(360.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            ui.heading("Fasteners");
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
        });
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
