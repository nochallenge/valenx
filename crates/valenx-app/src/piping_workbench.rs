//! The right-side **Piping Workbench** panel — native pipe-section sizing
//! over `valenx-piping`.
//!
//! Mirrors the springs / gears / geomatics workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_piping_workbench`,
//! toggled from the View menu. The form picks an NPS designation, a wall
//! schedule, a material, and a length; the "Compute" button reports the
//! outer / inner diameters and the flow / metal cross-section areas, the
//! wetted perimeter, and the external surface area, as a monospace readout.

use eframe::egui;
use nalgebra::Vector3;

use valenx_piping::{Material, PipeSection, PipingError, Schedule};
use valenx_viz::{project_point, OrbitCamera, ViewDirection};

use crate::ValenxApp;

/// The NPS designations the OD lookup table knows (`valenx_piping::dims`).
const NPS_OPTIONS: &[&str] = &[
    "1/8", "1/4", "3/8", "1/2", "3/4", "1", "1-1/4", "1-1/2", "2", "2-1/2", "3", "3-1/2", "4", "5",
    "6", "8", "10", "12", "14", "16", "18", "20", "24", "30", "36",
];

/// Persistent form + result state for the Piping Workbench.
pub struct PipingWorkbenchState {
    /// NPS designation (e.g. `"2"`), keyed into the OD table.
    nominal_size: String,
    /// Wall schedule (Sch 40 / 80 / 160).
    schedule: Schedule,
    /// Pipe material (BOM metadata; does not change the geometry).
    material: Material,
    /// Section length (mm).
    length_mm: f64,
    /// Formatted section readout (empty until the first compute).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
}

impl Default for PipingWorkbenchState {
    fn default() -> Self {
        Self {
            nominal_size: "2".to_string(),
            schedule: Schedule::Sch40,
            material: Material::CarbonSteel,
            length_mm: 1000.0,
            result: String::new(),
            error: None,
        }
    }
}

/// Draw the Piping Workbench right-side panel. A no-op when the
/// `show_piping_workbench` toggle is off.
pub fn draw_piping_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_piping_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_piping_workbench",
        "Piping",
        |app, ui| {
            ui.label(
                egui::RichText::new("native pipe-section sizing · valenx-piping")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.piping;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Nominal pipe size").strong());
                    egui::ComboBox::from_id_source("valenx_piping_nps")
                        .selected_text(format!("NPS {}", s.nominal_size))
                        .show_ui(ui, |ui| {
                            for nps in NPS_OPTIONS {
                                ui.selectable_value(&mut s.nominal_size, (*nps).to_string(), *nps);
                            }
                        });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Schedule").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.schedule, Schedule::Sch40, "Sch 40");
                        ui.radio_value(&mut s.schedule, Schedule::Sch80, "Sch 80");
                        ui.radio_value(&mut s.schedule, Schedule::Sch160, "Sch 160");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Material").strong());
                    ui.horizontal_wrapped(|ui| {
                        ui.radio_value(&mut s.material, Material::CarbonSteel, "Carbon steel");
                        ui.radio_value(&mut s.material, Material::StainlessSteel, "Stainless");
                        ui.radio_value(&mut s.material, Material::Copper, "Copper");
                        ui.radio_value(&mut s.material, Material::Pvc, "PVC");
                        ui.radio_value(&mut s.material, Material::Pex, "PEX");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Geometry").strong());
                    ui.horizontal(|ui| {
                        ui.label("length (mm)");
                        ui.add(egui::DragValue::new(&mut s.length_mm).speed(10.0));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("\u{25B6} Compute").strong())
                        .clicked()
                    {
                        run_piping(s);
                    }

                    if let Some(e) = &s.error {
                        ui.add_space(4.0);
                        ui.colored_label(egui::Color32::from_rgb(220, 90, 90), e);
                    }

                    if !s.result.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new("Section").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }

                    // Live cross-section annulus (OD cyan · ID orange, face-on).
                    if let Some(edges) = preview_rings(s) {
                        ui.add_space(6.0);
                        ui.label(
                            egui::RichText::new("Cross-section preview (OD \u{00B7} ID)").strong(),
                        );
                        draw_rings_preview(ui, &edges);
                    }
                });
        },
    );
    if close {
        app.show_piping_workbench = false;
    }
}

/// `n`-segment closed circle of the given `radius` in the XY plane (z = 0),
/// returned as `n + 1` points with the last equal to the first (so a
/// `windows(2)` walk draws a closed ring with no floating-point gap at θ=2π).
fn circle(radius: f64, n: usize) -> Vec<Vector3<f64>> {
    let n = n.max(1);
    let mut pts: Vec<Vector3<f64>> = (0..n)
        .map(|i| {
            let th = (i as f64 / n as f64) * std::f64::consts::TAU;
            Vector3::new(radius * th.cos(), radius * th.sin(), 0.0)
        })
        .collect();
    let first = pts[0];
    pts.push(first);
    pts
}

/// Build the cross-section preview — the OD ring (cyan), plus the bore/ID
/// ring (orange) when a wall thickness is tabulated, as coloured XY segments.
/// `None` only when the NPS has no tabulated OD (unknown size). The rings
/// depend only on NPS + schedule, so this queries a unit-length section (the
/// form's length feeds only the external-surface readout, not the section).
fn preview_rings(s: &PipingWorkbenchState) -> Option<Vec<([Vector3<f64>; 2], egui::Color32)>> {
    let section = PipeSection::new(
        Vector3::zeros(),
        Vector3::new(1.0, 0.0, 0.0),
        s.nominal_size.clone(),
        s.schedule,
        s.material,
    );
    let od = section.outer_diameter_mm().ok()?;
    if !(od.is_finite() && od > 0.0) {
        return None;
    }
    let cyan = egui::Color32::from_rgb(120, 200, 255);
    let orange = egui::Color32::from_rgb(255, 180, 90);
    let mut edges: Vec<([Vector3<f64>; 2], egui::Color32)> = Vec::new();
    for w in circle(od / 2.0, 64).windows(2) {
        edges.push(([w[0], w[1]], cyan));
    }
    // The bore ring is optional: several NPS/schedule pairs in the dropdown
    // have a tabulated OD but no wall thickness, so the ID is unavailable.
    // Draw the OD ring regardless rather than discarding a valid section.
    let id = section
        .inner_diameter_mm()
        .ok()
        .filter(|&v| v.is_finite() && v > 0.0);
    if let Some(id) = id {
        for w in circle(id / 2.0, 64).windows(2) {
            edges.push(([w[0], w[1]], orange));
        }
    }
    Some(edges)
}

/// Draw coloured XY ring segments as a face-on (front-view, looking down −Z)
/// wireframe in a fixed-height canvas, the camera auto-framed to the rings.
/// A segment is painted only when both endpoints project in front of the
/// camera, so the render path never panics.
fn draw_rings_preview(ui: &mut egui::Ui, edges: &[([Vector3<f64>; 2], egui::Color32)]) {
    let (response, painter) = ui.allocate_painter(
        egui::vec2(ui.available_width(), 200.0),
        egui::Sense::hover(),
    );
    let rect = response.rect;

    let mut min = [f32::INFINITY; 3];
    let mut max = [f32::NEG_INFINITY; 3];
    for (edge, _) in edges {
        for p in edge {
            for k in 0..3 {
                let v = p[k] as f32;
                min[k] = min[k].min(v);
                max[k] = max[k].max(v);
            }
        }
    }

    let mut cam = OrbitCamera::default();
    cam.set_view(ViewDirection::Front);
    cam.frame_bounds(min, max);

    let (w, h) = (rect.width(), rect.height());
    for (edge, color) in edges {
        let a = project_point(
            &cam,
            w,
            h,
            [edge[0].x as f32, edge[0].y as f32, edge[0].z as f32],
        );
        let b = project_point(
            &cam,
            w,
            h,
            [edge[1].x as f32, edge[1].y as f32, edge[1].z as f32],
        );
        if let (Some(a), Some(b)) = (a, b) {
            painter.line_segment(
                [
                    egui::pos2(rect.min.x + a.x, rect.min.y + a.y),
                    egui::pos2(rect.min.x + b.x, rect.min.y + b.y),
                ],
                egui::Stroke::new(1.5, *color),
            );
        }
    }
}

/// Build a [`PipeSection`] from the form, run the section calculations,
/// and format the readout. Extracted from the draw closure so it is
/// unit-testable. Surfaces any backend [`PipingError`] (unknown NPS, or a
/// schedule with no tabulated wall thickness) as the error string.
fn run_piping(s: &mut PipingWorkbenchState) {
    s.error = None;

    if !(s.length_mm.is_finite() && s.length_mm > 0.0) {
        s.error = Some("length must be positive".into());
        return;
    }

    // A horizontal section from the origin; only its length feeds the
    // external-surface calculation, so the orientation is irrelevant.
    let section = PipeSection::new(
        Vector3::zeros(),
        Vector3::new(s.length_mm, 0.0, 0.0),
        s.nominal_size.clone(),
        s.schedule,
        s.material,
    );

    let compute = || -> Result<(f64, f64, f64, f64, f64, f64), PipingError> {
        Ok((
            section.outer_diameter_mm()?,
            section.inner_diameter_mm()?,
            section.flow_area_mm2()?,
            section.metal_cross_section_mm2()?,
            section.wetted_perimeter_mm()?,
            section.external_surface_area_mm2()?,
        ))
    };
    let (od, id, flow, metal, wetted, external) = match compute() {
        Ok(t) => t,
        Err(e) => {
            s.error = Some(e.to_string());
            return;
        }
    };

    s.result = format!(
        "NPS            : {}  ({:?})\n\
         material       : {:?}\n\
         length         : {:.1} mm\n\n\
         outer diameter : {:.3} mm\n\
         inner diameter : {:.3} mm   (OD \u{2212} 2\u{00B7}wall)\n\
         flow area      : {:.2} mm\u{00B2}   (\u{03C0}\u{00B7}ID\u{00B2}/4)\n\
         metal section  : {:.2} mm\u{00B2}   (\u{03C0}\u{00B7}(OD\u{00B2}\u{2212}ID\u{00B2})/4)\n\
         wetted perim.  : {:.3} mm   (\u{03C0}\u{00B7}ID)\n\
         external surf. : {:.0} mm\u{00B2}   (\u{03C0}\u{00B7}OD\u{00B7}L)",
        s.nominal_size, s.schedule, s.material, s.length_mm, od, id, flow, metal, wetted, external,
    );
}

/// Build the **Piping (NPS)** result card for the Workbench+Agent bridge — a
/// DATA-ONLY [`crate::WorkspaceProduct`] (`mesh: None`) whose `lines` are the
/// genuine pipe geometry ([`run_piping`]) for the canonical default selection
/// (NPS 2, Sch 40, carbon steel: OD / ID / flow area / metal section). Registered
/// as the `"piping"` producer in [`crate::products_registry::lookup`]; the tile
/// renders it as a text card, not a 3-D view.
pub(crate) fn piping_product() -> crate::WorkspaceProduct {
    let mut s = PipingWorkbenchState::default();
    run_piping(&mut s);
    crate::WorkspaceProduct {
        title: "Piping (NPS)".into(),
        lines: crate::products_registry::lines_from_readout(&s.result),
        mesh: None,
        vertex_colors: None,
        camera: Default::default(),
        kind2d: None,
        last_export: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = PipingWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn compute_default_nps2_sch40() {
        let mut s = PipingWorkbenchState::default();
        run_piping(&mut s);
        assert!(s.error.is_none());
        assert!(!s.result.is_empty());
        assert!(s.result.contains("outer diameter"));
        assert!(s.result.contains("inner diameter"));
        assert!(s.result.contains("flow area"));
        assert!(s.result.contains("wetted perim"));
        // NPS 2 outer diameter = 2.375 in × 25.4 = 60.325 mm; the bore is
        // strictly smaller. Recompute via the backend to confirm the
        // default is the intended worked example.
        let section = PipeSection::new(
            Vector3::zeros(),
            Vector3::new(1000.0, 0.0, 0.0),
            "2",
            Schedule::Sch40,
            Material::CarbonSteel,
        );
        assert!((section.outer_diameter_mm().unwrap() - 60.325).abs() < 1e-6);
        assert!(section.inner_diameter_mm().unwrap() < section.outer_diameter_mm().unwrap());
    }

    #[test]
    fn compute_rejects_bad_length_and_unknown_nps() {
        // Non-positive length → validation error before any backend call.
        let mut s = PipingWorkbenchState {
            length_mm: 0.0,
            ..Default::default()
        };
        run_piping(&mut s);
        assert!(s.error.is_some());
        assert!(s.result.is_empty());
        // Unknown NPS → the backend's UnknownNps error is surfaced.
        let mut s2 = PipingWorkbenchState {
            nominal_size: "999".to_string(),
            ..Default::default()
        };
        run_piping(&mut s2);
        assert!(s2.error.is_some());
        assert!(s2.result.is_empty());
    }

    #[test]
    fn circle_is_closed_and_on_radius() {
        let pts = circle(5.0, 64);
        assert_eq!(pts.len(), 65); // 64 segments → 65 points
        assert_eq!(pts.first(), pts.last()); // closed exactly (last == first)
                                             // Every vertex lies on the circle of the given radius, in the z=0 plane.
        assert!(pts
            .iter()
            .all(|p| ((p.x * p.x + p.y * p.y).sqrt() - 5.0).abs() < 1e-9));
        assert!(pts.iter().all(|p| p.z == 0.0));
    }

    #[test]
    fn preview_rings_default_has_od_and_id_rings() {
        let s = PipingWorkbenchState::default();
        let edges = preview_rings(&s).expect("NPS 2 Sch40 has a tabulated section");
        let cyan = egui::Color32::from_rgb(120, 200, 255);
        let orange = egui::Color32::from_rgb(255, 180, 90);
        // Both rings present, 64 segments each, all planar (z = 0).
        assert_eq!(edges.iter().filter(|(_, c)| *c == cyan).count(), 64);
        assert_eq!(edges.iter().filter(|(_, c)| *c == orange).count(), 64);
        assert!(edges.iter().all(|(e, _)| e[0].z == 0.0 && e[1].z == 0.0));
        // Cyan ring radius = backend OD/2, orange = ID/2, and ID < OD.
        let section = PipeSection::new(
            Vector3::zeros(),
            Vector3::new(1.0, 0.0, 0.0),
            "2",
            Schedule::Sch40,
            Material::CarbonSteel,
        );
        let od = section.outer_diameter_mm().unwrap();
        let id = section.inner_diameter_mm().unwrap();
        let ring_radius = |want: egui::Color32| {
            edges
                .iter()
                .filter(|(_, c)| *c == want)
                .flat_map(|(e, _)| [e[0], e[1]])
                .map(|p| (p.x * p.x + p.y * p.y).sqrt())
                .fold(0.0_f64, f64::max)
        };
        assert!((ring_radius(cyan) - od / 2.0).abs() < 1e-6);
        assert!((ring_radius(orange) - id / 2.0).abs() < 1e-6);
        assert!(id < od);
    }

    #[test]
    fn preview_rings_none_for_unknown_nps() {
        let s = PipingWorkbenchState {
            nominal_size: "999".to_string(),
            ..Default::default()
        };
        assert!(preview_rings(&s).is_none());
    }

    #[test]
    fn preview_rings_draws_od_only_when_wall_untabulated() {
        // NPS 5 has a tabulated OD but no Sch40 wall, so the backend can give
        // an OD but not an ID. The preview must still render the OD ring
        // (cyan) rather than discarding the whole section.
        let s = PipingWorkbenchState {
            nominal_size: "5".to_string(),
            schedule: Schedule::Sch40,
            ..Default::default()
        };
        // Confirm the precondition this test depends on: OD ok, ID err.
        let section = PipeSection::new(
            Vector3::zeros(),
            Vector3::new(1.0, 0.0, 0.0),
            "5",
            Schedule::Sch40,
            Material::CarbonSteel,
        );
        assert!(section.outer_diameter_mm().is_ok());
        assert!(section.inner_diameter_mm().is_err());

        let edges = preview_rings(&s).expect("OD ring renders even without an ID");
        let cyan = egui::Color32::from_rgb(120, 200, 255);
        let orange = egui::Color32::from_rgb(255, 180, 90);
        assert_eq!(edges.iter().filter(|(_, c)| *c == cyan).count(), 64);
        assert_eq!(edges.iter().filter(|(_, c)| *c == orange).count(), 0);
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
            draw_piping_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_piping_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_piping_workbench = true;
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_a_result_and_error_without_panic() {
        let mut app = ValenxApp::default();
        app.show_piping_workbench = true;
        run_piping(&mut app.piping);
        app.piping.error = Some("invalid pipe parameters".to_string());
        draw_workbench(&mut app);
    }
}
