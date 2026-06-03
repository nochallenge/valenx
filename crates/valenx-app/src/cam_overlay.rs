//! CAM toolpath simulation overlay rendered atop the viewport.
//!
//! Phase 10E (Task 45): draw the active [`valenx_cam::Toolpath`] as
//! coloured polylines on top of the 3D viewport. Same egui-painter
//! pattern as [`crate::draft_overlay`] and [`crate::sketch_overlay`].
//!
//! ## Visual key
//!
//! - Gray — [`valenx_cam::MoveKind::Rapid`] (non-cutting traversal).
//! - Cyan — [`valenx_cam::MoveKind::Cut`] (XY cut along part).
//! - Red — [`valenx_cam::MoveKind::Plunge`] (vertical plunge).
//!
//! Behind-camera segments are silently culled by
//! [`valenx_viz::project_point`].

use eframe::egui;
use valenx_cam::{simulate::toolpath_polylines, MoveKind, Toolpath};
use valenx_viz::{project_point, OrbitCamera};

/// Draw the toolpath into `rect` using `painter`. Additive — only
/// lines on top of whatever the viewport already rendered.
pub fn draw(painter: &egui::Painter, rect: egui::Rect, camera: &OrbitCamera, toolpath: &Toolpath) {
    let w = rect.width();
    let h = rect.height();
    let origin = rect.min;
    let polylines = toolpath_polylines(toolpath);
    for (kind, points) in polylines {
        let stroke = stroke_for(kind);
        for pair in points.windows(2) {
            let a = match project_point(
                camera,
                w,
                h,
                [pair[0].x as f32, pair[0].y as f32, pair[0].z as f32],
            ) {
                Some(p) => p,
                None => continue,
            };
            let b = match project_point(
                camera,
                w,
                h,
                [pair[1].x as f32, pair[1].y as f32, pair[1].z as f32],
            ) {
                Some(p) => p,
                None => continue,
            };
            painter.line_segment(
                [
                    egui::pos2(origin.x + a.x, origin.y + a.y),
                    egui::pos2(origin.x + b.x, origin.y + b.y),
                ],
                stroke,
            );
        }
    }
}

fn stroke_for(kind: MoveKind) -> egui::Stroke {
    match kind {
        MoveKind::Rapid => egui::Stroke::new(1.0, egui::Color32::from_rgb(140, 140, 140)),
        MoveKind::Cut => egui::Stroke::new(1.5, egui::Color32::from_rgb(80, 220, 220)),
        MoveKind::Plunge => egui::Stroke::new(2.0, egui::Color32::from_rgb(240, 80, 80)),
        // Arc moves render with the cut colour (they're cut moves
        // optimized to G2/G3 by the arc-fitter pass).
        MoveKind::Arc { .. } => egui::Stroke::new(1.5, egui::Color32::from_rgb(120, 220, 180)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;
    use valenx_cam::toolpath::Move;

    #[test]
    fn empty_toolpath_does_not_panic() {
        // Smoke-only: with no rect+painter, we just construct a
        // toolpath and ensure stroke_for handles every variant.
        let _ = stroke_for(MoveKind::Rapid);
        let _ = stroke_for(MoveKind::Cut);
        let _ = stroke_for(MoveKind::Plunge);
        let _ = Toolpath::new();
    }

    #[test]
    fn stroke_colors_differ_by_kind() {
        let r = stroke_for(MoveKind::Rapid);
        let c = stroke_for(MoveKind::Cut);
        let p = stroke_for(MoveKind::Plunge);
        assert_ne!(r.color, c.color);
        assert_ne!(c.color, p.color);
    }

    #[test]
    fn polyline_grouping_via_simulate_matches() {
        // Sanity that the renderer's data source returns the right
        // shape — three polylines for rapid + plunge + cut.
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Rapid, Vector3::new(0.0, 0.0, 5.0), 0.0));
        tp.push(Move::new(
            MoveKind::Plunge,
            Vector3::new(0.0, 0.0, 0.0),
            200.0,
        ));
        tp.push(Move::new(
            MoveKind::Cut,
            Vector3::new(10.0, 0.0, 0.0),
            500.0,
        ));
        let polys = toolpath_polylines(&tp);
        assert_eq!(polys.len(), 3);
        assert_eq!(polys[0].0, MoveKind::Rapid);
        assert_eq!(polys[1].0, MoveKind::Plunge);
        assert_eq!(polys[2].0, MoveKind::Cut);
    }
}
