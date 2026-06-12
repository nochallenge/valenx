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
use valenx_cam::{
    simulate::{estimated_time, toolpath_polylines},
    MoveKind, Toolpath,
};
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

    // Toolpath-analytics readout — a monospace stats box at the top-left
    // of the viewport, drawn only when a toolpath is active.
    let stats = toolpath_stats(toolpath);
    if stats.moves > 0 {
        let bbox = match stats.extents {
            Some((dx, dy, dz)) => format!("{dx:.1}\u{00D7}{dy:.1}\u{00D7}{dz:.1} mm"),
            None => "\u{2014}".to_string(),
        };
        let text = format!(
            "CAM  {} moves \u{00B7} {:.1} mm total  (rapid {:.1} / cut {:.1})\n\
             ~{:.2} min \u{00B7} bbox {}",
            stats.moves, stats.total_mm, stats.rapid_mm, stats.cut_mm, stats.est_time_min, bbox,
        );
        painter.text(
            rect.left_top() + egui::vec2(8.0, 8.0),
            egui::Align2::LEFT_TOP,
            text,
            egui::FontId::monospace(12.0),
            egui::Color32::from_rgb(200, 220, 220),
        );
    }
}

/// Aggregate toolpath analytics — move count, total / rapid / cut path
/// length (mm), estimated cycle time (minutes), and the bounding-box
/// extents. A pure function over the public [`Toolpath`] API.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct CamStats {
    moves: usize,
    total_mm: f64,
    rapid_mm: f64,
    cut_mm: f64,
    est_time_min: f64,
    extents: Option<(f64, f64, f64)>,
}

fn toolpath_stats(tp: &Toolpath) -> CamStats {
    // Split the per-segment chord length (as `total_distance` measures it)
    // by the destination move's kind — the same classification
    // `estimated_time` uses, so rapid + cut == total.
    let mut rapid_mm = 0.0;
    let mut cut_mm = 0.0;
    for w in tp.moves.windows(2) {
        let d = (w[1].position - w[0].position).norm();
        if w[1].kind == MoveKind::Rapid {
            rapid_mm += d;
        } else {
            cut_mm += d;
        }
    }
    let extents = tp
        .bounding_box()
        .map(|(lo, hi)| (hi.x - lo.x, hi.y - lo.y, hi.z - lo.z));
    CamStats {
        moves: tp.len(),
        total_mm: tp.total_distance(),
        rapid_mm,
        cut_mm,
        est_time_min: estimated_time(tp),
        extents,
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

    #[test]
    fn toolpath_stats_splits_rapid_and_cut() {
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Rapid, Vector3::new(0.0, 0.0, 5.0), 0.0));
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(10.0, 0.0, 5.0),
            0.0,
        ));
        tp.push(Move::new(
            MoveKind::Plunge,
            Vector3::new(10.0, 0.0, 0.0),
            200.0,
        ));
        tp.push(Move::new(
            MoveKind::Cut,
            Vector3::new(20.0, 0.0, 0.0),
            500.0,
        ));
        let s = toolpath_stats(&tp);
        assert_eq!(s.moves, 4);
        // One rapid segment (10 mm); plunge 5 + cut 10 = 15 mm cutting.
        assert!((s.rapid_mm - 10.0).abs() < 1e-9);
        assert!((s.cut_mm - 15.0).abs() < 1e-9);
        // Every segment is classified, so rapid + cut == total_distance.
        assert!((s.rapid_mm + s.cut_mm - s.total_mm).abs() < 1e-9);
        assert!((s.total_mm - tp.total_distance()).abs() < 1e-9);
        assert_eq!(s.extents, Some((20.0, 0.0, 5.0)));
        // The cycle-time estimate matches the backend estimator (minutes).
        assert!((s.est_time_min - estimated_time(&tp)).abs() < 1e-12);
    }

    #[test]
    fn toolpath_stats_empty_is_zero_safe() {
        let s = toolpath_stats(&Toolpath::new());
        assert_eq!(s.moves, 0);
        assert_eq!(s.total_mm, 0.0);
        assert_eq!(s.rapid_mm, 0.0);
        assert_eq!(s.cut_mm, 0.0);
        assert_eq!(s.extents, None);
    }

    #[test]
    fn draw_with_and_without_toolpath_does_not_panic() {
        let cam = OrbitCamera::default();
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Rapid, Vector3::new(0.0, 0.0, 5.0), 0.0));
        tp.push(Move::new(
            MoveKind::Cut,
            Vector3::new(10.0, 0.0, 0.0),
            500.0,
        ));
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                let rect = ui.max_rect();
                let painter = ui.painter();
                draw(painter, rect, &cam, &tp);
                draw(painter, rect, &cam, &Toolpath::new());
            });
        });
    }
}
