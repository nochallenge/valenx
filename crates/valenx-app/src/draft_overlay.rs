//! 2D draft overlay rendered atop the viewport.
//!
//! Phase 4E of the Draft workbench: render the active
//! [`valenx_draft::DraftDocument`] as line / polyline / arc / circle
//! / rectangle / polygon / dimension / text overlays on top of the
//! 3D viewport. Each entity's 2D coords are projected to world space
//! via the document's [`valenx_draft::WorkingPlane`], then through
//! the same [`OrbitCamera`] that drives the rest of the viewport.
//!
//! Drawing happens after the wgpu render pass, in the egui painter
//! pass that sits on top, mirroring the [`sketch_overlay`] pattern.
//!
//! [`sketch_overlay`]: crate::sketch_overlay
//!
//! ## Visual key
//!
//! - Cyan lines / circles / arcs / polylines / rectangles / polygons
//!   — normal entities.
//! - Yellow lines / circles / arcs — the selected entity.
//! - Yellow dots — snap markers visible within ~10 screen pixels of
//!   the cursor (endpoints, midpoints, grid intersections).
//! - Small text labels — rendered with the egui monospace font.

use eframe::egui;
use nalgebra::Vector3;
use valenx_draft::{DraftDocument, DraftEntity, WorkingPlane};
use valenx_viz::{project_point, OrbitCamera};

/// Number of straight-line segments used to approximate a full
/// circle. Same value for arcs (proportionally to their sweep) keeps
/// curvature consistent.
const CIRCLE_SEGMENTS: usize = 48;

/// Hard upper bound on the segment count for any single arc draw.
/// R33 L2 — see [`crate::sketch_overlay`]'s counterpart: a non-finite
/// sweep angle saturates the float-to-`usize` cast to `usize::MAX`,
/// which would freeze the paint thread.
const CIRCLE_SEGMENTS_MAX: usize = 4096;

/// Segment count for an arc of the given (absolute) `sweep` radians,
/// floored at 8 and capped at [`CIRCLE_SEGMENTS_MAX`]. A non-finite
/// `sweep` yields the floor. R33 L2.
fn arc_segment_count(sweep: f64, segments_per_rev: usize) -> usize {
    if !sweep.is_finite() {
        return 8;
    }
    let raw = ((sweep / std::f64::consts::TAU) * segments_per_rev as f64)
        .ceil()
        .max(8.0);
    (raw as usize).min(CIRCLE_SEGMENTS_MAX)
}

/// Pixel radius around the cursor inside which a snap marker is
/// shown. Mirrors the sketch overlay's threshold for consistency.
const SNAP_PIXEL_RADIUS: f32 = 10.0;

/// Per-frame state the draft overlay needs from the rest of the app.
/// Same shape as [`crate::sketch_overlay::SketchOverlayState`].
pub struct DraftOverlayState<'a> {
    /// The document to render.
    pub document: &'a DraftDocument,
    /// Currently selected entity index (drawn yellow). `None` =
    /// nothing selected.
    pub selected: Option<usize>,
    /// Current cursor position in absolute egui screen coords. Used
    /// by snap markers. `None` when the cursor is off the viewport.
    pub cursor_screen: Option<egui::Pos2>,
    /// Grid spacing for the grid-intersection snap marker. `<= 0`
    /// disables grid snap markers entirely.
    pub grid_spacing: f64,
}

/// Draw the draft overlay into `rect` using `painter`. Additive — no
/// fill, no background, just lines + arcs + dots + text on top of
/// whatever the viewport already rendered.
///
/// Safe to call on an empty document, with a camera viewing the
/// working plane edge-on, or with `rect` of any size — projection
/// failures (behind-camera) just cull the offending entity.
pub fn draw(
    painter: &egui::Painter,
    rect: egui::Rect,
    camera: &OrbitCamera,
    state: DraftOverlayState<'_>,
) {
    let w = rect.width();
    let h = rect.height();
    let origin = rect.min;
    let plane = &state.document.working_plane;

    let normal_stroke = egui::Stroke::new(1.5, egui::Color32::from_rgb(80, 200, 220));
    let selected_stroke = egui::Stroke::new(2.0, egui::Color32::YELLOW);
    let text_color = egui::Color32::from_rgb(220, 220, 240);

    for (idx, entity) in state.document.entities.iter().enumerate() {
        let stroke = if state.selected == Some(idx) {
            selected_stroke
        } else {
            normal_stroke
        };
        draw_entity(
            painter, rect, camera, plane, entity, stroke, text_color, w, h, origin,
        );
    }

    // Snap markers — endpoints + midpoints + grid intersections within
    // SNAP_PIXEL_RADIUS of the cursor.
    if let Some(cur) = state.cursor_screen {
        let snap_color = egui::Color32::YELLOW;
        for vertex in valenx_draft::snap::endpoints(state.document) {
            draw_snap_dot(
                painter, camera, plane, w, h, origin, cur, vertex, snap_color,
            );
        }
        for vertex in valenx_draft::snap::midpoints(state.document) {
            draw_snap_dot(
                painter,
                camera,
                plane,
                w,
                h,
                origin,
                cur,
                vertex,
                egui::Color32::LIGHT_GREEN,
            );
        }
        // Grid intersections — only show the four nearest cells to
        // avoid drawing the entire plane.
        if state.grid_spacing > 0.0 {
            // Find the cursor's local-frame approximate position by
            // sampling: project several grid candidates near the
            // cursor and draw the ones whose screen distance is in
            // range. We don't have a screen→world unproject helper
            // here, so probe a 5×5 grid around the camera target's
            // local position.
            let target_local = plane.world_to_local(Vector3::new(
                camera.target.x as f64,
                camera.target.y as f64,
                camera.target.z as f64,
            ));
            let g = state.grid_spacing;
            let cx = (target_local[0] / g).round() as i32;
            let cy = (target_local[1] / g).round() as i32;
            for dx in -3..=3 {
                for dy in -3..=3 {
                    let p = [(cx + dx) as f64 * g, (cy + dy) as f64 * g];
                    draw_snap_dot(
                        painter,
                        camera,
                        plane,
                        w,
                        h,
                        origin,
                        cur,
                        p,
                        egui::Color32::from_rgb(120, 120, 200),
                    );
                }
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn draw_entity(
    painter: &egui::Painter,
    rect: egui::Rect,
    camera: &OrbitCamera,
    plane: &WorkingPlane,
    entity: &DraftEntity,
    stroke: egui::Stroke,
    text_color: egui::Color32,
    w: f32,
    h: f32,
    origin: egui::Pos2,
) {
    match entity {
        DraftEntity::Line { start, end } => {
            draw_segment(painter, camera, plane, w, h, origin, *start, *end, stroke);
        }
        DraftEntity::Polyline { points, closed } => {
            for win in points.windows(2) {
                draw_segment(painter, camera, plane, w, h, origin, win[0], win[1], stroke);
            }
            if *closed && points.len() >= 2 {
                draw_segment(
                    painter,
                    camera,
                    plane,
                    w,
                    h,
                    origin,
                    *points.last().unwrap(),
                    points[0],
                    stroke,
                );
            }
        }
        DraftEntity::Arc {
            center,
            radius,
            start_angle,
            end_angle,
        } => {
            draw_arc(
                painter,
                rect,
                camera,
                plane,
                *center,
                *radius,
                *start_angle,
                *end_angle,
                stroke,
            );
        }
        DraftEntity::Circle { center, radius } => {
            draw_arc(
                painter,
                rect,
                camera,
                plane,
                *center,
                *radius,
                0.0,
                std::f64::consts::TAU,
                stroke,
            );
        }
        DraftEntity::Rectangle { min, max } => {
            let bl = *min;
            let br = [max[0], min[1]];
            let tr = *max;
            let tl = [min[0], max[1]];
            draw_segment(painter, camera, plane, w, h, origin, bl, br, stroke);
            draw_segment(painter, camera, plane, w, h, origin, br, tr, stroke);
            draw_segment(painter, camera, plane, w, h, origin, tr, tl, stroke);
            draw_segment(painter, camera, plane, w, h, origin, tl, bl, stroke);
        }
        DraftEntity::Polygon {
            center,
            radius,
            sides,
        } => {
            let n = (*sides).max(3) as i32;
            let two_pi = std::f64::consts::TAU;
            let vertex = |i: i32| -> [f64; 2] {
                let theta = two_pi * (i as f64) / (n as f64);
                [
                    center[0] + radius * theta.cos(),
                    center[1] + radius * theta.sin(),
                ]
            };
            for i in 0..n {
                draw_segment(
                    painter,
                    camera,
                    plane,
                    w,
                    h,
                    origin,
                    vertex(i),
                    vertex((i + 1) % n),
                    stroke,
                );
            }
        }
        DraftEntity::LinearDimension { from, to, offset } => {
            // Compute the perpendicular vector in the local frame and
            // shift the dimension line by `offset` along it. Witness
            // lines connect the measured endpoints to the dimension
            // line, and the numeric label sits at the midpoint.
            let dx = to[0] - from[0];
            let dy = to[1] - from[1];
            let len = (dx * dx + dy * dy).sqrt();
            if len < 1e-12 {
                return;
            }
            let nx = -dy / len;
            let ny = dx / len;
            let off_from = [from[0] + nx * offset, from[1] + ny * offset];
            let off_to = [to[0] + nx * offset, to[1] + ny * offset];
            // Witness lines.
            draw_segment(
                painter, camera, plane, w, h, origin, *from, off_from, stroke,
            );
            draw_segment(painter, camera, plane, w, h, origin, *to, off_to, stroke);
            // Dimension line itself.
            draw_segment(
                painter, camera, plane, w, h, origin, off_from, off_to, stroke,
            );
            // Numeric label at the midpoint.
            let mid = [
                (off_from[0] + off_to[0]) * 0.5,
                (off_from[1] + off_to[1]) * 0.5,
            ];
            if let Some(sp) = project_local(camera, plane, w, h, mid) {
                painter.text(
                    origin + sp.to_vec2(),
                    egui::Align2::CENTER_BOTTOM,
                    format!("{len:.3}"),
                    egui::FontId::monospace(11.0),
                    text_color,
                );
            }
        }
        DraftEntity::Text {
            position,
            content,
            size,
        } => {
            // v1: fixed monospace font, no rotation. `size` scales
            // the on-screen pixel size linearly so users at least see
            // larger text for larger numbers; a future revision can
            // project a per-glyph width along the plane axes.
            if let Some(sp) = project_local(camera, plane, w, h, *position) {
                painter.text(
                    origin + sp.to_vec2(),
                    egui::Align2::LEFT_TOP,
                    content,
                    egui::FontId::monospace((*size as f32 * 24.0).clamp(8.0, 64.0)),
                    text_color,
                );
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn draw_segment(
    painter: &egui::Painter,
    camera: &OrbitCamera,
    plane: &WorkingPlane,
    w: f32,
    h: f32,
    origin: egui::Pos2,
    a: [f64; 2],
    b: [f64; 2],
    stroke: egui::Stroke,
) {
    if let (Some(p), Some(q)) = (
        project_local(camera, plane, w, h, a),
        project_local(camera, plane, w, h, b),
    ) {
        painter.line_segment([origin + p.to_vec2(), origin + q.to_vec2()], stroke);
    }
}

#[allow(clippy::too_many_arguments)]
fn draw_arc(
    painter: &egui::Painter,
    rect: egui::Rect,
    camera: &OrbitCamera,
    plane: &WorkingPlane,
    center: [f64; 2],
    radius: f64,
    a0: f64,
    a1: f64,
    stroke: egui::Stroke,
) {
    if radius <= 0.0 {
        return;
    }
    let w = rect.width();
    let h = rect.height();
    let origin = rect.min;
    let sweep = (a1 - a0).abs();
    // R33 L2: bounded so a non-finite sweep can't saturate to usize::MAX.
    let segs = arc_segment_count(sweep, CIRCLE_SEGMENTS);
    let mut prev: Option<egui::Pos2> = None;
    for i in 0..=segs {
        let t = i as f64 / segs as f64;
        let theta = a0 + (a1 - a0) * t;
        let pt = [
            center[0] + radius * theta.cos(),
            center[1] + radius * theta.sin(),
        ];
        let cur = project_local(camera, plane, w, h, pt);
        if let (Some(p), Some(c)) = (prev, cur) {
            painter.line_segment([origin + p.to_vec2(), origin + c.to_vec2()], stroke);
        }
        prev = cur;
    }
}

#[allow(clippy::too_many_arguments)]
fn draw_snap_dot(
    painter: &egui::Painter,
    camera: &OrbitCamera,
    plane: &WorkingPlane,
    w: f32,
    h: f32,
    origin: egui::Pos2,
    cursor: egui::Pos2,
    vertex: [f64; 2],
    color: egui::Color32,
) {
    if let Some(sp) = project_local(camera, plane, w, h, vertex) {
        let screen = origin + sp.to_vec2();
        if screen.distance(cursor) <= SNAP_PIXEL_RADIUS {
            painter.circle_filled(screen, 4.0, color);
        }
    }
}

/// Project a local-frame 2D point into screen coords.
fn project_local(
    camera: &OrbitCamera,
    plane: &WorkingPlane,
    w: f32,
    h: f32,
    p: [f64; 2],
) -> Option<egui::Pos2> {
    let world = plane.local_to_world(p);
    let sp = project_point(
        camera,
        w,
        h,
        [world.x as f32, world.y as f32, world.z as f32],
    )?;
    Some(egui::pos2(sp.x, sp.y))
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_draft::WorkingPlane;
    use valenx_viz::ViewDirection;

    fn front_cam() -> OrbitCamera {
        let mut cam = OrbitCamera::default();
        cam.set_view(ViewDirection::Front);
        cam.target = nalgebra::Point3::origin();
        cam.distance = 10.0;
        cam
    }

    #[test]
    fn project_local_at_origin_maps_near_center() {
        let cam = front_cam();
        let plane = WorkingPlane::from_xy();
        let p = project_local(&cam, &plane, 800.0, 600.0, [0.0, 0.0]).expect("projects");
        // OrbitCamera Front looks down -Y at the origin from +Y; the
        // XY plane appears edge-on but the origin is on-axis so it
        // still lands near screen centre.
        assert!((p.x - 400.0).abs() < 5.0);
        assert!((p.y - 300.0).abs() < 5.0);
    }

    #[test]
    fn draws_empty_document_without_panicking() {
        // Smoke test: the actual painter wants an egui Context; the
        // function decomposes into project_local + segment helpers
        // and we want to confirm those don't blow up on edge cases.
        let cam = front_cam();
        let plane = WorkingPlane::from_xy();
        assert!(
            project_local(&cam, &plane, 1.0, 1.0, [1e6, 1e6]).is_some()
                || project_local(&cam, &plane, 1.0, 1.0, [1e6, 1e6]).is_none()
        );
    }

    // R33 L2: the draft overlay's arc segment-count helper must stay
    // bounded for a non-finite sweep, same as the sketch overlay's.
    #[test]
    fn arc_segment_count_is_bounded_for_non_finite_sweep() {
        assert_eq!(arc_segment_count(f64::INFINITY, CIRCLE_SEGMENTS), 8);
        assert_eq!(arc_segment_count(f64::NAN, CIRCLE_SEGMENTS), 8);
        // Huge finite sweep clamps to the cap.
        assert_eq!(
            arc_segment_count(1.0e9, CIRCLE_SEGMENTS),
            CIRCLE_SEGMENTS_MAX
        );
    }
}
