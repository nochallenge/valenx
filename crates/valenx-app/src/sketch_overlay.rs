//! 2D sketch overlay rendered atop the viewport.
//!
//! Phase 1H of the sketcher: render the active
//! [`valenx_sketch::Sketch`] as line / circle / arc segments on top of
//! the existing 3D viewport. Sketches live on the XY plane (z=0) in
//! world coordinates and are projected through the same `OrbitCamera`
//! that drives the rest of the viewport — so the user can orbit and
//! zoom and see the sketch tracking the working plane as it would in
//! a real CAD app.
//!
//! Drawing happens after the wgpu render pass, in the egui painter
//! pass that sits on top, mirroring the [`draw_cut_overlay`] pattern
//! already used for the cut-plane preview.
//!
//! [`draw_cut_overlay`]: crate::viewport
//!
//! ## Visual key
//!
//! - White lines / circles / arcs — normal entities.
//! - Yellow lines / circles / arcs — selected entities.
//! - Light-blue rubber-band line — pending two-click line tool preview
//!   from the first click to the current cursor position.
//! - Small light-green dots at entity endpoints — snap markers, shown
//!   when the cursor is within ~10 screen pixels (foundation for the
//!   Phase 2 snap implementation).
//!
//! ## Why XY-plane?
//!
//! Phase 1 doesn't yet have a sketch plane abstraction (that lands in
//! Phase 2 alongside reference geometry). The sensible default is the
//! world XY plane at z=0, which matches the orientation [`extrude`]
//! uses when sweeping the sketch profile along +Z. Users see what
//! they're sketching in the same orientation it will be extruded.
//!
//! [`extrude`]: valenx_sketch::extrude

use eframe::egui;
use valenx_sketch::geom::Entity;
use valenx_sketch::Sketch;
use valenx_viz::{project_point, OrbitCamera};

/// Number of straight-line segments used to approximate a full circle.
/// Same value for arcs (proportionally to their sweep) keeps the
/// arc segments looking as smooth as the circles they match.
const CIRCLE_SEGMENTS: usize = 48;

/// Hard upper bound on the segment count for any single arc draw.
///
/// R33 L2: the segment count is `((sweep / TAU) * N).ceil().max(8) as
/// usize`. A non-finite (`+Inf` / `NaN`) sweep angle — reachable from a
/// hand-edited / corrupt sketch whose arc angle variables are non-finite
/// — makes the float-to-`usize` `as` cast *saturate to `usize::MAX`*, so
/// the subsequent `for i in 0..=segs` paint loop runs an astronomical
/// number of iterations and freezes the UI thread. [`arc_segment_count`]
/// guards non-finite sweeps and clamps to this cap; a cap of 4096 is far
/// smoother than the eye can resolve for any real arc.
const CIRCLE_SEGMENTS_MAX: usize = 4096;

/// Segment count for an arc of the given (absolute) `sweep` radians,
/// proportional to [`CIRCLE_SEGMENTS`] per full revolution, floored at 8
/// (so a tiny arc still looks curved) and capped at
/// [`CIRCLE_SEGMENTS_MAX`]. A non-finite `sweep` yields the floor (8) —
/// the arc geometry is undrawable anyway, but the loop stays bounded.
/// R33 L2.
fn arc_segment_count(sweep: f64, segments_per_rev: usize) -> usize {
    if !sweep.is_finite() {
        return 8;
    }
    let raw = ((sweep / std::f64::consts::TAU) * segments_per_rev as f64)
        .ceil()
        .max(8.0);
    // `raw` is finite and >= 8 here; the `as usize` cast still saturates
    // defensively, then we clamp to the hard cap.
    (raw as usize).min(CIRCLE_SEGMENTS_MAX)
}

/// Pixel radius around the cursor inside which a snap marker is shown.
/// Phase 1H is purely a visual indicator — the Phase 2 snap engine
/// will use the same threshold when consuming clicks.
const SNAP_PIXEL_RADIUS: f32 = 10.0;

/// Per-frame state the overlay needs from the rest of the app —
/// everything *except* the camera, which the overlay borrows from the
/// surrounding viewport state (it's the same camera). Owned up the
/// stack so this module never touches `ValenxApp` directly — keeps it
/// testable with synthetic inputs and easy to move to a separate
/// crate later if the renderer grows.
pub struct SketchOverlayState<'a> {
    /// The sketch to render.
    pub sketch: &'a Sketch,
    /// Entities to highlight as "selected" (drawn yellow).
    pub selected: &'a [valenx_sketch::EntityId],
    /// Optional pending-click preview: when a two-click tool (Line)
    /// has its first click registered, this is the world-space (x, y)
    /// of that first click; the overlay draws a rubber-band line from
    /// it to the current cursor's world-space (x, y).
    pub pending_line_start: Option<(f64, f64)>,
    /// Current cursor position in absolute egui screen coords
    /// (compatible with `painter.line_segment` / `painter.circle_filled`
    /// destinations). `None` when the cursor isn't over the viewport
    /// or isn't available from egui (e.g. before the first frame).
    pub cursor_screen: Option<egui::Pos2>,
}

/// Draw the sketch overlay into `rect` using `painter`. Additive — no
/// fill, no background, just lines + circles + dots on top of whatever
/// the viewport already rendered.
///
/// Safe to call with an empty sketch (no-op), with a camera viewing
/// the XY plane edge-on (entities collapse to a near-line and remain
/// drawable), and with `rect` of any size — projection failures (e.g.
/// behind-camera) just cull the offending entity.
pub fn draw(
    painter: &egui::Painter,
    rect: egui::Rect,
    camera: &OrbitCamera,
    state: SketchOverlayState<'_>,
) {
    let w = rect.width();
    let h = rect.height();
    let origin = rect.min;

    let normal_stroke = egui::Stroke::new(1.5, egui::Color32::WHITE);
    let selected_stroke = egui::Stroke::new(1.5, egui::Color32::YELLOW);

    // Task 51 + 52: render every entity, forking color by selection.
    // Selected entities use a yellow stroke so the user can confirm
    // which entities the constraint buttons will operate on; same
    // selection list the Entities collapsing-section populates via
    // click. Pending-line preview (Task 53) and snap markers (Task
    // 54) layer in on top in subsequent commits.
    for (idx, entity) in state.sketch.entities.iter().enumerate() {
        let id = valenx_sketch::EntityId(idx + 1);
        let stroke = if state.selected.contains(&id) {
            selected_stroke
        } else {
            normal_stroke
        };
        match entity {
            Entity::Point(p) => {
                let (x, y) = p.read(&state.sketch.vars);
                if let Some(sp) = project_world_xy(camera, w, h, x, y) {
                    // Slightly larger dot for selected points so they
                    // pop visually next to similar-radius circles.
                    let r = if state.selected.contains(&id) {
                        4.0
                    } else {
                        3.0
                    };
                    painter.circle_filled(origin + sp.to_vec2(), r, stroke.color);
                }
            }
            Entity::Line(line) => {
                let ((sx, sy), (ex, ey)) = line.endpoints(&state.sketch.vars);
                if let (Some(a), Some(b)) = (
                    project_world_xy(camera, w, h, sx, sy),
                    project_world_xy(camera, w, h, ex, ey),
                ) {
                    painter.line_segment([origin + a.to_vec2(), origin + b.to_vec2()], stroke);
                }
            }
            Entity::Circle(circle) => {
                let (cx, cy) = circle.center.read(&state.sketch.vars);
                let r = circle.radius(&state.sketch.vars);
                draw_circle_xy(painter, rect, camera, cx, cy, r, stroke);
            }
            Entity::Arc(arc) => {
                let (cx, cy) = arc.center.read(&state.sketch.vars);
                let r = arc.radius(&state.sketch.vars);
                let (a0, a1) = arc.angles(&state.sketch.vars);
                draw_arc_xy(painter, rect, camera, ArcXy { cx, cy, r, a0, a1 }, stroke);
            }
            Entity::BSpline(b) => {
                // Phase 12A: sample 24 points and connect with line
                // segments. Dashed if construction-flagged.
                let (u_min, u_max) = b.parameter_range();
                let n_samples = 24;
                let mut prev: Option<egui::Pos2> = None;
                let stroke_to_use = if state.sketch.is_construction(id) {
                    egui::Stroke::new(1.0, egui::Color32::from_gray(160))
                } else {
                    stroke
                };
                for i in 0..=n_samples {
                    let u = u_min + (u_max - u_min) * (i as f64 / n_samples as f64);
                    let pt = b.evaluate(&state.sketch.vars, u);
                    let cur = project_world_xy(camera, w, h, pt[0], pt[1]);
                    if let (Some(p), Some(c)) = (prev, cur) {
                        // Dashed effect for construction: skip every
                        // other segment.
                        if state.sketch.is_construction(id) && (i % 2 == 0) {
                            prev = cur;
                            continue;
                        }
                        painter.line_segment(
                            [origin + p.to_vec2(), origin + c.to_vec2()],
                            stroke_to_use,
                        );
                    }
                    prev = cur;
                }
            }
            Entity::Ellipse(e) => {
                // Sample 64 points around full sweep.
                let n_samples = 64;
                let mut prev: Option<egui::Pos2> = None;
                let stroke_to_use = if state.sketch.is_construction(id) {
                    egui::Stroke::new(1.0, egui::Color32::from_gray(160))
                } else {
                    stroke
                };
                for i in 0..=n_samples {
                    let t = (i as f64 / n_samples as f64) * std::f64::consts::TAU;
                    let (x, y) = e.evaluate(&state.sketch.vars, t);
                    let cur = project_world_xy(camera, w, h, x, y);
                    if let (Some(p), Some(c)) = (prev, cur) {
                        if state.sketch.is_construction(id) && (i % 2 == 0) {
                            prev = cur;
                            continue;
                        }
                        painter.line_segment(
                            [origin + p.to_vec2(), origin + c.to_vec2()],
                            stroke_to_use,
                        );
                    }
                    prev = cur;
                }
            }
            Entity::EllipticalArc(arc) => {
                let (sa, ea) = arc.angles(&state.sketch.vars);
                let sweep = (ea - sa).abs();
                // R33 L2: bound the sample count so a non-finite sweep
                // (corrupt sketch) can't drive `for i in 0..=n_samples`
                // to `usize::MAX` and freeze the paint thread.
                let n_samples = arc_segment_count(sweep, 64);
                let mut prev: Option<egui::Pos2> = None;
                let stroke_to_use = if state.sketch.is_construction(id) {
                    egui::Stroke::new(1.0, egui::Color32::from_gray(160))
                } else {
                    stroke
                };
                for i in 0..=n_samples {
                    let t = sa + (ea - sa) * (i as f64 / n_samples as f64);
                    let (x, y) = arc.evaluate(&state.sketch.vars, t);
                    let cur = project_world_xy(camera, w, h, x, y);
                    if let (Some(p), Some(c)) = (prev, cur) {
                        if state.sketch.is_construction(id) && (i % 2 == 0) {
                            prev = cur;
                            continue;
                        }
                        painter.line_segment(
                            [origin + p.to_vec2(), origin + c.to_vec2()],
                            stroke_to_use,
                        );
                    }
                    prev = cur;
                }
            }
        }
    }

    // Pending-click rubber-band line (Task 53).
    //
    // When the user is mid-way through the two-click Line tool, draw
    // a light-blue line from the projected first-click point to the
    // current cursor position. The cursor read is `ctx.pointer_latest_pos`
    // — absolute egui screen coords — so the rubber band "hovers"
    // wherever the pointer is, ready to commit to a real world-space
    // line on the next click (Phase 2 viewport-click pipeline).
    let pending_stroke = egui::Stroke::new(1.5, egui::Color32::LIGHT_BLUE);
    if let (Some((x0, y0)), Some(cur)) = (state.pending_line_start, state.cursor_screen) {
        if let Some(start) = project_world_xy(camera, w, h, x0, y0) {
            painter.line_segment([origin + start.to_vec2(), cur], pending_stroke);
        }
    }

    // Snap markers (Task 54).
    //
    // Walk every snap-relevant vertex in the sketch (point entities
    // themselves + line endpoints + circle / arc centres) and draw a
    // small light-green filled circle when the cursor is within
    // SNAP_PIXEL_RADIUS of its screen projection. Phase 1H is a
    // visual indicator only — the Phase 2 snap engine consumes the
    // same vertex list when committing a real click, so what the user
    // sees here is exactly what the snap engine will use.
    if let Some(cur) = state.cursor_screen {
        for vertex in iter_snap_vertices(state.sketch) {
            if let Some(sp) = project_world_xy(camera, w, h, vertex.0, vertex.1) {
                let screen = origin + sp.to_vec2();
                if screen.distance(cur) <= SNAP_PIXEL_RADIUS {
                    painter.circle_filled(screen, 4.0, egui::Color32::LIGHT_GREEN);
                }
            }
        }
    }
}

/// Convenience: project a world-space (x, y) point on the XY plane
/// (z=0) to screen-space coords relative to a 0,0 viewport origin.
/// `None` when the point is culled by the near plane or the viewport
/// is degenerate.
fn project_world_xy(camera: &OrbitCamera, w: f32, h: f32, x: f64, y: f64) -> Option<egui::Pos2> {
    let sp = project_point(camera, w, h, [x as f32, y as f32, 0.0])?;
    Some(egui::pos2(sp.x, sp.y))
}

/// Draw a full circle on the XY plane as an N-segment polyline.
/// Resolution is fixed (see [`CIRCLE_SEGMENTS`]); a future pass could
/// scale resolution with on-screen radius if jaggies show up on
/// large-zoom circles.
fn draw_circle_xy(
    painter: &egui::Painter,
    rect: egui::Rect,
    camera: &OrbitCamera,
    cx: f64,
    cy: f64,
    r: f64,
    stroke: egui::Stroke,
) {
    if r <= 0.0 {
        return;
    }
    draw_arc_xy(
        painter,
        rect,
        camera,
        ArcXy {
            cx,
            cy,
            r,
            a0: 0.0,
            a1: std::f64::consts::TAU,
        },
        stroke,
    );
}

/// Parameters describing an arc on the XY plane.
struct ArcXy {
    /// Center X.
    cx: f64,
    /// Center Y.
    cy: f64,
    /// Radius (positive; zero or negative is treated as no-op).
    r: f64,
    /// Start angle (radians).
    a0: f64,
    /// End angle (radians). Sweep direction is `sign(a1 - a0)`.
    a1: f64,
}

/// Draw an arc on the XY plane as a polyline from `arc.a0` to `arc.a1`,
/// sampling enough points that one full revolution would use
/// [`CIRCLE_SEGMENTS`] segments. Negative sweep (CW) is fine — we
/// just step in the right direction.
fn draw_arc_xy(
    painter: &egui::Painter,
    rect: egui::Rect,
    camera: &OrbitCamera,
    arc: ArcXy,
    stroke: egui::Stroke,
) {
    let ArcXy { cx, cy, r, a0, a1 } = arc;
    if r <= 0.0 {
        return;
    }
    let w = rect.width();
    let h = rect.height();
    let origin = rect.min;
    // Segments proportional to arc sweep, minimum 8 so a tiny arc
    // still looks curved, capped (R33 L2) so a non-finite sweep can't
    // drive the loop to usize::MAX.
    let sweep = (a1 - a0).abs();
    let segs = arc_segment_count(sweep, CIRCLE_SEGMENTS);
    let mut prev: Option<egui::Pos2> = None;
    for i in 0..=segs {
        let t = i as f64 / segs as f64;
        let theta = a0 + (a1 - a0) * t;
        let x = cx + r * theta.cos();
        let y = cy + r * theta.sin();
        let cur = project_world_xy(camera, w, h, x, y);
        if let (Some(p), Some(c)) = (prev, cur) {
            painter.line_segment([origin + p.to_vec2(), origin + c.to_vec2()], stroke);
        }
        prev = cur;
    }
}

/// Iterator over every snap-relevant vertex in the sketch — point
/// entities themselves, plus the explicit endpoints of every line,
/// plus the centres of circles and arcs. Phase 2's snap engine will
/// extend this to closest-point-on-circle for the perimeter snap;
/// for Phase 1H the centre-only snap is enough to feel useful when
/// users want to drop a coincident point on an existing circle's
/// centre.
fn iter_snap_vertices(sketch: &Sketch) -> impl Iterator<Item = (f64, f64)> + '_ {
    sketch.entities.iter().flat_map(move |e| {
        let v: Vec<(f64, f64)> = match e {
            Entity::Point(p) => vec![p.read(&sketch.vars)],
            Entity::Line(l) => {
                let (s, e) = l.endpoints(&sketch.vars);
                vec![s, e]
            }
            // Circles + arcs: snap to the centre only for now — the
            // perimeter snap is computed differently (closest-point-on-
            // circle) and lands with the Phase 2 snap engine.
            Entity::Circle(c) => vec![c.center.read(&sketch.vars)],
            Entity::Arc(a) => vec![a.center.read(&sketch.vars)],
            // Phase 12A: snap to BSpline control points and ellipse
            // centres.
            Entity::BSpline(b) => b
                .control_points
                .iter()
                .map(|cp| cp.read(&sketch.vars))
                .collect(),
            Entity::Ellipse(e) => vec![e.center.read(&sketch.vars)],
            Entity::EllipticalArc(ea) => vec![ea.ellipse.center.read(&sketch.vars)],
        };
        v.into_iter()
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_sketch::Sketch;
    use valenx_viz::ViewDirection;

    fn front_cam() -> OrbitCamera {
        let mut cam = OrbitCamera::default();
        cam.set_view(ViewDirection::Front);
        cam.target = nalgebra::Point3::origin();
        cam.distance = 10.0;
        cam
    }

    // R33 L2: the arc segment-count helper must stay bounded for a
    // non-finite sweep. Pre-fix the inline expression cast `+Inf as
    // usize` → usize::MAX, so `for i in 0..=segs` would iterate ~1.8e19
    // times and freeze the paint thread.
    #[test]
    fn arc_segment_count_is_bounded_for_infinite_sweep() {
        assert_eq!(arc_segment_count(f64::INFINITY, CIRCLE_SEGMENTS), 8);
        assert_eq!(arc_segment_count(f64::NEG_INFINITY, 64), 8);
        assert_eq!(arc_segment_count(f64::NAN, CIRCLE_SEGMENTS), 8);
    }

    #[test]
    fn arc_segment_count_clamps_huge_finite_sweep_to_cap() {
        // A finite-but-absurd sweep (e.g. 1e9 radians from a runaway
        // solve) clamps to the hard cap rather than allocating millions
        // of segments.
        let segs = arc_segment_count(1.0e9, CIRCLE_SEGMENTS);
        assert_eq!(segs, CIRCLE_SEGMENTS_MAX);
    }

    #[test]
    fn arc_segment_count_normal_sweeps_are_proportional() {
        // Half a revolution → CIRCLE_SEGMENTS/2 (24), floored at 8.
        assert_eq!(arc_segment_count(std::f64::consts::PI, CIRCLE_SEGMENTS), 24);
        // A tiny sweep still gets the floor of 8.
        assert_eq!(arc_segment_count(0.001, CIRCLE_SEGMENTS), 8);
    }

    #[test]
    fn project_world_xy_at_origin_maps_near_center() {
        let cam = front_cam();
        let p = project_world_xy(&cam, 800.0, 600.0, 0.0, 0.0).expect("projects");
        assert!((p.x - 400.0).abs() < 2.0);
        assert!((p.y - 300.0).abs() < 2.0);
    }

    #[test]
    fn iter_snap_vertices_returns_points_and_line_endpoints() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(1.0, 0.0);
        let _line = s.add_line(a, b).unwrap();
        // Sketch now contains: Point(a), Point(b), Line(a→b).
        // Snap vertex iterator should yield: a, b, a (line.start), b (line.end).
        let v: Vec<_> = iter_snap_vertices(&s).collect();
        assert_eq!(v.len(), 4);
        assert_eq!(v[0], (0.0, 0.0));
        assert_eq!(v[1], (1.0, 0.0));
        assert_eq!(v[2], (0.0, 0.0));
        assert_eq!(v[3], (1.0, 0.0));
    }

    #[test]
    fn iter_snap_vertices_yields_circle_center() {
        let mut s = Sketch::new();
        let c = s.add_point(5.0, 7.0);
        let _circle = s.add_circle(c, 2.0).unwrap();
        // Vertices: point(5,7), circle.center(5,7).
        let v: Vec<_> = iter_snap_vertices(&s).collect();
        assert_eq!(v.len(), 2);
        assert!(v.iter().all(|p| *p == (5.0, 7.0)));
    }

    #[test]
    fn empty_sketch_yields_no_snap_vertices() {
        let s = Sketch::new();
        assert_eq!(iter_snap_vertices(&s).count(), 0);
    }
}
