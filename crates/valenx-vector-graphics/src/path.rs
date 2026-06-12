//! Path utilities — bounding box + arc-length.

use nalgebra::Vector2;

use crate::entity::PathSegment;

/// Axis-aligned bounding box of a path. Returns `((min_x, min_y),
/// (max_x, max_y))`. Returns the origin pair for an empty path.
pub fn bbox(path: &[PathSegment]) -> (Vector2<f64>, Vector2<f64>) {
    let mut lo = Vector2::new(f64::INFINITY, f64::INFINITY);
    let mut hi = Vector2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut cur = Vector2::zeros();
    let mut start_of_subpath = Vector2::zeros();
    let add = |p: Vector2<f64>, lo: &mut Vector2<f64>, hi: &mut Vector2<f64>| {
        if p.x < lo.x {
            lo.x = p.x;
        }
        if p.y < lo.y {
            lo.y = p.y;
        }
        if p.x > hi.x {
            hi.x = p.x;
        }
        if p.y > hi.y {
            hi.y = p.y;
        }
    };
    for seg in path {
        match seg {
            PathSegment::MoveTo(p) => {
                cur = *p;
                start_of_subpath = *p;
                add(*p, &mut lo, &mut hi);
            }
            PathSegment::LineTo(p) => {
                add(*p, &mut lo, &mut hi);
                cur = *p;
            }
            PathSegment::CurveTo { c1, c2, end } => {
                // Bezier convex hull bbox.
                add(*c1, &mut lo, &mut hi);
                add(*c2, &mut lo, &mut hi);
                add(*end, &mut lo, &mut hi);
                cur = *end;
            }
            PathSegment::QuadTo { c, end } => {
                add(*c, &mut lo, &mut hi);
                add(*end, &mut lo, &mut hi);
                cur = *end;
            }
            PathSegment::Arc { rx, ry, end, .. } => {
                add(*end, &mut lo, &mut hi);
                // Conservative: expand by (rx, ry) around the
                // current point.
                add(cur - Vector2::new(*rx, *ry), &mut lo, &mut hi);
                add(cur + Vector2::new(*rx, *ry), &mut lo, &mut hi);
                cur = *end;
            }
            PathSegment::Close => {
                cur = start_of_subpath;
            }
        }
    }
    if !lo.x.is_finite() {
        return (Vector2::zeros(), Vector2::zeros());
    }
    (lo, hi)
}

/// Arc-length of a path. Curves are subdivided into `steps` linear
/// segments; arcs use the angular radius bound `pi * (rx + ry) / 2`.
pub fn length(path: &[PathSegment]) -> f64 {
    const STEPS: usize = 32;
    let mut total = 0.0;
    let mut cur = Vector2::zeros();
    let mut start_of_subpath = Vector2::zeros();
    for seg in path {
        match seg {
            PathSegment::MoveTo(p) => {
                cur = *p;
                start_of_subpath = *p;
            }
            PathSegment::LineTo(p) => {
                total += (p - cur).norm();
                cur = *p;
            }
            PathSegment::CurveTo { c1, c2, end } => {
                let mut prev = cur;
                for k in 1..=STEPS {
                    let t = k as f64 / STEPS as f64;
                    let p = cubic_bezier(cur, *c1, *c2, *end, t);
                    total += (p - prev).norm();
                    prev = p;
                }
                cur = *end;
            }
            PathSegment::QuadTo { c, end } => {
                let mut prev = cur;
                for k in 1..=STEPS {
                    let t = k as f64 / STEPS as f64;
                    let p = quad_bezier(cur, *c, *end, t);
                    total += (p - prev).norm();
                    prev = p;
                }
                cur = *end;
            }
            PathSegment::Arc { rx, ry, end, .. } => {
                // Conservative quarter-circumference bound.
                total += std::f64::consts::PI * (rx + ry) * 0.5;
                cur = *end;
            }
            PathSegment::Close => {
                total += (start_of_subpath - cur).norm();
                cur = start_of_subpath;
            }
        }
    }
    total
}

fn cubic_bezier(
    a: Vector2<f64>,
    b: Vector2<f64>,
    c: Vector2<f64>,
    d: Vector2<f64>,
    t: f64,
) -> Vector2<f64> {
    let u = 1.0 - t;
    a * (u * u * u) + b * (3.0 * u * u * t) + c * (3.0 * u * t * t) + d * (t * t * t)
}

fn quad_bezier(a: Vector2<f64>, b: Vector2<f64>, c: Vector2<f64>, t: f64) -> Vector2<f64> {
    let u = 1.0 - t;
    a * (u * u) + b * (2.0 * u * t) + c * (t * t)
}
