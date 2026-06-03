//! Raster scanning of a closed XY polygon.
//!
//! Three fill strategies for the Pocket + Face operations:
//!
//! - [`zigzag`] — back-and-forth lines spaced by `step_over`. Each
//!   line is reversed from the previous one (no lift-and-return),
//!   so the tool sweeps the area in a serpentine.
//! - [`parallel`] — one-way lines (lift-and-return between passes).
//!   Useful when climb-vs-conventional matters and the tool must
//!   only cut in one direction.
//! - [`spiral`] — concentric inward offset polygons, joined into
//!   one polyline. Best for round / convex pockets.
//!
//! ## Algorithm (zigzag + parallel)
//!
//! 1. Compute the polygon's XY axis-aligned bounding box.
//! 2. Pick a step direction `d̂` (unit vector at `angle_deg`) and a
//!    perpendicular step direction `p̂` (90° CCW of `d̂`).
//! 3. Project all polygon vertices onto `p̂` → find the min/max
//!    extent along the perpendicular.
//! 4. For each scan-line position `t = min_p + step_over * k`,
//!    intersect the infinite line `r0 = (cog) + t * p̂ + s * d̂`
//!    with every polygon edge. Sort hits by `s`.
//! 5. Pair up consecutive intersections — odd indices are "leaving
//!    the polygon", even are "entering". Emit one segment per pair.
//! 6. For zigzag, reverse every other scan line.
//!
//! ## v1 simplifications (documented)
//!
//! - **Convex / mostly-convex polygons only.** For heavily-concave
//!   polygons with multiple intersection segments per line, the
//!   "pair up consecutive intersections" rule still works (the
//!   even-odd parity rule for point-in-polygon), but the result is
//!   bumpy and may leave thin slivers unfilled.
//! - **No retract between zigzag lines.** Each segment is its own
//!   polyline; the caller (op generator) is responsible for adding
//!   the "rapid up, traverse, plunge" wrapper when emitting moves
//!   for `parallel`.
//! - **Spiral fill** uses [`crate::offset::polygon`] inward — same
//!   limitations apply (no self-intersection handling).

use nalgebra::{Vector2, Vector3};

use crate::offset;

const EPS: f64 = 1e-9;

/// Zig-zag raster fill — alternating-direction parallel lines
/// clipped to the polygon interior.
///
/// Returns a list of polylines, one per scan line. Each polyline
/// has exactly two vertices (start and end of the in-polygon
/// segment). The caller's op generator wraps these in rapid+plunge
/// motions.
pub fn zigzag(polygon: &[Vector3<f64>], step_over: f64, angle_deg: f64) -> Vec<Vec<Vector3<f64>>> {
    let mut lines = scan_lines(polygon, step_over, angle_deg);
    // Reverse every other line — that's the zigzag part.
    for (i, line) in lines.iter_mut().enumerate() {
        if i % 2 == 1 {
            line.reverse();
        }
    }
    lines
}

/// One-way parallel raster fill — every line in the same direction.
pub fn parallel(
    polygon: &[Vector3<f64>],
    step_over: f64,
    angle_deg: f64,
) -> Vec<Vec<Vector3<f64>>> {
    scan_lines(polygon, step_over, angle_deg)
}

/// Concentric inward offset spiral. Returns a single polyline from
/// the outermost ring inward.
///
/// Each ring is computed by inward-offsetting the previous one by
/// `step_over`. Iteration stops when the offset polygon collapses
/// (see [`crate::offset::polygon`]).
pub fn spiral(polygon: &[Vector3<f64>], step_over: f64) -> Vec<Vector3<f64>> {
    if polygon.len() < 3 || !(step_over > 0.0) {
        return Vec::new();
    }
    let mut out: Vec<Vector3<f64>> = Vec::new();
    let mut current = polygon.to_vec();
    // Cap iterations to avoid pathological infinite loops on
    // weird polygons.
    for _ in 0..4096 {
        // Close the ring by repeating the first vertex.
        out.extend_from_slice(&current);
        if let Some(first) = current.first() {
            out.push(*first);
        }
        let next = offset::polygon(&current, -step_over);
        if next.is_empty() {
            break;
        }
        current = next.into_iter().next().unwrap();
    }
    out
}

/// Build the parallel scan-line segments without any direction
/// alternation. Shared by [`zigzag`] and [`parallel`].
fn scan_lines(polygon: &[Vector3<f64>], step_over: f64, angle_deg: f64) -> Vec<Vec<Vector3<f64>>> {
    if polygon.len() < 3 || !(step_over > 0.0) {
        return Vec::new();
    }
    let theta = angle_deg.to_radians();
    let dir = Vector2::new(theta.cos(), theta.sin());
    // Perpendicular (90° CCW of dir).
    let perp = Vector2::new(-dir.y, dir.x);
    // Compute extents along dir and perp (so we can place scan
    // lines outside the polygon's bbox and clip inward).
    let mut s_min = f64::INFINITY;
    let mut s_max = f64::NEG_INFINITY;
    let mut t_min = f64::INFINITY;
    let mut t_max = f64::NEG_INFINITY;
    for v in polygon {
        let p2 = Vector2::new(v.x, v.y);
        let s = p2.dot(&dir);
        let t = p2.dot(&perp);
        if s < s_min {
            s_min = s;
        }
        if s > s_max {
            s_max = s;
        }
        if t < t_min {
            t_min = t;
        }
        if t > t_max {
            t_max = t;
        }
    }
    let z = polygon[0].z;
    let mut lines: Vec<Vec<Vector3<f64>>> = Vec::new();
    // Start half a step in so the first line is not exactly on the
    // boundary (which causes degenerate single-point intersections).
    let mut t = t_min + step_over * 0.5;
    while t < t_max - EPS {
        let segs = intersect_polygon_with_line(polygon, &dir, &perp, t, z, s_min, s_max);
        lines.extend(segs);
        t += step_over;
    }
    lines
}

/// Intersect an infinite line (defined by `perp · point = t` and
/// direction `dir`) with the closed polygon. Returns a list of
/// in-polygon segments (each as two vertices).
fn intersect_polygon_with_line(
    polygon: &[Vector3<f64>],
    dir: &Vector2<f64>,
    perp: &Vector2<f64>,
    t: f64,
    z: f64,
    _s_min: f64,
    _s_max: f64,
) -> Vec<Vec<Vector3<f64>>> {
    // Collect signed perpendicular distances of each vertex to the line.
    let n = polygon.len();
    let mut crossings: Vec<f64> = Vec::with_capacity(8);
    for i in 0..n {
        let a = polygon[i];
        let b = polygon[(i + 1) % n];
        let a2 = Vector2::new(a.x, a.y);
        let b2 = Vector2::new(b.x, b.y);
        let da = a2.dot(perp) - t;
        let db = b2.dot(perp) - t;
        // Skip edges entirely on one side.
        if (da > EPS && db > EPS) || (da < -EPS && db < -EPS) {
            continue;
        }
        // Both vertices on the line — skip (edge is collinear,
        // covered by neighbouring edges' crossings).
        if da.abs() < EPS && db.abs() < EPS {
            continue;
        }
        // Compute intersection point parameter and project on dir.
        let denom = da - db;
        if denom.abs() < EPS {
            continue;
        }
        let alpha = da / denom;
        let ip = a2 + (b2 - a2) * alpha;
        crossings.push(ip.dot(dir));
    }
    crossings.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    // Pair up: 0-1 = first in-polygon segment, 2-3 = second, etc.
    let mut segs = Vec::new();
    let mut i = 0;
    while i + 1 < crossings.len() {
        let s0 = crossings[i];
        let s1 = crossings[i + 1];
        // Avoid emitting zero-length segments from degenerate
        // double crossings at the same parameter.
        if (s1 - s0).abs() > EPS {
            let p0 = Vector2::new(0.0, 0.0) + *dir * s0 + *perp * t;
            let p1 = Vector2::new(0.0, 0.0) + *dir * s1 + *perp * t;
            segs.push(vec![
                Vector3::new(p0.x, p0.y, z),
                Vector3::new(p1.x, p1.y, z),
            ]);
        }
        i += 2;
    }
    segs
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: f64, y: f64) -> Vector3<f64> {
        Vector3::new(x, y, 0.0)
    }

    #[test]
    fn unit_square_zigzag_horizontal() {
        let sq = vec![p(0.0, 0.0), p(10.0, 0.0), p(10.0, 10.0), p(0.0, 10.0)];
        let lines = zigzag(&sq, 2.0, 0.0);
        // First line at y=1.0 (half a step in), then y=3.0, 5.0, 7.0, 9.0 → 5 lines.
        assert_eq!(lines.len(), 5, "got {} lines", lines.len());
        // Every line spans x=0..10 along y.
        for (i, line) in lines.iter().enumerate() {
            assert_eq!(line.len(), 2, "line {i} has {} vertices", line.len());
            let xs: Vec<f64> = line.iter().map(|v| v.x).collect();
            let min_x = xs.iter().cloned().fold(f64::INFINITY, f64::min);
            let max_x = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            assert!((min_x - 0.0).abs() < 1e-6, "min_x = {min_x}");
            assert!((max_x - 10.0).abs() < 1e-6, "max_x = {max_x}");
        }
    }

    #[test]
    fn zigzag_alternates_direction() {
        let sq = vec![p(0.0, 0.0), p(10.0, 0.0), p(10.0, 10.0), p(0.0, 10.0)];
        let lines = zigzag(&sq, 2.0, 0.0);
        // Line 0: increasing x. Line 1: decreasing x.
        assert!(lines[0][0].x < lines[0][1].x);
        assert!(lines[1][0].x > lines[1][1].x);
        assert!(lines[2][0].x < lines[2][1].x);
    }

    #[test]
    fn parallel_does_not_alternate() {
        let sq = vec![p(0.0, 0.0), p(10.0, 0.0), p(10.0, 10.0), p(0.0, 10.0)];
        let lines = parallel(&sq, 2.0, 0.0);
        assert_eq!(lines.len(), 5);
        // Every line goes in the same direction.
        for line in &lines {
            assert!(line[0].x < line[1].x);
        }
    }

    #[test]
    fn spiral_fills_with_concentric_squares() {
        let sq = vec![p(0.0, 0.0), p(10.0, 0.0), p(10.0, 10.0), p(0.0, 10.0)];
        let path = spiral(&sq, 1.0);
        // Outer ring is 5 verts (square + close), then inner rings shrink.
        // For a 10x10 square, inward by 1 collapses around step 5
        // (half the smaller side). So we get ≈4-5 rings × 5 verts.
        assert!(path.len() >= 15, "spiral too short: {} verts", path.len());
        assert!(path.len() <= 40, "spiral too long: {} verts", path.len());
        // First vertex should be on the input.
        assert_eq!(path[0], p(0.0, 0.0));
    }

    #[test]
    fn zigzag_45deg_works() {
        let sq = vec![p(0.0, 0.0), p(10.0, 0.0), p(10.0, 10.0), p(0.0, 10.0)];
        let lines = zigzag(&sq, 2.0, 45.0);
        // 45° lines through a 10x10 square — should still produce
        // several scan lines.
        assert!(
            lines.len() >= 3,
            "expected ≥3 lines at 45°, got {}",
            lines.len()
        );
    }

    #[test]
    fn empty_polygon_yields_empty() {
        assert!(zigzag(&[], 1.0, 0.0).is_empty());
        assert!(parallel(&[], 1.0, 0.0).is_empty());
        assert!(spiral(&[], 1.0).is_empty());
    }

    #[test]
    fn nonpositive_step_yields_empty() {
        let sq = vec![p(0.0, 0.0), p(10.0, 0.0), p(10.0, 10.0), p(0.0, 10.0)];
        assert!(zigzag(&sq, 0.0, 0.0).is_empty());
        assert!(spiral(&sq, -1.0).is_empty());
    }
}
