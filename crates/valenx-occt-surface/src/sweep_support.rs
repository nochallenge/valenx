//! Shared frame-transport machinery for the sweep family (Phases 71,
//! 89, 91).
//!
//! [`sweep_api_pipe`](fn@crate::sweep_api_pipe) already shipped a
//! Bishop / rotation-minimising frame for the basic single-profile
//! pipe. The multi-guide / auxiliary-spine / evolved sweeps in
//! Phases 71, 89 and 91 reuse exactly the same transport — this
//! module hoists the primitives so all four modules share one
//! implementation:
//!
//! - [`vertex_tangents`] — per-vertex unit tangents of a polyline.
//! - [`perp_basis`] — a seed cross-section frame perpendicular to a
//!   tangent.
//! - [`rotate_frame`] — Rodrigues parallel transport of a frame
//!   across one spine segment.
//! - [`arc_length_param`] / [`sample_polyline_at`] — normalised
//!   arc-length matching, used to pair a primary spine with an
//!   auxiliary spine / guide rail of a different point count.
//!
//! `sweep_api_pipe` keeps its own private copies (it predates this
//! module); the new sweeps route through here.

use crate::error::OcctSurfaceError;

/// A 3D point / vector as a plain array — the lingua franca of the
/// sweep modules' public signatures.
pub type Vec3 = [f64; 3];

/// `a - b`.
pub fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// `a + b`.
pub fn add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// `a · s`.
pub fn scale(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Dot product.
pub fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product.
pub fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Euclidean length.
pub fn norm(a: Vec3) -> f64 {
    dot(a, a).sqrt()
}

/// Unit vector — returns the input unchanged if it is (numerically)
/// zero-length.
pub fn normalize(a: Vec3) -> Vec3 {
    let l = norm(a);
    if l < 1e-20 {
        a
    } else {
        scale(a, 1.0 / l)
    }
}

/// Linear interpolation between `a` and `b` at parameter `t`.
pub fn lerp(a: Vec3, b: Vec3, t: f64) -> Vec3 {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
    ]
}

/// Per-vertex unit tangents of a polyline.
///
/// Interior vertices use the normalised average of the incoming and
/// outgoing segment directions; the end vertices use the single
/// adjacent segment.
///
/// # Errors
///
/// [`OcctSurfaceError::BadInput`] if the polyline has fewer than two
/// points or contains a zero-length segment (a vertex with no
/// well-defined tangent).
pub fn vertex_tangents(polyline: &[Vec3]) -> Result<Vec<Vec3>, OcctSurfaceError> {
    let n = polyline.len();
    if n < 2 {
        return Err(OcctSurfaceError::bad_input(
            "polyline",
            "need at least two points to define tangents",
        ));
    }
    let mut tangents = Vec::with_capacity(n);
    for i in 0..n {
        let incoming = if i > 0 {
            sub(polyline[i], polyline[i - 1])
        } else {
            [0.0; 3]
        };
        let outgoing = if i + 1 < n {
            sub(polyline[i + 1], polyline[i])
        } else {
            [0.0; 3]
        };
        let t = match (norm(incoming) > 1e-12, norm(outgoing) > 1e-12) {
            (true, true) => normalize(add(normalize(incoming), normalize(outgoing))),
            (true, false) => normalize(incoming),
            (false, true) => normalize(outgoing),
            (false, false) => {
                return Err(OcctSurfaceError::bad_input(
                    "polyline",
                    "polyline has a zero-length segment",
                ));
            }
        };
        tangents.push(t);
    }
    Ok(tangents)
}

/// Two orthonormal vectors spanning the plane perpendicular to the
/// unit tangent `t`. Seeds the first cross-section frame of a sweep.
pub fn perp_basis(t: Vec3) -> (Vec3, Vec3) {
    let seed = if t[0].abs() <= t[1].abs() && t[0].abs() <= t[2].abs() {
        [1.0, 0.0, 0.0]
    } else if t[1].abs() <= t[2].abs() {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let u = normalize(cross(t, seed));
    let v = cross(t, u);
    (u, v)
}

/// Parallel-transport (Bishop) step: rotate the orthonormal frame
/// `(u, v)` by the minimal rotation that carries the unit tangent
/// `t0` onto the unit tangent `t1`.
///
/// Parallel (no rotation) and antiparallel (180°, degenerate) tangent
/// pairs leave the frame unchanged — a 180° spine kink is unusual and
/// the visual error is local.
pub fn rotate_frame(u: Vec3, v: Vec3, t0: Vec3, t1: Vec3) -> (Vec3, Vec3) {
    let axis = cross(t0, t1);
    let sin_a = norm(axis);
    let cos_a = dot(t0, t1).clamp(-1.0, 1.0);
    if sin_a < 1e-12 {
        return (u, v);
    }
    let k = scale(axis, 1.0 / sin_a);
    (rodrigues(u, k, cos_a, sin_a), rodrigues(v, k, cos_a, sin_a))
}

/// Rodrigues rotation of `x` about the unit axis `k` by an angle whose
/// cosine / sine are `c` / `s`.
pub fn rodrigues(x: Vec3, k: Vec3, c: f64, s: f64) -> Vec3 {
    let kx = cross(k, x);
    let kdotx = dot(k, x);
    [
        x[0] * c + kx[0] * s + k[0] * kdotx * (1.0 - c),
        x[1] * c + kx[1] * s + k[1] * kdotx * (1.0 - c),
        x[2] * c + kx[2] * s + k[2] * kdotx * (1.0 - c),
    ]
}

/// Normalised arc-length parameter (0..=1) of vertex `i` along a
/// polyline — vertex 0 maps to 0, the last vertex to 1, interior
/// vertices to their cumulative-length fraction.
///
/// A zero-length polyline maps every vertex to 0.
pub fn arc_length_param(polyline: &[Vec3], i: usize) -> f64 {
    let n = polyline.len();
    if n < 2 || i == 0 {
        return 0.0;
    }
    let mut cum = vec![0.0_f64; n];
    for k in 1..n {
        cum[k] = cum[k - 1] + norm(sub(polyline[k], polyline[k - 1]));
    }
    let total = cum[n - 1];
    if total < 1e-12 {
        0.0
    } else {
        (cum[i.min(n - 1)] / total).clamp(0.0, 1.0)
    }
}

/// Point on a polyline at normalised arc-length parameter `s` (0..=1).
///
/// `s = 0` is the first vertex, `s = 1` the last; interior values
/// interpolate linearly along the segment that contains that arc
/// length. Used to pair a primary spine with an auxiliary spine /
/// guide rail of a different resolution.
pub fn sample_polyline_at(polyline: &[Vec3], s: f64) -> Vec3 {
    let n = polyline.len();
    if n == 0 {
        return [0.0; 3];
    }
    if n == 1 {
        return polyline[0];
    }
    let s = s.clamp(0.0, 1.0);
    // Cumulative arc length.
    let mut cum = vec![0.0_f64; n];
    for k in 1..n {
        cum[k] = cum[k - 1] + norm(sub(polyline[k], polyline[k - 1]));
    }
    let total = cum[n - 1];
    if total < 1e-12 {
        return polyline[0];
    }
    let target = s * total;
    // Find the segment containing `target`.
    for k in 1..n {
        if cum[k] >= target - 1e-12 {
            let seg = cum[k] - cum[k - 1];
            let local = if seg > 1e-12 {
                (target - cum[k - 1]) / seg
            } else {
                0.0
            };
            return lerp(polyline[k - 1], polyline[k], local.clamp(0.0, 1.0));
        }
    }
    polyline[n - 1]
}

/// Resample a closed polygon to exactly `n` vertices spaced by equal
/// arc length around its perimeter.
///
/// Lets profiles with different point counts be stitched into a
/// common ring — the same trick [`sweep_api_thru_sections`](fn@crate::sweep_api_thru_sections)
/// uses for lofting. A degenerate (zero-perimeter) polygon collapses
/// to `n` copies of its first vertex.
pub fn resample_closed_polygon(poly: &[Vec3], n: usize) -> Vec<Vec3> {
    let m = poly.len();
    if m == 0 || n == 0 {
        return Vec::new();
    }
    let mut seg_len = Vec::with_capacity(m);
    let mut total = 0.0;
    for i in 0..m {
        let d = norm(sub(poly[(i + 1) % m], poly[i]));
        seg_len.push(d);
        total += d;
    }
    if total < 1e-12 {
        return vec![poly[0]; n];
    }
    let step = total / n as f64;
    let mut out = Vec::with_capacity(n);
    let mut seg = 0usize;
    let mut seg_start = 0.0;
    for k in 0..n {
        let target = k as f64 * step;
        while seg + 1 < m && seg_start + seg_len[seg] < target {
            seg_start += seg_len[seg];
            seg += 1;
        }
        let local = if seg_len[seg] > 1e-12 {
            (target - seg_start) / seg_len[seg]
        } else {
            0.0
        };
        out.push(lerp(poly[seg], poly[(seg + 1) % m], local.clamp(0.0, 1.0)));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tangents_of_a_straight_line_all_point_the_same_way() {
        let line = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let t = vertex_tangents(&line).unwrap();
        for tan in &t {
            assert!((tan[0] - 1.0).abs() < 1e-9);
        }
    }

    #[test]
    fn tangents_reject_degenerate_segment() {
        let bad = [[0.0; 3], [0.0; 3], [1.0, 0.0, 0.0]];
        assert!(vertex_tangents(&bad).is_err());
    }

    #[test]
    fn perp_basis_is_orthonormal_to_the_tangent() {
        let t = normalize([1.0, 2.0, 3.0]);
        let (u, v) = perp_basis(t);
        assert!(dot(u, t).abs() < 1e-9);
        assert!(dot(v, t).abs() < 1e-9);
        assert!(dot(u, v).abs() < 1e-9);
        assert!((norm(u) - 1.0).abs() < 1e-9);
        assert!((norm(v) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn rotate_frame_carries_the_frame_with_the_tangent() {
        // A 90° tangent turn: +X → +Y. The frame must stay orthonormal
        // and perpendicular to the new tangent.
        let t0 = [1.0, 0.0, 0.0];
        let t1 = [0.0, 1.0, 0.0];
        let (u0, v0) = perp_basis(t0);
        let (u1, v1) = rotate_frame(u0, v0, t0, t1);
        assert!(dot(u1, t1).abs() < 1e-9, "u must stay ⟂ new tangent");
        assert!(dot(v1, t1).abs() < 1e-9, "v must stay ⟂ new tangent");
        assert!((norm(u1) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn arc_length_param_endpoints() {
        let poly = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [3.0, 4.0, 0.0]];
        assert!(arc_length_param(&poly, 0).abs() < 1e-12);
        assert!((arc_length_param(&poly, 2) - 1.0).abs() < 1e-12);
        // The polyline has two segments: 3 units (0→1) and 4 units
        // (1→2), total length 7. Vertex 1 sits 3 units along → 3/7.
        assert!((arc_length_param(&poly, 1) - 3.0 / 7.0).abs() < 1e-9);
    }

    #[test]
    fn sample_polyline_midpoint() {
        let poly = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let mid = sample_polyline_at(&poly, 0.5);
        assert!((mid[0] - 5.0).abs() < 1e-9);
        let start = sample_polyline_at(&poly, 0.0);
        let end = sample_polyline_at(&poly, 1.0);
        assert!((start[0] - 0.0).abs() < 1e-9);
        assert!((end[0] - 10.0).abs() < 1e-9);
    }

    #[test]
    fn sample_polyline_matches_at_quarter() {
        // An L-shape: 0→(4,0,0)→(4,3,0), total length 7. At s=4/7 we
        // are exactly at the corner.
        let poly = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [4.0, 3.0, 0.0]];
        let corner = sample_polyline_at(&poly, 4.0 / 7.0);
        assert!((corner[0] - 4.0).abs() < 1e-6);
        assert!(corner[1].abs() < 1e-6);
    }

    #[test]
    fn resample_gives_the_requested_count_on_the_perimeter() {
        // A unit square — every resampled point still lies on an edge.
        let square = [
            [-1.0, -1.0, 0.0],
            [1.0, -1.0, 0.0],
            [1.0, 1.0, 0.0],
            [-1.0, 1.0, 0.0],
        ];
        let rs = resample_closed_polygon(&square, 16);
        assert_eq!(rs.len(), 16);
        for p in &rs {
            let on_edge = (p[0].abs() - 1.0).abs() < 1e-9 || (p[1].abs() - 1.0).abs() < 1e-9;
            assert!(on_edge, "resampled point off the perimeter: {p:?}");
        }
    }
}
