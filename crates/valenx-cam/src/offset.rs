//! 2D polygon offset — per-vertex bisector translation.
//!
//! Offsets a closed XY polygon by a signed distance:
//!
//! - **`distance > 0`** — outward offset (polygon grows).
//! - **`distance < 0`** — inward offset (polygon shrinks). Used by
//!   the Pocket operation to derive the boundary-from-tool-radius
//!   keep-out polygon.
//!
//! ## Algorithm
//!
//! For each vertex `Vᵢ` with incoming edge `Vᵢ₋₁→Vᵢ` and outgoing
//! edge `Vᵢ→Vᵢ₊₁`, we:
//!
//! 1. Compute the left-hand normals `n_in`, `n_out` of the two
//!    edges (for a CCW polygon, these point *outward*).
//! 2. Average them and renormalise → the angle bisector direction.
//! 3. Scale by `distance / sin(α/2)` where α is the turn angle,
//!    so the offset polygon stays parallel to the originals.
//! 4. Translate `Vᵢ` along the bisector.
//!
//! ## v1 simplifications (documented)
//!
//! - **No self-intersection detection.** Heavily-concave polygons
//!   or extreme inward offsets that would collapse a region can
//!   self-cross. The caller (Pocket op) is expected to use sensible
//!   `step_over` values that keep the offset polygon non-self-
//!   intersecting.
//! - **Single-polygon output only.** Returns `Vec<Vec<Vector3<f64>>>`
//!   for API stability with later versions that split disjoint
//!   regions, but the current implementation always returns at most
//!   one polygon (or zero if degenerate).
//! - **180° turns** (cusps) translate along the edge normal — not a
//!   true offset but avoids `NaN` from a degenerate bisector.
//!
//! See [`polygon`] for the entry point.

use nalgebra::{Vector2, Vector3};

/// Offset a closed polygon by `distance` mm.
///
/// Input: a list of vertices in XY (the Z coordinate is preserved
/// from each input vertex). The polygon is treated as closed —
/// `polygon[last]` connects back to `polygon[0]`.
///
/// Output: a list of polygons (each as a `Vec<Vector3<f64>>`). v1
/// always returns either 0 or 1 polygons. A 0-result means the
/// polygon collapsed (e.g. inward offset by more than the
/// inscribed radius).
///
/// See the module docs for the full algorithm + limitations.
pub fn polygon(polygon: &[Vector3<f64>], distance: f64) -> Vec<Vec<Vector3<f64>>> {
    if polygon.len() < 3 {
        return Vec::new();
    }
    if distance.abs() < 1e-12 {
        return vec![polygon.to_vec()];
    }
    let n = polygon.len();
    let mut out = Vec::with_capacity(n);
    // Detect winding to flip the sign convention if needed: we want
    // distance > 0 to push the polygon outward for the caller's
    // sake regardless of input winding.
    let area_2 = signed_area_2(polygon);
    let sign = if area_2 >= 0.0 { 1.0 } else { -1.0 };
    for i in 0..n {
        let prev = polygon[(i + n - 1) % n];
        let curr = polygon[i];
        let next = polygon[(i + 1) % n];
        let e_in = Vector2::new(curr.x - prev.x, curr.y - prev.y);
        let e_out = Vector2::new(next.x - curr.x, next.y - curr.y);
        let e_in_n = e_in.norm();
        let e_out_n = e_out.norm();
        if e_in_n < 1e-12 || e_out_n < 1e-12 {
            // Degenerate edge — keep the vertex in place.
            out.push(curr);
            continue;
        }
        let e_in_u = e_in / e_in_n;
        let e_out_u = e_out / e_out_n;
        // Right-hand normals (perpendicular, -90° rotation): for a
        // CCW polygon these point outward. `sign` flips this for CW.
        let n_in = Vector2::new(e_in_u.y, -e_in_u.x);
        let n_out = Vector2::new(e_out_u.y, -e_out_u.x);
        let bisect = n_in + n_out;
        let bisect_norm = bisect.norm();
        let translate = if bisect_norm < 1e-9 {
            // Near-180° turn (cusp): bisector degenerate. Translate
            // along one edge normal (caller documented in module
            // rustdoc — this is a v1 limitation).
            n_in * distance * sign
        } else {
            let b = bisect / bisect_norm;
            // sin(half-turn) = |bisect| / 2 (for unit normals).
            let sin_half = (bisect_norm * 0.5).clamp(1e-9, 1.0);
            b * (distance * sign / sin_half)
        };
        out.push(Vector3::new(
            curr.x + translate.x,
            curr.y + translate.y,
            curr.z,
        ));
    }
    // Cheap collapse check: if any offset edge points in the
    // opposite direction from the input edge, the offset over-shrank
    // past zero (vertices crossed and the polygon would be inverted
    // or self-overlapping). Also catch the near-zero-area
    // degenerate case.
    if collapsed(&out, area_2, polygon) {
        return Vec::new();
    }
    vec![out]
}

/// Twice the signed XY area of `polygon`. Positive ⇒ CCW; negative ⇒
/// CW. Wraps from `polygon[last]` back to `polygon[0]`.
pub fn signed_area_2(polygon: &[Vector3<f64>]) -> f64 {
    let n = polygon.len();
    if n < 3 {
        return 0.0;
    }
    let mut s = 0.0;
    for i in 0..n {
        let a = polygon[i];
        let b = polygon[(i + 1) % n];
        s += a.x * b.y - b.x * a.y;
    }
    s
}

/// `true` if the offset polygon has degenerate area, its winding
/// flipped, or any edge reversed direction relative to the input —
/// all signal that the offset over-shrank.
fn collapsed(offset: &[Vector3<f64>], original_area_2: f64, original: &[Vector3<f64>]) -> bool {
    let n = offset.len();
    if n < 3 || original.len() != n {
        return true;
    }
    let a = signed_area_2(offset);
    if a.abs() < 1e-6 {
        return true;
    }
    if (a > 0.0) != (original_area_2 > 0.0) {
        return true;
    }
    // Edge-direction check: every offset edge should still point in
    // roughly the same direction as its source edge (positive dot
    // product). A flipped edge means vertices crossed.
    for i in 0..n {
        let a0 = original[i];
        let a1 = original[(i + 1) % n];
        let b0 = offset[i];
        let b1 = offset[(i + 1) % n];
        let e_orig = Vector2::new(a1.x - a0.x, a1.y - a0.y);
        let e_off = Vector2::new(b1.x - b0.x, b1.y - b0.y);
        if e_orig.dot(&e_off) < 0.0 {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: f64, y: f64) -> Vector3<f64> {
        Vector3::new(x, y, 0.0)
    }

    #[test]
    fn unit_square_offset_outward_grows() {
        let sq = vec![p(0.0, 0.0), p(1.0, 0.0), p(1.0, 1.0), p(0.0, 1.0)];
        let out = polygon(&sq, 1.0);
        assert_eq!(out.len(), 1);
        let off = &out[0];
        // Each vertex pushed outward along the bisector by 1.0 *
        // sqrt(2) for a square's 90° corners (since
        // distance / sin(45°) = 1 / (sqrt(2)/2) = sqrt(2)).
        // So a corner at (0,0) goes to (-1, -1).
        assert!((off[0] - p(-1.0, -1.0)).norm() < 1e-9, "got {:?}", off[0]);
        assert!((off[1] - p(2.0, -1.0)).norm() < 1e-9);
        assert!((off[2] - p(2.0, 2.0)).norm() < 1e-9);
        assert!((off[3] - p(-1.0, 2.0)).norm() < 1e-9);
    }

    #[test]
    fn unit_square_offset_inward_shrinks() {
        let sq = vec![p(0.0, 0.0), p(2.0, 0.0), p(2.0, 2.0), p(0.0, 2.0)];
        let out = polygon(&sq, -0.5);
        assert_eq!(out.len(), 1);
        let off = &out[0];
        assert!((off[0] - p(0.5, 0.5)).norm() < 1e-9, "got {:?}", off[0]);
        assert!((off[1] - p(1.5, 0.5)).norm() < 1e-9);
        assert!((off[2] - p(1.5, 1.5)).norm() < 1e-9);
        assert!((off[3] - p(0.5, 1.5)).norm() < 1e-9);
    }

    #[test]
    fn over_shrink_collapses() {
        let sq = vec![p(0.0, 0.0), p(1.0, 0.0), p(1.0, 1.0), p(0.0, 1.0)];
        // Inward by 1.0 > inscribed-half-extent 0.5 ⇒ collapse.
        let out = polygon(&sq, -1.0);
        assert!(out.is_empty(), "expected collapse, got {out:?}");
    }

    #[test]
    fn pentagon_offset_inward_preserves_count() {
        // Regular pentagon at radius 5.
        let mut poly = Vec::with_capacity(5);
        for k in 0..5 {
            let t = (k as f64) * std::f64::consts::TAU / 5.0;
            poly.push(Vector3::new(5.0 * t.cos(), 5.0 * t.sin(), 0.0));
        }
        let out = polygon(&poly, -0.5);
        assert_eq!(out.len(), 1);
        let off = &out[0];
        assert_eq!(off.len(), 5);
        // Every offset vertex should sit at radius ≈ 5 - 0.5 / cos(36°)
        // (since the pentagon's interior bisector forms 36° with the
        // edge normal). Use a loose tolerance.
        let expected = 5.0 - 0.5 / (std::f64::consts::PI / 5.0).cos();
        for v in off {
            let r = (v.x * v.x + v.y * v.y).sqrt();
            assert!((r - expected).abs() < 1e-6, "radius {r} ≠ {expected}");
        }
    }

    #[test]
    fn zero_distance_returns_input() {
        let sq = vec![p(0.0, 0.0), p(1.0, 0.0), p(1.0, 1.0), p(0.0, 1.0)];
        let out = polygon(&sq, 0.0);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0], sq);
    }

    #[test]
    fn cw_polygon_offset_outward_is_consistent() {
        // CW square — offset by +0.5 should still grow.
        let sq = vec![p(0.0, 0.0), p(0.0, 1.0), p(1.0, 1.0), p(1.0, 0.0)];
        let out = polygon(&sq, 0.5);
        assert_eq!(out.len(), 1);
        // Expanded square should have AABB extending below 0 and
        // beyond 1.
        let bb_min_x = out[0].iter().map(|v| v.x).fold(f64::INFINITY, f64::min);
        let bb_max_x = out[0].iter().map(|v| v.x).fold(f64::NEG_INFINITY, f64::max);
        assert!(
            bb_min_x < -0.4 && bb_max_x > 1.4,
            "AABB {bb_min_x}..{bb_max_x}"
        );
    }

    #[test]
    fn signed_area_2_sign() {
        let ccw = vec![p(0.0, 0.0), p(1.0, 0.0), p(1.0, 1.0), p(0.0, 1.0)];
        assert!(signed_area_2(&ccw) > 0.0);
        let cw = vec![p(0.0, 0.0), p(0.0, 1.0), p(1.0, 1.0), p(1.0, 0.0)];
        assert!(signed_area_2(&cw) < 0.0);
    }
}
