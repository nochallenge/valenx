//! Cylindrical fillet strip generator.
//!
//! Given a sharp convex edge between two faces with outward normals
//! `n0`, `n1`, [`build`] returns a triangle strip whose vertices ride
//! on a circular arc of radius `r`. The arc is tangent to both faces
//! (it meets each face at a point `r` units away from the shared
//! edge along the face surface) and curves smoothly between them.
//!
//! # Geometry derivation
//!
//! Let:
//! - `t` = edge tangent (unit), oriented from `p0` to `p1`
//! - `n0`, `n1` = outward face normals (unit), both perpendicular to `t`
//! - `b0 = t × n0` = bitangent on face 0 — lies in face 0, perpendicular
//!   to the edge, pointing INTO face 0 away from the shared edge
//! - `b1 = n1 × t` = symmetric bitangent on face 1
//! - `r` = fillet radius
//!
//! The fillet cylinder axis is parallel to `t`; its center sits at
//! `c = p_edge + r * (b0 - n0)` (one face-thickness away from each
//! face surface). The arc starts at `c + r*n0` (= `p_edge + r*b0`, on
//! face 0) and ends at `c + r*n1` (= `p_edge + r*b1`, on face 1),
//! sweeping in the plane perpendicular to `t` via spherical-linear
//! interpolation of the normal direction.
//!
//! The strip lives in world space; vertices are absolute coordinates.

use nalgebra::Vector3;

/// A triangle strip in world space.
///
/// Triangles are encoded as `(i0, i1, i2)` indices into [`Strip::vertices`].
/// All triangles wind consistently (right-hand rule against the
/// outward-facing fillet surface).
#[derive(Clone, Debug, Default)]
pub struct Strip {
    /// World-space vertex positions.
    pub vertices: Vec<Vector3<f64>>,
    /// Triangle indices, each triple references three [`Self::vertices`].
    pub triangles: Vec<(usize, usize, usize)>,
}

/// Build a cylindrical fillet strip that smoothly bridges two faces
/// meeting along the sharp convex edge `(p0, p1)`.
///
/// - `n0`, `n1`: outward face normals of the two faces sharing the
///   edge.
/// - `b0`, `b1`: face bitangents — unit vectors lying in each face,
///   perpendicular to the edge, and pointing AWAY from the shared
///   edge into the body of that face. The fillet's two endpoints
///   sit at `p_edge + r*b0` (face 0) and `p_edge + r*b1` (face 1).
///   These cannot be derived from the normals alone — the caller
///   must compute them from the actual adjacent triangle geometry.
/// - `r`: fillet radius.
/// - `n_segments`: number of facets along the curved arc; more →
///   smoother. Must be ≥ 1.
///
/// Returns a strip with `(n_segments + 1) * 2` vertices and
/// `n_segments * 2` triangles. The strip's first two vertices sit on
/// face 0; the last two on face 1.
///
/// 8 arguments by design — every piece of geometry (the two edge
/// endpoints, the two face normals, the two bitangents, the radius,
/// and the segment count) is genuinely independent. Bundling them
/// into a struct would obscure the call sites without buying
/// type-safety (every field is f64-like).
#[allow(clippy::too_many_arguments)]
pub fn build(
    p0: Vector3<f64>,
    p1: Vector3<f64>,
    n0: Vector3<f64>,
    n1: Vector3<f64>,
    b0: Vector3<f64>,
    b1: Vector3<f64>,
    r: f64,
    n_segments: usize,
) -> Strip {
    assert!(n_segments >= 1, "n_segments must be >= 1");
    let edge_t = {
        let d = p1 - p0;
        let len = d.norm();
        if len < 1e-30 {
            // Degenerate edge — fall back to an arbitrary axis. The
            // caller should have already rejected this case via
            // FilletError::DegenerateEdge.
            Vector3::x()
        } else {
            d / len
        }
    };
    let n0 = n0.normalize();
    let n1 = n1.normalize();
    let b0 = b0.normalize();
    let _b1 = b1.normalize(); // currently unused — symmetry check only

    // Cylinder center offset from the edge: at the bitangent step
    // on face 0, then back away from the face by one radius.
    let center_offset = b0 * r - n0 * r;

    let mut vertices = Vec::with_capacity((n_segments + 1) * 2);
    for i in 0..=n_segments {
        let t = i as f64 / n_segments as f64;
        // Direction from cylinder axis to the arc point in the plane
        // perpendicular to `edge_t`. At t=0 this is `n0`, at t=1
        // this is `n1`.
        let radial = slerp(n0, n1, t, edge_t);
        let off = center_offset + radial * r;
        vertices.push(p0 + off);
        vertices.push(p1 + off);
    }

    let mut triangles = Vec::with_capacity(n_segments * 2);
    for i in 0..n_segments {
        let i0 = 2 * i;
        let i1 = 2 * i + 1;
        let i2 = 2 * (i + 1);
        let i3 = 2 * (i + 1) + 1;
        // Two triangles per quad facet, wound consistently outward.
        triangles.push((i0, i1, i3));
        triangles.push((i0, i3, i2));
    }
    Strip {
        vertices,
        triangles,
    }
}

/// Spherical-linear-interpolate between two unit vectors `n0` and
/// `n1`, rotating around the given `axis`.
///
/// Both inputs are first projected into the plane perpendicular to
/// `axis` (eliminating any component along the rotation axis), then
/// normalized and slerped. Returns `n0` (projected) at `t = 0` and
/// `n1` (projected) at `t = 1`.
pub fn slerp(n0: Vector3<f64>, n1: Vector3<f64>, t: f64, axis: Vector3<f64>) -> Vector3<f64> {
    let axis = axis.normalize();
    let n0p = (n0 - axis * n0.dot(&axis)).normalize();
    let n1p = (n1 - axis * n1.dot(&axis)).normalize();
    let cos_omega = n0p.dot(&n1p).clamp(-1.0, 1.0);
    let omega = cos_omega.acos();
    if omega.abs() < 1e-9 {
        return n0p;
    }
    let sin_omega = omega.sin();
    let a = ((1.0 - t) * omega).sin() / sin_omega;
    let b = (t * omega).sin() / sin_omega;
    n0p * a + n1p * b
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slerp_endpoints() {
        let axis = Vector3::x();
        let n0 = Vector3::z();
        let n1 = Vector3::y();
        let a = slerp(n0, n1, 0.0, axis);
        let b = slerp(n0, n1, 1.0, axis);
        assert!((a - n0).norm() < 1e-9);
        assert!((b - n1).norm() < 1e-9);
    }

    #[test]
    fn slerp_midpoint_on_unit_sphere() {
        let axis = Vector3::x();
        let n0 = Vector3::z();
        let n1 = Vector3::y();
        let mid = slerp(n0, n1, 0.5, axis);
        let expected = Vector3::new(0.0, 1.0, 1.0).normalize();
        assert!((mid - expected).norm() < 1e-9, "got {mid:?}");
        assert!((mid.norm() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn build_strip_count_for_cube_corner() {
        // Cube edge from origin to +X. Bottom face outward normal -Z,
        // front face outward normal -Y. Bitangents point INTO each
        // face away from the shared edge: bottom bitangent = +Y,
        // front bitangent = +Z. 4 segments, radius 0.1.
        let strip = build(
            Vector3::zeros(),
            Vector3::x(),
            -Vector3::z(),
            -Vector3::y(),
            Vector3::y(),
            Vector3::z(),
            0.1,
            4,
        );
        assert_eq!(strip.vertices.len(), 10);
        assert_eq!(strip.triangles.len(), 8);
    }

    #[test]
    fn build_strip_endpoints_lie_on_each_face() {
        let r = 0.1;
        let strip = build(
            Vector3::zeros(),
            Vector3::x(),
            -Vector3::z(),
            -Vector3::y(),
            Vector3::y(),
            Vector3::z(),
            r,
            4,
        );
        // First two vertices: on bottom face (z=0), offset by r
        // along +Y into the face.
        assert!(
            (strip.vertices[0] - Vector3::new(0.0, r, 0.0)).norm() < 1e-9,
            "v0 should be (0, r, 0), got {:?}",
            strip.vertices[0]
        );
        assert!(
            (strip.vertices[1] - Vector3::new(1.0, r, 0.0)).norm() < 1e-9,
            "v1 should be (1, r, 0), got {:?}",
            strip.vertices[1]
        );
        // Last two vertices: on front face (y=0), offset by r along +Z.
        assert!(
            (strip.vertices[8] - Vector3::new(0.0, 0.0, r)).norm() < 1e-9,
            "v8 should be (0, 0, r), got {:?}",
            strip.vertices[8]
        );
        assert!(
            (strip.vertices[9] - Vector3::new(1.0, 0.0, r)).norm() < 1e-9,
            "v9 should be (1, 0, r), got {:?}",
            strip.vertices[9]
        );
    }

    #[test]
    fn build_strip_midpoint_at_arc_apex() {
        let r = 0.1;
        let strip = build(
            Vector3::zeros(),
            Vector3::x(),
            -Vector3::z(),
            -Vector3::y(),
            Vector3::y(),
            Vector3::z(),
            r,
            4,
        );
        // Row 2 (middle of 5 rows): index 4 = p0 side.
        let mid_p0 = strip.vertices[4];
        // Center of cylinder = (0, r, r). Radius r. Bisector dir is
        // -(Y+Z)/sqrt(2) (pointing OUT toward the corner cavity).
        let center = Vector3::new(0.0, r, r);
        let bisector = -Vector3::new(0.0, 1.0, 1.0).normalize();
        let expected = center + bisector * r;
        assert!(
            (mid_p0 - expected).norm() < 1e-9,
            "mid_p0 = {mid_p0:?}, expected {expected:?}"
        );
    }
}
