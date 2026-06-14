//! **Class-A continuity inspection** — measure how well two surfaces meet
//! along a shared edge, the way a surfacing tool's "connect checker" does.
//!
//! Where [`crate::sew`] *enforces* continuity (stitching control points to a
//! target class), this module *measures* it: sample the shared boundary and
//! report the deviation at each continuity level —
//!
//! - **G0 (position)**: the positional gap `‖Pₐ − P_b‖` between the two edges;
//! - **G1 (tangent)**: the angle between the surfaces' tangent planes (their
//!   unit normals) across the seam, in degrees;
//! - **G2 (curvature)**: the difference in surface mean curvature across the
//!   seam (a curvature-continuity indicator).
//!
//! That is exactly the diagnostic Class-A surfacing depends on: a body looks
//! "Class-A" when adjacent patches are at least G2 (curvature-continuous) so
//! reflections flow unbroken across the seams.
//!
//! Honest scope: a sampled, research-grade connect-checker. It assumes the two
//! edges are parametrised over matching arcs (with an explicit `reverse` flag
//! for the common opposite-direction case); it does not reparametrise mismatched
//! edges, and the G2 metric is the mean-curvature difference (a scalar
//! indicator), not the full directional second-fundamental-form comparison a
//! production checker reports. It is a measurement tool toward — not equal to —
//! CATIA-class surface analysis.

use crate::nurbs_surface::NurbsSurface;
use crate::sew::Edge;

/// `(u, v)` on `surf`'s `edge` at fractional position `t ∈ [0, 1]`.
fn edge_uv(surf: &NurbsSurface, edge: Edge, t: f64) -> (f64, f64) {
    let (u0, u1) = surf.u_range();
    let (v0, v1) = surf.v_range();
    let lerp = |a: f64, b: f64| a + t * (b - a);
    match edge {
        Edge::UMin => (u0, lerp(v0, v1)),
        Edge::UMax => (u1, lerp(v0, v1)),
        Edge::VMin => (lerp(u0, u1), v0),
        Edge::VMax => (lerp(u0, u1), v1),
    }
}

/// The deviation between two surfaces along a shared edge, at each continuity
/// level.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContinuityReport {
    /// Number of samples taken along the edge.
    pub samples: usize,
    /// Maximum positional gap `‖Pₐ − P_b‖` (model units) — the G0 deviation.
    pub g0_max_gap: f64,
    /// Mean positional gap.
    pub g0_mean_gap: f64,
    /// Maximum tangent-plane (normal) angle across the seam (degrees) — the G1
    /// deviation.
    pub g1_max_angle_deg: f64,
    /// Mean tangent-plane angle (degrees).
    pub g1_mean_angle_deg: f64,
    /// Maximum mean-curvature difference across the seam — the G2 deviation
    /// indicator.
    pub g2_max_curvature_diff: f64,
}

impl ContinuityReport {
    /// Best continuity class the seam achieves at the given tolerances:
    /// `2` (G2), `1` (G1), `0` (G0), or `-1` (the edges do not even meet).
    pub fn continuity_class(
        &self,
        position_tol: f64,
        angle_tol_deg: f64,
        curvature_tol: f64,
    ) -> i32 {
        if self.g0_max_gap > position_tol {
            return -1;
        }
        if self.g1_max_angle_deg > angle_tol_deg {
            return 0;
        }
        if self.g2_max_curvature_diff > curvature_tol {
            return 1;
        }
        2
    }
}

/// Measure the continuity between surface `a`'s `a_edge` and surface `b`'s
/// `b_edge`, sampling `samples + 1` points along the shared boundary. Set
/// `b_reversed` when `b`'s edge runs in the opposite parametric direction to
/// `a`'s (so `a` at `t` is paired with `b` at `1 − t`).
pub fn measure_edge_continuity(
    a: &NurbsSurface,
    a_edge: Edge,
    b: &NurbsSurface,
    b_edge: Edge,
    b_reversed: bool,
    samples: usize,
) -> ContinuityReport {
    let n = samples.max(1);
    let (mut g0_max, mut g0_sum) = (0.0_f64, 0.0_f64);
    let (mut g1_max, mut g1_sum) = (0.0_f64, 0.0_f64);
    let mut g2_max = 0.0_f64;
    let count = n + 1;
    for i in 0..=n {
        let ta = i as f64 / n as f64;
        let tb = if b_reversed { 1.0 - ta } else { ta };
        let (ua, va) = edge_uv(a, a_edge, ta);
        let (ub, vb) = edge_uv(b, b_edge, tb);

        let gap = (a.evaluate(ua, va) - b.evaluate(ub, vb)).norm();
        g0_max = g0_max.max(gap);
        g0_sum += gap;

        let na = a.normal(ua, va);
        let nb = b.normal(ub, vb);
        // Tangent-plane deviation: angle between the unit normals, taking the
        // acute angle so an opposite-facing-but-coplanar pair still reads 0.
        let cos = na.dot(&nb).abs().clamp(0.0, 1.0);
        let ang = cos.acos().to_degrees();
        g1_max = g1_max.max(ang);
        g1_sum += ang;

        let dh = (a.mean_curvature(ua, va) - b.mean_curvature(ub, vb)).abs();
        g2_max = g2_max.max(dh);
    }
    ContinuityReport {
        samples: count,
        g0_max_gap: g0_max,
        g0_mean_gap: g0_sum / count as f64,
        g1_max_angle_deg: g1_max,
        g1_mean_angle_deg: g1_sum / count as f64,
        g2_max_curvature_diff: g2_max,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    /// A flat bilinear patch from its four corners (u0v0, u1v0, u0v1, u1v1).
    fn plane(
        p_u0v0: Vector3<f64>,
        p_u1v0: Vector3<f64>,
        p_u0v1: Vector3<f64>,
        p_u1v1: Vector3<f64>,
    ) -> NurbsSurface {
        NurbsSurface::new(
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            vec![vec![p_u0v0, p_u0v1], vec![p_u1v0, p_u1v1]],
            vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        )
        .expect("valid bilinear patch")
    }

    /// A degree-2 patch that arches in +z away from its `v = 0` edge (so it has
    /// nonzero mean curvature), spanning x∈[0,1], y∈[1,2].
    fn arched_patch() -> NurbsSurface {
        let z_mid = 0.45;
        let cp = |x: f64, y: f64, z: f64| Vector3::new(x, y, z);
        let xs = [0.0, 0.5, 1.0];
        let rows: Vec<Vec<Vector3<f64>>> = xs
            .iter()
            .map(|&x| {
                vec![
                    cp(x, 1.0, 0.0),   // v=0 seam (y=1, z=0)
                    cp(x, 1.5, z_mid), // v mid, arched up
                    cp(x, 2.0, 0.0),   // v=1
                ]
            })
            .collect();
        NurbsSurface::new(
            2,
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            rows,
            vec![vec![1.0; 3]; 3],
        )
        .expect("valid arched patch")
    }

    #[test]
    fn coplanar_planes_are_g2() {
        // Two unit squares in z=0, adjacent along y=1. A's VMax meets B's VMin.
        let a = plane(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
        );
        let b = plane(
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 2.0, 0.0),
            Vector3::new(1.0, 2.0, 0.0),
        );
        let r = measure_edge_continuity(&a, Edge::VMax, &b, Edge::VMin, false, 24);
        assert!(r.g0_max_gap < 1e-9, "g0 {}", r.g0_max_gap);
        assert!(r.g1_max_angle_deg < 1e-6, "g1 {}", r.g1_max_angle_deg);
        assert!(
            r.g2_max_curvature_diff < 1e-6,
            "g2 {}",
            r.g2_max_curvature_diff
        );
        assert_eq!(r.continuity_class(1e-6, 0.5, 1e-4), 2);
    }

    #[test]
    fn dihedral_angle_is_a_g1_break() {
        let theta = 30.0_f64.to_radians();
        let a = plane(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
        );
        // B shares the y=1 edge but tilts up out of plane by `theta`.
        let (c, s) = (theta.cos(), theta.sin());
        let b = plane(
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0 + c, s),
            Vector3::new(1.0, 1.0 + c, s),
        );
        let r = measure_edge_continuity(&a, Edge::VMax, &b, Edge::VMin, false, 24);
        assert!(r.g0_max_gap < 1e-9, "edges still meet, g0 {}", r.g0_max_gap);
        assert!(
            (r.g1_max_angle_deg - 30.0).abs() < 0.5,
            "g1 should equal the 30° dihedral, got {}",
            r.g1_max_angle_deg
        );
        // Both patches are planar → no curvature break; class is G0 only.
        assert_eq!(r.continuity_class(1e-6, 1.0, 1e-4), 0);
    }

    #[test]
    fn curvature_break_shows_in_g2() {
        // Flat patch meeting a curved (arched) patch along y=1: the arched
        // patch has nonzero mean curvature, so the seam is not G2.
        let flat = plane(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
        );
        let curved = arched_patch();
        let r = measure_edge_continuity(&flat, Edge::VMax, &curved, Edge::VMin, false, 24);
        assert!(r.g0_max_gap < 1e-9, "edges meet, g0 {}", r.g0_max_gap);
        assert!(
            r.g2_max_curvature_diff > 1e-2,
            "curvature break should register in g2, got {}",
            r.g2_max_curvature_diff
        );
    }

    #[test]
    fn a_surface_is_g2_continuous_with_itself() {
        let a = arched_patch();
        let r = measure_edge_continuity(&a, Edge::VMax, &a, Edge::VMax, false, 16);
        assert!(r.g0_max_gap < 1e-9);
        assert!(r.g1_max_angle_deg < 1e-6);
        assert!(r.g2_max_curvature_diff < 1e-6);
        assert_eq!(r.continuity_class(1e-6, 0.5, 1e-4), 2);
    }

    #[test]
    fn is_deterministic() {
        let a = arched_patch();
        let b = arched_patch();
        let r1 = measure_edge_continuity(&a, Edge::VMin, &b, Edge::VMin, false, 20);
        let r2 = measure_edge_continuity(&a, Edge::VMin, &b, Edge::VMin, false, 20);
        assert_eq!(r1, r2);
    }
}
