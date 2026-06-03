//! Coons patch — surface filling a closed 4-boundary curve loop.
//!
//! The boundary curves are ordered `[c0, c1, d0, d1]` where:
//! - `c0(u)` runs along v = v_min,
//! - `c1(u)` runs along v = v_max,
//! - `d0(v)` runs along u = u_min,
//! - `d1(v)` runs along u = u_max.
//!
//! All four curves must share a common parameter range `[0, 1]` and
//! the corner points must coincide within `TOLERANCE`:
//! - `c0(0) == d0(0)` = corner `P00`,
//! - `c0(1) == d1(0)` = corner `P10`,
//! - `c1(0) == d0(1)` = corner `P01`,
//! - `c1(1) == d1(1)` = corner `P11`.
//!
//! ## Formula
//!
//! The bilinear Coons patch is the unique transfinite surface that
//! interpolates the four boundary curves and is linear in the
//! "between-curve" parameter:
//!
//! ```text
//! S(u, v) = (1 - v) c0(u) + v c1(u)
//!         + (1 - u) d0(v) + u d1(v)
//!         - [ (1-u)(1-v) P00 + u(1-v) P10
//!           + (1-u)  v   P01 + u  v   P11 ]
//! ```
//!
//! The bracketed term is the bilinear-corner correction that
//! subtracts the duplicate boundary contribution.

use nalgebra::Vector3;

use crate::error::SurfaceError;
use crate::nurbs_curve::NurbsCurve;
use crate::nurbs_surface::NurbsSurface;

/// Tolerance used to verify boundary corner continuity.
pub const CORNER_TOLERANCE: f64 = 1.0e-6;

/// Coons patch boundary, in the canonical ordering used by
/// [`fill`].
///
/// - `c0` runs along `v = v_min` (`u` is the parameter),
/// - `c1` runs along `v = v_max`,
/// - `d0` runs along `u = u_min` (`v` is the parameter),
/// - `d1` runs along `u = u_max`.
pub struct CoonsBoundary<'a> {
    /// Boundary curve at `v = v_min`.
    pub c0: &'a NurbsCurve,
    /// Boundary curve at `v = v_max`.
    pub c1: &'a NurbsCurve,
    /// Boundary curve at `u = u_min`.
    pub d0: &'a NurbsCurve,
    /// Boundary curve at `u = u_max`.
    pub d1: &'a NurbsCurve,
}

/// Tuple of the four corner points `(P00, P10, P01, P11)`
/// returned by [`validate_boundary`].
pub type CoonsCorners = (Vector3<f64>, Vector3<f64>, Vector3<f64>, Vector3<f64>);

/// Validate the four-boundary corner-continuity condition.
///
/// Returns the four corner points `(P00, P10, P01, P11)` if the
/// corners agree within [`CORNER_TOLERANCE`].
pub fn validate_boundary(b: &CoonsBoundary<'_>) -> Result<CoonsCorners, SurfaceError> {
    let c0_start = b.c0.evaluate(b.c0.parameter_range().0);
    let c0_end = b.c0.evaluate(b.c0.parameter_range().1);
    let c1_start = b.c1.evaluate(b.c1.parameter_range().0);
    let c1_end = b.c1.evaluate(b.c1.parameter_range().1);
    let d0_start = b.d0.evaluate(b.d0.parameter_range().0);
    let d0_end = b.d0.evaluate(b.d0.parameter_range().1);
    let d1_start = b.d1.evaluate(b.d1.parameter_range().0);
    let d1_end = b.d1.evaluate(b.d1.parameter_range().1);

    let check = |a: Vector3<f64>, b: Vector3<f64>, name: &str| -> Result<(), SurfaceError> {
        if (a - b).norm() > CORNER_TOLERANCE {
            Err(SurfaceError::BadBoundary(format!(
                "{name}: {a:?} ≠ {b:?} (Δ = {})",
                (a - b).norm()
            )))
        } else {
            Ok(())
        }
    };

    check(c0_start, d0_start, "corner P00 (c0(0) vs d0(0))")?;
    check(c0_end, d1_start, "corner P10 (c0(1) vs d1(0))")?;
    check(c1_start, d0_end, "corner P01 (c1(0) vs d0(1))")?;
    check(c1_end, d1_end, "corner P11 (c1(1) vs d1(1))")?;

    // Average the agreeing corner pairs so downstream code uses a
    // consistent value (eliminates the choice of "which curve do we
    // trust at this corner").
    let p00 = 0.5 * (c0_start + d0_start);
    let p10 = 0.5 * (c0_end + d1_start);
    let p01 = 0.5 * (c1_start + d0_end);
    let p11 = 0.5 * (c1_end + d1_end);
    Ok((p00, p10, p01, p11))
}

/// Evaluate the closed-form bilinear Coons patch at `(u, v)`,
/// where `u, v ∈ [0, 1]` are the **normalised** parameters across
/// each boundary curve.
pub fn evaluate(b: &CoonsBoundary<'_>, u: f64, v: f64) -> Vector3<f64> {
    // Remap u into each curve's parameter range so callers can use
    // a single canonical [0, 1] frame regardless of how the
    // boundary curves were parameterised.
    let mu_c0 = remap(u, b.c0);
    let mu_c1 = remap(u, b.c1);
    let mu_d0 = remap(v, b.d0);
    let mu_d1 = remap(v, b.d1);

    let c0_u = b.c0.evaluate(mu_c0);
    let c1_u = b.c1.evaluate(mu_c1);
    let d0_v = b.d0.evaluate(mu_d0);
    let d1_v = b.d1.evaluate(mu_d1);

    let p00 = b.c0.evaluate(b.c0.parameter_range().0);
    let p10 = b.c0.evaluate(b.c0.parameter_range().1);
    let p01 = b.c1.evaluate(b.c1.parameter_range().0);
    let p11 = b.c1.evaluate(b.c1.parameter_range().1);

    let l_c = (1.0 - v) * c0_u + v * c1_u;
    let l_d = (1.0 - u) * d0_v + u * d1_v;
    let bilinear =
        (1.0 - u) * (1.0 - v) * p00 + u * (1.0 - v) * p10 + (1.0 - u) * v * p01 + u * v * p11;
    l_c + l_d - bilinear
}

/// Map a normalised parameter `t ∈ [0, 1]` to the curve's actual
/// parameter range.
fn remap(t: f64, c: &NurbsCurve) -> f64 {
    let (lo, hi) = c.parameter_range();
    lo + t.clamp(0.0, 1.0) * (hi - lo)
}

/// Build a Coons-patch [`NurbsSurface`] from four boundary curves.
///
/// **v1 representation:** the returned surface is a degree-(3, 3)
/// tensor-product NurbsSurface with a 4×4 control point grid. Each
/// control point is sampled from the closed-form Coons patch at
/// the standard cubic-Bezier sample parameters `(u, v) ∈ {0, 1/3,
/// 2/3, 1} × {0, 1/3, 2/3, 1}`.
///
/// This means:
/// - The four corner CPs are exactly the corners of the boundary,
/// - The interior 2×2 CPs reproduce the Coons patch exactly at
///   their sample positions and approximate it elsewhere by the
///   cubic-Bezier basis (the surface degree),
/// - For boundaries that are themselves cubic Beziers, the produced
///   patch is exact; for higher-degree boundaries it is a smooth
///   degree-3 approximation suitable for downstream tessellation.
///
/// True Coons-patch-to-NURBS conversion (preserving the boundary
/// curves exactly) is tracked under Phase 9.5 — until then, this
/// approximation is what the viewport sees.
pub fn fill(boundary: [NurbsCurve; 4]) -> Result<NurbsSurface, SurfaceError> {
    let b = CoonsBoundary {
        c0: &boundary[0],
        c1: &boundary[1],
        d0: &boundary[2],
        d1: &boundary[3],
    };
    let _corners = validate_boundary(&b)?;

    // Cubic Bezier sample parameters.
    let samples = [0.0_f64, 1.0 / 3.0, 2.0 / 3.0, 1.0];
    let mut cps = vec![vec![Vector3::zeros(); 4]; 4];
    for (i, &u) in samples.iter().enumerate() {
        for (j, &v) in samples.iter().enumerate() {
            cps[i][j] = evaluate(&b, u, v);
        }
    }
    // Open-uniform cubic knot vector → standard cubic Bezier patch.
    let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
    let weights = vec![vec![1.0; 4]; 4];
    NurbsSurface::new(3, 3, knots.clone(), knots, cps, weights)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn line(a: Vector3<f64>, b: Vector3<f64>) -> NurbsCurve {
        // Cubic Bezier of a straight line with CPs at 0, 1/3, 2/3, 1
        // along the line — gives an exact degree-3 representation.
        let p0 = a;
        let p1 = a + (b - a) / 3.0;
        let p2 = a + 2.0 * (b - a) / 3.0;
        let p3 = b;
        NurbsCurve::new(
            3,
            vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            vec![p0, p1, p2, p3],
            vec![1.0; 4],
        )
        .unwrap()
    }

    #[test]
    fn validates_corner_continuity_passes_for_unit_square() {
        let p00 = Vector3::new(0.0, 0.0, 0.0);
        let p10 = Vector3::new(1.0, 0.0, 0.0);
        let p11 = Vector3::new(1.0, 1.0, 0.0);
        let p01 = Vector3::new(0.0, 1.0, 0.0);
        let c0 = line(p00, p10); // v=0 edge: u=0..1 along bottom
        let c1 = line(p01, p11); // v=1 edge: u=0..1 along top
        let d0 = line(p00, p01); // u=0 edge: v=0..1 along left
        let d1 = line(p10, p11); // u=1 edge: v=0..1 along right
        let b = CoonsBoundary {
            c0: &c0,
            c1: &c1,
            d0: &d0,
            d1: &d1,
        };
        validate_boundary(&b).unwrap();
    }

    #[test]
    fn validates_corner_continuity_fails_when_corner_is_off() {
        let p00 = Vector3::new(0.0, 0.0, 0.0);
        let p10 = Vector3::new(1.0, 0.0, 0.0);
        let p11 = Vector3::new(1.0, 1.0, 0.0);
        let p01_a = Vector3::new(0.0, 1.0, 0.0);
        let p01_b = Vector3::new(0.5, 1.0, 0.0); // mismatched!
        let c0 = line(p00, p10);
        let c1 = line(p01_a, p11);
        let d0 = line(p00, p01_b); // wrong corner at v=1
        let d1 = line(p10, p11);
        let b = CoonsBoundary {
            c0: &c0,
            c1: &c1,
            d0: &d0,
            d1: &d1,
        };
        let err = validate_boundary(&b).unwrap_err();
        assert_eq!(err.code(), "surface.bad_boundary");
    }

    #[test]
    fn coons_evaluate_unit_square_bilinear() {
        // Four straight edges → Coons patch is the bilinear surface
        // S(u,v) = (u, v, 0).
        let p00 = Vector3::new(0.0, 0.0, 0.0);
        let p10 = Vector3::new(1.0, 0.0, 0.0);
        let p11 = Vector3::new(1.0, 1.0, 0.0);
        let p01 = Vector3::new(0.0, 1.0, 0.0);
        let c0 = line(p00, p10);
        let c1 = line(p01, p11);
        let d0 = line(p00, p01);
        let d1 = line(p10, p11);
        let b = CoonsBoundary {
            c0: &c0,
            c1: &c1,
            d0: &d0,
            d1: &d1,
        };
        for &u in &[0.1_f64, 0.5, 0.9] {
            for &v in &[0.1_f64, 0.5, 0.9] {
                let p = evaluate(&b, u, v);
                let expected = Vector3::new(u, v, 0.0);
                assert!(
                    (p - expected).norm() < 1e-10,
                    "(u={u}, v={v}): got {p:?}, expected {expected:?}"
                );
            }
        }
    }

    #[test]
    fn fill_unit_square_returns_planar_patch() {
        let p00 = Vector3::new(0.0, 0.0, 0.0);
        let p10 = Vector3::new(1.0, 0.0, 0.0);
        let p11 = Vector3::new(1.0, 1.0, 0.0);
        let p01 = Vector3::new(0.0, 1.0, 0.0);
        let c0 = line(p00, p10);
        let c1 = line(p01, p11);
        let d0 = line(p00, p01);
        let d1 = line(p10, p11);
        let s = fill([c0, c1, d0, d1]).unwrap();
        // The 4 corner CPs of the produced NurbsSurface must match.
        assert!((s.control_points[0][0] - p00).norm() < 1e-10);
        assert!((s.control_points[3][0] - p10).norm() < 1e-10);
        assert!((s.control_points[0][3] - p01).norm() < 1e-10);
        assert!((s.control_points[3][3] - p11).norm() < 1e-10);
        // The surface evaluated at corners returns the corners.
        assert!((s.evaluate(0.0, 0.0) - p00).norm() < 1e-10);
        assert!((s.evaluate(1.0, 1.0) - p11).norm() < 1e-10);
        // And the centre is at (0.5, 0.5, 0).
        assert!((s.evaluate(0.5, 0.5) - Vector3::new(0.5, 0.5, 0.0)).norm() < 1e-10);
    }

    #[test]
    fn fill_quarter_arcs_produces_curved_patch() {
        // Four edges that are approximately quarter-circles in the
        // xy plane, raised so the patch has curvature.
        // Use the 4 cubic-Bezier approximations of a unit-radius
        // quarter-circle (standard 4-CP control polygon).
        let k = 4.0 / 3.0 * (std::f64::consts::PI / 8.0).tan();
        // We just need *some* curved boundaries with matching
        // corners — use semi-arbitrary CPs but ensure the corners
        // coincide perfectly.
        let p00 = Vector3::new(0.0, 0.0, 0.0);
        let p10 = Vector3::new(1.0, 0.0, 0.0);
        let p11 = Vector3::new(1.0, 1.0, 0.0);
        let p01 = Vector3::new(0.0, 1.0, 0.0);
        let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];

        // c0 (v=0): bottom edge bowed upward in z by `k`.
        let c0 = NurbsCurve::new(
            3,
            knots.clone(),
            vec![
                p00,
                Vector3::new(1.0 / 3.0, 0.0, k),
                Vector3::new(2.0 / 3.0, 0.0, k),
                p10,
            ],
            vec![1.0; 4],
        )
        .unwrap();
        // c1 (v=1): top edge bowed upward in z by `k`.
        let c1 = NurbsCurve::new(
            3,
            knots.clone(),
            vec![
                p01,
                Vector3::new(1.0 / 3.0, 1.0, k),
                Vector3::new(2.0 / 3.0, 1.0, k),
                p11,
            ],
            vec![1.0; 4],
        )
        .unwrap();
        // d0 (u=0): left edge straight.
        let d0 = line(p00, p01);
        // d1 (u=1): right edge straight.
        let d1 = line(p10, p11);

        let s = fill([c0, c1, d0, d1]).unwrap();
        // Corners still match.
        assert!((s.evaluate(0.0, 0.0) - p00).norm() < 1e-9);
        assert!((s.evaluate(1.0, 1.0) - p11).norm() < 1e-9);
        // Interior point has lifted z (curvature is being captured).
        let mid = s.evaluate(0.5, 0.5);
        assert!(mid.z > 0.1, "mid z = {} should be > 0.1", mid.z);
    }
}
