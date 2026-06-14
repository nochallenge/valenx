//! Frenet-Serret frame and torsion of a NURBS curve.
//!
//! For a regular 3D curve `C(u)` the Frenet-Serret apparatus attaches
//! a right-handed orthonormal frame to every non-singular point. The
//! unit tangent `T = C'/|C'|` gives the direction of travel; the
//! principal normal `N` is the unit vector toward the centre of
//! curvature, lying in the osculating plane and orthogonal to `T`; and
//! the binormal `B = T x N` is normal to the osculating plane.
//!
//! Two scalar invariants accompany the triad: the curvature `kappa`
//! (how fast `T` turns) and the torsion `tau` (how fast `B` turns,
//! i.e. how fast the curve twists out of its osculating plane).
//!
//! This module derives the frame from the curve's first three
//! derivatives `d1 = C'`, `d2 = C''`, `d3 = C'''` using the standard
//! cross-product identities, which are valid for an arbitrary
//! (not necessarily arc-length) parameterisation:
//!
//! ```text
//! T     = d1 / |d1|
//! B     = (d1 x d2) / |d1 x d2|
//! N     = B x T
//! kappa = |d1 x d2| / |d1|^3
//! tau   = ((d1 x d2) . d3) / |d1 x d2|^2
//! ```
//!
//! The frame is undefined at two kinds of degenerate point and
//! [`frame_at`] returns an error there. The first is a singular point
//! where `|d1| ~= 0` (the tangent vanishes). The second is an
//! inflection or locally-straight point where `|d1 x d2| ~= 0`: there
//! `d1` and `d2` are parallel, so the osculating plane — and hence `N`,
//! `B`, and `tau` — is not determined.
//!
//! ## Scope caveat
//!
//! This is a **research / preliminary-design-grade** differential-
//! geometry utility. The derivatives come from [`NurbsCurve::derivative`],
//! which uses central finite differences, so the tangent / normal /
//! binormal and curvature (first/second derivative) are accurate to
//! roughly 6-7 significant figures, while the torsion (which needs the
//! third derivative, formed by nested differencing) is markedly noisier
//! near degeneracies and should be read as indicative rather than
//! certified. It is **not** a substitute for the validated geometric-
//! modelling kernels in CATIA, NX, Ansys, or Adams; do not use it for
//! certification, contact analysis, or any safety-critical workflow.

use nalgebra::Vector3;

use crate::error::SurfaceError;
use crate::nurbs_curve::NurbsCurve;

/// Absolute tolerance below which the tangent magnitude `|C'|` is
/// treated as zero (a singular point). The first derivative is the
/// cleanest finite-difference term, so a tight absolute floor is safe.
const SPEED_EPS: f64 = 1e-9;

/// Relative tolerance for the osculating-plane test. The point is
/// treated as an inflection / locally-straight point when
/// `|C' x C''| < CROSS_EPS_REL * |C'|`, i.e. the component of `C''`
/// perpendicular to `C'` is negligible compared with the speed. Using a
/// ratio makes the test scale-invariant and keeps it well above the
/// nested-finite-difference noise of [`NurbsCurve::derivative`] (which,
/// for an exactly straight segment, leaves only ~1e-10-relative residue)
/// while staying far below the O(1) ratio of any genuinely curved span.
const CROSS_EPS_REL: f64 = 1e-6;

/// The Frenet-Serret frame of a curve at one parameter value.
///
/// `tangent`, `normal`, and `binormal` form a right-handed orthonormal
/// triad (`binormal = tangent x normal` up to floating-point error).
/// All three are unit vectors. `curvature` is non-negative; `torsion`
/// may be negative (a left-handed twist) or positive (right-handed).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrenetFrame {
    /// Curve parameter `u` at which the frame was evaluated.
    pub u: f64,
    /// Point on the curve, `C(u)`.
    pub point: Vector3<f64>,
    /// Unit tangent `T = C'/|C'|`.
    pub tangent: Vector3<f64>,
    /// Unit principal normal `N` (points toward the centre of
    /// curvature; orthogonal to `tangent`).
    pub normal: Vector3<f64>,
    /// Unit binormal `B = tangent x normal` (normal to the osculating
    /// plane).
    pub binormal: Vector3<f64>,
    /// Curvature `kappa = |C' x C''| / |C'|^3` (>= 0). The reciprocal
    /// `1/kappa` is the radius of the osculating circle.
    pub curvature: f64,
    /// Torsion `tau = ((C' x C'') . C''') / |C' x C''|^2`. Zero for a
    /// planar curve; its sign encodes the handedness of the twist.
    pub torsion: f64,
}

/// Compute the Frenet-Serret frame and torsion of `curve` at parameter
/// `u`.
///
/// `u` is clamped-evaluated through [`NurbsCurve::derivative`]; values
/// outside the curve's [`NurbsCurve::parameter_range`] yield
/// [`SurfaceError::EvaluationOutOfRange`].
///
/// # Errors
///
/// Returns [`SurfaceError::EvaluationOutOfRange`] if `u` is outside the
/// valid knot range, and [`SurfaceError::IntersectionFailed`] (the
/// crate's `"geometry"`-category variant) if the point is degenerate —
/// either a singular point (`|C'| ~= 0`) or an inflection / locally-
/// straight point (`|C' x C''| ~= 0`) where the principal normal,
/// binormal, and torsion are not defined.
pub fn frame_at(curve: &NurbsCurve, u: f64) -> Result<FrenetFrame, SurfaceError> {
    let (u_min, u_max) = curve.parameter_range();
    // Allow a hair of slack so callers evaluating exactly at the
    // endpoints (subject to rounding) are not rejected.
    let slack = (u_max - u_min).abs() * 1e-9 + 1e-12;
    if u < u_min - slack || u > u_max + slack {
        return Err(SurfaceError::EvaluationOutOfRange { u });
    }

    let point = curve.evaluate(u);
    let d1 = curve.derivative(u, 1);
    let d2 = curve.derivative(u, 2);
    let d3 = curve.derivative(u, 3);

    let speed = d1.norm();
    if speed < SPEED_EPS {
        return Err(SurfaceError::IntersectionFailed(format!(
            "Frenet frame undefined at u={u}: singular point (|C'|={speed:.3e} ~= 0)"
        )));
    }

    let cross = d1.cross(&d2);
    let cross_norm = cross.norm();
    if cross_norm < CROSS_EPS_REL * speed {
        return Err(SurfaceError::IntersectionFailed(format!(
            "Frenet frame undefined at u={u}: inflection / locally-straight \
             point (|C' x C''|={cross_norm:.3e} negligible vs |C'|={speed:.3e}, \
             normal & torsion undetermined)"
        )));
    }

    let tangent = d1 / speed;
    let binormal = cross / cross_norm;
    // N = B x T is already unit-length (B _|_ T, both unit) and points
    // toward the centre of curvature.
    let normal = binormal.cross(&tangent);

    let curvature = cross_norm / speed.powi(3);
    let torsion = cross.dot(&d3) / (cross_norm * cross_norm);

    Ok(FrenetFrame {
        u,
        point,
        tangent,
        normal,
        binormal,
        curvature,
        torsion,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Standard cubic Bezier knot vector — four control points, one
    /// segment. `[0,0,0,0,1,1,1,1]`.
    fn bezier_knots() -> Vec<f64> {
        vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
    }

    fn cubic_bezier(cps: [Vector3<f64>; 4]) -> NurbsCurve {
        NurbsCurve::new(3, bezier_knots(), cps.to_vec(), vec![1.0; 4]).unwrap()
    }

    /// Assert a triad is orthonormal: unit norms and zero pairwise dot
    /// products, both to a loose finite-difference-aware tolerance.
    fn assert_orthonormal(f: &FrenetFrame) {
        let tol = 1e-6;
        assert!(
            (f.tangent.norm() - 1.0).abs() < tol,
            "|T|={}",
            f.tangent.norm()
        );
        assert!(
            (f.normal.norm() - 1.0).abs() < tol,
            "|N|={}",
            f.normal.norm()
        );
        assert!(
            (f.binormal.norm() - 1.0).abs() < tol,
            "|B|={}",
            f.binormal.norm()
        );
        assert!(
            f.tangent.dot(&f.normal).abs() < tol,
            "T.N={}",
            f.tangent.dot(&f.normal)
        );
        assert!(
            f.tangent.dot(&f.binormal).abs() < tol,
            "T.B={}",
            f.tangent.dot(&f.binormal)
        );
        assert!(
            f.normal.dot(&f.binormal).abs() < tol,
            "N.B={}",
            f.normal.dot(&f.binormal)
        );
        // Right-handedness: B should equal T x N.
        let txn = f.tangent.cross(&f.normal);
        assert!((txn - f.binormal).norm() < tol, "B != T x N: {txn:?}");
    }

    #[test]
    fn planar_curve_has_zero_torsion_everywhere() {
        // A non-degenerate cubic whose control points all lie in the
        // z = 0 plane: the whole curve is planar, so torsion == 0 at
        // every interior parameter (the curve never leaves its
        // osculating plane). Control points are arranged so the curve
        // genuinely bends (non-zero curvature) and never inflects on
        // the sampled interior.
        let c = cubic_bezier([
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 2.0, 0.0),
            Vector3::new(3.0, 2.0, 0.0),
            Vector3::new(4.0, 0.0, 0.0),
        ]);
        for &u in &[0.15_f64, 0.3, 0.5, 0.7, 0.85] {
            let f = frame_at(&c, u).unwrap();
            // Third-derivative finite differencing is the noisy term;
            // a planar curve nonetheless pins torsion to ~0.
            assert!(
                f.torsion.abs() < 1e-4,
                "u={u}: torsion={} should be ~0 for a planar curve",
                f.torsion
            );
            // For a planar curve the binormal is the constant plane
            // normal (+z here, up to sign).
            assert!(
                f.binormal.x.abs() < 1e-6 && f.binormal.y.abs() < 1e-6,
                "u={u}: binormal {:?} should be the +/-z plane normal",
                f.binormal
            );
        }
    }

    #[test]
    fn frame_is_orthonormal_on_planar_curve() {
        let c = cubic_bezier([
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 2.0, 0.0),
            Vector3::new(3.0, 2.0, 0.0),
            Vector3::new(4.0, 0.0, 0.0),
        ]);
        for &u in &[0.2_f64, 0.5, 0.8] {
            let f = frame_at(&c, u).unwrap();
            assert_orthonormal(&f);
        }
    }

    #[test]
    fn frame_is_orthonormal_on_nonplanar_curve() {
        // A genuinely 3D (twisting) cubic: control points span all
        // three axes so the curve has non-zero torsion. Orthonormality
        // of the triad must still hold exactly.
        let c = cubic_bezier([
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 2.0, 0.0),
            Vector3::new(2.0, 0.0, 3.0),
            Vector3::new(4.0, 1.0, 1.0),
        ]);
        for &u in &[0.25_f64, 0.5, 0.75] {
            let f = frame_at(&c, u).unwrap();
            assert_orthonormal(&f);
        }
    }

    #[test]
    fn point_matches_curve_evaluation() {
        let c = cubic_bezier([
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 2.0, 0.0),
            Vector3::new(3.0, 2.0, 0.0),
            Vector3::new(4.0, 0.0, 0.0),
        ]);
        let u = 0.4;
        let f = frame_at(&c, u).unwrap();
        assert!((f.point - c.evaluate(u)).norm() < 1e-12);
        // Tangent agrees in direction with the raw first derivative.
        let d1 = c.derivative(u, 1);
        let tdir = d1 / d1.norm();
        assert!((f.tangent - tdir).norm() < 1e-9);
    }

    #[test]
    fn curvature_of_circular_arc_matches_radius() {
        // A quarter circle of radius R in the xy-plane modelled as a
        // rational quadratic Bezier (a NURBS exactly represents conics)
        // has constant curvature kappa = 1/R. Control points:
        //   P0=(R,0), P1=(R,R), P2=(0,R), weights (1, 1/sqrt2, 1),
        //   knots [0,0,0,1,1,1].
        let r = 2.0_f64;
        let w = 1.0 / 2.0_f64.sqrt();
        let c = NurbsCurve::new(
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![
                Vector3::new(r, 0.0, 0.0),
                Vector3::new(r, r, 0.0),
                Vector3::new(0.0, r, 0.0),
            ],
            vec![1.0, w, 1.0],
        )
        .unwrap();
        // Check curvature at the symmetric midpoint, where finite-
        // difference noise is smallest.
        let f = frame_at(&c, 0.5).unwrap();
        let expected = 1.0 / r;
        assert!(
            (f.curvature - expected).abs() < 1e-3,
            "kappa={} expected {expected} (1/R)",
            f.curvature
        );
        // A planar arc still has zero torsion.
        assert!(f.torsion.abs() < 1e-3, "arc torsion={}", f.torsion);
        assert_orthonormal(&f);
    }

    #[test]
    fn straight_segment_is_rejected_as_degenerate() {
        // Collinear control points => |C' x C''| == 0 everywhere, so
        // the principal normal / binormal / torsion are undefined and
        // frame_at must report a geometry error rather than return NaNs.
        let c = cubic_bezier([
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
        ]);
        let err = frame_at(&c, 0.5).unwrap_err();
        assert_eq!(err.code(), "surface.intersection_failed");
        assert_eq!(err.category(), "geometry");
    }

    #[test]
    fn out_of_range_parameter_is_rejected() {
        let c = cubic_bezier([
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 2.0, 0.0),
            Vector3::new(3.0, 2.0, 0.0),
            Vector3::new(4.0, 0.0, 0.0),
        ]);
        let err = frame_at(&c, 2.0).unwrap_err();
        assert_eq!(err.code(), "surface.evaluation_out_of_range");
    }
}
