//! Phase 156 — `GeomLProp_CLProps::Curvature` — curve curvature at
//! parameter `t`.
//!
//! ## What OCCT does
//!
//! For a curve `C(t)`, curvature is `|C' × C''| / |C'|^3`. OCCT's
//! `GeomLProp_CLProps` packages this together with the tangent and
//! normal vectors as a one-shot "curve-local-properties" probe at a
//! given parameter.
//!
//! Used for offsetting (the offset radius must exceed the curve's
//! minimum radius of curvature; offsetting through the curvature
//! limit produces cusps), for variable-radius fillets, and for
//! generating curvature plots in surfacing UIs.
//!
//! ## v1 status
//!
//! **Honest v1.** Computes both derivatives via
//! `NurbsCurve::derivative` and assembles the cross-product formula
//! directly. Surface-curvature (`GeomLProp_SLProps`) is the
//! parametric analog and is a separate v1.5 — needs analytic second
//! partials of the surface, which we haven't wired through yet.

use nalgebra::Vector3;
use valenx_surface::NurbsCurve;

use crate::error::OcctAdvancedError;

/// Curvature report from [`geom_lib_curvature_at_point`].
#[derive(Clone, Debug, PartialEq)]
pub struct CurvatureReport {
    /// Scalar curvature `κ = |C' × C''| / |C'|^3`. Units are
    /// inverse length.
    pub curvature: f64,
    /// Radius of curvature `1/κ`, set to `f64::INFINITY` when `κ` is
    /// below `1e-12` (locally straight).
    pub radius: f64,
    /// First derivative (tangent direction, not normalised — caller
    /// can normalise via [`crate::geom_lib_tangent_at_point()`] if
    /// desired).
    pub first_derivative: Vector3<f64>,
    /// Second derivative (curvature direction).
    pub second_derivative: Vector3<f64>,
}

/// Compute curvature at parameter `t` on `curve`.
///
/// # Errors
///
/// - [`OcctAdvancedError::BadInput`] when `t` falls outside the
///   curve's valid parameter range.
/// - [`OcctAdvancedError::Defect`] when the first derivative
///   vanishes (cusp / degenerate — curvature undefined).
pub fn geom_lib_curvature_at_point(
    curve: &NurbsCurve,
    t: f64,
) -> Result<CurvatureReport, OcctAdvancedError> {
    let (t_min, t_max) = curve.parameter_range();
    if !t.is_finite() {
        return Err(OcctAdvancedError::bad_input(
            "t",
            "parameter must be finite",
        ));
    }
    if t < t_min || t > t_max {
        return Err(OcctAdvancedError::bad_input(
            "t",
            format!("{t} outside [{t_min}, {t_max}]"),
        ));
    }
    let d1 = curve.derivative(t, 1);
    let d2 = curve.derivative(t, 2);
    let d1_norm = d1.norm();
    if d1_norm < 1e-12 {
        return Err(OcctAdvancedError::defect(
            format!("t={t}"),
            "first derivative vanishes; curvature undefined",
        ));
    }
    let cross = d1.cross(&d2);
    let kappa = cross.norm() / d1_norm.powi(3);
    let radius = if kappa > 1e-12 { 1.0 / kappa } else { f64::INFINITY };
    Ok(CurvatureReport {
        curvature: kappa,
        radius,
        first_derivative: d1,
        second_derivative: d2,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn line(a: Vector3<f64>, b: Vector3<f64>) -> NurbsCurve {
        let p1 = a + (b - a) / 3.0;
        let p2 = a + 2.0 * (b - a) / 3.0;
        NurbsCurve::new(
            3,
            vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            vec![a, p1, p2, b],
            vec![1.0; 4],
        )
        .unwrap()
    }

    #[test]
    fn rejects_out_of_range() {
        let c = line(Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0));
        let err = geom_lib_curvature_at_point(&c, 2.0).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn straight_line_has_zero_curvature() {
        let c = line(Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0));
        let r = geom_lib_curvature_at_point(&c, 0.5).unwrap();
        // A straight cubic Bezier (collinear control points) has
        // curvature ≈ 0 everywhere — within finite-difference noise.
        assert!(r.curvature < 1e-3, "expected ~0, got {}", r.curvature);
        assert!(
            r.radius.is_infinite() || r.radius > 1e3,
            "expected ~∞, got {}",
            r.radius
        );
    }

    #[test]
    fn cusp_curve_flagged() {
        // Construct a degenerate curve: all four control points are
        // (0,0,0) — derivative collapses to zero everywhere.
        let cps = vec![Vector3::zeros(); 4];
        let c = NurbsCurve::new(
            3,
            vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            cps,
            vec![1.0; 4],
        )
        .unwrap();
        let err = geom_lib_curvature_at_point(&c, 0.5).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.defect");
    }
}
