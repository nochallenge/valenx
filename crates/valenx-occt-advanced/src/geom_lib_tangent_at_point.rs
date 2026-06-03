//! Phase 155 — `GeomLib::Tangent` for curves — return the unit
//! tangent at parameter `t`.
//!
//! ## What OCCT does
//!
//! `Geom_Curve::D1(t, P, dT)` returns the curve point and its first
//! derivative `dT`. The unit tangent is `dT / |dT|`. Most callers
//! want the unit form (frame transport, draw-arrow visualisations,
//! offset-curve construction) — that's what this op returns.
//!
//! ## v1 status
//!
//! **Honest v1.** Wraps `NurbsCurve::derivative(t, 1)` with
//! validation, normalization, and a `Defect` for the degenerate
//! case where the first derivative vanishes (cusp / collapsed
//! segment).

use nalgebra::Vector3;
use valenx_surface::NurbsCurve;

use crate::error::OcctAdvancedError;

/// Compute the unit tangent at parameter `t` on `curve`.
///
/// # Errors
///
/// - [`OcctAdvancedError::BadInput`] when `t` falls outside the
///   curve's valid parameter range.
/// - [`OcctAdvancedError::Defect`] when the first derivative is
///   near-zero (cusp / degenerate).
pub fn geom_lib_tangent_at_point(
    curve: &NurbsCurve,
    t: f64,
) -> Result<Vector3<f64>, OcctAdvancedError> {
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
    let d = curve.derivative(t, 1);
    let d_norm = d.norm();
    if d_norm < 1e-12 {
        return Err(OcctAdvancedError::defect(
            format!("t={t}"),
            format!("first-derivative magnitude {d_norm:.3e} (cusp / degenerate)"),
        ));
    }
    Ok(d / d_norm)
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
        let err = geom_lib_tangent_at_point(&c, -0.1).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_non_finite() {
        let c = line(Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0));
        let err = geom_lib_tangent_at_point(&c, f64::INFINITY).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn straight_line_tangent_is_unit_direction() {
        // x-axis line from (0,0,0) to (1,0,0) — tangent is +x.
        let c = line(Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0));
        let t = geom_lib_tangent_at_point(&c, 0.5).unwrap();
        assert!((t.x - 1.0).abs() < 1e-6, "tx should be 1, got {}", t.x);
        assert!(t.y.abs() < 1e-6, "ty should be 0, got {}", t.y);
        assert!(t.z.abs() < 1e-6, "tz should be 0, got {}", t.z);
    }

    #[test]
    fn diagonal_line_tangent_is_unit_diagonal() {
        // Diagonal from (0,0,0) to (1,1,1) — tangent is (1,1,1)/sqrt(3).
        let c = line(Vector3::zeros(), Vector3::new(1.0, 1.0, 1.0));
        let t = geom_lib_tangent_at_point(&c, 0.5).unwrap();
        let expected = 1.0 / 3.0_f64.sqrt();
        assert!((t.x - expected).abs() < 1e-6);
        assert!((t.y - expected).abs() < 1e-6);
        assert!((t.z - expected).abs() < 1e-6);
        // Verify unit norm.
        assert!((t.norm() - 1.0).abs() < 1e-9);
    }
}
