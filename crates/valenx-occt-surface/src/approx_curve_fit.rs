//! Phase 78 — `Approx_Curve3d`: fit a B-spline curve to scattered 3D
//! points within a tolerance.
//!
//! ## What OCCT does
//!
//! `Approx_Curve3d` is OCCT's least-squares curve approximator. The
//! caller supplies an ordered sequence of `gp_Pnt`s plus the desired
//! degree and tolerance; the algorithm picks an optimal knot vector
//! and CP count, solves the LSQ system, and returns a
//! `Geom_BSplineCurve`. The most common downstream uses are:
//!
//! - Refit intersection polylines produced by `BRepAlgoAPI_Section`.
//! - Smooth out laser-scanner traces before sweeping them.
//! - Recover an analytic curve from a polyline stored in DXF/DWG.
//!
//! ## v1 status
//!
//! **Honest implementation.** `valenx-surface::fit::nurbs_curve_through_points`
//! solves the same LSQ problem. We delegate, returning the achieved
//! RMS error and verifying the tolerance budget the caller supplied.

use nalgebra::Vector3;
use valenx_surface::{fit, NurbsCurve};

use crate::error::OcctSurfaceError;

/// Output of a curve approximation.
#[derive(Clone, Debug)]
pub struct ApproxCurve {
    /// Fitted NURBS curve.
    pub curve: NurbsCurve,
    /// RMS error in the input points' metric.
    pub rms_error: f64,
}

/// Fit a NURBS curve through `points` using `n_cps` control points
/// and `degree`, capped at the caller-supplied `tolerance`.
///
/// # Errors
///
/// [`OcctSurfaceError::BadInput`] for malformed inputs or when the
/// achieved RMS exceeds `tolerance`;
/// [`OcctSurfaceError::TruckLimit`] when the underlying
/// `valenx-surface` fitter rejects the request.
pub fn approx_curve_fit(
    points: &[Vector3<f64>],
    degree: usize,
    n_cps: usize,
    tolerance: f64,
) -> Result<ApproxCurve, OcctSurfaceError> {
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(OcctSurfaceError::bad_input(
            "tolerance",
            format!("must be a positive finite number, got {tolerance}"),
        ));
    }
    let result = fit::nurbs_curve_through_points(points, degree, n_cps).map_err(|e| {
        OcctSurfaceError::TruckLimit(format!("fit::nurbs_curve_through_points: {e:?}"))
    })?;
    if result.rms_error > tolerance {
        return Err(OcctSurfaceError::bad_input(
            "tolerance",
            format!(
                "achieved RMS error {:.6} > requested tolerance {:.6}; try more control points",
                result.rms_error, tolerance
            ),
        ));
    }
    Ok(ApproxCurve {
        curve: result.curve,
        rms_error: result.rms_error,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fit_collinear_points_returns_zero_error() {
        // Five collinear points along the X axis — a degree-1 line
        // fits them exactly.
        let points: Vec<Vector3<f64>> = (0..5)
            .map(|i| Vector3::new(i as f64, 0.0, 0.0))
            .collect();
        let result = approx_curve_fit(&points, 1, 2, 1e-9)
            .expect("collinear points fit a line at machine precision");
        assert!(result.rms_error < 1e-9);
    }

    #[test]
    fn fit_rejects_bad_tolerance() {
        let points: Vec<Vector3<f64>> = (0..3)
            .map(|i| Vector3::new(i as f64, 0.0, 0.0))
            .collect();
        let err = approx_curve_fit(&points, 1, 2, 0.0).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }
}
