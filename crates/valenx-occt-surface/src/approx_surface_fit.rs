//! Phase 77 — `BRepApprox_Approx`: fit a B-spline surface to a
//! gridded point set within a tolerance.
//!
//! ## What OCCT does
//!
//! `BRepApprox_Approx` is the OCCT approximation algorithm used after
//! IGES/STEP import to refit imported surfaces with cleaner knot
//! vectors, or after `BRepFill_OffsetSurface` to convert the implicit
//! offset to an explicit B-spline. The caller supplies:
//!
//! - The data points (gridded `TColgp_Array2OfPnt` of size m×n).
//! - Target degrees in u and v (typically 3 each).
//! - Caller-supplied tolerance in model units — the fit iterates,
//!   adding knot spans until the maximum deviation falls below this
//!   value.
//!
//! Output is a `Geom_BSplineSurface` plus the achieved max deviation
//! and a per-knot diagnostic vector.
//!
//! ## v1 status
//!
//! **Honest implementation.** `valenx-surface::fit::nurbs_surface_through_grid`
//! provides a least-squares grid fit with a caller-chosen
//! `(n_cps_u, n_cps_v)`. We expose it through OCCT's calling
//! convention and convert the tolerance argument into the RMS-error
//! cap (callers get `BadInput` if the achieved error exceeds it).

use nalgebra::Vector3;
use valenx_surface::{fit, NurbsSurface};

use crate::error::OcctSurfaceError;

/// Output of an approximation fit.
#[derive(Clone, Debug)]
pub struct ApproxSurface {
    /// Fitted NURBS surface.
    pub surface: NurbsSurface,
    /// RMS error achieved across all input data points.
    pub rms_error: f64,
}

/// Fit a B-spline surface through a structured `(m x n)` grid of 3D
/// points, with the given degrees and target control-point counts.
///
/// `tolerance` is the maximum acceptable RMS error in model units.
/// If the fit's RMS exceeds it, the caller gets `BadInput` so they
/// can retry with more control points.
///
/// # Errors
///
/// [`OcctSurfaceError::BadInput`] for malformed inputs or when the
/// achieved RMS exceeds `tolerance`;
/// [`OcctSurfaceError::TruckLimit`] when the underlying
/// `valenx-surface` fitter rejects the request.
pub fn approx_surface_fit(
    points_uv: &[Vec<Vector3<f64>>],
    degree_u: usize,
    degree_v: usize,
    n_cps_u: usize,
    n_cps_v: usize,
    tolerance: f64,
) -> Result<ApproxSurface, OcctSurfaceError> {
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(OcctSurfaceError::bad_input(
            "tolerance",
            format!("must be a positive finite number, got {tolerance}"),
        ));
    }
    let result = fit::nurbs_surface_through_grid(
        points_uv, degree_u, degree_v, n_cps_u, n_cps_v,
    )
    .map_err(|e| OcctSurfaceError::TruckLimit(format!("fit::nurbs_surface_through_grid: {e:?}")))?;

    if result.rms_error > tolerance {
        return Err(OcctSurfaceError::bad_input(
            "tolerance",
            format!(
                "achieved RMS error {:.6} > requested tolerance {:.6}; try more control points",
                result.rms_error, tolerance
            ),
        ));
    }
    Ok(ApproxSurface {
        surface: result.surface,
        rms_error: result.rms_error,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flat_grid_3x3() -> Vec<Vec<Vector3<f64>>> {
        let mut grid = Vec::with_capacity(3);
        for i in 0..3 {
            let mut row = Vec::with_capacity(3);
            for j in 0..3 {
                row.push(Vector3::new(i as f64, j as f64, 0.0));
            }
            grid.push(row);
        }
        grid
    }

    #[test]
    fn fit_flat_3x3_grid_returns_low_error() {
        let grid = flat_grid_3x3();
        let result = approx_surface_fit(&grid, 2, 2, 3, 3, 1e-6)
            .expect("flat grid should fit a degree-(2,2) bicubic at machine precision");
        assert!(result.rms_error < 1e-6);
    }

    #[test]
    fn fit_rejects_bad_tolerance() {
        let grid = flat_grid_3x3();
        let err = approx_surface_fit(&grid, 2, 2, 3, 3, -0.1).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }
}
