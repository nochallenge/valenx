//! Phase 157 — `GeomLib::BuildCurveFromEvolute` — compute the evolute
//! (curve of centres of curvature) of an input curve.
//!
//! ## What OCCT does
//!
//! For a curve `C(t)` with curvature `κ(t)` and principal normal
//! direction `N(t)`, the evolute is the locus of curvature centres:
//!
//! ```text
//!   E(t) = C(t) + (1 / κ(t)) * N(t)
//! ```
//!
//! Geometric interpretation: it's the curve traced by the centres of
//! the osculating circles. Used in classical differential geometry,
//! involute-gear design (the gear-tooth profile is the involute of
//! the base circle — its inverse is the evolute), and offsetting
//! (the parallel-curve cusps when the offset radius hits the
//! evolute).
//!
//! ## v1 status
//!
//! **Honest v1.** Samples the input curve at `n_samples` parameter
//! values, computes `E(t)` at each, and returns the resulting
//! polyline as a sequence of 3D points. The closed-form evolute as a
//! `NurbsCurve` is Phase 157.5 — depends on a "fit-NURBS-to-polyline"
//! op (which valenx-surface's `approx_curve_fit` provides; the
//! integration is just plumbing).

use nalgebra::Vector3;
use valenx_surface::NurbsCurve;

use crate::error::OcctAdvancedError;
use crate::geom_lib_curvature_at_point::geom_lib_curvature_at_point;

/// Evolute polyline returned by [`geom_lib_evolute_curve`].
#[derive(Clone, Debug, PartialEq)]
pub struct EvolutePoint {
    /// Parameter on the source curve.
    pub t: f64,
    /// Centre of curvature in 3D.
    pub center: Vector3<f64>,
    /// Radius of curvature at this parameter (`1 / κ`).
    pub radius: f64,
}

/// Default sample count. 64 gives a smooth evolute polyline for
/// cubic NURBS at viewport-scale tolerance.
pub const DEFAULT_SAMPLES: usize = 64;

/// Compute the evolute of `curve` sampled at `n_samples` parameters.
///
/// Skips samples where the local curvature is below `min_curvature`
/// (those would project to infinity along the normal — typically
/// straight stretches). Returns the surviving points in parameter
/// order.
///
/// # Errors
///
/// - [`OcctAdvancedError::BadInput`] for `n_samples < 2` or
///   non-positive `min_curvature`.
pub fn geom_lib_evolute_curve(
    curve: &NurbsCurve,
    n_samples: usize,
    min_curvature: f64,
) -> Result<Vec<EvolutePoint>, OcctAdvancedError> {
    if n_samples < 2 {
        return Err(OcctAdvancedError::bad_input(
            "n_samples",
            "need ≥2 sample points",
        ));
    }
    if !min_curvature.is_finite() || min_curvature <= 0.0 {
        return Err(OcctAdvancedError::bad_input(
            "min_curvature",
            "must be positive finite",
        ));
    }

    let (t_min, t_max) = curve.parameter_range();
    let mut out = Vec::new();
    for i in 0..n_samples {
        let t = t_min + (t_max - t_min) * (i as f64 / (n_samples - 1) as f64);
        // Curvature probe may legitimately Defect at endpoints with
        // zero derivative — skip those samples rather than propagate
        // the failure.
        let cur = match geom_lib_curvature_at_point(curve, t) {
            Ok(c) => c,
            Err(_) => continue,
        };
        if cur.curvature < min_curvature {
            continue;
        }
        let pt = curve.evaluate(t);
        // Principal normal direction: ((C' × C'') × C') normalized,
        // which is the in-plane component of C'' perpendicular to C'.
        let d1 = cur.first_derivative;
        let d2 = cur.second_derivative;
        let binormal = d1.cross(&d2);
        let normal_direction = binormal.cross(&d1);
        let nd_norm = normal_direction.norm();
        if nd_norm < 1e-12 {
            continue;
        }
        let n = normal_direction / nd_norm;
        let center = pt + n * cur.radius;
        out.push(EvolutePoint {
            t,
            center,
            radius: cur.radius,
        });
    }
    Ok(out)
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
    fn rejects_one_sample() {
        let c = line(Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0));
        let err = geom_lib_evolute_curve(&c, 1, 1e-6).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_zero_min_curvature() {
        let c = line(Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0));
        let err = geom_lib_evolute_curve(&c, DEFAULT_SAMPLES, 0.0).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn straight_line_evolute_is_empty() {
        // A straight line has zero curvature everywhere — no evolute
        // points survive the min-curvature filter.
        let c = line(Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0));
        let pts = geom_lib_evolute_curve(&c, DEFAULT_SAMPLES, 1e-6).unwrap();
        assert!(pts.is_empty(), "expected empty, got {} points", pts.len());
    }

    #[test]
    fn curved_bezier_yields_points() {
        // Cubic Bezier with non-collinear control points — has real
        // curvature in the interior.
        let cps = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.5, 1.0, 0.0),
            Vector3::new(1.5, 1.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
        ];
        let c = NurbsCurve::new(
            3,
            vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            cps,
            vec![1.0; 4],
        )
        .unwrap();
        let pts = geom_lib_evolute_curve(&c, DEFAULT_SAMPLES, 1e-6).unwrap();
        assert!(!pts.is_empty(), "expected ≥1 evolute point");
        // Every survivor has finite radius.
        for p in &pts {
            assert!(p.radius.is_finite(), "radius should be finite, got {}", p.radius);
        }
    }
}
