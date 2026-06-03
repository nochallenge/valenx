//! Phase 153 — `GeomAPI_ProjectPointOnCurve` — find the parameter on
//! a curve closest to a query point.
//!
//! ## What OCCT does
//!
//! `GeomAPI_ProjectPointOnCurve(P, curve, u_min, u_max)` solves the
//! orthogonal projection problem: find `u*` such that the line from
//! `P` to `curve.evaluate(u*)` is perpendicular to the curve tangent
//! at `u*`. Returns all extrema (multiple if the curve self-
//! intersects or has multiple equally-distant minima).
//!
//! OCCT solves this by:
//! 1. Sample the curve at `n_samples` parameter values, find the
//!    sample with minimum distance to `P`.
//! 2. Apply Newton's method from that seed to refine `u*` to the
//!    true minimum (`f(u) = (P - C(u)) · C'(u)`, find `f(u) = 0`).
//! 3. Optionally repeat from secondary minima for additional extrema.
//!
//! ## v1 status
//!
//! **Honest v1.** Implements the OCCT algorithm: dense sampling +
//! Newton refinement. Returns the single closest point per call —
//! the multi-extrema variant is Phase 153.5 (needs a peak-detection
//! pass over the sample distances).

use nalgebra::Vector3;
use valenx_surface::NurbsCurve;

use crate::error::OcctAdvancedError;

/// Result of the project-point-on-curve op.
#[derive(Clone, Debug, PartialEq)]
pub struct ProjectionResult {
    /// Parameter `u*` on the curve.
    pub u: f64,
    /// Foot point: `curve.evaluate(u*)`.
    pub foot: Vector3<f64>,
    /// Euclidean distance from query point to foot.
    pub distance: f64,
}

/// Default sample count for the seeding pass. 64 is a good balance
/// between speed and seed quality for cubic NURBS — tighten for
/// higher-degree curves with many wiggles.
pub const DEFAULT_SAMPLES: usize = 64;
/// Default Newton iteration count. 8 is plenty for cubic Newton
/// convergence from a good seed.
pub const DEFAULT_NEWTON_ITERS: usize = 8;

/// Project `point` onto `curve` and return the closest point.
///
/// # Errors
///
/// - [`OcctAdvancedError::BadInput`] for `n_samples < 2` or
///   non-finite `point` coordinates.
pub fn geom_lib_intersect_point_curve(
    point: Vector3<f64>,
    curve: &NurbsCurve,
    n_samples: usize,
) -> Result<ProjectionResult, OcctAdvancedError> {
    if n_samples < 2 {
        return Err(OcctAdvancedError::bad_input(
            "n_samples",
            "need ≥2 seed samples",
        ));
    }
    if !point.iter().all(|c| c.is_finite()) {
        return Err(OcctAdvancedError::bad_input(
            "point",
            "coordinates must be finite",
        ));
    }

    let (u_min, u_max) = curve.parameter_range();
    if u_max <= u_min {
        return Err(OcctAdvancedError::bad_input(
            "curve",
            format!("parameter range [{u_min},{u_max}] is empty"),
        ));
    }

    // Pass 1: sample and find best seed.
    let mut best_u = u_min;
    let mut best_d2 = f64::INFINITY;
    let mut best_foot = curve.evaluate(u_min);
    for i in 0..n_samples {
        let t = u_min + (u_max - u_min) * (i as f64 / (n_samples - 1) as f64);
        let p = curve.evaluate(t);
        let d2 = (p - point).norm_squared();
        if d2 < best_d2 {
            best_d2 = d2;
            best_u = t;
            best_foot = p;
        }
    }

    // Pass 2: Newton refinement. f(u) = (C(u) - P) · C'(u);
    // f'(u) = C'(u) · C'(u) + (C(u) - P) · C''(u). Step: u -= f/f'.
    let mut u = best_u;
    for _ in 0..DEFAULT_NEWTON_ITERS {
        let c = curve.evaluate(u);
        let c_prime = curve.derivative(u, 1);
        let c_dprime = curve.derivative(u, 2);
        let diff = c - point;
        let f = diff.dot(&c_prime);
        let fp = c_prime.dot(&c_prime) + diff.dot(&c_dprime);
        if fp.abs() < 1e-12 {
            // Singular Jacobian — bail with current best.
            break;
        }
        let step = f / fp;
        u -= step;
        u = u.clamp(u_min, u_max);
        if step.abs() < 1e-10 {
            break;
        }
    }

    let foot = curve.evaluate(u);
    let distance = (foot - point).norm();
    // Take the better of the refined point and the seed (Newton can
    // diverge near the endpoints for degenerate seeds).
    if (foot - point).norm_squared() <= best_d2 {
        Ok(ProjectionResult { u, foot, distance })
    } else {
        Ok(ProjectionResult {
            u: best_u,
            foot: best_foot,
            distance: best_d2.sqrt(),
        })
    }
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
        let err = geom_lib_intersect_point_curve(Vector3::zeros(), &c, 1).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_non_finite_point() {
        let c = line(Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0));
        let err = geom_lib_intersect_point_curve(
            Vector3::new(f64::NAN, 0.0, 0.0),
            &c,
            DEFAULT_SAMPLES,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn point_on_curve_projects_to_itself() {
        // Point at (0.5, 0, 0) lies on the unit-x straight line; should
        // project to u ≈ 0.5 with near-zero distance.
        let c = line(Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0));
        let p = Vector3::new(0.5, 0.0, 0.0);
        let r = geom_lib_intersect_point_curve(p, &c, DEFAULT_SAMPLES).unwrap();
        assert!(r.distance < 1e-9, "expected near-zero, got {}", r.distance);
        assert!((r.u - 0.5).abs() < 1e-4, "expected u≈0.5, got {}", r.u);
    }

    #[test]
    fn perpendicular_point_projects_orthogonally() {
        // Point at (0.3, 1.0, 0): closest point on the x-axis line
        // should be (0.3, 0, 0) at u = 0.3, distance 1.0.
        let c = line(Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0));
        let p = Vector3::new(0.3, 1.0, 0.0);
        let r = geom_lib_intersect_point_curve(p, &c, DEFAULT_SAMPLES).unwrap();
        assert!(
            (r.distance - 1.0).abs() < 1e-6,
            "expected ~1.0, got {}",
            r.distance
        );
        assert!((r.u - 0.3).abs() < 1e-3, "expected u≈0.3, got {}", r.u);
    }

    #[test]
    fn endpoint_projection() {
        // Point well past the curve end: closest is curve(1.0).
        let c = line(Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0));
        let p = Vector3::new(2.0, 0.0, 0.0);
        let r = geom_lib_intersect_point_curve(p, &c, DEFAULT_SAMPLES).unwrap();
        assert!((r.u - 1.0).abs() < 1e-6, "expected u≈1.0, got {}", r.u);
    }
}
