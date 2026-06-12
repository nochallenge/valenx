//! NURBS curve — Non-Uniform Rational B-Spline.
//!
//! A NURBS curve is defined by:
//! - **degree** `p` (e.g. 3 for a cubic),
//! - **knot vector** `U` of length `m = n + p + 2` where `n + 1` is
//!   the number of control points,
//! - **control points** `P_i` in 3D, and
//! - **weights** `w_i` (one per control point; `1.0` reduces to a
//!   plain B-spline).
//!
//! Evaluation uses standard Cox-de Boor basis-function recursion +
//! the rational denominator. See [`NurbsCurve::evaluate`].

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::SurfaceError;

/// A 3D NURBS curve.
///
/// Built via [`NurbsCurve::new`] (validated) or
/// [`NurbsCurve::new_unchecked`] (no validation — only call when the
/// inputs were produced by another validated NURBS).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NurbsCurve {
    /// Polynomial degree `p` (1 = polyline, 3 = cubic, etc.).
    pub degree: usize,
    /// Knot vector — non-decreasing, length `n + p + 2`.
    pub knots: Vec<f64>,
    /// Control points — n+1 vectors in 3D.
    pub control_points: Vec<Vector3<f64>>,
    /// Per-control-point weights. Same length as `control_points`.
    /// All-ones reduces to a non-rational B-spline.
    pub weights: Vec<f64>,
}

impl NurbsCurve {
    /// Construct a validated NURBS curve.
    ///
    /// Validates:
    /// - `degree >= 1` and `<= 9` (practical upper bound),
    /// - `control_points.len() >= degree + 1`,
    /// - `weights.len() == control_points.len()`,
    /// - `knots.len() == control_points.len() + degree + 1`,
    /// - knots are non-decreasing.
    pub fn new(
        degree: usize,
        knots: Vec<f64>,
        control_points: Vec<Vector3<f64>>,
        weights: Vec<f64>,
    ) -> Result<Self, SurfaceError> {
        let curve = Self::new_unchecked(degree, knots, control_points, weights);
        curve.validate()?;
        Ok(curve)
    }

    /// Check the NURBS invariant on an already-built curve:
    /// - `degree >= 1` and `<= 9`,
    /// - `control_points.len() >= degree + 1`,
    /// - `weights.len() == control_points.len()`,
    /// - `knots.len() == control_points.len() + degree + 1`,
    /// - knots are non-decreasing.
    ///
    /// [`Self::new`] runs this at construction. A curve obtained another
    /// way — deserialised from an untrusted file, or via
    /// [`Self::new_unchecked`] — can re-check itself with this before the
    /// evaluation methods (which index `knots`/`control_points`) are
    /// called.
    pub fn validate(&self) -> Result<(), SurfaceError> {
        if !(1..=9).contains(&self.degree) {
            return Err(SurfaceError::BadDegree(self.degree));
        }
        let n_cp = self.control_points.len();
        if n_cp < self.degree + 1 {
            return Err(SurfaceError::BadKnotVector {
                reason: format!(
                    "need at least {} control points for degree {}, got {}",
                    self.degree + 1,
                    self.degree,
                    n_cp
                ),
            });
        }
        if self.weights.len() != n_cp {
            return Err(SurfaceError::BadKnotVector {
                reason: format!(
                    "weights len {} ≠ control_points len {}",
                    self.weights.len(),
                    n_cp
                ),
            });
        }
        let expected = n_cp + self.degree + 1;
        if self.knots.len() != expected {
            return Err(SurfaceError::BadKnotVector {
                reason: format!(
                    "expected {expected} knots (n_cp + degree + 1), got {}",
                    self.knots.len()
                ),
            });
        }
        for w in self.knots.windows(2) {
            if w[1] < w[0] {
                return Err(SurfaceError::BadKnotVector {
                    reason: "knots must be non-decreasing".into(),
                });
            }
        }
        Ok(())
    }

    /// Skip validation — caller asserts the inputs are well-formed.
    /// Used by surface tessellation that's already produced valid
    /// intermediate curves.
    pub fn new_unchecked(
        degree: usize,
        knots: Vec<f64>,
        control_points: Vec<Vector3<f64>>,
        weights: Vec<f64>,
    ) -> Self {
        Self {
            degree,
            knots,
            control_points,
            weights,
        }
    }

    /// Number of control points.
    pub fn n_control_points(&self) -> usize {
        self.control_points.len()
    }

    /// Valid parameter range: `[knots[degree], knots[n]]` where
    /// `n = n_cp`. The clamped-endpoint convention means
    /// `evaluate(u_min)` returns the first CP and `evaluate(u_max)`
    /// returns the last CP whenever the knot vector is clamped.
    pub fn parameter_range(&self) -> (f64, f64) {
        let n = self.control_points.len();
        (self.knots[self.degree], self.knots[n])
    }

    /// Find the knot span `k` such that
    /// `knots[k] <= u < knots[k+1]`, clamped to the valid range so
    /// `u == u_max` returns `n - 1` (avoids edge-case crashes when
    /// callers evaluate at the upper endpoint).
    ///
    /// Linear scan — fine for typical knot vectors. Binary search is
    /// a v1.5 optimisation if it ever shows up in a profile.
    pub fn find_knot_span(&self, u: f64) -> usize {
        find_knot_span(u, &self.knots, self.degree, self.control_points.len())
    }

    /// Basis-function values for span `k`. Returns `degree + 1`
    /// floats: `N_{k-p,p}(u), N_{k-p+1,p}(u), ..., N_{k,p}(u)`.
    ///
    /// Standard Cox-de Boor triangular recursion — see e.g. *The
    /// NURBS Book* (Piegl & Tiller), Algorithm A2.2.
    pub fn basis_functions(&self, span: usize, u: f64) -> Vec<f64> {
        basis_functions(span, u, self.degree, &self.knots)
    }

    /// Evaluate the NURBS curve at parameter `u`.
    ///
    /// Returns the weighted sum of control points × basis functions,
    /// divided by the weighted basis sum.
    ///
    /// For a clamped knot vector:
    /// - `evaluate(u_min)` returns `control_points[0]`,
    /// - `evaluate(u_max)` returns the last control point,
    /// - intermediate values are convex combinations of CPs in the
    ///   local span (so the curve always lies inside the convex hull
    ///   of the control polygon).
    pub fn evaluate(&self, u: f64) -> Vector3<f64> {
        let span = self.find_knot_span(u);
        let basis = self.basis_functions(span, u);
        let mut num = Vector3::zeros();
        let mut den = 0.0_f64;
        for (i, b) in basis.iter().enumerate() {
            let cp_idx = span - self.degree + i;
            let w = self.weights[cp_idx];
            let wb = w * b;
            num += self.control_points[cp_idx] * wb;
            den += wb;
        }
        if den.abs() < 1e-30 {
            // Pathological: all-zero weights at this u. Return the
            // un-normalised numerator so the caller still gets a
            // finite point near the control polygon.
            num
        } else {
            num / den
        }
    }

    /// k-th derivative at parameter `u`.
    ///
    /// v1 uses central finite differences with step `h = 1e-4` over
    /// the valid parameter range. Analytic derivatives via shifted
    /// control points are a v1.5 optimisation; for tessellation +
    /// tangent-vector callers (this crate's only consumers), finite
    /// differences give 6-7 digits of accuracy which is plenty.
    ///
    /// Returns the zero vector if `k == 0`-th derivative is requested
    /// after the first via repeated convolution — callers normally
    /// just want the first or second derivative.
    pub fn derivative(&self, u: f64, k: usize) -> Vector3<f64> {
        if k == 0 {
            return self.evaluate(u);
        }
        let (u_min, u_max) = self.parameter_range();
        let h = 1e-4_f64.max((u_max - u_min) * 1e-5);
        // Centre the stencil but clamp to the valid range so we
        // don't hit EvaluationOutOfRange near the endpoints.
        let u_lo = (u - h).max(u_min);
        let u_hi = (u + h).min(u_max);
        let denom = u_hi - u_lo;
        if denom.abs() < 1e-30 {
            return Vector3::zeros();
        }
        if k == 1 {
            (self.evaluate(u_hi) - self.evaluate(u_lo)) / denom
        } else {
            // Recursive central difference for higher orders.
            (self.derivative(u_hi, k - 1) - self.derivative(u_lo, k - 1)) / denom
        }
    }
}

/// Free helper — exposed so `nurbs_surface` can use it without
/// going through a curve instance.
pub fn find_knot_span(u: f64, knots: &[f64], degree: usize, n_cp: usize) -> usize {
    // n_cp is the number of control points; the canonical
    // "highest valid span index" is n_cp - 1 (when u == u_max).
    let n = n_cp;
    if u >= knots[n] {
        return n - 1;
    }
    if u <= knots[degree] {
        return degree;
    }
    // Linear scan from `degree` upward. `knots.len() == n + degree + 1`,
    // so the highest index we can probe safely is `n - 1`.
    for k in degree..n {
        if knots[k] <= u && u < knots[k + 1] {
            return k;
        }
    }
    n - 1
}

/// Cox-de Boor basis functions for span `k`. Returns `degree + 1`
/// values. Standard Piegl & Tiller A2.2.
pub fn basis_functions(span: usize, u: f64, degree: usize, knots: &[f64]) -> Vec<f64> {
    let p = degree;
    let mut n = vec![0.0_f64; p + 1];
    let mut left = vec![0.0_f64; p + 1];
    let mut right = vec![0.0_f64; p + 1];
    n[0] = 1.0;
    for j in 1..=p {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            let denom = right[r + 1] + left[j - r];
            // Guard the zero/zero case (repeated knots): contribution
            // is just dropped.
            let temp = if denom.abs() < 1e-30 {
                0.0
            } else {
                n[r] / denom
            };
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }
    n
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

    #[test]
    fn rejects_bad_degree() {
        let err = NurbsCurve::new(0, bezier_knots(), vec![Vector3::zeros(); 4], vec![1.0; 4])
            .unwrap_err();
        assert_eq!(err.code(), "surface.bad_degree");

        let err = NurbsCurve::new(99, bezier_knots(), vec![Vector3::zeros(); 4], vec![1.0; 4])
            .unwrap_err();
        assert_eq!(err.code(), "surface.bad_degree");
    }

    #[test]
    fn rejects_short_knot_vector() {
        // 4 CPs + degree 3 → expect 8 knots; provide 7.
        let err = NurbsCurve::new(
            3,
            vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![Vector3::zeros(); 4],
            vec![1.0; 4],
        )
        .unwrap_err();
        assert_eq!(err.code(), "surface.bad_knot_vector");
    }

    #[test]
    fn rejects_non_decreasing_knots() {
        let err = NurbsCurve::new(
            3,
            vec![0.0, 0.0, 0.0, 0.5, 0.3, 1.0, 1.0, 1.0],
            vec![Vector3::zeros(); 4],
            vec![1.0; 4],
        )
        .unwrap_err();
        assert_eq!(err.code(), "surface.bad_knot_vector");
    }

    #[test]
    fn accepts_well_formed_cubic_bezier() {
        let c = cubic_bezier([
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
        ]);
        assert_eq!(c.degree, 3);
        assert_eq!(c.n_control_points(), 4);
        assert_eq!(c.parameter_range(), (0.0, 1.0));
    }

    // ===== Phase 9B — evaluation tests =====

    #[test]
    fn find_knot_span_basic() {
        // Open uniform knot vector for 5 CPs, degree 2 →
        // [0,0,0,0.5,1,1,1] of length 7.
        let knots = vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0];
        // span at u=0.25 should be 2 (knots[2]=0 <= 0.25 < knots[3]=0.5)
        assert_eq!(find_knot_span(0.25, &knots, 2, 4), 2);
        assert_eq!(find_knot_span(0.5, &knots, 2, 4), 3);
        assert_eq!(find_knot_span(0.75, &knots, 2, 4), 3);
        // upper endpoint clamps to n - 1
        assert_eq!(find_knot_span(1.0, &knots, 2, 4), 3);
    }

    #[test]
    fn basis_functions_partition_of_unity() {
        // Cox-de Boor basis functions sum to 1 at every interior u.
        let knots = bezier_knots();
        for &u in &[0.0_f64, 0.1, 0.25, 0.5, 0.7, 0.9, 1.0] {
            let span = find_knot_span(u, &knots, 3, 4);
            let n = basis_functions(span, u, 3, &knots);
            let sum: f64 = n.iter().sum();
            assert!((sum - 1.0).abs() < 1e-10, "u={u}: sum={sum} != 1");
        }
    }

    #[test]
    fn cubic_bezier_endpoint_clamps() {
        // At u=0 the cubic Bezier returns the first CP exactly.
        // At u=1 it returns the last CP exactly.
        let cps = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 2.0, 3.0),
            Vector3::new(2.0, 5.0, 1.0),
            Vector3::new(4.0, 0.0, 0.0),
        ];
        let c = cubic_bezier(cps);
        let start = c.evaluate(0.0);
        let end = c.evaluate(1.0);
        assert!((start - cps[0]).norm() < 1e-10, "start={start:?}");
        assert!((end - cps[3]).norm() < 1e-10, "end={end:?}");
    }

    #[test]
    fn cubic_bezier_midpoint_matches_decasteljau() {
        // Canonical reference: cubic Bezier with knots [0,0,0,0,1,1,1,1]
        // at u=0.5 gives the de Casteljau midpoint:
        //   m_01 = (P0 + P1) / 2
        //   m_12 = (P1 + P2) / 2
        //   m_23 = (P2 + P3) / 2
        //   m_012 = (m_01 + m_12) / 2
        //   m_123 = (m_12 + m_23) / 2
        //   B(0.5) = (m_012 + m_123) / 2
        //         = (P0 + 3 P1 + 3 P2 + P3) / 8
        let p0 = Vector3::new(0.0, 0.0, 0.0);
        let p1 = Vector3::new(1.0, 2.0, 0.0);
        let p2 = Vector3::new(3.0, 2.0, 0.0);
        let p3 = Vector3::new(4.0, 0.0, 0.0);
        let c = cubic_bezier([p0, p1, p2, p3]);
        let got = c.evaluate(0.5);
        let expected = (p0 + 3.0 * p1 + 3.0 * p2 + p3) / 8.0;
        assert!(
            (got - expected).norm() < 1e-10,
            "got {got:?}, expected {expected:?}"
        );
    }

    #[test]
    fn rational_weight_one_matches_polynomial() {
        // Weights of 1 → polynomial Bezier (the basis denominator
        // collapses to the partition-of-unity sum = 1).
        let cps = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(2.0, 1.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
        ];
        let c = cubic_bezier(cps);
        for &u in &[0.25_f64, 0.5, 0.75] {
            let polyn = (1.0 - u).powi(3) * cps[0]
                + 3.0 * (1.0 - u).powi(2) * u * cps[1]
                + 3.0 * (1.0 - u) * u.powi(2) * cps[2]
                + u.powi(3) * cps[3];
            let got = c.evaluate(u);
            assert!(
                (got - polyn).norm() < 1e-10,
                "u={u}: got {got:?}, expected {polyn:?}"
            );
        }
    }

    #[test]
    fn weighted_curve_passes_through_endpoints() {
        // Even with non-uniform weights, the endpoints are returned
        // exactly because basis at u=u_min/u_max is concentrated on
        // the corresponding control point.
        let cps = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(2.0, 1.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
        ];
        let weights = vec![1.0, 2.0, 0.5, 1.0];
        let c = NurbsCurve::new(3, bezier_knots(), cps.clone(), weights).unwrap();
        let start = c.evaluate(0.0);
        let end = c.evaluate(1.0);
        assert!((start - cps[0]).norm() < 1e-10);
        assert!((end - cps[3]).norm() < 1e-10);
    }

    #[test]
    fn straight_line_derivative_is_constant() {
        // 4 collinear CPs along x → tangent at any u is the +x axis.
        let c = cubic_bezier([
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
        ]);
        let d_mid = c.derivative(0.5, 1);
        // Derivative of a degree-3 Bezier with collinear CPs along x
        // from x=0 to x=3 over u∈[0,1] is the constant 3.0 (the
        // total length scaled by the unit parameter).
        assert!((d_mid.x - 3.0).abs() < 1e-3, "d_mid.x = {}", d_mid.x);
        assert!(d_mid.y.abs() < 1e-3);
        assert!(d_mid.z.abs() < 1e-3);
    }
}
