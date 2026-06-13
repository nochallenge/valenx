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

    /// The **unit tangent** vector at parameter `u` — the normalised first
    /// derivative `C'(u)/|C'(u)|`, the direction of travel along the curve.
    ///
    /// Returns the zero vector at a singular point where `|C'(u)| ≈ 0`
    /// (a cusp, or a momentarily stationary parameterisation), where the
    /// tangent direction is undefined.
    pub fn unit_tangent(&self, u: f64) -> Vector3<f64> {
        let d = self.derivative(u, 1);
        let n = d.norm();
        if n < 1e-12 {
            Vector3::zeros()
        } else {
            d / n
        }
    }

    /// The **curvature** `κ(u)` at parameter `u` — the reciprocal of the
    /// osculating-circle radius, `κ = |C'(u) × C''(u)| / |C'(u)|³`.
    ///
    /// This is the parameterisation-independent bending rate: a straight
    /// segment has `κ = 0` everywhere, and a circle of radius `r` has the
    /// constant `κ = 1/r` at every point (whatever its NURBS
    /// parameterisation). `C'` and `C''` are taken from a single centred
    /// 3-point stencil `C(u−h), C(u), C(u+h)` so the two are mutually
    /// consistent; the stencil centre is clamped to keep all three samples
    /// inside the valid [`parameter_range`](Self::parameter_range).
    ///
    /// Returns `0.0` at a singular point (`|C'(u)| ≈ 0`), where the osculating
    /// circle — and hence the curvature — is undefined, and for a domain too
    /// small to form the stencil.
    pub fn curvature(&self, u: f64) -> f64 {
        let (u_min, u_max) = self.parameter_range();
        let span = u_max - u_min;
        let h = (span * 1e-3).max(1e-6);
        if span < 4.0 * h {
            return 0.0;
        }
        let uc = u.max(u_min + h).min(u_max - h);
        let c_minus = self.evaluate(uc - h);
        let c_0 = self.evaluate(uc);
        let c_plus = self.evaluate(uc + h);
        let d1 = (c_plus - c_minus) / (2.0 * h);
        let d2 = (c_plus - 2.0 * c_0 + c_minus) / (h * h);
        let speed = d1.norm();
        if speed < 1e-12 {
            return 0.0;
        }
        d1.cross(&d2).norm() / speed.powi(3)
    }

    /// The **arc length** of the whole curve — the geometric length of the
    /// traced path over its full [`parameter_range`](Self::parameter_range),
    /// `∫ |C'(u)| du`. A convenience for
    /// [`arc_length_between`](Self::arc_length_between)`(u_min, u_max)`.
    pub fn arc_length(&self) -> f64 {
        let (u_min, u_max) = self.parameter_range();
        self.arc_length_between(u_min, u_max)
    }

    /// The **arc length between two parameters** — `∫_{u0}^{u1} |C'(u)| du`,
    /// the geometric length of the curve segment from `u0` to `u1`.
    ///
    /// Computed by composite **Simpson's rule** on the speed `|C'(u)|`: the
    /// interval is first split at every interior knot (so no panel straddles a
    /// knot, where the speed is only `C^{p-1}`-continuous), then each smooth
    /// span is integrated with a 64-panel composite Simpson rule (4th-order
    /// accurate). The result is parameterisation-independent — it is the true
    /// geometric length, not the parameter span.
    ///
    /// `u0`/`u1` are clamped to the valid range; a reversed or degenerate
    /// interval (`u1 ≤ u0`) returns `0.0`.
    pub fn arc_length_between(&self, u0: f64, u1: f64) -> f64 {
        const PANELS_PER_SPAN: usize = 64; // even, for composite Simpson

        let (u_min, u_max) = self.parameter_range();
        let a = u0.clamp(u_min, u_max);
        let b = u1.clamp(u_min, u_max);
        if b <= a {
            return 0.0;
        }
        // Breakpoints: a, the distinct interior knots in (a, b), then b.
        let eps = 1e-12;
        let mut breaks = vec![a];
        let mut prev = a;
        for &k in &self.knots {
            if k > a + eps && k < b - eps && k > prev + eps {
                breaks.push(k);
                prev = k;
            }
        }
        breaks.push(b);

        let speed = |u: f64| self.derivative(u, 1).norm();
        let mut total = 0.0;
        for seg in breaks.windows(2) {
            let (sa, sb) = (seg[0], seg[1]);
            let n = PANELS_PER_SPAN;
            let h = (sb - sa) / n as f64;
            let mut sum = speed(sa) + speed(sb);
            for i in 1..n {
                let u = sa + h * i as f64;
                sum += if i % 2 == 1 { 4.0 } else { 2.0 } * speed(u);
            }
            total += h / 3.0 * sum;
        }
        total
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
    fn rational_quadratic_traces_exact_circle() {
        // GROUND TRUTH: a rational quadratic NURBS with the canonical
        // weights represents a conic EXACTLY. The standard 90° arc uses
        // control points P0=(r,0), P1=(r,r), P2=(0,r) with weights
        // (1, √2/2, 1) and clamped knots [0,0,0,1,1,1]; the resulting
        // curve is the EXACT quarter circle x² + y² = r², not an
        // approximation. (Middle weight w = cos(45°) = √2/2 for a 90°
        // sweep; the general rule is w = cos(half-angle).)
        //
        // This pins the exact-conic property that partition-of-unity and
        // endpoint-clamp tests do not cover. Tolerance is 1e-12 because
        // the identity is algebraically exact in f64 — the only error is
        // floating-point round-off in the Bernstein/denominator sums.
        let r = 2.5_f64;
        let w = std::f64::consts::FRAC_1_SQRT_2; // √2/2 = cos45°
        let cps = vec![
            Vector3::new(r, 0.0, 0.0),
            Vector3::new(r, r, 0.0),
            Vector3::new(0.0, r, 0.0),
        ];
        let weights = vec![1.0, w, 1.0];
        let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let c = NurbsCurve::new(2, knots, cps, weights).unwrap();
        // Sample across the whole domain, including the interior where a
        // mere polygon/parabola would deviate from the true circle.
        for &u in &[0.0_f64, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
            let p = c.evaluate(u);
            let r2 = p.x * p.x + p.y * p.y;
            assert!(
                (r2 - r * r).abs() < 1e-12,
                "u={u}: point {p:?} has r²={r2} ≠ {} (off circle)",
                r * r
            );
            assert!(p.z.abs() < 1e-12, "u={u}: arc left the z=0 plane");
        }
        // Spot-check the geometric midpoint: a 90° rational-quadratic arc
        // evaluated at u=0.5 lands at (r/√2, r/√2) — the 45° point —
        // which a non-rational parabola through the same CPs would miss.
        let mid = c.evaluate(0.5);
        let half = r * std::f64::consts::FRAC_1_SQRT_2;
        assert!(
            (mid.x - half).abs() < 1e-12 && (mid.y - half).abs() < 1e-12,
            "midpoint {mid:?} != (r/√2, r/√2) = ({half}, {half})"
        );
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

    // ===== differential geometry: tangent / curvature / arc length =====

    /// The canonical exact-circle rational quadratic: a 90° arc of radius `r`,
    /// CPs (r,0),(r,r),(0,r), weights (1, √2/2, 1), clamped knots
    /// [0,0,0,1,1,1] — the EXACT quarter circle (see
    /// `rational_quadratic_traces_exact_circle`).
    fn quarter_circle(r: f64) -> NurbsCurve {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        NurbsCurve::new(
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![
                Vector3::new(r, 0.0, 0.0),
                Vector3::new(r, r, 0.0),
                Vector3::new(0.0, r, 0.0),
            ],
            vec![1.0, w, 1.0],
        )
        .unwrap()
    }

    fn x_line() -> NurbsCurve {
        cubic_bezier([
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
        ])
    }

    #[test]
    fn arc_length_of_quarter_circle_is_quarter_circumference() {
        // GROUND TRUTH: the rational-quadratic quarter circle is the EXACT
        // circle x²+y²=r², so its arc length is exactly a quarter of the
        // circumference, (π/2)·r — independent of the NURBS parameterization.
        // For r=2.5 that is π·2.5/2 = 3.9269908169872414. The composite-Simpson
        // integral of |C'(u)| recovers it to the finite-difference floor.
        let r = 2.5_f64;
        let c = quarter_circle(r);
        let expected = std::f64::consts::FRAC_PI_2 * r; // (π/2)·r
        let len = c.arc_length();
        assert!(
            (len - expected).abs() < 1e-5,
            "quarter-circle arc length {len} != {expected} = (π/2)·r"
        );
    }

    #[test]
    fn arc_length_of_straight_line_is_endpoint_distance() {
        // Collinear cubic Bezier from x=0 to x=3 → C(u)=(3u,0,0), a straight
        // segment of length exactly 3 (the integrand |C'|=3 is constant, so
        // Simpson is exact to round-off).
        let c = x_line();
        assert!(
            (c.arc_length() - 3.0).abs() < 1e-9,
            "line length {}",
            c.arc_length()
        );
        // Partial length over the first parameter-half is 1.5 (uniform here).
        let (u0, u1) = c.parameter_range();
        let half = c.arc_length_between(u0, 0.5 * (u0 + u1));
        assert!((half - 1.5).abs() < 1e-9, "half-line length {half}");
    }

    #[test]
    fn arc_length_between_is_additive_and_clamped() {
        let c = quarter_circle(1.0);
        let (u0, u1) = c.parameter_range();
        let mid = 0.5 * (u0 + u1);
        let whole = c.arc_length();
        let part_a = c.arc_length_between(u0, mid);
        let part_b = c.arc_length_between(mid, u1);
        // Additivity: the two halves sum to the whole — to the integration
        // floor, since the whole and the halves use different Simpson panel
        // subdivisions and finite-difference sample points (~1e-7), not
        // bit-exactly. (The straight-line case above, with constant |C'|, is
        // exact and holds at 1e-9.)
        assert!(
            (part_a + part_b - whole).abs() < 1e-6,
            "{part_a} + {part_b} != {whole}"
        );
        // Degenerate and reversed intervals return 0.
        assert_eq!(c.arc_length_between(mid, mid), 0.0);
        assert_eq!(c.arc_length_between(u1, u0), 0.0);
    }

    #[test]
    fn curvature_of_circle_is_inverse_radius() {
        // GROUND TRUTH: a circle of radius r has constant curvature κ = 1/r at
        // every point, independent of parameterization. The exact rational-
        // quadratic quarter circle of radius r=2 must read κ = 0.5 everywhere.
        let r = 2.0_f64;
        let c = quarter_circle(r);
        for &u in &[0.2_f64, 0.4, 0.5, 0.6, 0.8] {
            let k = c.curvature(u);
            assert!(
                (k - 1.0 / r).abs() < 1e-4,
                "u={u}: curvature {k} != 1/r = {}",
                1.0 / r
            );
        }
    }

    #[test]
    fn curvature_of_straight_line_is_zero() {
        // A straight segment has zero curvature everywhere (C'' = 0).
        let c = x_line();
        for &u in &[0.25_f64, 0.5, 0.75] {
            assert!(
                c.curvature(u).abs() < 1e-6,
                "line curvature {} at u={u}",
                c.curvature(u)
            );
        }
    }

    #[test]
    fn unit_tangent_is_normalized_and_directionally_correct() {
        // Straight line along +x → unit tangent (1,0,0) everywhere.
        let line = x_line();
        let t = line.unit_tangent(0.5);
        assert!(
            (t - Vector3::new(1.0, 0.0, 0.0)).norm() < 1e-6,
            "line tangent {t:?}"
        );
        assert!((t.norm() - 1.0).abs() < 1e-9);
        // The quarter circle starts at (r,0) sweeping CCW toward (0,r), so the
        // tangent at the start points in +y.
        let circ = quarter_circle(2.0);
        let t0 = circ.unit_tangent(0.0);
        assert!(
            (t0 - Vector3::new(0.0, 1.0, 0.0)).norm() < 1e-3,
            "circle start tangent {t0:?}"
        );
        assert!((t0.norm() - 1.0).abs() < 1e-9);
    }
}
