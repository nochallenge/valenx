//! Natural cubic spline interpolation.
//!
//! Given knots `(x_0, y_0) .. (x_n, y_n)` with strictly ascending
//! abscissae, the natural cubic spline `S(x)` is the unique function
//! that
//!
//! 1. is a cubic polynomial on each segment `[x_i, x_{i+1}]`,
//! 2. interpolates the data: `S(x_i) = y_i`,
//! 3. is twice continuously differentiable (`C2`) across every
//!    interior knot, and
//! 4. has zero second derivative at both ends
//!    (`S''(x_0) = S''(x_n) = 0`, the *natural* boundary condition).
//!
//! ## Solving for the moments
//!
//! With `M_i = S''(x_i)` and segment widths `h_i = x_{i+1} - x_i`,
//! `C2` continuity gives, for `i = 1 .. n-1`,
//!
//! `h_{i-1} M_{i-1} + 2 (h_{i-1} + h_i) M_i + h_i M_{i+1}`
//! `  = 6 ( (y_{i+1} - y_i)/h_i - (y_i - y_{i-1})/h_{i-1} )`.
//!
//! Together with `M_0 = M_n = 0` this is a symmetric, diagonally
//! dominant tridiagonal system, solved here with the Thomas algorithm
//! (an `O(n)` specialisation of Gaussian elimination). On each
//! segment,
//!
//! `S(x) = M_i (x_{i+1} - x)^3 / (6 h_i)`
//! `     + M_{i+1} (x - x_i)^3 / (6 h_i)`
//! `     + (y_i/h_i - M_i h_i/6) (x_{i+1} - x)`
//! `     + (y_{i+1}/h_i - M_{i+1} h_i/6) (x - x_i)`.

use nalgebra::DVector;

use crate::error::InterpolationError;
use crate::points::DataPoints;

/// A natural cubic spline through a set of knots.
///
/// Build with [`CubicSpline::new`]; query the value, first, or second
/// derivative with [`CubicSpline::value_at`],
/// [`CubicSpline::derivative_at`], and
/// [`CubicSpline::second_derivative_at`].
#[derive(Clone, Debug, PartialEq)]
pub struct CubicSpline {
    data: DataPoints,
    /// Second derivatives `M_i = S''(x_i)` at each knot.
    moments: Vec<f64>,
}

impl CubicSpline {
    /// Build a natural cubic spline from validated [`DataPoints`].
    ///
    /// # Errors
    ///
    /// Returns [`InterpolationError::TooFewPoints`] if fewer than two
    /// knots are supplied (a spline needs at least one segment).
    pub fn new(data: DataPoints) -> Result<Self, InterpolationError> {
        let n = data.len();
        if n < 2 {
            return Err(InterpolationError::TooFewPoints { needed: 2, got: n });
        }
        let moments = solve_moments(&data);
        Ok(Self { data, moments })
    }

    /// Build directly from `x` / `y` slices.
    ///
    /// # Errors
    ///
    /// Propagates any [`DataPoints::new`] validation error, then any
    /// [`CubicSpline::new`] error.
    pub fn from_points(xs: &[f64], ys: &[f64]) -> Result<Self, InterpolationError> {
        Self::new(DataPoints::new(xs, ys)?)
    }

    /// The underlying knots.
    pub fn data(&self) -> &DataPoints {
        &self.data
    }

    /// Second derivatives `S''(x_i)` at the knots, also called the
    /// spline moments. The first and last are zero for the natural
    /// boundary condition.
    pub fn moments(&self) -> &[f64] {
        &self.moments
    }

    /// Evaluate `S(x)` at `at`.
    ///
    /// # Errors
    ///
    /// Returns [`InterpolationError::OutOfDomain`] if `at` lies
    /// outside `[x_first, x_last]`.
    pub fn value_at(&self, at: f64) -> Result<f64, InterpolationError> {
        self.data.check_domain(at)?;
        let (i, h, a, b) = self.locate(at);
        let ys = self.data.ys();
        let (mi, mj) = (self.moments[i], self.moments[i + 1]);
        // a = x_{i+1} - x, b = x - x_i (both >= 0).
        let term_curv = (mi * a.powi(3) + mj * b.powi(3)) / (6.0 * h);
        let term_lin = (ys[i] / h - mi * h / 6.0) * a + (ys[i + 1] / h - mj * h / 6.0) * b;
        Ok(term_curv + term_lin)
    }

    /// Evaluate the first derivative `S'(x)` at `at`.
    ///
    /// # Errors
    ///
    /// Returns [`InterpolationError::OutOfDomain`] if `at` lies
    /// outside `[x_first, x_last]`.
    pub fn derivative_at(&self, at: f64) -> Result<f64, InterpolationError> {
        self.data.check_domain(at)?;
        let (i, h, a, b) = self.locate(at);
        let ys = self.data.ys();
        let (mi, mj) = (self.moments[i], self.moments[i + 1]);
        // d/dx of the segment formula. d a/dx = -1, d b/dx = +1.
        let d_curv = (-mi * a.powi(2) + mj * b.powi(2)) / (2.0 * h);
        let d_lin = -(ys[i] / h - mi * h / 6.0) + (ys[i + 1] / h - mj * h / 6.0);
        Ok(d_curv + d_lin)
    }

    /// Evaluate the second derivative `S''(x)` at `at`.
    ///
    /// `S''` is linear on each segment, interpolating the knot moments
    /// `M_i` and `M_{i+1}`; this is what makes the spline `C2`.
    ///
    /// # Errors
    ///
    /// Returns [`InterpolationError::OutOfDomain`] if `at` lies
    /// outside `[x_first, x_last]`.
    pub fn second_derivative_at(&self, at: f64) -> Result<f64, InterpolationError> {
        self.data.check_domain(at)?;
        let (i, h, a, b) = self.locate(at);
        let (mi, mj) = (self.moments[i], self.moments[i + 1]);
        Ok((mi * a + mj * b) / h)
    }

    /// Locate the bracketing segment and return `(i, h_i, a, b)` with
    /// `a = x_{i+1} - at` and `b = at - x_i`.
    fn locate(&self, at: f64) -> (usize, f64, f64, f64) {
        let xs = self.data.xs();
        let i = self.data.segment_index(at);
        let h = xs[i + 1] - xs[i];
        let a = xs[i + 1] - at;
        let b = at - xs[i];
        (i, h, a, b)
    }
}

/// Solve the tridiagonal moment system for the natural cubic spline.
///
/// Returns the vector of knot second derivatives `M_0 .. M_n` with
/// `M_0 = M_n = 0`. For the interior unknowns the symmetric
/// tridiagonal system is solved with the Thomas algorithm; the small
/// dense right-hand side is assembled through an [`nalgebra::DVector`]
/// so the linear-algebra dependency carries the data.
fn solve_moments(data: &DataPoints) -> Vec<f64> {
    let xs = data.xs();
    let ys = data.ys();
    let n = xs.len();
    let mut moments = vec![0.0_f64; n];
    // Only interior knots 1 .. n-1 are unknown; with < 3 knots there
    // are no interior unknowns and the spline is the connecting line
    // (both moments already zero).
    if n < 3 {
        return moments;
    }
    let m = n - 2; // number of interior unknowns
    let h: Vec<f64> = (0..n - 1).map(|i| xs[i + 1] - xs[i]).collect();

    // Right-hand side d_k = 6 * (slope_{k} - slope_{k-1}) for interior knot k.
    let mut rhs = DVector::<f64>::zeros(m);
    for k in 0..m {
        let i = k + 1; // interior knot index in 1 .. n-1
        let slope_right = (ys[i + 1] - ys[i]) / h[i];
        let slope_left = (ys[i] - ys[i - 1]) / h[i - 1];
        rhs[k] = 6.0 * (slope_right - slope_left);
    }

    // Tridiagonal coefficients for interior rows (row k, interior knot i = k+1):
    //   sub_k = h[i-1]              (below diagonal)
    //   diag_k = 2 (h[i-1] + h[i])
    //   sup_k = h[i]               (above diagonal)
    let sub: Vec<f64> = (0..m).map(|k| h[k]).collect(); // h[i-1] with i=k+1 -> h[k]
    let diag: Vec<f64> = (0..m).map(|k| 2.0 * (h[k] + h[k + 1])).collect();
    let sup: Vec<f64> = (0..m).map(|k| h[k + 1]).collect(); // h[i] with i=k+1 -> h[k+1]

    // Thomas algorithm (forward sweep then back substitution).
    let mut c_prime = vec![0.0_f64; m];
    let mut d_prime = vec![0.0_f64; m];
    c_prime[0] = sup[0] / diag[0];
    d_prime[0] = rhs[0] / diag[0];
    for k in 1..m {
        let denom = diag[k] - sub[k] * c_prime[k - 1];
        c_prime[k] = sup[k] / denom;
        d_prime[k] = (rhs[k] - sub[k] * d_prime[k - 1]) / denom;
    }
    let mut sol = vec![0.0_f64; m];
    sol[m - 1] = d_prime[m - 1];
    for k in (0..m - 1).rev() {
        sol[k] = d_prime[k] - c_prime[k] * sol[k + 1];
    }
    for (k, &mk) in sol.iter().enumerate() {
        moments[k + 1] = mk;
    }
    moments
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    fn close(a: f64, b: f64) {
        assert!((a - b).abs() < EPS, "expected {a} ~= {b}");
    }

    #[test]
    fn passes_through_every_knot() {
        let xs = [0.0, 1.0, 2.0, 4.0, 5.0];
        let ys = [0.0, 1.0, 0.0, 2.0, -1.0];
        let s = CubicSpline::from_points(&xs, &ys).unwrap();
        for (&x, &y) in xs.iter().zip(ys.iter()) {
            close(s.value_at(x).unwrap(), y);
        }
    }

    #[test]
    fn natural_boundary_second_derivative_is_zero() {
        let xs = [0.0, 1.0, 2.0, 3.0];
        let ys = [1.0, 3.0, 2.0, 5.0];
        let s = CubicSpline::from_points(&xs, &ys).unwrap();
        close(s.second_derivative_at(0.0).unwrap(), 0.0);
        close(s.second_derivative_at(3.0).unwrap(), 0.0);
        // moments at the ends are exactly zero.
        close(s.moments()[0], 0.0);
        close(s.moments()[3], 0.0);
    }

    #[test]
    fn c2_continuity_at_interior_knots() {
        // The value and both derivatives must have no jump across each
        // interior knot. A genuine discontinuity would leave a jump of
        // O(1) independent of the probe offset `d`; a continuous spline
        // leaves only an O(d) difference between the one-sided limits
        // sampled at `xk - d` and `xk + d`, which vanishes as d -> 0.
        let xs = [0.0, 1.0, 2.0, 3.0, 4.0];
        let ys = [0.0, 2.0, 1.0, 3.0, 0.0];
        let s = CubicSpline::from_points(&xs, &ys).unwrap();
        let d = 1e-7;
        for (k, &xk) in xs.iter().enumerate().take(xs.len() - 1).skip(1) {
            // C0: the value at the knot equals its ordinate exactly,
            // and the one-sided values straddle it within O(d).
            close(s.value_at(xk).unwrap(), ys[k]);
            let vl = s.value_at(xk - d).unwrap();
            let vr = s.value_at(xk + d).unwrap();
            assert!((vl - vr).abs() < 1e-4, "C0 break at {xk}: {vl} vs {vr}");
            // C1: one-sided first derivatives agree across the knot.
            let dl = s.derivative_at(xk - d).unwrap();
            let dr = s.derivative_at(xk + d).unwrap();
            assert!((dl - dr).abs() < 1e-4, "C1 break at {xk}: {dl} vs {dr}");
            // C2: one-sided second derivatives agree across the knot.
            let sl = s.second_derivative_at(xk - d).unwrap();
            let sr = s.second_derivative_at(xk + d).unwrap();
            assert!((sl - sr).abs() < 1e-4, "C2 break at {xk}: {sl} vs {sr}");
        }
    }

    #[test]
    fn first_derivative_matches_finite_difference() {
        let xs = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let ys = [0.0, 1.0, 4.0, 9.0, 7.0, 2.0];
        let s = CubicSpline::from_points(&xs, &ys).unwrap();
        let d = 1e-6;
        for k in 1..=9 {
            let x = 0.5 * k as f64;
            let analytic = s.derivative_at(x).unwrap();
            let fd = (s.value_at(x + d).unwrap() - s.value_at(x - d).unwrap()) / (2.0 * d);
            assert!((analytic - fd).abs() < 1e-4, "deriv mismatch at {x}");
        }
    }

    #[test]
    fn second_derivative_matches_finite_difference() {
        let xs = [0.0, 1.0, 2.0, 3.0, 4.0];
        let ys = [1.0, 0.0, 2.0, -1.0, 3.0];
        let s = CubicSpline::from_points(&xs, &ys).unwrap();
        let d = 1e-4;
        for k in 1..=7 {
            let x = 0.5 * k as f64;
            let analytic = s.second_derivative_at(x).unwrap();
            let fd = (s.value_at(x + d).unwrap() - 2.0 * s.value_at(x).unwrap()
                + s.value_at(x - d).unwrap())
                / (d * d);
            assert!((analytic - fd).abs() < 1e-2, "2nd deriv mismatch at {x}");
        }
    }

    #[test]
    fn straight_line_data_is_reproduced_exactly() {
        // A natural cubic spline through collinear points is the line
        // itself (all moments vanish).
        let f = |x: f64| 3.0 * x - 2.0;
        let xs = [0.0, 1.0, 2.0, 3.0, 4.0];
        let ys: Vec<f64> = xs.iter().map(|&x| f(x)).collect();
        let s = CubicSpline::from_points(&xs, &ys).unwrap();
        for &m in s.moments() {
            close(m, 0.0);
        }
        for k in 0..=16 {
            let x = 4.0 * (k as f64 / 16.0);
            close(s.value_at(x).unwrap(), f(x));
        }
    }

    #[test]
    fn constant_data_gives_constant_output() {
        let s = CubicSpline::from_points(&[0.0, 1.0, 2.0, 3.0, 4.0], &[6.0; 5]).unwrap();
        for k in 0..=16 {
            let x = 4.0 * (k as f64 / 16.0);
            close(s.value_at(x).unwrap(), 6.0);
            close(s.derivative_at(x).unwrap(), 0.0);
            close(s.second_derivative_at(x).unwrap(), 0.0);
        }
    }

    #[test]
    fn two_points_is_the_connecting_line() {
        let s = CubicSpline::from_points(&[0.0, 2.0], &[1.0, 5.0]).unwrap();
        close(s.value_at(1.0).unwrap(), 3.0);
        close(s.value_at(0.5).unwrap(), 2.0);
        close(s.second_derivative_at(1.0).unwrap(), 0.0);
    }

    #[test]
    fn rejects_single_point() {
        let err = CubicSpline::from_points(&[0.0], &[1.0]).unwrap_err();
        assert_eq!(err.code(), "interpolation.too_few_points");
    }

    #[test]
    fn rejects_extrapolation() {
        let s = CubicSpline::from_points(&[0.0, 1.0, 2.0], &[0.0, 1.0, 0.0]).unwrap();
        assert_eq!(
            s.value_at(2.5).unwrap_err().code(),
            "interpolation.out_of_domain"
        );
    }
}
