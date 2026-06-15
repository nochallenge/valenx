//! Lagrange polynomial interpolation.
//!
//! Through `n + 1` distinct knots there is a unique polynomial of
//! degree at most `n` that passes through all of them. This module
//! evaluates that polynomial using the numerically stable
//! barycentric form (Berrut and Trefethen, 2004) rather than the
//! naive product-of-basis-functions form.
//!
//! ## Barycentric weights
//!
//! For knots `x_0 .. x_n` the weights are
//! `w_j = 1 / prod_{k != j} (x_j - x_k)`. Then for any query `x` not
//! equal to a knot,
//! `p(x) = (sum_j w_j y_j / (x - x_j)) / (sum_j w_j / (x - x_j))`,
//! and `p(x_j) = y_j` exactly at the knots.

use crate::error::InterpolationError;
use crate::points::DataPoints;

/// A Lagrange interpolating polynomial in barycentric form.
///
/// Build with [`LagrangePolynomial::new`]; query with
/// [`LagrangePolynomial::value_at`]. [`LagrangePolynomial::degree`]
/// reports the polynomial degree, which equals `n` for `n + 1` knots.
#[derive(Clone, Debug, PartialEq)]
pub struct LagrangePolynomial {
    data: DataPoints,
    weights: Vec<f64>,
}

impl LagrangePolynomial {
    /// Build the interpolating polynomial from validated
    /// [`DataPoints`].
    ///
    /// The barycentric weights are precomputed once here so each
    /// later evaluation is `O(n)`.
    pub fn new(data: DataPoints) -> Self {
        let xs = data.xs();
        let n = xs.len();
        let mut weights = vec![1.0_f64; n];
        for j in 0..n {
            let mut prod = 1.0_f64;
            for (k, &xk) in xs.iter().enumerate() {
                if k != j {
                    prod *= xs[j] - xk;
                }
            }
            weights[j] = 1.0 / prod;
        }
        Self { data, weights }
    }

    /// Build directly from `x` / `y` slices.
    ///
    /// # Errors
    ///
    /// Propagates any [`DataPoints::new`] validation error (distinct,
    /// sorted, finite abscissae are required, exactly as the
    /// polynomial's existence demands).
    pub fn from_points(xs: &[f64], ys: &[f64]) -> Result<Self, InterpolationError> {
        Ok(Self::new(DataPoints::new(xs, ys)?))
    }

    /// The underlying knots.
    pub fn data(&self) -> &DataPoints {
        &self.data
    }

    /// Degree of the interpolating polynomial: `n` for `n + 1` knots.
    pub fn degree(&self) -> usize {
        self.data.len() - 1
    }

    /// Evaluate the polynomial at `at`.
    ///
    /// Unlike the piecewise methods, the Lagrange polynomial is a
    /// single global polynomial and is mathematically defined for all
    /// real `at`; for consistency with the other interpolants this
    /// crate still restricts queries to `[x_first, x_last]` and
    /// rejects out-of-domain `at`.
    ///
    /// # Errors
    ///
    /// Returns [`InterpolationError::OutOfDomain`] if `at` lies
    /// outside `[x_first, x_last]`.
    pub fn value_at(&self, at: f64) -> Result<f64, InterpolationError> {
        self.data.check_domain(at)?;
        Ok(self.eval_unchecked(at))
    }

    /// Evaluate the polynomial at any real `at`, including outside the
    /// knot span.
    ///
    /// The Lagrange polynomial is globally defined, so this method
    /// performs no domain check. Use it deliberately when polynomial
    /// extrapolation is wanted; prefer [`LagrangePolynomial::value_at`]
    /// for the domain-checked behaviour shared with the other
    /// interpolants.
    pub fn value_at_unchecked(&self, at: f64) -> f64 {
        self.eval_unchecked(at)
    }

    /// Core barycentric evaluation, with the exact-knot special case.
    fn eval_unchecked(&self, at: f64) -> f64 {
        let xs = self.data.xs();
        let ys = self.data.ys();
        // Exact hit on a knot: return its ordinate (and avoid a
        // division by zero in the barycentric quotient).
        for (j, &xj) in xs.iter().enumerate() {
            if at == xj {
                return ys[j];
            }
        }
        let mut numer = 0.0_f64;
        let mut denom = 0.0_f64;
        for j in 0..xs.len() {
            let term = self.weights[j] / (at - xs[j]);
            numer += term * ys[j];
            denom += term;
        }
        numer / denom
    }
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
        let xs = [-1.0, 0.0, 2.0, 3.0, 5.0];
        let ys = [4.0, 1.0, -2.0, 6.0, 0.5];
        let p = LagrangePolynomial::from_points(&xs, &ys).unwrap();
        for (&x, &y) in xs.iter().zip(ys.iter()) {
            close(p.value_at(x).unwrap(), y);
        }
    }

    #[test]
    fn degree_is_n_for_n_plus_one_points() {
        let p = LagrangePolynomial::from_points(&[0.0, 1.0, 2.0], &[0.0, 1.0, 4.0]).unwrap();
        assert_eq!(p.degree(), 2);
        let q = LagrangePolynomial::from_points(&[0.0, 1.0, 2.0, 3.0, 4.0], &[1.0; 5]).unwrap();
        assert_eq!(q.degree(), 4);
    }

    #[test]
    fn recovers_the_generating_quadratic() {
        // Sample y = 2x^2 - 3x + 1 at three points; the unique degree-2
        // interpolant must reproduce it everywhere on the span.
        let f = |x: f64| 2.0 * x * x - 3.0 * x + 1.0;
        let xs = [0.0, 1.0, 2.0];
        let ys = [f(0.0), f(1.0), f(2.0)];
        let p = LagrangePolynomial::from_points(&xs, &ys).unwrap();
        for k in 0..=20 {
            let x = 2.0 * (k as f64 / 20.0);
            close(p.value_at(x).unwrap(), f(x));
        }
    }

    #[test]
    fn recovers_a_generating_cubic() {
        // y = x^3 - x sampled at four points -> exact cubic recovery.
        let f = |x: f64| x * x * x - x;
        let xs = [-2.0, -1.0, 1.0, 2.0];
        let ys = [f(-2.0), f(-1.0), f(1.0), f(2.0)];
        let p = LagrangePolynomial::from_points(&xs, &ys).unwrap();
        for k in 0..=20 {
            let x = -2.0 + 4.0 * (k as f64 / 20.0);
            close(p.value_at(x).unwrap(), f(x));
        }
    }

    #[test]
    fn two_points_give_the_straight_line() {
        // Degree 1 Lagrange == linear interpolation.
        let p = LagrangePolynomial::from_points(&[0.0, 4.0], &[1.0, 9.0]).unwrap();
        assert_eq!(p.degree(), 1);
        close(p.value_at(2.0).unwrap(), 5.0);
        close(p.value_at(1.0).unwrap(), 3.0);
    }

    #[test]
    fn constant_data_gives_constant_output() {
        let p = LagrangePolynomial::from_points(&[0.0, 1.0, 2.0, 3.0], &[7.0; 4]).unwrap();
        for k in 0..=15 {
            let x = 3.0 * (k as f64 / 15.0);
            close(p.value_at(x).unwrap(), 7.0);
        }
    }

    #[test]
    fn unchecked_extrapolation_extends_the_polynomial() {
        // y = x^2 at three knots; unchecked eval at x=3 (outside span
        // [0,2]) must still give 9.
        let p = LagrangePolynomial::from_points(&[0.0, 1.0, 2.0], &[0.0, 1.0, 4.0]).unwrap();
        close(p.value_at_unchecked(3.0), 9.0);
        // The domain-checked entry point rejects the same query.
        assert!(p.value_at(3.0).is_err());
    }
}
