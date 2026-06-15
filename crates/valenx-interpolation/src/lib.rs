//! # valenx-interpolation
//!
//! One-dimensional interpolation through scattered `(x, y)` samples.
//!
//! ## What
//!
//! Three classic interpolants, each built from a validated set of
//! knots and queried by abscissa:
//!
//! - [`LinearInterpolator`] — piecewise-linear (`C0`); on each
//!   segment the interpolant is the straight line joining the two
//!   bracketing knots.
//! - [`LagrangePolynomial`] — the unique global polynomial of degree
//!   at most `n` through `n + 1` knots, evaluated in the
//!   numerically stable barycentric form.
//! - [`CubicSpline`] — the natural cubic spline: cubic on each
//!   segment, twice continuously differentiable (`C2`), with zero
//!   second derivative at both ends.
//!
//! All three share the validated [`DataPoints`] container, so the
//! "strictly ascending, finite, distinct abscissae" precondition is
//! enforced once at construction.
//!
//! ## Model
//!
//! - Linear: `S(x) = y_i + (x - x_i)/(x_{i+1} - x_i) * (y_{i+1} - y_i)`
//!   on the segment containing `x`.
//! - Lagrange (barycentric weights `w_j = 1 / prod_{k != j}(x_j - x_k)`):
//!   `p(x) = (sum_j w_j y_j/(x - x_j)) / (sum_j w_j/(x - x_j))`, with the
//!   knots returned exactly.
//! - Cubic spline: the knot second derivatives (moments) `M_i` solve a
//!   symmetric, diagonally dominant tridiagonal system (Thomas
//!   algorithm) derived from `C2` continuity plus the natural end
//!   conditions `M_0 = M_n = 0`. Value, first, and second derivatives
//!   then come from the closed-form per-segment cubic.
//!
//! Every interpolant passes through all data points exactly; constant
//! data yields a constant interpolant; and queries are restricted to
//! the closed data domain `[x_first, x_last]` (extrapolation is
//! rejected rather than silently returned).
//!
//! ## Honest scope
//!
//! Research/educational grade. These are textbook closed-form and
//! numerical models (piecewise-linear, barycentric Lagrange, natural
//! cubic spline via a tridiagonal solve) intended for learning,
//! prototyping, and well-conditioned data. This crate is NOT a
//! clinical, medical, or production engineering tool, and makes no
//! safety or fitness-for-purpose guarantee. In particular, high-degree
//! Lagrange interpolation of equispaced data is subject to Runge
//! oscillation, and no interpolant here extrapolates beyond the data
//! span. Validate against your own ground truth before relying on any
//! result.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod lagrange;
pub mod linear;
pub mod points;
pub mod spline;

pub use error::{ErrorCategory, InterpolationError};
pub use lagrange::LagrangePolynomial;
pub use linear::LinearInterpolator;
pub use points::DataPoints;
pub use spline::CubicSpline;

#[cfg(test)]
mod integration_tests {
    //! Cross-method ground-truth checks exercising the public surface.

    use super::*;

    const EPS: f64 = 1e-9;

    fn close(a: f64, b: f64) {
        assert!((a - b).abs() < EPS, "expected {a} ~= {b}");
    }

    #[test]
    fn all_three_agree_at_the_knots() {
        let xs = [0.0, 1.0, 2.0, 3.0];
        let ys = [1.0, 4.0, 9.0, 16.0];
        let lin = LinearInterpolator::from_points(&xs, &ys).unwrap();
        let lag = LagrangePolynomial::from_points(&xs, &ys).unwrap();
        let spl = CubicSpline::from_points(&xs, &ys).unwrap();
        for (&x, &y) in xs.iter().zip(ys.iter()) {
            close(lin.value_at(x).unwrap(), y);
            close(lag.value_at(x).unwrap(), y);
            close(spl.value_at(x).unwrap(), y);
        }
    }

    #[test]
    fn all_three_reproduce_a_straight_line() {
        // y = -2x + 5 sampled at four knots: every method is exact.
        let f = |x: f64| -2.0 * x + 5.0;
        let xs = [0.0, 1.5, 2.0, 4.0];
        let ys: Vec<f64> = xs.iter().map(|&x| f(x)).collect();
        let lin = LinearInterpolator::from_points(&xs, &ys).unwrap();
        let lag = LagrangePolynomial::from_points(&xs, &ys).unwrap();
        let spl = CubicSpline::from_points(&xs, &ys).unwrap();
        for k in 0..=12 {
            let x = 4.0 * (k as f64 / 12.0);
            close(lin.value_at(x).unwrap(), f(x));
            close(lag.value_at(x).unwrap(), f(x));
            close(spl.value_at(x).unwrap(), f(x));
        }
    }

    #[test]
    fn shared_data_points_rejects_unsorted_for_every_method() {
        let bad_x = [0.0, 3.0, 1.0];
        let y = [0.0, 0.0, 0.0];
        assert!(LinearInterpolator::from_points(&bad_x, &y).is_err());
        assert!(LagrangePolynomial::from_points(&bad_x, &y).is_err());
        assert!(CubicSpline::from_points(&bad_x, &y).is_err());
    }
}
