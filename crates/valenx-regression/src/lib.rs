//! # valenx-regression
//!
//! Least-squares regression primitives for Valenx.
//!
//! ## What
//!
//! Closed-form and numerical curve fitting over 1-D data:
//!
//! - [`linear::fit`] — simple ordinary-least-squares line `y = slope * x +
//!   intercept`, returning a [`linear::LinearFit`] with `slope`,
//!   `intercept`, and the coefficient of determination `R^2`. The fit
//!   also yields per-point [`linear::LinearFit::residuals`] and a
//!   [`linear::LinearFit::predict`] evaluator.
//! - [`polynomial::fit`] — polynomial regression of any `degree` via the
//!   normal equations, returning a [`polynomial::PolynomialFit`] with
//!   ascending-power `coefficients` and `R^2`.
//! - [`stats::correlation`] — the Pearson product-moment correlation
//!   coefficient `r`, plus the supporting [`stats::mean`],
//!   [`stats::variance`], [`stats::covariance`], and [`stats::std_dev`]
//!   moments.
//!
//! ## Model
//!
//! The simple-linear fit uses the textbook closed form
//!
//! ```text
//! slope     = cov(x, y) / var(x)
//! intercept = mean_y - slope * mean_x
//! ```
//!
//! so the line always passes through the centroid `(mean_x, mean_y)`.
//! Polynomial fits assemble the Vandermonde design matrix `X` and solve
//! the normal equations `(X^T X) c = X^T y` with an LU factorisation from
//! `nalgebra`. Goodness of fit is `R^2 = 1 - SS_res / SS_tot`, clamped to
//! `[0, 1]`; for the one-predictor line `R^2` equals `r^2`, the square of
//! the Pearson correlation coefficient.
//!
//! ## Honest scope
//!
//! Research/educational grade. These are textbook closed-form / numerical
//! least-squares models — ordinary least squares with the monomial basis
//! and the population (`1/n`) moments. They are **not** a clinical,
//! medical, or production engineering tool, and deliberately omit the
//! machinery a statistical-inference package would provide: no weighted or
//! robust regression, no standard errors, confidence or prediction
//! intervals, hypothesis tests, or p-values; no regularisation; and no
//! orthogonal-polynomial basis, so the monomial normal equations become
//! ill-conditioned for high polynomial degrees over wide `x` ranges. Do
//! not use the outputs for decisions where those guarantees matter.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod linear;
pub mod polynomial;
pub mod stats;

pub use error::{ErrorCategory, RegressionError, Result};
pub use linear::{fit as fit_linear, LinearFit};
pub use polynomial::{fit as fit_polynomial, PolynomialFit};
pub use stats::{correlation, covariance, mean, std_dev, variance};

#[cfg(test)]
mod integration_tests {
    //! Cross-module ground-truth checks that exercise the public API the
    //! way a caller would, validating the relationships the individual
    //! module tests assume in isolation.

    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn perfect_line_end_to_end() {
        // y = 2x + 1: line recovered exactly, R^2 = 1, residuals 0,
        // r^2 = R^2.
        let xs = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let ys = [1.0, 3.0, 5.0, 7.0, 9.0, 11.0];

        let fit = fit_linear(&xs, &ys).unwrap();
        assert!((fit.slope - 2.0).abs() < EPS, "slope {}", fit.slope);
        assert!(
            (fit.intercept - 1.0).abs() < EPS,
            "intercept {}",
            fit.intercept
        );
        assert!((fit.r_squared - 1.0).abs() < EPS, "r2 {}", fit.r_squared);

        for r in fit.residuals(&xs, &ys).unwrap() {
            assert!(r.abs() < EPS, "residual {r}");
        }

        let r = correlation(&xs, &ys).unwrap();
        assert!((fit.r_squared - r * r).abs() < EPS, "R^2 != r^2");
    }

    #[test]
    fn polynomial_degree_one_agrees_with_linear() {
        let xs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        let ys = [2.2, 3.8, 6.1, 7.9, 10.2, 11.8, 14.1];

        let lin = fit_linear(&xs, &ys).unwrap();
        let poly = fit_polynomial(&xs, &ys, 1).unwrap();

        assert!((poly.coefficients[0] - lin.intercept).abs() < EPS);
        assert!((poly.coefficients[1] - lin.slope).abs() < EPS);
        assert!((poly.r_squared - lin.r_squared).abs() < EPS);
    }

    #[test]
    fn r_squared_always_in_unit_interval() {
        let xs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
        let ys = [3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0, 5.0, 3.0];
        let fit = fit_linear(&xs, &ys).unwrap();
        assert!((0.0..=1.0).contains(&fit.r_squared), "r2 {}", fit.r_squared);

        let poly = fit_polynomial(&xs, &ys, 3).unwrap();
        assert!(
            (0.0..=1.0).contains(&poly.r_squared),
            "r2 {}",
            poly.r_squared
        );
    }

    /// Compile-time check that the public fit structs implement serde's
    /// derive traits (the crate enables `serde`'s `derive` feature). No
    /// concrete format crate is pulled in, so this only asserts the bound
    /// is satisfied rather than performing an actual round-trip.
    #[test]
    fn fit_structs_are_serde() {
        fn assert_serde<T: serde::Serialize + serde::de::DeserializeOwned>() {}
        assert_serde::<LinearFit>();
        assert_serde::<PolynomialFit>();
    }
}
