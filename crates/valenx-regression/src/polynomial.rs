//! Polynomial least-squares regression via the normal equations.
//!
//! ## Model
//!
//! Fits `y = c0 + c1*x + c2*x^2 + ... + c_d*x^d` of a chosen `degree`
//! `d`, minimising the sum of squared vertical residuals. With the
//! Vandermonde design matrix `X` (column `j` is `x^j`), the
//! least-squares coefficients solve the normal equations
//!
//! ```text
//! (X^T X) c = X^T y
//! ```
//!
//! The `(d + 1) x (d + 1)` symmetric system `X^T X` is assembled directly
//! and solved with an LU factorisation from `nalgebra`. Degree `1`
//! reproduces the simple-linear fit; the coefficient vector is then
//! `[intercept, slope]`.
//!
//! Conditioning note: the monomial (Vandermonde) basis becomes
//! ill-conditioned for high degrees over wide `x` ranges. This is a
//! research/educational implementation — for large degrees prefer an
//! orthogonal-polynomial basis. A rank-deficient system surfaces as
//! [`crate::error::RegressionError::Degenerate`].

use crate::error::{RegressionError, Result};
use nalgebra::{DMatrix, DVector};
use serde::{Deserialize, Serialize};

/// Result of a polynomial least-squares fit.
///
/// `coefficients[j]` multiplies `x^j`, so element `0` is the constant
/// term and the last element multiplies `x^degree`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PolynomialFit {
    /// Coefficients in ascending power order: `[c0, c1, ..., c_degree]`.
    pub coefficients: Vec<f64>,
    /// Coefficient of determination `R^2` in `[0, 1]`.
    pub r_squared: f64,
}

impl PolynomialFit {
    /// Degree of the fitted polynomial (`coefficients.len() - 1`).
    pub fn degree(&self) -> usize {
        self.coefficients.len() - 1
    }

    /// Evaluate the fitted polynomial at `x` using Horner's scheme.
    pub fn predict(&self, x: f64) -> f64 {
        let mut acc = 0.0;
        for &c in self.coefficients.iter().rev() {
            acc = acc * x + c;
        }
        acc
    }

    /// Residuals `y_i - predict(x_i)` for every observation, in input
    /// order.
    ///
    /// # Errors
    ///
    /// Returns [`RegressionError::LengthMismatch`] if `xs` and `ys` differ
    /// in length.
    pub fn residuals(&self, xs: &[f64], ys: &[f64]) -> Result<Vec<f64>> {
        if xs.len() != ys.len() {
            return Err(RegressionError::length_mismatch(xs.len(), ys.len()));
        }
        Ok(xs
            .iter()
            .zip(ys.iter())
            .map(|(&x, &y)| y - self.predict(x))
            .collect())
    }
}

/// Fit a degree-`degree` polynomial by least squares through the normal
/// equations.
///
/// `degree` `1` is a straight line; the resulting `coefficients` are then
/// `[intercept, slope]`, matching [`crate::linear::fit`].
///
/// # Errors
///
/// - [`RegressionError::LengthMismatch`] if `xs` and `ys` differ in
///   length.
/// - [`RegressionError::TooFewPoints`] if fewer than `degree + 1` points
///   are supplied (a degree-`degree` polynomial has `degree + 1`
///   coefficients to determine).
/// - [`RegressionError::Degenerate`] if the normal matrix `X^T X` is
///   rank-deficient and the LU solve cannot factor it (for example, every
///   `x` identical, or duplicated abscissae that leave the system
///   singular).
pub fn fit(xs: &[f64], ys: &[f64], degree: usize) -> Result<PolynomialFit> {
    if xs.len() != ys.len() {
        return Err(RegressionError::length_mismatch(xs.len(), ys.len()));
    }
    let n_coeffs = degree + 1;
    if xs.len() < n_coeffs {
        return Err(RegressionError::too_few_points(n_coeffs, xs.len()));
    }

    // Vandermonde design matrix X: rows = observations, cols = powers
    // 0..=degree.
    let n_rows = xs.len();
    let mut design = DMatrix::<f64>::zeros(n_rows, n_coeffs);
    for (i, &x) in xs.iter().enumerate() {
        let mut power = 1.0;
        for j in 0..n_coeffs {
            design[(i, j)] = power;
            power *= x;
        }
    }

    let y = DVector::<f64>::from_row_slice(ys);

    // Normal equations: (X^T X) c = X^T y.
    let xt = design.transpose();
    let normal = &xt * &design;
    let rhs = &xt * &y;

    let solution = normal
        .lu()
        .solve(&rhs)
        .ok_or_else(|| RegressionError::degenerate("normal matrix X^T X is singular (rank-deficient design); reduce the degree or supply more distinct x values"))?;

    let coefficients: Vec<f64> = solution.iter().copied().collect();

    let r_squared = r_squared_for(xs, ys, &coefficients);

    Ok(PolynomialFit {
        coefficients,
        r_squared,
    })
}

/// Coefficient of determination for the given ascending-power
/// `coefficients` against the data.
///
/// `R^2 = 1 - SS_res / SS_tot`, clamped to `[0, 1]`. A constant response
/// (`SS_tot == 0`) is treated as a perfect fit (`1.0`).
fn r_squared_for(xs: &[f64], ys: &[f64], coefficients: &[f64]) -> f64 {
    let n = ys.len() as f64;
    let mean_y: f64 = ys.iter().sum::<f64>() / n;
    let mut ss_res = 0.0;
    let mut ss_tot = 0.0;
    for (&x, &y) in xs.iter().zip(ys.iter()) {
        // Horner evaluation of the polynomial at x.
        let mut pred = 0.0;
        for &c in coefficients.iter().rev() {
            pred = pred * x + c;
        }
        let res = y - pred;
        ss_res += res * res;
        let dev = y - mean_y;
        ss_tot += dev * dev;
    }
    if ss_tot == 0.0 {
        return 1.0;
    }
    (1.0 - ss_res / ss_tot).clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linear;

    const EPS: f64 = 1e-7;

    #[test]
    fn exact_quadratic_is_recovered() {
        // y = 1 + 2x + 3x^2, exactly.
        let xs = [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0];
        let ys: Vec<f64> = xs.iter().map(|&x| 1.0 + 2.0 * x + 3.0 * x * x).collect();
        let fit = fit(&xs, &ys, 2).unwrap();
        assert!(
            (fit.coefficients[0] - 1.0).abs() < EPS,
            "c0 {:?}",
            fit.coefficients
        );
        assert!(
            (fit.coefficients[1] - 2.0).abs() < EPS,
            "c1 {:?}",
            fit.coefficients
        );
        assert!(
            (fit.coefficients[2] - 3.0).abs() < EPS,
            "c2 {:?}",
            fit.coefficients
        );
    }

    #[test]
    fn exact_quadratic_has_unit_r_squared_and_zero_residuals() {
        let xs = [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0];
        let ys: Vec<f64> = xs.iter().map(|&x| 1.0 + 2.0 * x + 3.0 * x * x).collect();
        let fit = fit(&xs, &ys, 2).unwrap();
        assert!((fit.r_squared - 1.0).abs() < EPS, "r2 {}", fit.r_squared);
        for r in fit.residuals(&xs, &ys).unwrap() {
            assert!(r.abs() < EPS, "residual {r}");
        }
    }

    #[test]
    fn exact_cubic_is_recovered() {
        // y = -1 + 0.5x - 2x^2 + x^3.
        let xs = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0];
        let ys: Vec<f64> = xs
            .iter()
            .map(|&x| -1.0 + 0.5 * x - 2.0 * x * x + x * x * x)
            .collect();
        let fit = fit(&xs, &ys, 3).unwrap();
        let expected = [-1.0, 0.5, -2.0, 1.0];
        for (got, exp) in fit.coefficients.iter().zip(expected.iter()) {
            assert!((got - exp).abs() < 1e-6, "coeffs {:?}", fit.coefficients);
        }
    }

    #[test]
    fn degree_one_matches_simple_linear() {
        let xs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let ys = [2.1, 3.9, 6.2, 7.8, 10.1, 11.7];
        let poly = fit(&xs, &ys, 1).unwrap();
        let lin = linear::fit(&xs, &ys).unwrap();
        // ascending-power: coefficients[0] = intercept, [1] = slope.
        assert!(
            (poly.coefficients[0] - lin.intercept).abs() < 1e-9,
            "intercept poly {} lin {}",
            poly.coefficients[0],
            lin.intercept
        );
        assert!(
            (poly.coefficients[1] - lin.slope).abs() < 1e-9,
            "slope poly {} lin {}",
            poly.coefficients[1],
            lin.slope
        );
        assert!(
            (poly.r_squared - lin.r_squared).abs() < 1e-9,
            "r2 poly {} lin {}",
            poly.r_squared,
            lin.r_squared
        );
    }

    #[test]
    fn predict_uses_horner_correctly() {
        // coefficients [1, 2, 3] -> 1 + 2x + 3x^2; at x = 2 -> 17.
        let fit = PolynomialFit {
            coefficients: vec![1.0, 2.0, 3.0],
            r_squared: 1.0,
        };
        assert!(
            (fit.predict(2.0) - 17.0).abs() < 1e-12,
            "{}",
            fit.predict(2.0)
        );
        assert_eq!(fit.degree(), 2);
    }

    #[test]
    fn r_squared_in_unit_interval_for_noisy_quadratic() {
        let xs = [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0];
        let ys = [5.1, 0.2, 0.9, 4.3, 9.1, 16.2, 24.8];
        let fit = fit(&xs, &ys, 2).unwrap();
        assert!((0.0..=1.0).contains(&fit.r_squared), "r2 {}", fit.r_squared);
    }

    #[test]
    fn too_few_points_for_degree_errors() {
        // A cubic needs 4 points; supply 3.
        let xs = [1.0, 2.0, 3.0];
        let ys = [1.0, 4.0, 9.0];
        assert!(fit(&xs, &ys, 3).is_err());
    }

    #[test]
    fn zero_variance_predictor_is_degenerate() {
        // All x identical -> Vandermonde columns collinear -> singular.
        let xs = [2.0, 2.0, 2.0, 2.0];
        let ys = [1.0, 2.0, 3.0, 4.0];
        assert!(fit(&xs, &ys, 2).is_err());
    }

    #[test]
    fn length_mismatch_errors() {
        assert!(fit(&[1.0, 2.0, 3.0], &[1.0, 2.0], 1).is_err());
    }

    #[test]
    fn degree_zero_fits_the_mean() {
        // A degree-0 polynomial is the best-fit constant = mean(y).
        let xs = [1.0, 2.0, 3.0, 4.0];
        let ys = [2.0, 4.0, 6.0, 8.0];
        let fit = fit(&xs, &ys, 0).unwrap();
        assert!(
            (fit.coefficients[0] - 5.0).abs() < 1e-9,
            "{:?}",
            fit.coefficients
        );
    }
}
