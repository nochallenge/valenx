//! Simple (one-predictor) ordinary-least-squares linear regression.
//!
//! ## Model
//!
//! Fits `y = slope * x + intercept` by minimising the sum of squared
//! vertical residuals. The closed-form ordinary-least-squares solution is
//!
//! ```text
//! slope     = cov(x, y) / var(x)
//! intercept = mean_y - slope * mean_x
//! ```
//!
//! Because the slope is a ratio of a covariance to a variance, the `1/n`
//! population moments from [`crate::stats`] give exactly the same answer
//! as the `1/(n - 1)` sample moments. The fitted line therefore always
//! passes through the centroid `(mean_x, mean_y)`.
//!
//! The goodness-of-fit `R^2 = 1 - SS_res / SS_tot` is the fraction of the
//! response variance the line explains: `1.0` for a perfect fit (all
//! residuals zero), `0.0` when the line does no better than the
//! horizontal mean line. For this one-predictor model `R^2` equals the
//! square of the Pearson correlation coefficient `r` (see
//! [`crate::stats::correlation`]).

use crate::error::{RegressionError, Result};
use crate::stats::{covariance, mean, variance};
use serde::{Deserialize, Serialize};

/// Result of a simple least-squares linear fit `y = slope * x +
/// intercept`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LinearFit {
    /// Fitted slope `cov(x, y) / var(x)`.
    pub slope: f64,
    /// Fitted intercept `mean_y - slope * mean_x`.
    pub intercept: f64,
    /// Coefficient of determination `R^2` in `[0, 1]`.
    pub r_squared: f64,
}

impl LinearFit {
    /// Predict the response `y` at a predictor value `x` from the fitted
    /// line.
    pub fn predict(&self, x: f64) -> f64 {
        self.slope * x + self.intercept
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

/// Fit `y = slope * x + intercept` by ordinary least squares.
///
/// The slope is the closed-form `cov(x, y) / var(x)`; the intercept then
/// forces the line through the centroid `(mean_x, mean_y)`. `R^2` is
/// computed from the resulting residuals and clamped to `[0, 1]` to
/// absorb floating-point overshoot at a perfect fit.
///
/// # Errors
///
/// - [`RegressionError::LengthMismatch`] if `xs` and `ys` differ in
///   length.
/// - [`RegressionError::TooFewPoints`] if fewer than two points are
///   supplied (a line is underdetermined by one point).
/// - [`RegressionError::Degenerate`] if the predictor has zero variance
///   (every `x` identical), making the slope infinite.
pub fn fit(xs: &[f64], ys: &[f64]) -> Result<LinearFit> {
    if xs.len() != ys.len() {
        return Err(RegressionError::length_mismatch(xs.len(), ys.len()));
    }
    if xs.len() < 2 {
        return Err(RegressionError::too_few_points(2, xs.len()));
    }

    let var_x = variance(xs)?;
    if var_x == 0.0 {
        return Err(RegressionError::degenerate(
            "predictor has zero variance: every x is identical, slope is undefined",
        ));
    }

    let slope = covariance(xs, ys)? / var_x;
    let mean_x = mean(xs)?;
    let mean_y = mean(ys)?;
    let intercept = mean_y - slope * mean_x;

    let r_squared = r_squared_for(xs, ys, slope, intercept, mean_y);

    Ok(LinearFit {
        slope,
        intercept,
        r_squared,
    })
}

/// Coefficient of determination for an explicit `slope` / `intercept`
/// line against the data, given the precomputed `mean_y`.
///
/// `R^2 = 1 - SS_res / SS_tot`, clamped to `[0, 1]`. When the response is
/// constant (`SS_tot == 0`) the residual sum is also zero for any
/// horizontal line through it, so the fit is treated as perfect and `1.0`
/// is returned.
fn r_squared_for(xs: &[f64], ys: &[f64], slope: f64, intercept: f64, mean_y: f64) -> f64 {
    let mut ss_res = 0.0;
    let mut ss_tot = 0.0;
    for (&x, &y) in xs.iter().zip(ys.iter()) {
        let pred = slope * x + intercept;
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
    use crate::stats::correlation;

    const EPS: f64 = 1e-9;

    #[test]
    fn perfect_line_recovers_slope_intercept() {
        // y = 2x + 1, exactly.
        let xs = [0.0, 1.0, 2.0, 3.0, 4.0];
        let ys = [1.0, 3.0, 5.0, 7.0, 9.0];
        let fit = fit(&xs, &ys).unwrap();
        assert!((fit.slope - 2.0).abs() < EPS, "slope {}", fit.slope);
        assert!(
            (fit.intercept - 1.0).abs() < EPS,
            "intercept {}",
            fit.intercept
        );
    }

    #[test]
    fn perfect_line_r_squared_is_one() {
        let xs = [0.0, 1.0, 2.0, 3.0, 4.0];
        let ys = [1.0, 3.0, 5.0, 7.0, 9.0];
        let fit = fit(&xs, &ys).unwrap();
        assert!((fit.r_squared - 1.0).abs() < EPS, "r2 {}", fit.r_squared);
    }

    #[test]
    fn perfect_line_residuals_are_zero() {
        let xs = [0.0, 1.0, 2.0, 3.0, 4.0];
        let ys = [1.0, 3.0, 5.0, 7.0, 9.0];
        let fit = fit(&xs, &ys).unwrap();
        let res = fit.residuals(&xs, &ys).unwrap();
        for r in res {
            assert!(r.abs() < EPS, "residual {r}");
        }
    }

    #[test]
    fn slope_equals_cov_over_var() {
        let xs = [1.0, 2.0, 4.0, 5.0, 7.0, 9.0];
        let ys = [2.0, 3.0, 3.0, 6.0, 5.0, 8.0];
        let fit = fit(&xs, &ys).unwrap();
        let expected = covariance(&xs, &ys).unwrap() / variance(&xs).unwrap();
        assert!(
            (fit.slope - expected).abs() < EPS,
            "slope {} expected {expected}",
            fit.slope
        );
    }

    #[test]
    fn line_passes_through_centroid() {
        let xs = [1.0, 2.0, 4.0, 5.0, 7.0, 9.0];
        let ys = [2.0, 3.0, 3.0, 6.0, 5.0, 8.0];
        let fit = fit(&xs, &ys).unwrap();
        let mx = mean(&xs).unwrap();
        let my = mean(&ys).unwrap();
        assert!((fit.predict(mx) - my).abs() < EPS, "centroid not on line");
    }

    #[test]
    fn r_squared_in_unit_interval_for_noisy_data() {
        let xs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let ys = [2.3, 2.9, 4.1, 3.8, 5.2, 6.1, 6.8, 8.4];
        let fit = fit(&xs, &ys).unwrap();
        assert!((0.0..=1.0).contains(&fit.r_squared), "r2 {}", fit.r_squared);
    }

    #[test]
    fn r_squared_equals_r_squared_for_simple_linear() {
        let xs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let ys = [2.3, 2.9, 4.1, 3.8, 5.2, 6.1, 6.8, 8.4];
        let fit = fit(&xs, &ys).unwrap();
        let r = correlation(&xs, &ys).unwrap();
        assert!(
            (fit.r_squared - r * r).abs() < EPS,
            "r2 {} r^2 {}",
            fit.r_squared,
            r * r
        );
    }

    #[test]
    fn horizontal_data_has_zero_slope() {
        // Constant y: best fit is a flat line, slope 0.
        let xs = [1.0, 2.0, 3.0, 4.0, 5.0];
        let ys = [7.0, 7.0, 7.0, 7.0, 7.0];
        let fit = fit(&xs, &ys).unwrap();
        assert!(fit.slope.abs() < EPS, "slope {}", fit.slope);
        assert!(
            (fit.intercept - 7.0).abs() < EPS,
            "intercept {}",
            fit.intercept
        );
    }

    #[test]
    fn residuals_sum_to_zero_for_ols() {
        // A property of OLS with an intercept: residuals sum to zero.
        let xs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let ys = [2.1, 3.9, 6.2, 7.8, 10.1, 11.7];
        let fit = fit(&xs, &ys).unwrap();
        let sum: f64 = fit.residuals(&xs, &ys).unwrap().iter().sum();
        assert!(sum.abs() < 1e-9, "residual sum {sum}");
    }

    #[test]
    fn decreasing_line_has_negative_slope() {
        let xs = [0.0, 1.0, 2.0, 3.0];
        let ys = [10.0, 7.0, 4.0, 1.0]; // y = -3x + 10
        let fit = fit(&xs, &ys).unwrap();
        assert!((fit.slope + 3.0).abs() < EPS, "slope {}", fit.slope);
        assert!(
            (fit.intercept - 10.0).abs() < EPS,
            "intercept {}",
            fit.intercept
        );
    }

    #[test]
    fn zero_variance_predictor_errors() {
        let xs = [2.0, 2.0, 2.0, 2.0];
        let ys = [1.0, 2.0, 3.0, 4.0];
        assert!(fit(&xs, &ys).is_err());
    }

    #[test]
    fn too_few_points_errors() {
        assert!(fit(&[1.0], &[2.0]).is_err());
    }

    #[test]
    fn length_mismatch_errors() {
        assert!(fit(&[1.0, 2.0, 3.0], &[1.0, 2.0]).is_err());
    }

    #[test]
    fn residuals_length_mismatch_errors() {
        let fit = LinearFit {
            slope: 1.0,
            intercept: 0.0,
            r_squared: 1.0,
        };
        assert!(fit.residuals(&[1.0, 2.0], &[1.0]).is_err());
    }
}
