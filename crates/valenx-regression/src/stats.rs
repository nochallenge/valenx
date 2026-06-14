//! Descriptive statistics shared by the regression fits.
//!
//! All routines here are the textbook population (biased, `1/n`)
//! estimators. They are the building blocks the simple-linear fit uses to
//! express the slope as `cov(x, y) / var(x)`, and they back the Pearson
//! correlation coefficient. Because the simple-linear slope is a *ratio*
//! of a covariance to a variance, the choice of `1/n` versus `1/(n - 1)`
//! cancels — so the ordinary-least-squares slope is identical either way.
//!
//! The Pearson [`correlation`] coefficient is likewise scale- and
//! normalisation-invariant, so these biased moments give exactly the same
//! `r` as the sample (`1/(n - 1)`) moments would.

use crate::error::{RegressionError, Result};

/// Arithmetic mean of a non-empty slice.
///
/// # Errors
///
/// Returns [`RegressionError::TooFewPoints`] if `values` is empty.
pub fn mean(values: &[f64]) -> Result<f64> {
    if values.is_empty() {
        return Err(RegressionError::too_few_points(1, 0));
    }
    let sum: f64 = values.iter().sum();
    Ok(sum / values.len() as f64)
}

/// Population variance `var(x) = mean((x - mean_x)^2)` (the `1/n`
/// estimator).
///
/// A return value of exactly `0.0` means every element is identical.
///
/// # Errors
///
/// Returns [`RegressionError::TooFewPoints`] if `values` is empty.
pub fn variance(values: &[f64]) -> Result<f64> {
    let m = mean(values)?;
    let n = values.len() as f64;
    let ss: f64 = values.iter().map(|&v| (v - m) * (v - m)).sum();
    Ok(ss / n)
}

/// Population covariance `cov(x, y) = mean((x - mean_x)(y - mean_y))` (the
/// `1/n` estimator).
///
/// # Errors
///
/// Returns [`RegressionError::LengthMismatch`] if the slices differ in
/// length, or [`RegressionError::TooFewPoints`] if they are empty.
pub fn covariance(xs: &[f64], ys: &[f64]) -> Result<f64> {
    if xs.len() != ys.len() {
        return Err(RegressionError::length_mismatch(xs.len(), ys.len()));
    }
    if xs.is_empty() {
        return Err(RegressionError::too_few_points(1, 0));
    }
    let mx = mean(xs)?;
    let my = mean(ys)?;
    let n = xs.len() as f64;
    let acc: f64 = xs
        .iter()
        .zip(ys.iter())
        .map(|(&x, &y)| (x - mx) * (y - my))
        .sum();
    Ok(acc / n)
}

/// Population standard deviation `sqrt(var(x))`.
///
/// # Errors
///
/// Returns [`RegressionError::TooFewPoints`] if `values` is empty.
pub fn std_dev(values: &[f64]) -> Result<f64> {
    Ok(variance(values)?.sqrt())
}

/// Pearson product-moment correlation coefficient
/// `r = cov(x, y) / (sd_x * sd_y)`, clamped to the closed interval
/// `[-1, 1]` to absorb floating-point overshoot.
///
/// `r` is `+1` for a perfectly increasing line, `-1` for a perfectly
/// decreasing line, and `0` when the predictor and response are linearly
/// unrelated.
///
/// # Errors
///
/// Returns [`RegressionError::LengthMismatch`] if the slices differ in
/// length, [`RegressionError::TooFewPoints`] if fewer than two points are
/// supplied, or [`RegressionError::Degenerate`] if either variable has
/// zero variance (a constant series, for which correlation is undefined).
pub fn correlation(xs: &[f64], ys: &[f64]) -> Result<f64> {
    if xs.len() != ys.len() {
        return Err(RegressionError::length_mismatch(xs.len(), ys.len()));
    }
    if xs.len() < 2 {
        return Err(RegressionError::too_few_points(2, xs.len()));
    }
    let sx = std_dev(xs)?;
    let sy = std_dev(ys)?;
    if sx == 0.0 || sy == 0.0 {
        return Err(RegressionError::degenerate(
            "correlation undefined: a variable has zero variance",
        ));
    }
    let r = covariance(xs, ys)? / (sx * sy);
    Ok(r.clamp(-1.0, 1.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    #[test]
    fn mean_of_known_values() {
        let got = mean(&[1.0, 2.0, 3.0, 4.0]).unwrap();
        assert!((got - 2.5).abs() < EPS, "got {got}");
    }

    #[test]
    fn mean_empty_errors() {
        assert!(mean(&[]).is_err());
    }

    #[test]
    fn variance_of_known_values() {
        // population variance of {2, 4, 4, 4, 5, 5, 7, 9} = 4.
        let got = variance(&[2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]).unwrap();
        assert!((got - 4.0).abs() < EPS, "got {got}");
    }

    #[test]
    fn variance_constant_is_zero() {
        let got = variance(&[3.0, 3.0, 3.0]).unwrap();
        assert!(got.abs() < EPS, "got {got}");
    }

    #[test]
    fn std_dev_matches_sqrt_variance() {
        let data = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let got = std_dev(&data).unwrap();
        assert!((got - 2.0).abs() < EPS, "got {got}");
    }

    #[test]
    fn covariance_with_self_equals_variance() {
        let data = [1.0, 2.0, 3.0, 7.0, 11.0];
        let cov = covariance(&data, &data).unwrap();
        let var = variance(&data).unwrap();
        assert!((cov - var).abs() < EPS, "cov {cov} var {var}");
    }

    #[test]
    fn correlation_perfect_positive_line() {
        let xs = [1.0, 2.0, 3.0, 4.0, 5.0];
        let ys = [3.0, 5.0, 7.0, 9.0, 11.0]; // y = 2x + 1
        let r = correlation(&xs, &ys).unwrap();
        assert!((r - 1.0).abs() < 1e-9, "got {r}");
    }

    #[test]
    fn correlation_perfect_negative_line() {
        let xs = [1.0, 2.0, 3.0, 4.0, 5.0];
        let ys = [10.0, 8.0, 6.0, 4.0, 2.0]; // y = -2x + 12
        let r = correlation(&xs, &ys).unwrap();
        assert!((r + 1.0).abs() < 1e-9, "got {r}");
    }

    #[test]
    fn correlation_uncorrelated_is_zero() {
        // Symmetric V: x and y have zero covariance by construction.
        let xs = [-2.0, -1.0, 0.0, 1.0, 2.0];
        let ys = [4.0, 1.0, 0.0, 1.0, 4.0];
        let r = correlation(&xs, &ys).unwrap();
        assert!(r.abs() < 1e-12, "got {r}");
    }

    #[test]
    fn correlation_in_unit_interval() {
        let xs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let ys = [2.1, 3.9, 6.2, 7.8, 10.1, 11.7];
        let r = correlation(&xs, &ys).unwrap();
        assert!((-1.0..=1.0).contains(&r), "got {r}");
    }

    #[test]
    fn correlation_zero_variance_errors() {
        let xs = [1.0, 1.0, 1.0, 1.0];
        let ys = [1.0, 2.0, 3.0, 4.0];
        assert!(correlation(&xs, &ys).is_err());
    }

    #[test]
    fn correlation_length_mismatch_errors() {
        assert!(correlation(&[1.0, 2.0], &[1.0]).is_err());
    }
}
