//! Summary statistics over a sample of scalar outputs.
//!
//! These are the forward-propagation read-outs: given the model outputs for a
//! batch of sampled inputs, summarise the resulting distribution. Every
//! routine guards an **empty sample** — it returns `None` / [`UqError`] rather
//! than dividing by zero and producing `NaN`.
//!
//! Variance and standard deviation use the **unbiased (sample)** estimator
//! with denominator `n - 1` (Bessel's correction); they therefore require at
//! least two samples.

use crate::error::UqError;

/// Arithmetic mean of `data`. Returns `None` for an empty sample.
#[must_use]
pub fn mean(data: &[f64]) -> Option<f64> {
    if data.is_empty() {
        return None;
    }
    let sum: f64 = data.iter().sum();
    Some(sum / data.len() as f64)
}

/// Unbiased (sample) variance of `data`, denominator `n - 1`.
///
/// Returns `None` if there are fewer than two samples (the `n - 1` estimator
/// is undefined for `n < 2`).
#[must_use]
pub fn variance(data: &[f64]) -> Option<f64> {
    let n = data.len();
    if n < 2 {
        return None;
    }
    let m = mean(data)?;
    let ss: f64 = data.iter().map(|&x| (x - m) * (x - m)).sum();
    Some(ss / (n as f64 - 1.0))
}

/// Sample standard deviation (square root of [`variance`]).
///
/// Returns `None` for fewer than two samples.
#[must_use]
pub fn std(data: &[f64]) -> Option<f64> {
    variance(data).map(f64::sqrt)
}

/// The `p`-th percentile of `data`, `p ∈ [0, 100]`, by linear interpolation
/// between the closest ranks on the sorted sample.
///
/// `p = 0` returns the minimum and `p = 100` the maximum. The interpolation
/// convention matches NumPy's default (`linear`) `percentile`.
///
/// # Errors
/// * [`UqError::EmptyInput`] if `data` is empty.
/// * [`UqError::OutOfRange`] if `p` is outside `[0, 100]` or non-finite.
pub fn percentile(data: &[f64], p: f64) -> Result<f64, UqError> {
    if data.is_empty() {
        return Err(UqError::EmptyInput("percentile of an empty sample".into()));
    }
    if !p.is_finite() || !(0.0..=100.0).contains(&p) {
        return Err(UqError::OutOfRange(format!(
            "percentile p must be in [0, 100] (got {p})"
        )));
    }

    let mut sorted: Vec<f64> = data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let n = sorted.len();
    if n == 1 {
        return Ok(sorted[0]);
    }

    // Rank in [0, n-1]; interpolate between the two bracketing order
    // statistics.
    let rank = (p / 100.0) * (n as f64 - 1.0);
    let lo = rank.floor() as usize;
    let hi = rank.ceil() as usize;
    if lo == hi {
        return Ok(sorted[lo]);
    }
    let frac = rank - lo as f64;
    Ok(sorted[lo] * (1.0 - frac) + sorted[hi] * frac)
}

/// The median (50th percentile).
///
/// # Errors
/// [`UqError::EmptyInput`] if `data` is empty.
pub fn median(data: &[f64]) -> Result<f64, UqError> {
    percentile(data, 50.0)
}

/// A two-sided empirical confidence interval at confidence `level`
/// (e.g. `0.95` for 95 %), returned as `(lower, upper)`.
///
/// This is the **percentile** interval: the `(1 - level)/2` and
/// `(1 + level)/2` quantiles of the sample. It is non-parametric — it makes no
/// Gaussian assumption — and is the natural interval for a Monte-Carlo output
/// distribution.
///
/// # Errors
/// * [`UqError::EmptyInput`] if `data` is empty.
/// * [`UqError::OutOfRange`] if `level` is not in the open interval `(0, 1)`.
pub fn confidence_interval(data: &[f64], level: f64) -> Result<(f64, f64), UqError> {
    if data.is_empty() {
        return Err(UqError::EmptyInput(
            "confidence interval of an empty sample".into(),
        ));
    }
    if !level.is_finite() || level <= 0.0 || level >= 1.0 {
        return Err(UqError::OutOfRange(format!(
            "confidence level must be in (0, 1) (got {level})"
        )));
    }
    let tail = (1.0 - level) / 2.0;
    let lower = percentile(data, tail * 100.0)?;
    let upper = percentile(data, (1.0 - tail) * 100.0)?;
    Ok((lower, upper))
}
