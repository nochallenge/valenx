//! Shared input validation — the validated constructors behind every
//! estimator.
//!
//! Each estimator funnels its preconditions through these helpers so the
//! same checks (non-empty, minimum count, all-finite, in-range, positive
//! scale) are enforced uniformly and produce the same [`StatsError`]
//! variants everywhere. Keeping them in one place is what makes the error
//! taxonomy in [`super::error`] meaningful: a `mean` of an empty slice and
//! a `median` of an empty slice raise the identical, well-typed error.

use crate::error::{Result, StatsError};

/// Ensure a sample slice is non-empty.
///
/// Returns [`StatsError::EmptySample`] tagged with `estimator` when `data`
/// is empty.
pub fn non_empty(data: &[f64], estimator: &'static str) -> Result<()> {
    if data.is_empty() {
        return Err(StatsError::EmptySample { estimator });
    }
    Ok(())
}

/// Ensure a sample slice holds at least `needed` observations.
///
/// Returns [`StatsError::TooFewObservations`] when `data.len() < needed`.
pub fn at_least(data: &[f64], needed: usize, estimator: &'static str) -> Result<()> {
    let got = data.len();
    if got < needed {
        return Err(StatsError::TooFewObservations {
            estimator,
            needed,
            got,
        });
    }
    Ok(())
}

/// Ensure a single scalar is finite.
///
/// Returns [`StatsError::NonFinite`] tagged with `name` when `value` is
/// `NaN`, `+∞` or `-∞`.
pub fn finite(value: f64, name: &'static str) -> Result<()> {
    if !value.is_finite() {
        return Err(StatsError::NonFinite { name });
    }
    Ok(())
}

/// Ensure every element of a slice is finite.
///
/// Returns [`StatsError::NonFinite`] tagged with `name` at the first
/// non-finite element encountered.
pub fn all_finite(data: &[f64], name: &'static str) -> Result<()> {
    for &x in data {
        finite(x, name)?;
    }
    Ok(())
}

/// Ensure a quantile / probability lies in the closed unit interval.
///
/// Returns [`StatsError::OutOfRange`] when `q` is outside `[0, 1]`, after
/// first checking finiteness via [`finite`].
pub fn unit_closed(q: f64, name: &'static str) -> Result<()> {
    finite(q, name)?;
    if !(0.0..=1.0).contains(&q) {
        return Err(StatsError::OutOfRange {
            name,
            value: q,
            expected: "[0, 1]",
        });
    }
    Ok(())
}

/// Ensure a confidence level lies in the open unit interval `(0, 1)`.
///
/// A confidence level of exactly `0` or `1` is rejected: zero gives a
/// degenerate zero-width interval and one needs an infinite critical value.
/// Returns [`StatsError::OutOfRange`] otherwise, after a finiteness check.
pub fn unit_open(c: f64, name: &'static str) -> Result<()> {
    finite(c, name)?;
    if !(c > 0.0 && c < 1.0) {
        return Err(StatsError::OutOfRange {
            name,
            value: c,
            expected: "(0, 1)",
        });
    }
    Ok(())
}

/// Ensure a scale parameter (a standard deviation / sigma) is finite and
/// strictly positive.
///
/// Returns [`StatsError::NonFinite`] for a non-finite value, otherwise
/// [`StatsError::NonPositiveScale`] when `value <= 0`.
pub fn positive_scale(value: f64, name: &'static str) -> Result<()> {
    finite(value, name)?;
    if value <= 0.0 {
        return Err(StatsError::NonPositiveScale { name, value });
    }
    Ok(())
}
