//! Error taxonomy for `valenx-statistics`.
//!
//! Every fallible public function returns [`Result<_, StatsError>`]. The
//! variants are deliberately coarse: a caller usually only needs to know
//! whether a sample was empty, whether it was too small for the estimator
//! (sample variance needs at least two observations), whether a parameter
//! was out of range (a probability / quantile outside `[0, 1]`, a
//! non-positive confidence level), or whether a non-finite value (`NaN` /
//! `±∞`) slipped into the data.
//!
//! Construct errors through the validated helpers in [`super::validate`]
//! rather than building the variants by hand — the helpers centralise the
//! finiteness and range checks so every estimator enforces the same
//! preconditions.

use thiserror::Error;

/// Errors produced by `valenx-statistics`.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum StatsError {
    /// The sample slice was empty. The mean, median and quartiles all need
    /// at least one observation; there is no meaningful value for none.
    #[error("empty sample: `{estimator}` needs at least one observation")]
    EmptySample {
        /// The estimator that was asked for, e.g. `"mean"`.
        estimator: &'static str,
    },

    /// The sample was non-empty but smaller than the estimator's minimum.
    /// The sample variance / standard deviation divide by `n - 1`, so they
    /// require `n >= 2`; with a single point the dispersion is undefined.
    #[error("sample too small: `{estimator}` needs at least {needed} observations, got {got}")]
    TooFewObservations {
        /// The estimator that was asked for, e.g. `"sample_variance"`.
        estimator: &'static str,
        /// The minimum number of observations the estimator requires.
        needed: usize,
        /// The number of observations actually supplied.
        got: usize,
    },

    /// A parameter fell outside its valid range — a probability or quantile
    /// `q` not in `[0, 1]`, or a confidence level not in the open `(0, 1)`.
    #[error("parameter `{name}` out of range: {value} is not in {expected}")]
    OutOfRange {
        /// The parameter name, e.g. `"q"` or `"confidence"`.
        name: &'static str,
        /// The offending value.
        value: f64,
        /// A human-readable description of the valid range, e.g. `"[0, 1]"`.
        expected: &'static str,
    },

    /// A supplied value was not finite (`NaN`, `+∞` or `-∞`). Order
    /// statistics and the closed-form estimators are only meaningful over
    /// finite reals, so non-finite inputs are rejected up front rather than
    /// silently producing `NaN` results.
    #[error("non-finite value for `{name}`: every input must be finite (no NaN or infinity)")]
    NonFinite {
        /// What contained the non-finite value, e.g. `"sample"` or `"sigma"`.
        name: &'static str,
    },

    /// A scale parameter that must be strictly positive was zero or
    /// negative — a standard deviation / sigma used as the denominator of a
    /// z-score, or supplied as the known population spread of an interval.
    #[error("non-positive scale `{name}`: {value} must be strictly greater than zero")]
    NonPositiveScale {
        /// The parameter name, e.g. `"std"` or `"sigma"`.
        name: &'static str,
        /// The offending value.
        value: f64,
    },
}

/// Convenience alias for `Result<T, StatsError>`.
pub type Result<T> = std::result::Result<T, StatsError>;
