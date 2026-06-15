//! Error taxonomy for the spring-combination calculator.
//!
//! Every fallible entry point in this crate returns
//! [`Result<_, SpringError>`]. The error type carries a stable
//! kebab-cased [`code`](SpringError::code) and a coarse
//! [`category`](SpringError::category) so callers can branch on the
//! failure mode (telemetry, user-facing messages, retry logic) without
//! string-matching the `Display` text.

use thiserror::Error;

/// Failure modes of the spring-combination models.
///
/// Variants are deliberately specific: an out-of-range spring rate is a
/// different problem from an empty combination, which is different again
/// from a series combination that contains a zero-rate (infinitely soft)
/// spring. Each carries enough context to render a precise message.
#[derive(Debug, Clone, PartialEq, Error)]
#[non_exhaustive]
pub enum SpringError {
    /// A spring rate `k` was not a finite, strictly positive value.
    ///
    /// A physical linear spring has a stiffness `k > 0` measured in
    /// newtons per metre (N/m). Zero, negative, infinite and NaN rates
    /// are rejected here rather than silently producing nonsense
    /// downstream.
    #[error("spring rate must be finite and > 0, got {value} (in {context})")]
    NonPositiveRate {
        /// The offending rate value, in N/m.
        value: f64,
        /// Where the bad value was supplied (which constructor / call).
        context: &'static str,
    },

    /// A combination was requested over an empty set of springs.
    ///
    /// Both the parallel and series reductions are folds over the member
    /// rates; with no members there is no well-defined equivalent rate
    /// (the parallel sum would be `0`, the series reciprocal sum would
    /// divide by zero), so the empty case is rejected explicitly.
    #[error("a {combination} combination needs at least one spring, got none")]
    EmptyCombination {
        /// The combination kind that was requested ("parallel" /
        /// "series").
        combination: &'static str,
    },

    /// A displacement was not a finite number.
    ///
    /// Force and energy are evaluated at a deflection `x` (in metres).
    /// NaN and infinite deflections are rejected; any finite value
    /// (including negative, i.e. compression on the opposite sign
    /// convention, and zero) is accepted.
    #[error("displacement must be finite, got {value}")]
    NonFiniteDisplacement {
        /// The offending displacement value, in metres.
        value: f64,
    },
}

impl SpringError {
    /// A stable, kebab-cased identifier for this error.
    ///
    /// Unlike the human-readable [`Display`](std::fmt::Display) text the
    /// code is guaranteed not to change across releases, so it is safe to
    /// log, compare against, or surface in an API.
    ///
    /// ```
    /// use valenx_springcombination::SpringError;
    /// let e = SpringError::NonFiniteDisplacement { value: f64::NAN };
    /// assert_eq!(e.code(), "spring.non-finite-displacement");
    /// ```
    pub fn code(&self) -> &'static str {
        match self {
            SpringError::NonPositiveRate { .. } => "spring.non-positive-rate",
            SpringError::EmptyCombination { .. } => "spring.empty-combination",
            SpringError::NonFiniteDisplacement { .. } => "spring.non-finite-displacement",
        }
    }

    /// The coarse category this error belongs to.
    ///
    /// Useful for deciding, at a glance, whether the caller passed bad
    /// input ([`ErrorCategory::Input`]) versus asked for something the
    /// model cannot evaluate ([`ErrorCategory::Domain`]).
    pub fn category(&self) -> ErrorCategory {
        match self {
            SpringError::NonPositiveRate { .. } => ErrorCategory::Input,
            SpringError::EmptyCombination { .. } => ErrorCategory::Domain,
            SpringError::NonFiniteDisplacement { .. } => ErrorCategory::Input,
        }
    }
}

/// A coarse bucket for [`SpringError`] values.
///
/// This is intentionally tiny — it groups the specific variants into the
/// two failure classes a caller usually cares about.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
#[non_exhaustive]
pub enum ErrorCategory {
    /// The caller supplied a value outside the model's accepted range
    /// (a non-positive rate, a non-finite displacement).
    Input,
    /// The request itself is not evaluable by the model (an empty
    /// combination).
    Domain,
}
