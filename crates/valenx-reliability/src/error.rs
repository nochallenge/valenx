//! Error taxonomy for the reliability models.
//!
//! Every fallible constructor in this crate validates its arguments and
//! returns a [`ReliabilityError`] on bad input rather than panicking or
//! silently producing a nonsensical reliability value. The variants are
//! deliberately fine-grained so callers (and tests) can distinguish a
//! non-positive rate from a negative time from an out-of-range
//! probability.

use thiserror::Error;

/// Errors raised when constructing or evaluating a reliability model.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum ReliabilityError {
    /// A rate parameter (failure rate `lambda`) was not strictly
    /// positive. A constant failure rate must be `> 0` for the
    /// exponential model to be well defined (a zero rate is an
    /// immortal component; a negative rate is unphysical).
    #[error("failure rate lambda must be > 0, got {value}")]
    NonPositiveRate {
        /// The offending value supplied by the caller.
        value: f64,
    },

    /// A Weibull shape parameter `beta` was not strictly positive.
    #[error("Weibull shape beta must be > 0, got {value}")]
    NonPositiveShape {
        /// The offending value supplied by the caller.
        value: f64,
    },

    /// A Weibull scale (characteristic life) parameter `eta` was not
    /// strictly positive.
    #[error("Weibull scale eta must be > 0, got {value}")]
    NonPositiveScale {
        /// The offending value supplied by the caller.
        value: f64,
    },

    /// A time argument was negative. Reliability is defined for
    /// `t >= 0`; negative ages have no meaning.
    #[error("time t must be >= 0, got {value}")]
    NegativeTime {
        /// The offending value supplied by the caller.
        value: f64,
    },

    /// A supplied component reliability fell outside the unit interval
    /// `[0, 1]`. Reliabilities are probabilities.
    #[error("reliability must lie in [0, 1], got {value}")]
    ProbabilityOutOfRange {
        /// The offending value supplied by the caller.
        value: f64,
    },

    /// A non-finite (`NaN` or infinite) value was supplied where a
    /// finite number is required.
    #[error("value for `{what}` must be finite, got {value}")]
    NotFinite {
        /// Which parameter was non-finite.
        what: &'static str,
        /// The offending value supplied by the caller.
        value: f64,
    },

    /// A reliability block diagram was built with no components.
    #[error("a {kind} system needs at least one component, got none")]
    EmptySystem {
        /// The system kind that was empty (for example `"series"`).
        kind: &'static str,
    },

    /// A `k`-out-of-`n` structure was requested with `k` outside the
    /// valid range `1 ..= n`.
    #[error("k-out-of-n requires 1 <= k <= n, got k = {k}, n = {n}")]
    InvalidKofN {
        /// The required number of working components.
        k: usize,
        /// The total number of components.
        n: usize,
    },
}

/// Coarse classification of a [`ReliabilityError`].
///
/// Useful for callers that want to react to *kinds* of failure (for
/// example, surfacing input errors to a user versus logging an internal
/// numerical problem) without matching every individual variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The caller supplied an invalid value (bad parameter, out-of-range
    /// probability, negative time, empty system, ...).
    Input,
    /// A value was non-finite where a finite number is required.
    Numeric,
}

impl ReliabilityError {
    /// A stable, kebab-cased identifier for this error.
    ///
    /// The string is part of the crate's contract and is intended to be
    /// matched programmatically or logged; it will not change for a
    /// given variant.
    #[must_use]
    pub fn code(&self) -> &'static str {
        match self {
            ReliabilityError::NonPositiveRate { .. } => "reliability.non-positive-rate",
            ReliabilityError::NonPositiveShape { .. } => "reliability.non-positive-shape",
            ReliabilityError::NonPositiveScale { .. } => "reliability.non-positive-scale",
            ReliabilityError::NegativeTime { .. } => "reliability.negative-time",
            ReliabilityError::ProbabilityOutOfRange { .. } => {
                "reliability.probability-out-of-range"
            }
            ReliabilityError::NotFinite { .. } => "reliability.not-finite",
            ReliabilityError::EmptySystem { .. } => "reliability.empty-system",
            ReliabilityError::InvalidKofN { .. } => "reliability.invalid-k-of-n",
        }
    }

    /// The coarse [`ErrorCategory`] this error belongs to.
    #[must_use]
    pub fn category(&self) -> ErrorCategory {
        match self {
            ReliabilityError::NotFinite { .. } => ErrorCategory::Numeric,
            _ => ErrorCategory::Input,
        }
    }
}

/// Validate that `value` for parameter `what` is finite, returning a
/// [`ReliabilityError::NotFinite`] otherwise.
pub(crate) fn require_finite(what: &'static str, value: f64) -> Result<f64, ReliabilityError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ReliabilityError::NotFinite { what, value })
    }
}

/// Validate that `t` is a finite, non-negative time.
pub(crate) fn require_time(t: f64) -> Result<f64, ReliabilityError> {
    let t = require_finite("time", t)?;
    if t < 0.0 {
        return Err(ReliabilityError::NegativeTime { value: t });
    }
    Ok(t)
}

/// Validate that `r` is a finite probability in `[0, 1]`.
pub(crate) fn require_probability(r: f64) -> Result<f64, ReliabilityError> {
    let r = require_finite("reliability", r)?;
    if !(0.0..=1.0).contains(&r) {
        return Err(ReliabilityError::ProbabilityOutOfRange { value: r });
    }
    Ok(r)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_are_stable_and_distinct() {
        let errors = [
            ReliabilityError::NonPositiveRate { value: 0.0 },
            ReliabilityError::NonPositiveShape { value: 0.0 },
            ReliabilityError::NonPositiveScale { value: 0.0 },
            ReliabilityError::NegativeTime { value: -1.0 },
            ReliabilityError::ProbabilityOutOfRange { value: 2.0 },
            ReliabilityError::NotFinite {
                what: "x",
                value: f64::NAN,
            },
            ReliabilityError::EmptySystem { kind: "series" },
            ReliabilityError::InvalidKofN { k: 3, n: 2 },
        ];
        let mut codes: Vec<&str> = errors.iter().map(ReliabilityError::code).collect();
        let count = codes.len();
        codes.sort_unstable();
        codes.dedup();
        assert_eq!(codes.len(), count, "every variant must have a unique code");
    }

    #[test]
    fn categories_split_input_from_numeric() {
        assert_eq!(
            ReliabilityError::NegativeTime { value: -1.0 }.category(),
            ErrorCategory::Input
        );
        assert_eq!(
            ReliabilityError::NotFinite {
                what: "lambda",
                value: f64::INFINITY,
            }
            .category(),
            ErrorCategory::Numeric
        );
    }

    #[test]
    fn require_time_rejects_negative_and_nonfinite() {
        assert!(matches!(
            require_time(-0.5),
            Err(ReliabilityError::NegativeTime { .. })
        ));
        assert!(matches!(
            require_time(f64::NAN),
            Err(ReliabilityError::NotFinite { .. })
        ));
        assert_eq!(require_time(3.0).unwrap(), 3.0);
    }

    #[test]
    fn require_probability_enforces_unit_interval() {
        assert!(matches!(
            require_probability(-0.01),
            Err(ReliabilityError::ProbabilityOutOfRange { .. })
        ));
        assert!(matches!(
            require_probability(1.01),
            Err(ReliabilityError::ProbabilityOutOfRange { .. })
        ));
        assert_eq!(require_probability(0.0).unwrap(), 0.0);
        assert_eq!(require_probability(1.0).unwrap(), 1.0);
    }

    #[test]
    fn display_messages_include_offending_value() {
        let e = ReliabilityError::NonPositiveRate { value: -2.5 };
        let msg = format!("{e}");
        assert!(msg.contains("-2.5"), "message should echo the value: {msg}");
    }
}
