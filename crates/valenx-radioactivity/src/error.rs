//! Error taxonomy for `valenx-radioactivity`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, RadioactivityError>`]. The variants are deliberately
//! narrow — a decay calculation only fails in three ways:
//!
//! 1. A physical quantity that must be strictly positive (a decay
//!    constant, a half-life, an initial population) was given as zero,
//!    negative, or non-finite ([`RadioactivityError::NonPositive`]).
//! 2. A dimensionless fraction that must lie in a half-open or closed
//!    interval (a remaining fraction in `(0, 1]`) was out of range
//!    ([`RadioactivityError::OutOfRange`]).
//! 3. A two-member chain was built with a parent and daughter sharing
//!    the *same* decay constant, for which the Bateman solution has a
//!    removable `0 / 0` singularity that this crate does not special-case
//!    ([`RadioactivityError::DegenerateChain`]).
//!
//! Use [`RadioactivityError::code`] for stable log / telemetry tagging.
//! The shape mirrors `valenx-springs`'s `SpringsError`.

use thiserror::Error;

/// Errors raised by radioactive-decay calculations.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum RadioactivityError {
    /// A quantity that must be strictly positive and finite (decay
    /// constant `lambda`, half-life, mean life, initial population `N0`,
    /// elapsed time when a positive value is required) was not.
    #[error("`{what}` must be a strictly positive, finite number, got {value}")]
    NonPositive {
        /// Logical parameter name (e.g. `"lambda"`, `"half_life"`).
        what: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A dimensionless value that must lie inside a stated interval was
    /// out of range (e.g. a remaining fraction outside `(0, 1]`).
    #[error("`{what}` must lie in {interval}, got {value}")]
    OutOfRange {
        /// Logical parameter name (e.g. `"remaining_fraction"`).
        what: &'static str,
        /// The offending value.
        value: f64,
        /// Human-readable interval, e.g. `"(0, 1]"`.
        interval: &'static str,
    },

    /// A two-member [`crate::chain::DecayChain`] was constructed with the
    /// parent and daughter decay constants equal (within a tolerance).
    /// The Bateman daughter term carries a `1 / (lambda_d - lambda_p)`
    /// factor whose `lambda_d == lambda_p` limit is finite but is not
    /// special-cased here; supply distinct constants instead.
    #[error(
        "degenerate chain: parent and daughter decay constants are equal \
         ({lambda} per unit time); the lambda_d == lambda_p limit is not handled"
    )]
    DegenerateChain {
        /// The shared decay constant the two members were given.
        lambda: f64,
    },
}

/// Coarse error category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this rather
/// than on every error variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A caller-supplied argument is outside its physical domain.
    Input,
    /// A model would hit a singularity the crate does not handle.
    Model,
}

impl RadioactivityError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"radioactivity.<sub_id>"`. Codes never change
    /// across minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            RadioactivityError::NonPositive { .. } => "radioactivity.non_positive",
            RadioactivityError::OutOfRange { .. } => "radioactivity.out_of_range",
            RadioactivityError::DegenerateChain { .. } => "radioactivity.degenerate_chain",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            RadioactivityError::NonPositive { .. } | RadioactivityError::OutOfRange { .. } => {
                ErrorCategory::Input
            }
            RadioactivityError::DegenerateChain { .. } => ErrorCategory::Model,
        }
    }

    /// Validating constructor for [`RadioactivityError::NonPositive`].
    ///
    /// Returns `Ok(value)` when `value` is finite and strictly greater
    /// than zero, otherwise the error. This is the single chokepoint used
    /// by every public constructor in the crate to reject bad inputs.
    pub fn require_positive(what: &'static str, value: f64) -> Result<f64> {
        if value.is_finite() && value > 0.0 {
            Ok(value)
        } else {
            Err(RadioactivityError::NonPositive { what, value })
        }
    }

    /// Validating constructor that requires a finite, non-negative value
    /// (zero is allowed — used for elapsed times, where `t = 0` is the
    /// legitimate initial instant).
    pub fn require_non_negative(what: &'static str, value: f64) -> Result<f64> {
        if value.is_finite() && value >= 0.0 {
            Ok(value)
        } else {
            Err(RadioactivityError::NonPositive { what, value })
        }
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, RadioactivityError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn require_positive_accepts_positive_finite() {
        assert_eq!(RadioactivityError::require_positive("lambda", 2.5), Ok(2.5));
    }

    #[test]
    fn require_positive_rejects_zero_negative_and_nonfinite() {
        for &bad in &[0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let err = RadioactivityError::require_positive("lambda", bad)
                .expect_err("expected rejection");
            assert_eq!(err.code(), "radioactivity.non_positive");
            assert_eq!(err.category(), ErrorCategory::Input);
        }
    }

    #[test]
    fn require_non_negative_allows_zero_but_not_negative() {
        assert_eq!(RadioactivityError::require_non_negative("t", 0.0), Ok(0.0));
        assert!(RadioactivityError::require_non_negative("t", -0.1).is_err());
        assert!(RadioactivityError::require_non_negative("t", f64::NAN).is_err());
    }

    #[test]
    fn codes_and_categories_are_stable() {
        let e = RadioactivityError::OutOfRange {
            what: "remaining_fraction",
            value: 2.0,
            interval: "(0, 1]",
        };
        assert_eq!(e.code(), "radioactivity.out_of_range");
        assert_eq!(e.category(), ErrorCategory::Input);

        let e = RadioactivityError::DegenerateChain { lambda: 0.3 };
        assert_eq!(e.code(), "radioactivity.degenerate_chain");
        assert_eq!(e.category(), ErrorCategory::Model);
    }

    #[test]
    fn display_mentions_the_parameter_and_value() {
        let msg = RadioactivityError::NonPositive {
            what: "half_life",
            value: -3.0,
        }
        .to_string();
        assert!(msg.contains("half_life"), "got: {msg}");
        assert!(msg.contains("-3"), "got: {msg}");
    }

    #[test]
    fn error_is_a_std_error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(RadioactivityError::OutOfRange {
            what: "f",
            value: 9.0,
            interval: "(0, 1]",
        });
        assert!(err.to_string().contains('f'));
    }
}
