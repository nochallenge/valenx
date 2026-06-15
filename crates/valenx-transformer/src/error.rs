//! Error taxonomy for `valenx-transformer`.
//!
//! Every fallible constructor and function in this crate returns
//! [`Result<_, TransformerError>`](crate::Result). The variants are
//! deliberately coarse — a transformer-relations caller usually only
//! needs to distinguish two situations:
//!
//! 1. A scalar argument is outside its physical domain — a
//!    non-positive turns count, a zero or negative voltage where a ratio
//!    is about to be taken, an efficiency outside `(0, 1]`
//!    ([`TransformerError::Invalid`]).
//! 2. A computation would divide by a quantity that has been validated
//!    to be non-zero elsewhere but reached zero here — a guard that
//!    should never fire for well-constructed inputs
//!    ([`TransformerError::DivideByZero`]).
//!
//! Use [`TransformerError::code`] for stable log / telemetry tagging and
//! [`TransformerError::category`] to bucket failures without matching
//! every variant. The pattern mirrors the sibling physics crates'
//! `*Error` taxonomies.

use thiserror::Error;

/// Errors raised by the transformer relations.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum TransformerError {
    /// A caller-supplied scalar is outside the domain the relation can
    /// accept: a non-positive number of turns, a non-positive voltage or
    /// impedance used as a ratio denominator, an efficiency that is not
    /// in the half-open interval `(0, 1]`, and so on. A property of the
    /// *call's arguments*, not of an internal state.
    #[error("invalid `{what}`: {reason}")]
    Invalid {
        /// Logical parameter name (for example `"turns_primary"`,
        /// `"efficiency"`, `"voltage_secondary"`).
        what: &'static str,
        /// Human-readable reason the value was rejected.
        reason: String,
    },

    /// A relation reached a zero denominator at the point of division.
    /// Constructors validate their inputs up front, so this is a
    /// defensive guard that should not fire for well-formed inputs; it
    /// exists so the crate never panics on a stray zero.
    #[error("division by zero while computing `{what}`")]
    DivideByZero {
        /// Logical quantity whose computation hit a zero denominator.
        what: &'static str,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this rather
/// than on every error variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// Caller-supplied input is outside its physical domain.
    Input,
    /// A defensive numerical guard fired (zero denominator).
    Numeric,
}

impl TransformerError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"transformer.<sub_id>"`. Codes never change
    /// across minor versions.
    #[must_use]
    pub fn code(&self) -> &'static str {
        match self {
            TransformerError::Invalid { .. } => "transformer.invalid",
            TransformerError::DivideByZero { .. } => "transformer.divide_by_zero",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    #[must_use]
    pub fn category(&self) -> &'static str {
        match self {
            TransformerError::Invalid { .. } => "input",
            TransformerError::DivideByZero { .. } => "numeric",
        }
    }

    /// Typed category enum (for callers that want to `match` instead of
    /// comparing the [`category`](Self::category) string).
    #[must_use]
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            TransformerError::Invalid { .. } => ErrorCategory::Input,
            TransformerError::DivideByZero { .. } => ErrorCategory::Numeric,
        }
    }

    /// Convenience constructor for [`TransformerError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        TransformerError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`TransformerError::DivideByZero`].
    #[must_use]
    pub fn divide_by_zero(what: &'static str) -> Self {
        TransformerError::DivideByZero { what }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = TransformerError::invalid("turns_primary", "must be positive");
        assert_eq!(err.code(), "transformer.invalid");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = TransformerError::divide_by_zero("turns_ratio");
        assert_eq!(err.code(), "transformer.divide_by_zero");
        assert_eq!(err.category(), "numeric");
        assert_eq!(err.category_enum(), ErrorCategory::Numeric);
    }

    #[test]
    fn display_is_informative() {
        let msg = TransformerError::invalid("efficiency", "must be in (0, 1]").to_string();
        assert!(msg.contains("efficiency"), "got: {msg}");
        assert!(msg.contains("(0, 1]"), "got: {msg}");

        let msg = TransformerError::divide_by_zero("turns_ratio").to_string();
        assert!(msg.contains("turns_ratio"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(TransformerError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }
}
