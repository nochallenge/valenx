//! Error taxonomy for `valenx-pulley`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, PulleyError>`]. The variants are intentionally coarse — a
//! pulley-mechanics caller usually only cares about three things:
//!
//! 1. Did the caller pass a non-physical number — a negative load, a zero
//!    rope count, an efficiency outside `(0, 1]`
//!    ([`PulleyError::Invalid`])?
//! 2. Would the requested arrangement divide by zero or otherwise be
//!    geometrically degenerate — a block-and-tackle with no supporting
//!    rope segments ([`PulleyError::Degenerate`])?
//! 3. Do two inputs disagree on a quantity the model needs to reconcile —
//!    a measured effort that is below the friction-free ideal, which would
//!    imply efficiency above `1` ([`PulleyError::Inconsistent`])?
//!
//! Use [`PulleyError::code`] for stable log / telemetry tagging and
//! [`PulleyError::category`] (or [`PulleyError::category_enum`]) to bucket
//! failures without matching every variant. The pattern mirrors
//! `valenx-springs`'s `SpringsError` and `valenx-popgen`'s `PopgenError`.

use thiserror::Error;

/// Errors produced by `valenx-pulley`.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum PulleyError {
    /// Caller passed a value the model cannot accept: a negative or
    /// non-finite load, a rope count of zero, an efficiency outside the
    /// half-open interval `(0, 1]`, etc. A property of the *call*, not of
    /// any reconciliation between inputs.
    #[error("invalid `{what}`: {reason}")]
    Invalid {
        /// Logical parameter name (e.g. `"load"`, `"supporting_ropes"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// The requested arrangement is geometrically degenerate — it would
    /// divide by zero or has no rope segment supporting the movable block.
    #[error("degenerate pulley system: {reason}")]
    Degenerate {
        /// Human-readable reason.
        reason: String,
    },

    /// Two supplied quantities cannot be reconciled: most commonly a
    /// measured actual effort that is *below* the friction-free ideal
    /// effort, which would imply a (physically impossible) efficiency
    /// greater than `1`.
    #[error("inconsistent inputs: {reason}")]
    Inconsistent {
        /// Human-readable reason.
        reason: String,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this rather
/// than on the individual error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// User-supplied input is wrong (bad argument value).
    Input,
    /// The configured arrangement is geometrically degenerate.
    Geometry,
    /// Two inputs disagree and cannot be reconciled.
    Consistency,
}

impl PulleyError {
    /// Stable snake-cased error code suitable for log / telemetry tagging.
    /// Format: `"pulley.<sub_id>"`. Codes never change across minor
    /// versions.
    pub fn code(&self) -> &'static str {
        match self {
            PulleyError::Invalid { .. } => "pulley.invalid",
            PulleyError::Degenerate { .. } => "pulley.degenerate",
            PulleyError::Inconsistent { .. } => "pulley.inconsistent",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            PulleyError::Invalid { .. } => "input",
            PulleyError::Degenerate { .. } => "geometry",
            PulleyError::Inconsistent { .. } => "consistency",
        }
    }

    /// Typed category enum (for callers that want to `match` instead of
    /// comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            PulleyError::Invalid { .. } => ErrorCategory::Input,
            PulleyError::Degenerate { .. } => ErrorCategory::Geometry,
            PulleyError::Inconsistent { .. } => ErrorCategory::Consistency,
        }
    }

    /// Convenience constructor for [`PulleyError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        PulleyError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`PulleyError::Degenerate`].
    pub fn degenerate(reason: impl Into<String>) -> Self {
        PulleyError::Degenerate {
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`PulleyError::Inconsistent`].
    pub fn inconsistent(reason: impl Into<String>) -> Self {
        PulleyError::Inconsistent {
            reason: reason.into(),
        }
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, PulleyError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = PulleyError::invalid("load", "must be non-negative");
        assert_eq!(err.code(), "pulley.invalid");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = PulleyError::degenerate("no supporting ropes");
        assert_eq!(err.code(), "pulley.degenerate");
        assert_eq!(err.category(), "geometry");
        assert_eq!(err.category_enum(), ErrorCategory::Geometry);

        let err = PulleyError::inconsistent("effort below ideal");
        assert_eq!(err.code(), "pulley.inconsistent");
        assert_eq!(err.category(), "consistency");
        assert_eq!(err.category_enum(), ErrorCategory::Consistency);
    }

    #[test]
    fn display_is_informative() {
        let msg = PulleyError::invalid("supporting_ropes", "must be >= 1").to_string();
        assert!(msg.contains("supporting_ropes"), "got: {msg}");
        assert!(msg.contains("must be >= 1"), "got: {msg}");

        let msg = PulleyError::degenerate("divide by zero").to_string();
        assert!(msg.contains("divide by zero"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(PulleyError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }
}
