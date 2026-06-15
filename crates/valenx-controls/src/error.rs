//! Error taxonomy for `valenx-controls`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, ControlsError>`]. The variants are deliberately coarse — a
//! control-systems caller usually only cares about two things:
//!
//! 1. Did the caller pass a non-physical parameter — a negative natural
//!    frequency, a negative damping ratio, a non-finite gain, a
//!    non-positive sample time ([`ControlsError::InvalidParameter`])?
//! 2. Did a closed-form metric apply outside its domain — e.g. the
//!    underdamped peak-time / overshoot formulae requested for a
//!    critically- or over-damped system, where the damped frequency is
//!    zero and the expression is singular ([`ControlsError::DomainError`])?
//!
//! Use [`ControlsError::code`] for stable log / telemetry tagging and
//! [`ControlsError::category`] to bucket failures without matching every
//! variant. The shape mirrors `valenx-astro`'s `AstroError` and
//! `valenx-springs`'s `SpringsError`.

use thiserror::Error;

/// Shorthand for `Result<T, ControlsError>`.
pub type Result<T> = core::result::Result<T, ControlsError>;

/// Anything that can go wrong constructing a model or evaluating a metric.
///
/// This enum is `#[non_exhaustive]`: new variants may be added in future
/// releases without it being a breaking change, so downstream `match`
/// arms must include a wildcard.
#[derive(Debug, Error, Clone, PartialEq)]
#[non_exhaustive]
pub enum ControlsError {
    /// A scalar parameter was non-physical in a way that would otherwise
    /// feed a silent `NaN`/`Inf` into a `√`/`ln`/`exp` or a division —
    /// a non-finite or non-positive natural frequency, a non-finite or
    /// negative damping ratio, a non-finite gain, a non-positive sample
    /// time, etc. A property of the *argument*, not of a model's state.
    #[error("invalid parameter `{what}`: {reason}")]
    InvalidParameter {
        /// Logical parameter name (e.g. `"wn"`, `"zeta"`, `"dt"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// A closed-form metric was requested outside the domain on which it
    /// is defined — e.g. the underdamped peak-time or percent-overshoot
    /// formulae for a `ζ ≥ 1` (critically- or over-damped) system, where
    /// the damped natural frequency is zero and the expression is
    /// singular. Carries a short reason.
    #[error("domain error: {0}")]
    DomainError(&'static str),
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this rather
/// than on every error variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// User-supplied input is wrong (a bad argument value).
    Input,
    /// A closed-form was evaluated outside its mathematical domain.
    Domain,
}

impl ControlsError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"controls.<sub_id>"`. Codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            ControlsError::InvalidParameter { .. } => "controls.invalid_parameter",
            ControlsError::DomainError(_) => "controls.domain_error",
        }
    }

    /// Coarse [`ErrorCategory`] for routing / display.
    pub fn category(&self) -> ErrorCategory {
        match self {
            ControlsError::InvalidParameter { .. } => ErrorCategory::Input,
            ControlsError::DomainError(_) => ErrorCategory::Domain,
        }
    }

    /// Convenience constructor for [`ControlsError::InvalidParameter`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        ControlsError::InvalidParameter {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`ControlsError::DomainError`].
    pub fn domain(reason: &'static str) -> Self {
        ControlsError::DomainError(reason)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = ControlsError::invalid("wn", "must be > 0");
        assert_eq!(err.code(), "controls.invalid_parameter");
        assert_eq!(err.category(), ErrorCategory::Input);

        let err = ControlsError::domain("peak time undefined for zeta >= 1");
        assert_eq!(err.code(), "controls.domain_error");
        assert_eq!(err.category(), ErrorCategory::Domain);
    }

    #[test]
    fn display_is_informative() {
        let msg = ControlsError::invalid("zeta", "must be non-negative").to_string();
        assert!(msg.contains("zeta"), "got: {msg}");
        assert!(msg.contains("non-negative"), "got: {msg}");

        let msg = ControlsError::domain("singular").to_string();
        assert!(msg.contains("singular"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(ControlsError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }
}
