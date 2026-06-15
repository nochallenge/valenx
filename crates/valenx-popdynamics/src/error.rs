//! Error taxonomy for `valenx-popdynamics`.
//!
//! Every fallible public function returns
//! [`Result<_, PopError>`]. The variants are deliberately coarse — a
//! population-dynamics caller usually only cares about two things:
//!
//! 1. Did the caller pass a parameter the model cannot accept — a
//!    negative carrying capacity, a non-positive integration step, a
//!    `t_end` that is not strictly after `t_start`
//!    ([`PopError::Invalid`])?
//! 2. Did the integration ask for an unreasonable amount of work — a
//!    `(t_end - t_start) / dt` step count that would overflow the
//!    result buffer ([`PopError::TooManySteps`])?
//!
//! Use [`PopError::code`] for stable log / telemetry tagging and
//! [`PopError::category`] to bucket failures without matching every
//! variant. The pattern mirrors the `valenx-sysbio` and
//! `valenx-springs` error modules.

use thiserror::Error;

/// Errors produced by `valenx-popdynamics`.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum PopError {
    /// Caller passed an argument the model cannot accept: a negative
    /// rate or carrying capacity, a non-positive step size, a time
    /// window whose end does not strictly follow its start, or an
    /// initial state with a negative compartment. A property of the
    /// *call*, not of any numerical failure.
    #[error("invalid `{what}`: {reason}")]
    Invalid {
        /// Logical parameter name (e.g. `"dt"`, `"k"`, `"t_end"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// The requested time window and step size imply more integration
    /// steps than the configured ceiling. Guards against a tiny `dt`
    /// over a huge horizon silently allocating an enormous trajectory.
    #[error("too many steps: {requested} requested, ceiling is {ceiling}")]
    TooManySteps {
        /// Number of steps the `(t_end - t_start) / dt` window implies.
        requested: u64,
        /// Hard ceiling enforced by the integrator.
        ceiling: u64,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this
/// rather than on the individual error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// User-supplied input is wrong (bad parameter or bad time window).
    Input,
    /// A resource limit was hit (the step-count ceiling).
    Limit,
}

impl PopError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"popdynamics.<sub_id>"`. Codes never change
    /// across minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            PopError::Invalid { .. } => "popdynamics.invalid",
            PopError::TooManySteps { .. } => "popdynamics.too_many_steps",
        }
    }

    /// Coarse category — see [`ErrorCategory`].
    pub fn category(&self) -> ErrorCategory {
        match self {
            PopError::Invalid { .. } => ErrorCategory::Input,
            PopError::TooManySteps { .. } => ErrorCategory::Limit,
        }
    }

    /// Convenience constructor for [`PopError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        PopError::Invalid {
            what,
            reason: reason.into(),
        }
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, PopError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = PopError::invalid("dt", "must be positive");
        assert_eq!(err.code(), "popdynamics.invalid");
        assert_eq!(err.category(), ErrorCategory::Input);

        let err = PopError::TooManySteps {
            requested: 10,
            ceiling: 5,
        };
        assert_eq!(err.code(), "popdynamics.too_many_steps");
        assert_eq!(err.category(), ErrorCategory::Limit);
    }

    #[test]
    fn display_is_informative() {
        let msg = PopError::invalid("k", "must be positive").to_string();
        assert!(msg.contains('k'), "got: {msg}");
        assert!(msg.contains("positive"), "got: {msg}");

        let msg = PopError::TooManySteps {
            requested: 99,
            ceiling: 8,
        }
        .to_string();
        assert!(msg.contains("99"), "got: {msg}");
        assert!(msg.contains('8'), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(PopError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }
}
