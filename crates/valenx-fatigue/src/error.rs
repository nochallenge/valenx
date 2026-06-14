//! Error taxonomy for `valenx-fatigue`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, FatigueError>`]. The variants are intentionally coarse —
//! a fatigue caller usually only cares about three things:
//!
//! 1. Did the caller pass a parameter outside its physical domain — a
//!    non-positive stress, a negative cycle count, a Basquin exponent
//!    that is not negative ([`FatigueError::Invalid`])?
//! 2. Are two inputs inconsistent with each other — an alternating
//!    stress at or above the ultimate strength, a mean stress that
//!    consumes the whole Goodman line ([`FatigueError::Domain`])?
//! 3. Did a damage accumulation receive mismatched cycle / capacity
//!    lists ([`FatigueError::Dimension`])?
//!
//! Use [`FatigueError::code`] for stable log / telemetry tagging and
//! [`FatigueError::category`] to bucket failures without matching every
//! variant. The pattern mirrors `valenx-springs`' `SpringsError` and
//! `valenx-popgen`'s `PopgenError`.

use thiserror::Error;

/// Errors produced by `valenx-fatigue`.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum FatigueError {
    /// A single argument is outside the domain the model accepts: a
    /// non-positive stress amplitude, a negative cycle count, a Basquin
    /// strength coefficient that is not strictly positive, or a Basquin
    /// exponent that is not strictly negative. A property of one *value*.
    #[error("invalid `{what}`: {reason}")]
    Invalid {
        /// Logical parameter name (e.g. `"stress_amplitude"`, `"b"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// Two inputs are individually valid but jointly inconsistent: an
    /// alternating stress at or above the ultimate tensile strength, a
    /// mean stress that reaches the static strength so no alternating
    /// stress is admissible, or a requested design factor below one.
    #[error("out of domain: {reason}")]
    Domain {
        /// Human-readable reason.
        reason: String,
    },

    /// Two parallel lists disagree on length — a per-block applied-cycle
    /// vector and a per-block capacity vector of different sizes.
    #[error("dimension mismatch for {context}: expected {expected}, got {actual}")]
    Dimension {
        /// What was expected.
        expected: usize,
        /// What was actually supplied.
        actual: usize,
        /// Short context label (e.g. `"damage blocks"`).
        context: &'static str,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this rather
/// than on every error variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A supplied value is outside its physical domain.
    Input,
    /// Two otherwise-valid inputs are mutually inconsistent.
    Domain,
}

impl FatigueError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"fatigue.<sub_id>"`. Codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            FatigueError::Invalid { .. } => "fatigue.invalid",
            FatigueError::Domain { .. } => "fatigue.domain",
            FatigueError::Dimension { .. } => "fatigue.dimension",
        }
    }

    /// Coarse category — see [`ErrorCategory`].
    pub fn category(&self) -> ErrorCategory {
        match self {
            FatigueError::Invalid { .. } | FatigueError::Dimension { .. } => ErrorCategory::Input,
            FatigueError::Domain { .. } => ErrorCategory::Domain,
        }
    }

    /// Convenience constructor for [`FatigueError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        FatigueError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`FatigueError::Domain`].
    pub fn domain(reason: impl Into<String>) -> Self {
        FatigueError::Domain {
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`FatigueError::Dimension`].
    pub fn dimension(expected: usize, actual: usize, context: &'static str) -> Self {
        FatigueError::Dimension {
            expected,
            actual,
            context,
        }
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, FatigueError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = FatigueError::invalid("stress_amplitude", "must be positive");
        assert_eq!(err.code(), "fatigue.invalid");
        assert_eq!(err.category(), ErrorCategory::Input);

        let err = FatigueError::domain("mean stress reaches ultimate");
        assert_eq!(err.code(), "fatigue.domain");
        assert_eq!(err.category(), ErrorCategory::Domain);

        let err = FatigueError::dimension(3, 2, "damage blocks");
        assert_eq!(err.code(), "fatigue.dimension");
        assert_eq!(err.category(), ErrorCategory::Input);
    }

    #[test]
    fn display_is_informative() {
        let msg = FatigueError::dimension(4, 9, "damage blocks").to_string();
        assert!(msg.contains('4') && msg.contains('9'), "got: {msg}");
        assert!(msg.contains("damage blocks"), "got: {msg}");

        let msg = FatigueError::invalid("b", "must be negative").to_string();
        assert!(msg.contains('b'), "got: {msg}");
        assert!(msg.contains("must be negative"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(FatigueError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }
}
