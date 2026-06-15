//! Lever-mechanics error taxonomy.
//!
//! All fallible constructors and computations in this crate funnel
//! through [`LeverError`]. Each variant carries a stable kebab-cased
//! [`LeverError::code`] and a coarse [`ErrorCategory`] so callers can
//! branch on the failure class without string-matching the display
//! message.

use thiserror::Error;

/// Errors raised while constructing or evaluating a lever.
#[derive(Debug, Error)]
pub enum LeverError {
    /// A geometric arm length was not a finite, strictly positive value.
    ///
    /// Arm lengths are measured from the fulcrum to the line of action
    /// of the effort or load force; a zero, negative, NaN, or infinite
    /// arm has no physical meaning for an ideal rigid lever.
    #[error("non-positive arm `{name}`: expected a finite value > 0, got {value}")]
    NonPositiveArm {
        /// Which arm was rejected (`"effort_arm"` or `"load_arm"`).
        name: &'static str,
        /// The offending value, echoed for diagnostics.
        value: f64,
    },

    /// A force magnitude was not finite (NaN or infinite).
    ///
    /// Forces may be zero (a lever at rest under no effort), but they
    /// must be finite for the moment-balance algebra to be meaningful.
    #[error("non-finite force `{name}`: expected a finite value, got {value}")]
    NonFiniteForce {
        /// Which force was rejected (`"effort"` or `"load"`).
        name: &'static str,
        /// The offending value, echoed for diagnostics.
        value: f64,
    },
}

/// Coarse classification of a [`LeverError`].
///
/// Mirrors the category split used across the Valenx workspace so a
/// caller can present a uniform "bad input vs. internal" distinction.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The caller supplied an out-of-domain quantity (bad arm or force).
    Input,
}

impl LeverError {
    /// Stable, kebab-cased identifier for this error.
    ///
    /// Suitable for logs, metrics labels, and machine-readable error
    /// surfaces; stable across releases independent of the human-facing
    /// [`Display`](std::fmt::Display) text.
    pub fn code(&self) -> &'static str {
        match self {
            LeverError::NonPositiveArm { .. } => "leverage.non-positive-arm",
            LeverError::NonFiniteForce { .. } => "leverage.non-finite-force",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    ///
    /// Every variant currently originates from caller-supplied input, so
    /// this always returns [`ErrorCategory::Input`]; the method exists so
    /// downstream match arms stay exhaustive if richer categories are
    /// added later.
    pub fn category(&self) -> ErrorCategory {
        match self {
            LeverError::NonPositiveArm { .. } | LeverError::NonFiniteForce { .. } => {
                ErrorCategory::Input
            }
        }
    }
}

/// Validate that an arm length is finite and strictly positive.
///
/// Returns the value unchanged on success; otherwise a
/// [`LeverError::NonPositiveArm`] tagged with `name`.
pub(crate) fn validate_arm(name: &'static str, value: f64) -> Result<f64, LeverError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(LeverError::NonPositiveArm { name, value })
    }
}

/// Validate that a force magnitude is finite (zero is allowed).
///
/// Returns the value unchanged on success; otherwise a
/// [`LeverError::NonFiniteForce`] tagged with `name`.
pub(crate) fn validate_force(name: &'static str, value: f64) -> Result<f64, LeverError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(LeverError::NonFiniteForce { name, value })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arm_rejects_zero_negative_and_nonfinite() {
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let err = validate_arm("effort_arm", bad).unwrap_err();
            assert!(matches!(err, LeverError::NonPositiveArm { name, .. } if name == "effort_arm"));
            assert_eq!(err.code(), "leverage.non-positive-arm");
            assert_eq!(err.category(), ErrorCategory::Input);
        }
    }

    #[test]
    fn arm_accepts_positive_finite() {
        let got = validate_arm("load_arm", 2.5).unwrap();
        assert!((got - 2.5).abs() < 1e-12);
    }

    #[test]
    fn force_rejects_only_nonfinite() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let err = validate_force("load", bad).unwrap_err();
            assert!(matches!(err, LeverError::NonFiniteForce { name, .. } if name == "load"));
            assert_eq!(err.code(), "leverage.non-finite-force");
        }
    }

    #[test]
    fn force_accepts_zero_and_negative_finite() {
        // Zero effort (lever at rest) and a sign-bearing force are both fine.
        assert!((validate_force("effort", 0.0).unwrap() - 0.0).abs() < 1e-12);
        assert!((validate_force("effort", -3.0).unwrap() + 3.0).abs() < 1e-12);
    }

    #[test]
    fn display_text_mentions_the_named_quantity() {
        let arm = LeverError::NonPositiveArm {
            name: "effort_arm",
            value: 0.0,
        };
        assert!(arm.to_string().contains("effort_arm"));
        let force = LeverError::NonFiniteForce {
            name: "load",
            value: f64::NAN,
        };
        assert!(force.to_string().contains("load"));
    }
}
