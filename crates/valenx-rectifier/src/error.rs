//! Error taxonomy for the rectifier closed-form models.
//!
//! Every fallible entry point validates its inputs and returns a
//! [`RectifierError`] rather than panicking, so callers can surface a
//! stable [`RectifierError::code`] / [`RectifierError::category`] pair.

use thiserror::Error;

/// Errors raised when validating rectifier inputs.
///
/// The variants are constructed exclusively through the checked helpers
/// [`RectifierError::positive`] and [`RectifierError::non_negative`] so
/// the carried message is always consistent.
#[derive(Debug, Error)]
pub enum RectifierError {
    /// A quantity that must be strictly positive was zero or negative.
    #[error("parameter `{name}` must be > 0, got {value}")]
    NonPositive {
        /// Offending parameter name (stable, `snake_case`).
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A quantity that must be finite was `NaN` or infinite.
    #[error("parameter `{name}` must be finite, got {value}")]
    NotFinite {
        /// Offending parameter name (stable, `snake_case`).
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A quantity that must be zero-or-positive was negative.
    #[error("parameter `{name}` must be >= 0, got {value}")]
    Negative {
        /// Offending parameter name (stable, `snake_case`).
        name: &'static str,
        /// The rejected value.
        value: f64,
    },
}

/// Coarse classification of a [`RectifierError`], for UI grouping.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The supplied value is outside the model's valid domain.
    Domain,
}

impl RectifierError {
    /// Validate that `value` is finite and strictly positive.
    ///
    /// Returns `Ok(value)` on success so the check composes inline:
    /// `let v = RectifierError::positive("v_peak", v)?;`.
    ///
    /// # Errors
    ///
    /// Returns [`RectifierError::NotFinite`] when `value` is `NaN` or
    /// infinite, and [`RectifierError::NonPositive`] when `value <= 0`.
    pub fn positive(name: &'static str, value: f64) -> Result<f64, RectifierError> {
        if !value.is_finite() {
            return Err(RectifierError::NotFinite { name, value });
        }
        if value <= 0.0 {
            return Err(RectifierError::NonPositive { name, value });
        }
        Ok(value)
    }

    /// Validate that `value` is finite and zero-or-positive.
    ///
    /// Returns `Ok(value)` on success.
    ///
    /// # Errors
    ///
    /// Returns [`RectifierError::NotFinite`] when `value` is `NaN` or
    /// infinite, and [`RectifierError::Negative`] when `value < 0`.
    pub fn non_negative(name: &'static str, value: f64) -> Result<f64, RectifierError> {
        if !value.is_finite() {
            return Err(RectifierError::NotFinite { name, value });
        }
        if value < 0.0 {
            return Err(RectifierError::Negative { name, value });
        }
        Ok(value)
    }

    /// Stable, kebab-cased identifier for the variant.
    ///
    /// Suitable for logs, telemetry, and localization keys; the string is
    /// part of the crate's public contract.
    pub fn code(&self) -> &'static str {
        match self {
            RectifierError::NonPositive { .. } => "rectifier.non-positive",
            RectifierError::NotFinite { .. } => "rectifier.not-finite",
            RectifierError::Negative { .. } => "rectifier.negative",
        }
    }

    /// Coarse [`ErrorCategory`] for the variant.
    pub fn category(&self) -> ErrorCategory {
        match self {
            RectifierError::NonPositive { .. }
            | RectifierError::NotFinite { .. }
            | RectifierError::Negative { .. } => ErrorCategory::Domain,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn positive_accepts_positive() {
        let v = RectifierError::positive("x", 3.0).expect("3 is positive");
        assert!((v - 3.0).abs() < 1e-12);
    }

    #[test]
    fn positive_rejects_zero_and_negative() {
        let zero = RectifierError::positive("x", 0.0).unwrap_err();
        assert!(matches!(zero, RectifierError::NonPositive { .. }));
        let neg = RectifierError::positive("x", -1.0).unwrap_err();
        assert!(matches!(neg, RectifierError::NonPositive { .. }));
    }

    #[test]
    fn positive_rejects_non_finite() {
        let nan = RectifierError::positive("x", f64::NAN).unwrap_err();
        assert!(matches!(nan, RectifierError::NotFinite { .. }));
        let inf = RectifierError::positive("x", f64::INFINITY).unwrap_err();
        assert!(matches!(inf, RectifierError::NotFinite { .. }));
    }

    #[test]
    fn non_negative_accepts_zero() {
        let v = RectifierError::non_negative("x", 0.0).expect("0 is non-negative");
        assert!(v.abs() < 1e-12);
    }

    #[test]
    fn non_negative_rejects_negative() {
        let neg = RectifierError::non_negative("x", -0.5).unwrap_err();
        assert!(matches!(neg, RectifierError::Negative { .. }));
    }

    #[test]
    fn code_and_category_are_stable() {
        let e = RectifierError::NonPositive {
            name: "x",
            value: -1.0,
        };
        assert_eq!(e.code(), "rectifier.non-positive");
        assert_eq!(e.category(), ErrorCategory::Domain);

        let f = RectifierError::NotFinite {
            name: "x",
            value: f64::NAN,
        };
        assert_eq!(f.code(), "rectifier.not-finite");

        let g = RectifierError::Negative {
            name: "x",
            value: -1.0,
        };
        assert_eq!(g.code(), "rectifier.negative");
    }

    #[test]
    fn display_includes_name_and_value() {
        let e = RectifierError::NonPositive {
            name: "v_peak",
            value: -2.0,
        };
        let msg = format!("{e}");
        assert!(msg.contains("v_peak"), "message was: {msg}");
        assert!(msg.contains("-2"), "message was: {msg}");
    }
}
