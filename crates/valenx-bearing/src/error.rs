//! Error taxonomy for the bearing rating-life calculator.
//!
//! Every fallible constructor and calculation in this crate returns a
//! [`BearingError`]. The variants distinguish *invalid user input*
//! (a non-positive load, a negative speed) from a *domain violation*
//! that the formulae cannot evaluate. Each error carries a stable,
//! kebab-cased [`BearingError::code`] and a coarse
//! [`BearingError::category`] so callers can branch on the failure
//! class without string-matching the human-readable message.

use thiserror::Error;

/// Errors raised while validating inputs or evaluating bearing-life
/// formulae.
#[derive(Debug, Error)]
pub enum BearingError {
    /// A numeric parameter was outside its valid domain.
    ///
    /// `name` is the stable parameter identifier (for example
    /// `"dynamic_load_rating"`); `reason` explains the violation.
    #[error("invalid parameter `{name}`: {reason}")]
    InvalidParameter {
        /// Stable parameter identifier.
        name: &'static str,
        /// Human-readable explanation of what was wrong.
        reason: String,
    },

    /// A parameter was required to be a finite number but was `NaN`
    /// or infinite.
    #[error("parameter `{name}` must be finite, got {value}")]
    NotFinite {
        /// Stable parameter identifier.
        name: &'static str,
        /// The offending value (may be `NaN` or `±inf`).
        value: f64,
    },
}

/// Coarse classification of a [`BearingError`], for telemetry and for
/// callers that want to react to a *class* of failure rather than a
/// specific variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// The caller supplied an invalid value (non-positive load,
    /// negative speed, non-finite input).
    Input,
    /// The requested calculation is outside the model's analytic
    /// domain.
    Domain,
}

impl BearingError {
    /// Stable, kebab-cased identifier for this error, suitable for
    /// logging and equality checks across releases.
    ///
    /// ```
    /// use valenx_bearing::BearingError;
    /// let e = BearingError::InvalidParameter {
    ///     name: "rpm",
    ///     reason: "must be positive".to_string(),
    /// };
    /// assert_eq!(e.code(), "bearing.invalid-parameter");
    /// ```
    #[must_use]
    pub fn code(&self) -> &'static str {
        match self {
            BearingError::InvalidParameter { .. } => "bearing.invalid-parameter",
            BearingError::NotFinite { .. } => "bearing.not-finite",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    #[must_use]
    pub fn category(&self) -> ErrorCategory {
        match self {
            BearingError::InvalidParameter { .. } | BearingError::NotFinite { .. } => {
                ErrorCategory::Input
            }
        }
    }
}

/// Validate that `value` is strictly positive and finite, returning it
/// unchanged on success.
///
/// This is the workhorse guard behind the validated constructors and
/// free functions: loads, ratings, exponents and speeds are all
/// physically meaningless at zero or below.
///
/// # Errors
///
/// Returns [`BearingError::NotFinite`] when `value` is `NaN` or
/// infinite, and [`BearingError::InvalidParameter`] when it is finite
/// but not greater than zero.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64, BearingError> {
    if !value.is_finite() {
        return Err(BearingError::NotFinite { name, value });
    }
    if value <= 0.0 {
        return Err(BearingError::InvalidParameter {
            name,
            reason: format!("must be greater than zero, got {value}"),
        });
    }
    Ok(value)
}

/// Validate that `value` is non-negative and finite, returning it
/// unchanged on success.
///
/// Used for quantities that may legitimately be zero — for example a
/// purely radial load has zero axial component.
///
/// # Errors
///
/// Returns [`BearingError::NotFinite`] when `value` is `NaN` or
/// infinite, and [`BearingError::InvalidParameter`] when it is finite
/// but negative.
pub(crate) fn require_non_negative(name: &'static str, value: f64) -> Result<f64, BearingError> {
    if !value.is_finite() {
        return Err(BearingError::NotFinite { name, value });
    }
    if value < 0.0 {
        return Err(BearingError::InvalidParameter {
            name,
            reason: format!("must not be negative, got {value}"),
        });
    }
    Ok(value)
}
