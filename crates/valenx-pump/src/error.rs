//! Error taxonomy for centrifugal-pump calculations.
//!
//! Every fallible constructor and free function in the crate returns
//! [`PumpError`]. The variants distinguish a bad scalar input
//! ([`PumpError::BadParameter`]) from a request that is dimensionally or
//! physically inconsistent ([`PumpError::Inconsistent`]) — for example a
//! system whose static head already exceeds the pump's shut-off head, so
//! no operating point exists.

use thiserror::Error;

/// Errors raised by pump hydraulics.
#[derive(Debug, Error)]
pub enum PumpError {
    /// A scalar parameter was outside its valid domain (non-positive
    /// where a positive value is required, a negative resistance
    /// coefficient, an efficiency outside `(0, 1]`, and so on).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Name of the offending parameter.
        name: &'static str,
        /// Human-readable explanation of why it was rejected.
        reason: String,
    },

    /// The inputs were each individually valid but cannot be satisfied
    /// together — e.g. a pump and system whose head curves never cross
    /// at a non-negative flow.
    #[error("inconsistent inputs: {0}")]
    Inconsistent(String),
}

/// Coarse category for an error, for callers that bucket failures.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The user supplied an out-of-domain value.
    Input,
    /// The model has no physical solution for these inputs.
    Algorithm,
}

impl PumpError {
    /// A stable, kebab-cased identifier for this error, suitable for
    /// logging or matching in tests without depending on the display
    /// message.
    pub fn code(&self) -> &'static str {
        match self {
            PumpError::BadParameter { .. } => "pump.bad_parameter",
            PumpError::Inconsistent(_) => "pump.inconsistent",
        }
    }

    /// The coarse [`ErrorCategory`] this error belongs to.
    pub fn category(&self) -> ErrorCategory {
        match self {
            PumpError::BadParameter { .. } => ErrorCategory::Input,
            PumpError::Inconsistent(_) => ErrorCategory::Algorithm,
        }
    }

    /// Build a [`PumpError::BadParameter`] from a parameter name and a
    /// reason. A small internal convenience used by the validating
    /// constructors throughout the crate.
    pub(crate) fn bad(name: &'static str, reason: impl Into<String>) -> Self {
        PumpError::BadParameter {
            name,
            reason: reason.into(),
        }
    }
}

/// Internal helper: reject `value` unless it is finite and strictly
/// positive, attributing any failure to the parameter `name`.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64, PumpError> {
    if !value.is_finite() {
        return Err(PumpError::bad(name, format!("must be finite, got {value}")));
    }
    if value <= 0.0 {
        return Err(PumpError::bad(name, format!("must be > 0, got {value}")));
    }
    Ok(value)
}

/// Internal helper: reject `value` unless it is finite and non-negative,
/// attributing any failure to the parameter `name`.
pub(crate) fn require_non_negative(name: &'static str, value: f64) -> Result<f64, PumpError> {
    if !value.is_finite() {
        return Err(PumpError::bad(name, format!("must be finite, got {value}")));
    }
    if value < 0.0 {
        return Err(PumpError::bad(name, format!("must be >= 0, got {value}")));
    }
    Ok(value)
}

/// Internal helper: reject `value` unless it is finite (used for signed
/// quantities such as a static head that may legitimately be negative).
pub(crate) fn require_finite(name: &'static str, value: f64) -> Result<f64, PumpError> {
    if !value.is_finite() {
        return Err(PumpError::bad(name, format!("must be finite, got {value}")));
    }
    Ok(value)
}
