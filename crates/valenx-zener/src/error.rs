//! Zener-regulator error taxonomy.

use thiserror::Error;

/// Errors raised when constructing or analysing a zener regulator.
#[derive(Debug, Error)]
pub enum ZenerError {
    /// A parameter was non-finite, negative, or otherwise outside its
    /// admissible range.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Offending parameter name.
        name: &'static str,
        /// Human-readable reason the value was rejected.
        reason: String,
    },

    /// The diode cannot regulate at the requested operating point — for
    /// example the supply does not exceed the zener voltage, so no
    /// headroom remains across the series resistor.
    #[error("regulator cannot operate: {0}")]
    CannotRegulate(String),
}

/// Coarse category for routing/telemetry of a [`ZenerError`].
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// Caller-supplied input was invalid.
    Input,
    /// The requested operating point is physically unreachable.
    Operating,
}

impl ZenerError {
    /// Stable kebab-cased identifier for logs and tests.
    pub fn code(&self) -> &'static str {
        match self {
            ZenerError::BadParameter { .. } => "zener.bad_parameter",
            ZenerError::CannotRegulate(_) => "zener.cannot_regulate",
        }
    }

    /// Coarse category for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            ZenerError::BadParameter { .. } => ErrorCategory::Input,
            ZenerError::CannotRegulate(_) => ErrorCategory::Operating,
        }
    }
}

/// Internal helper: reject a value that is not finite and strictly
/// positive, attributing the failure to `name`.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64, ZenerError> {
    if !value.is_finite() {
        return Err(ZenerError::BadParameter {
            name,
            reason: format!("must be finite, got {value}"),
        });
    }
    if value <= 0.0 {
        return Err(ZenerError::BadParameter {
            name,
            reason: format!("must be > 0, got {value}"),
        });
    }
    Ok(value)
}

/// Internal helper: reject a value that is not finite and non-negative,
/// attributing the failure to `name`. Used for quantities (e.g. load
/// current) where exactly zero is physically valid.
pub(crate) fn require_non_negative(name: &'static str, value: f64) -> Result<f64, ZenerError> {
    if !value.is_finite() {
        return Err(ZenerError::BadParameter {
            name,
            reason: format!("must be finite, got {value}"),
        });
    }
    if value < 0.0 {
        return Err(ZenerError::BadParameter {
            name,
            reason: format!("must be >= 0, got {value}"),
        });
    }
    Ok(value)
}
