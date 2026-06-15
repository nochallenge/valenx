//! DC-motor model error taxonomy.

use thiserror::Error;

/// Errors raised when constructing or evaluating a [`crate::motor::DcMotor`].
#[derive(Debug, Error, Clone, PartialEq)]
pub enum DcMotorError {
    /// A scalar parameter was non-finite (`NaN` or infinite).
    #[error("parameter `{name}` must be finite, got {value}")]
    NotFinite {
        /// Parameter name.
        name: &'static str,
        /// Offending value.
        value: f64,
    },

    /// A scalar parameter was outside its allowed range.
    #[error("parameter `{name}` {reason}, got {value}")]
    OutOfRange {
        /// Parameter name.
        name: &'static str,
        /// Human-readable constraint, e.g. `"must be > 0"`.
        reason: &'static str,
        /// Offending value.
        value: f64,
    },
}

/// Coarse classification used by callers that group errors for the UI.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The value is structurally invalid (`NaN`/infinite).
    Domain,
    /// The value is finite but violates a physical constraint.
    Range,
}

impl DcMotorError {
    /// Stable, dot-separated identifier suitable for logs and tests.
    pub fn code(&self) -> &'static str {
        match self {
            DcMotorError::NotFinite { .. } => "dcmotor.not_finite",
            DcMotorError::OutOfRange { .. } => "dcmotor.out_of_range",
        }
    }

    /// Coarse category for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            DcMotorError::NotFinite { .. } => ErrorCategory::Domain,
            DcMotorError::OutOfRange { .. } => ErrorCategory::Range,
        }
    }

    /// Validate that `value` is finite, returning [`DcMotorError::NotFinite`]
    /// otherwise. Used internally by the validated constructors.
    pub(crate) fn require_finite(name: &'static str, value: f64) -> Result<(), DcMotorError> {
        if value.is_finite() {
            Ok(())
        } else {
            Err(DcMotorError::NotFinite { name, value })
        }
    }

    /// Validate that `value` is finite and strictly positive.
    pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<(), DcMotorError> {
        DcMotorError::require_finite(name, value)?;
        if value > 0.0 {
            Ok(())
        } else {
            Err(DcMotorError::OutOfRange {
                name,
                reason: "must be > 0",
                value,
            })
        }
    }
}
