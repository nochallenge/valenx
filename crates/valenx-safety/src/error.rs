//! Error taxonomy for safety consolidation.

use thiserror::Error;

/// Errors raised while building flags or a risk report.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum SafetyError {
    /// A required text field was empty.
    #[error("empty {what}")]
    Empty {
        /// What was empty.
        what: &'static str,
    },

    /// A value that must be finite was `NaN` or infinite.
    #[error("non-finite {what}")]
    NonFinite {
        /// What was non-finite.
        what: &'static str,
    },

    /// A calibrated confidence fell outside `[0, 1]`.
    #[error("confidence {value} is outside [0, 1]")]
    ConfidenceOutOfRange {
        /// The offending confidence.
        value: f64,
    },
}

impl SafetyError {
    /// A short, stable machine-readable code for this error.
    pub fn code(&self) -> &'static str {
        match self {
            SafetyError::Empty { .. } => "empty",
            SafetyError::NonFinite { .. } => "non_finite",
            SafetyError::ConfidenceOutOfRange { .. } => "confidence_out_of_range",
        }
    }
}

pub(crate) fn check_finite(what: &'static str, v: f64) -> Result<(), SafetyError> {
    if v.is_finite() {
        Ok(())
    } else {
        Err(SafetyError::NonFinite { what })
    }
}
