//! Springs workbench error taxonomy.

use thiserror::Error;

/// Errors raised by spring generation.
#[derive(Debug, Error)]
pub enum SpringsError {
    /// Bad parameter (non-positive, etc).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Geometric inversion (e.g. wire bigger than coil).
    #[error("degenerate spring: {0}")]
    Degenerate(String),
}

/// Coarse error category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// Tunable knob.
    Config,
    /// Algorithm domain.
    Algorithm,
}

impl SpringsError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            SpringsError::BadParameter { .. } => "springs.bad_parameter",
            SpringsError::Degenerate(_) => "springs.degenerate",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            SpringsError::BadParameter { .. } => ErrorCategory::Config,
            SpringsError::Degenerate(_) => ErrorCategory::Algorithm,
        }
    }
}
