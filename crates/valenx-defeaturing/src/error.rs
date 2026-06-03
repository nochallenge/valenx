//! Defeaturing workbench error taxonomy.

use thiserror::Error;

/// Errors raised by defeaturing.
#[derive(Debug, Error)]
pub enum DefeatureError {
    /// Bad parameter.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Underlying tessellation failed.
    #[error("tessellation: {0}")]
    Tessellation(String),
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

impl DefeatureError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            DefeatureError::BadParameter { .. } => "defeature.bad_parameter",
            DefeatureError::Tessellation(_) => "defeature.tessellation",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            DefeatureError::BadParameter { .. } => ErrorCategory::Config,
            DefeatureError::Tessellation(_) => ErrorCategory::Algorithm,
        }
    }
}
