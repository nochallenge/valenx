//! Threaded-profiles error taxonomy.

use thiserror::Error;

/// Errors raised by thread table lookups and helix sweeps.
#[derive(Debug, Error)]
pub enum ThreadsProError {
    /// Bad parameter.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Designation not found in any thread table.
    #[error("designation `{0}` not found in any thread table")]
    UnknownDesignation(String),
}

/// Coarse error category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// Algorithm domain.
    Algorithm,
}

impl ThreadsProError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            ThreadsProError::BadParameter { .. } => "threads_pro.bad_parameter",
            ThreadsProError::UnknownDesignation(_) => "threads_pro.unknown_designation",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            ThreadsProError::BadParameter { .. } => ErrorCategory::Input,
            ThreadsProError::UnknownDesignation(_) => ErrorCategory::Algorithm,
        }
    }
}
