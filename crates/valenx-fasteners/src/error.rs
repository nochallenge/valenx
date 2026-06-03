//! Fasteners workbench error taxonomy.

use thiserror::Error;

/// Errors raised by fastener generation.
#[derive(Debug, Error)]
pub enum FastenerError {
    /// Bad parameter (negative length, etc).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Unknown size lookup.
    #[error("unknown size `{0}`")]
    UnknownSize(String),
}

/// Coarse error category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// Tunable knob.
    Config,
}

impl FastenerError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            FastenerError::BadParameter { .. } => "fasteners.bad_parameter",
            FastenerError::UnknownSize(_) => "fasteners.unknown_size",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            FastenerError::BadParameter { .. } => ErrorCategory::Config,
            FastenerError::UnknownSize(_) => ErrorCategory::Input,
        }
    }
}
