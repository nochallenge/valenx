//! Piping error taxonomy.

use thiserror::Error;

/// Errors raised by the piping crate.
#[derive(Debug, Error)]
pub enum PipingError {
    /// Bad parameter.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// NPS designation not in the table.
    #[error("unknown NPS `{0}`")]
    UnknownNps(String),

    /// CAD kernel error wrapped.
    #[error("cad: {0}")]
    Cad(String),
}

/// Coarse error category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// Algorithm.
    Algorithm,
    /// Backend.
    Backend,
}

impl PipingError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            PipingError::BadParameter { .. } => "piping.bad_parameter",
            PipingError::UnknownNps(_) => "piping.unknown_nps",
            PipingError::Cad(_) => "piping.cad",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            PipingError::BadParameter { .. } => ErrorCategory::Input,
            PipingError::UnknownNps(_) => ErrorCategory::Algorithm,
            PipingError::Cad(_) => ErrorCategory::Backend,
        }
    }
}
