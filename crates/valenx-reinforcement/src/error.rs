//! Reinforcement workbench error taxonomy.

use thiserror::Error;

/// Errors raised by reinforcement cage generation.
#[derive(Debug, Error)]
pub enum ReinforcementError {
    /// Bad parameter (negative dimension, etc).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// IO error.
    #[error("io: {0}")]
    Io(#[from] std::io::Error),

    /// RON error.
    #[error("ron: {0}")]
    Ron(String),
}

/// Coarse error category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// Tunable knob.
    Config,
    /// Transient / IO / parse.
    Runtime,
}

impl ReinforcementError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            ReinforcementError::BadParameter { .. } => "reinforcement.bad_parameter",
            ReinforcementError::Io(_) => "reinforcement.io",
            ReinforcementError::Ron(_) => "reinforcement.ron",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            ReinforcementError::BadParameter { .. } => ErrorCategory::Config,
            ReinforcementError::Io(_) | ReinforcementError::Ron(_) => ErrorCategory::Runtime,
        }
    }
}
