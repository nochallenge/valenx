//! Collision detection error taxonomy.

use thiserror::Error;

/// Errors raised by collision checks.
#[derive(Debug, Error)]
pub enum CollisionError {
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

impl CollisionError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            CollisionError::BadParameter { .. } => "collision.bad_parameter",
            CollisionError::Tessellation(_) => "collision.tessellation",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            CollisionError::BadParameter { .. } => ErrorCategory::Config,
            CollisionError::Tessellation(_) => ErrorCategory::Algorithm,
        }
    }
}
