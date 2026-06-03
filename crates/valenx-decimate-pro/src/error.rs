//! Decimate-Pro error taxonomy.

use thiserror::Error;

/// Errors raised by curvature / UV / feature-aware decimation.
#[derive(Debug, Error)]
pub enum DecimateProError {
    /// Bad parameter.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Mismatch between mesh size and a supplied side channel
    /// (UVs / feature mask / curvature override).
    #[error("size mismatch on `{name}`: mesh has {mesh}, got {got}")]
    SizeMismatch {
        /// Channel name.
        name: &'static str,
        /// Mesh count.
        mesh: usize,
        /// Provided count.
        got: usize,
    },
}

/// Coarse error category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// Algorithm domain.
    Algorithm,
}

impl DecimateProError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            DecimateProError::BadParameter { .. } => "decimate_pro.bad_parameter",
            DecimateProError::SizeMismatch { .. } => "decimate_pro.size_mismatch",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            DecimateProError::BadParameter { .. } | DecimateProError::SizeMismatch { .. } => {
                ErrorCategory::Input
            }
        }
    }
}
