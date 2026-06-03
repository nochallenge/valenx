//! Frames workbench error taxonomy.

use thiserror::Error;

/// Errors raised by frame ops.
#[derive(Debug, Error)]
pub enum FramesError {
    /// Bad parameter (non-positive dimension, etc).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Path needs at least 2 vertices.
    #[error("degenerate path: {0}")]
    DegeneratePath(String),

    /// Index out of range.
    #[error("bad index {got} (have {n})")]
    BadIndex {
        /// Asked-for id.
        got: usize,
        /// Available entries.
        n: usize,
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
    /// Algorithm domain.
    Algorithm,
    /// Transient.
    Runtime,
}

impl FramesError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            FramesError::BadParameter { .. } => "frames.bad_parameter",
            FramesError::DegeneratePath(_) => "frames.degenerate_path",
            FramesError::BadIndex { .. } => "frames.bad_index",
            FramesError::Io(_) => "frames.io",
            FramesError::Ron(_) => "frames.ron",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            FramesError::BadParameter { .. } => ErrorCategory::Config,
            FramesError::DegeneratePath(_) | FramesError::BadIndex { .. } => {
                ErrorCategory::Algorithm
            }
            FramesError::Io(_) | FramesError::Ron(_) => ErrorCategory::Runtime,
        }
    }
}
