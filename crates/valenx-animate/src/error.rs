//! Animation workbench error taxonomy.

use thiserror::Error;

/// Errors raised by animation construction or playback.
#[derive(Debug, Error)]
pub enum AnimateError {
    /// Animation has no keyframes — can't sample.
    #[error("empty animation")]
    Empty,

    /// Keyframes are not in monotonically-increasing time order.
    #[error("keyframes must be in monotonic time order ({reason})")]
    NotMonotonic {
        /// Detail string.
        reason: String,
    },

    /// User parameter out of range.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Offending parameter.
        name: &'static str,
        /// Human-readable explanation.
        reason: String,
    },

    /// Wraps an upstream assembly error from `apply_all_joints`.
    #[error("assembly: {0}")]
    Assembly(String),

    /// IO error.
    #[error("io: {0}")]
    Io(#[from] std::io::Error),

    /// RON error.
    #[error("ron: {0}")]
    Ron(String),
}

impl From<valenx_assembly::AssemblyError> for AnimateError {
    fn from(e: valenx_assembly::AssemblyError) -> Self {
        AnimateError::Assembly(e.to_string())
    }
}

/// Coarse category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// Tunable knob.
    Config,
    /// Bug in upstream assembly.
    Internal,
    /// Transient IO / parse.
    Runtime,
}

impl AnimateError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            AnimateError::Empty => "animate.empty",
            AnimateError::NotMonotonic { .. } => "animate.not_monotonic",
            AnimateError::BadParameter { .. } => "animate.bad_parameter",
            AnimateError::Assembly(_) => "animate.assembly",
            AnimateError::Io(_) => "animate.io",
            AnimateError::Ron(_) => "animate.ron",
        }
    }

    /// High-level classification.
    pub fn category(&self) -> ErrorCategory {
        match self {
            AnimateError::Empty | AnimateError::NotMonotonic { .. } => ErrorCategory::Input,
            AnimateError::BadParameter { .. } => ErrorCategory::Config,
            AnimateError::Assembly(_) => ErrorCategory::Internal,
            AnimateError::Io(_) | AnimateError::Ron(_) => ErrorCategory::Runtime,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_and_cats() {
        assert_eq!(AnimateError::Empty.code(), "animate.empty");
        assert_eq!(
            AnimateError::BadParameter {
                name: "fps",
                reason: "0".into()
            }
            .category(),
            ErrorCategory::Config
        );
    }
}
