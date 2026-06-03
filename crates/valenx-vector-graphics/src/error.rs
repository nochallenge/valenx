//! Typed errors for `valenx-vector-graphics`.

use thiserror::Error;

/// Errors raised by the SVG-style primitive layer.
#[derive(Debug, Error)]
pub enum VectorError {
    /// Bad parameter.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// SVG parse error.
    #[error("svg parse error at byte {byte_offset}: {message}")]
    Parse {
        /// Byte offset in the input.
        byte_offset: usize,
        /// Diagnostic message.
        message: String,
    },

    /// Empty path / polygon / etc.
    #[error("empty: {0}")]
    Empty(&'static str),
}

/// Coarse category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// SVG format.
    Format,
}

impl VectorError {
    /// Stable kebab code.
    pub fn code(&self) -> &'static str {
        match self {
            VectorError::BadParameter { .. } => "vector.bad_parameter",
            VectorError::Parse { .. } => "vector.parse",
            VectorError::Empty(_) => "vector.empty",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            VectorError::BadParameter { .. } | VectorError::Empty(_) => ErrorCategory::Input,
            VectorError::Parse { .. } => ErrorCategory::Format,
        }
    }
}
