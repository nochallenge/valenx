//! Typed errors for `valenx-librecad-2d`.

use thiserror::Error;

/// Errors raised by the LibreCAD-style 2D workbench.
#[derive(Debug, Error)]
pub enum LibreCadError {
    /// Bad parameter.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// I/O error wrapping the underlying message (file paths handled
    /// by the caller; this crate never opens a [`std::fs`] handle in
    /// `lib.rs` directly, but `dxf::write_full` does).
    #[error("io: {0}")]
    Io(String),

    /// DXF parse error — usually a malformed group-code pair or
    /// unknown entity kind.
    #[error("dxf parse error at line {line}: {message}")]
    DxfParse {
        /// Line number (1-based) where the parser tripped.
        line: usize,
        /// Diagnostic message.
        message: String,
    },

    /// Round-trip semantic compare failed — entity count or per-entity
    /// data mismatch.
    #[error("round-trip mismatch: {0}")]
    RoundTripMismatch(String),
}

/// Coarse category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// I/O.
    Io,
    /// Parser / format.
    Format,
}

impl LibreCadError {
    /// Stable kebab code.
    pub fn code(&self) -> &'static str {
        match self {
            LibreCadError::BadParameter { .. } => "librecad.bad_parameter",
            LibreCadError::Io(_) => "librecad.io",
            LibreCadError::DxfParse { .. } => "librecad.dxf_parse",
            LibreCadError::RoundTripMismatch(_) => "librecad.round_trip",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            LibreCadError::BadParameter { .. } => ErrorCategory::Input,
            LibreCadError::Io(_) => ErrorCategory::Io,
            LibreCadError::DxfParse { .. } | LibreCadError::RoundTripMismatch(_) => {
                ErrorCategory::Format
            }
        }
    }
}
