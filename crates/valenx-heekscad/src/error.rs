//! Typed errors for `valenx-heekscad`.

use thiserror::Error;

/// Errors raised by the HeeksCAD primitives + CAM ops.
#[derive(Debug, Error)]
pub enum HeeksCadError {
    /// Bad parameter.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Layer / object lookup miss.
    #[error("unknown {kind} `{name}`")]
    Unknown {
        /// `"layer"` / `"object"` / `"tool"`.
        kind: &'static str,
        /// Name that didn't resolve.
        name: String,
    },

    /// IO.
    #[error("io: {0}")]
    Io(String),

    /// RON persistence.
    #[error("persist: {0}")]
    Persist(String),
}

/// Coarse category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// Reference resolution.
    Lookup,
    /// IO / persistence.
    Io,
}

impl HeeksCadError {
    /// Stable kebab code.
    pub fn code(&self) -> &'static str {
        match self {
            HeeksCadError::BadParameter { .. } => "heekscad.bad_parameter",
            HeeksCadError::Unknown { .. } => "heekscad.unknown",
            HeeksCadError::Io(_) => "heekscad.io",
            HeeksCadError::Persist(_) => "heekscad.persist",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            HeeksCadError::BadParameter { .. } => ErrorCategory::Input,
            HeeksCadError::Unknown { .. } => ErrorCategory::Lookup,
            HeeksCadError::Io(_) | HeeksCadError::Persist(_) => ErrorCategory::Io,
        }
    }
}
