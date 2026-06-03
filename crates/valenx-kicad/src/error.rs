//! KiCad workbench error taxonomy.

use thiserror::Error;

/// Errors raised by KiCad PCB I/O.
#[derive(Debug, Error)]
pub enum KicadError {
    /// Bad parameter.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// S-expression parser error.
    #[error("kicad_pcb parse: {0}")]
    Parse(String),

    /// IO error.
    #[error("io: {0}")]
    Io(#[from] std::io::Error),

    /// Required field missing.
    #[error("missing field: {0}")]
    MissingField(&'static str),

    /// Feature deferred to a later sub-phase.
    #[error("not yet implemented: {0}")]
    NotImplemented(&'static str),
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
    /// Algorithm domain.
    Algorithm,
}

impl KicadError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            KicadError::BadParameter { .. } => "kicad.bad_parameter",
            KicadError::Parse(_) => "kicad.parse",
            KicadError::Io(_) => "kicad.io",
            KicadError::MissingField(_) => "kicad.missing_field",
            KicadError::NotImplemented(_) => "kicad.not_implemented",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            KicadError::BadParameter { .. } => ErrorCategory::Config,
            KicadError::MissingField(_) => ErrorCategory::Input,
            KicadError::Parse(_) | KicadError::Io(_) => ErrorCategory::Runtime,
            KicadError::NotImplemented(_) => ErrorCategory::Algorithm,
        }
    }
}
