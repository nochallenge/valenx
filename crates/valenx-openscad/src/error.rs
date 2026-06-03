//! Typed error taxonomy for the OpenSCAD CSG engine.

use thiserror::Error;

use valenx_openscad_import::OpenScadError;

/// Errors raised by the engine.
#[derive(Debug, Error)]
pub enum OpenScadCsgError {
    /// Bad parameter.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Lex / parse / eval failure forwarded from the importer.
    #[error("openscad: {0}")]
    Inner(#[from] OpenScadError),

    /// CAD-kernel failure (e.g. boolean tolerance).
    #[error("cad: {op}: {reason}")]
    Cad {
        /// Op label.
        op: &'static str,
        /// Reason.
        reason: String,
    },

    /// Operation not implemented yet (e.g. `minkowski`).
    #[error("not implemented: {op}")]
    NotImplemented {
        /// Op label.
        op: &'static str,
    },
}

/// Coarse error category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input (parameters, source text).
    Input,
    /// Algorithm / kernel.
    Algorithm,
    /// Not yet implemented.
    NotImplemented,
}

impl OpenScadCsgError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            OpenScadCsgError::BadParameter { .. } => "openscad.bad_parameter",
            OpenScadCsgError::Inner(_) => "openscad.inner",
            OpenScadCsgError::Cad { .. } => "openscad.cad",
            OpenScadCsgError::NotImplemented { .. } => "openscad.not_implemented",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            OpenScadCsgError::BadParameter { .. } | OpenScadCsgError::Inner(_) => {
                ErrorCategory::Input
            }
            OpenScadCsgError::Cad { .. } => ErrorCategory::Algorithm,
            OpenScadCsgError::NotImplemented { .. } => ErrorCategory::NotImplemented,
        }
    }
}
