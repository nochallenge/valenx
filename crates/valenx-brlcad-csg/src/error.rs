//! Typed errors for `valenx-brlcad-csg`.

use thiserror::Error;

/// Errors raised by the BRL-CAD-style CSG evaluator + parser.
#[derive(Debug, Error)]
pub enum BrlCadError {
    /// Bad parameter on a primitive (e.g. negative radius).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// MGED parser error.
    #[error("parse error at line {line}: {message}")]
    Parse {
        /// 1-based line number.
        line: usize,
        /// Diagnostic message.
        message: String,
    },

    /// Empty tree — `evaluate` called on `()`.
    #[error("empty CSG tree")]
    EmptyTree,
}

/// Coarse category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input / parameter.
    Input,
    /// Parser / format.
    Format,
}

impl BrlCadError {
    /// Stable kebab code.
    pub fn code(&self) -> &'static str {
        match self {
            BrlCadError::BadParameter { .. } => "brlcad.bad_parameter",
            BrlCadError::Parse { .. } => "brlcad.parse",
            BrlCadError::EmptyTree => "brlcad.empty_tree",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            BrlCadError::BadParameter { .. } => ErrorCategory::Input,
            BrlCadError::Parse { .. } | BrlCadError::EmptyTree => ErrorCategory::Format,
        }
    }
}
