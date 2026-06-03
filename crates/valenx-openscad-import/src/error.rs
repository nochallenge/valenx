//! OpenSCAD importer error taxonomy.

use thiserror::Error;

/// Errors raised by the OpenSCAD importer.
#[derive(Debug, Error)]
pub enum OpenScadError {
    /// Lexer hit a character it doesn't recognise.
    #[error("lex error at byte {pos}: {reason}")]
    Lex {
        /// Byte offset.
        pos: usize,
        /// Human-readable reason.
        reason: String,
    },

    /// Parser hit an unexpected token.
    #[error("parse error: {reason}")]
    Parse {
        /// Human-readable reason.
        reason: String,
    },

    /// Evaluator hit a construct it can't handle.
    #[error("eval error: {reason}")]
    Eval {
        /// Human-readable reason.
        reason: String,
    },

    /// Underlying CAD kernel raised an error.
    #[error("cad: {0}")]
    Cad(String),
}

/// Coarse error category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// Input syntax.
    Input,
    /// Semantic / interpretation.
    Algorithm,
    /// CAD backend.
    Backend,
}

impl OpenScadError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            OpenScadError::Lex { .. } => "openscad.lex",
            OpenScadError::Parse { .. } => "openscad.parse",
            OpenScadError::Eval { .. } => "openscad.eval",
            OpenScadError::Cad(_) => "openscad.cad",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            OpenScadError::Lex { .. } | OpenScadError::Parse { .. } => ErrorCategory::Input,
            OpenScadError::Eval { .. } => ErrorCategory::Algorithm,
            OpenScadError::Cad(_) => ErrorCategory::Backend,
        }
    }
}
