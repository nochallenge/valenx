//! Symbols workbench error taxonomy.

use thiserror::Error;

/// Errors raised by schematic / symbol ops.
#[derive(Debug, Error)]
pub enum SymbolError {
    /// Bad parameter (negative rotation, etc).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Wire polyline has fewer than 2 vertices.
    #[error("degenerate wire: {0}")]
    DegenerateWire(String),

    /// Placement references a wire/symbol index out of range.
    #[error("bad index {got} (have {n})")]
    BadIndex {
        /// The id we asked for.
        got: usize,
        /// Available entries.
        n: usize,
    },

    /// IO error (placeholder for future load/save hooks).
    #[error("io: {0}")]
    Io(#[from] std::io::Error),

    /// RON serialization error.
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
    /// Transient / IO / parse.
    Runtime,
}

impl SymbolError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            SymbolError::BadParameter { .. } => "symbols.bad_parameter",
            SymbolError::DegenerateWire(_) => "symbols.degenerate_wire",
            SymbolError::BadIndex { .. } => "symbols.bad_index",
            SymbolError::Io(_) => "symbols.io",
            SymbolError::Ron(_) => "symbols.ron",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            SymbolError::BadParameter { .. } => ErrorCategory::Config,
            SymbolError::DegenerateWire(_) | SymbolError::BadIndex { .. } => {
                ErrorCategory::Algorithm
            }
            SymbolError::Io(_) | SymbolError::Ron(_) => ErrorCategory::Runtime,
        }
    }
}
