//! Geomatics workbench error taxonomy.

use thiserror::Error;

/// Errors raised by geomatics ops.
#[derive(Debug, Error)]
pub enum GeomaticsError {
    /// Bad parameter (negative grid size, lat out of range, etc).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// XYZ ASCII parser detected an irregular grid.
    #[error("irregular grid: {0}")]
    IrregularGrid(String),

    /// XYZ ASCII parser hit a malformed row.
    #[error("parse error at line {line}: {msg}")]
    Parse {
        /// 1-based line number.
        line: usize,
        /// Detail.
        msg: String,
    },

    /// IO error.
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
}

/// Coarse error category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// Tunable knob.
    Config,
    /// Parse error.
    Parse,
    /// Transient / IO.
    Runtime,
}

impl GeomaticsError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            GeomaticsError::BadParameter { .. } => "geomatics.bad_parameter",
            GeomaticsError::IrregularGrid(_) => "geomatics.irregular_grid",
            GeomaticsError::Parse { .. } => "geomatics.parse",
            GeomaticsError::Io(_) => "geomatics.io",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            GeomaticsError::BadParameter { .. } => ErrorCategory::Config,
            GeomaticsError::IrregularGrid(_) => ErrorCategory::Input,
            GeomaticsError::Parse { .. } => ErrorCategory::Parse,
            GeomaticsError::Io(_) => ErrorCategory::Runtime,
        }
    }
}
