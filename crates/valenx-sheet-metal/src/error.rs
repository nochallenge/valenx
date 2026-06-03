//! Sheet metal workbench error taxonomy.

use thiserror::Error;

/// Errors raised by sheet metal ops.
#[derive(Debug, Error)]
pub enum SheetMetalError {
    /// Bad parameter (negative thickness, etc).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Outline polygon is non-simple / has fewer than 3 vertices.
    #[error("bad polygon: {0}")]
    BadPolygon(String),

    /// Edge id out of range for the sheet outline.
    #[error("bad edge index {got} (outline has {n} edges)")]
    BadEdge {
        /// The id we asked for.
        got: usize,
        /// Number of outline edges available.
        n: usize,
    },

    /// IO error.
    #[error("io: {0}")]
    Io(#[from] std::io::Error),

    /// RON error.
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

impl SheetMetalError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            SheetMetalError::BadParameter { .. } => "sheet_metal.bad_parameter",
            SheetMetalError::BadPolygon(_) => "sheet_metal.bad_polygon",
            SheetMetalError::BadEdge { .. } => "sheet_metal.bad_edge",
            SheetMetalError::Io(_) => "sheet_metal.io",
            SheetMetalError::Ron(_) => "sheet_metal.ron",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            SheetMetalError::BadParameter { .. } => ErrorCategory::Config,
            SheetMetalError::BadPolygon(_) | SheetMetalError::BadEdge { .. } => {
                ErrorCategory::Algorithm
            }
            SheetMetalError::Io(_) | SheetMetalError::Ron(_) => ErrorCategory::Runtime,
        }
    }
}
