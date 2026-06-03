//! Manipulator workbench error taxonomy.

use thiserror::Error;

/// Errors raised by manipulator ops.
#[derive(Debug, Error)]
pub enum ManipulatorError {
    /// Bad parameter.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Index out of range.
    #[error("bad index {got} (have {n})")]
    BadIndex {
        /// Asked-for id.
        got: usize,
        /// Available entries.
        n: usize,
    },

    /// Underlying tessellation failed.
    #[error("tessellation: {0}")]
    Tessellation(String),

    /// The op needs BRep topology that isn't available on a mesh-
    /// backed solid.
    #[error("requires BRep input: {0}")]
    RequiresBrep(&'static str),
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
}

impl ManipulatorError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            ManipulatorError::BadParameter { .. } => "manipulator.bad_parameter",
            ManipulatorError::BadIndex { .. } => "manipulator.bad_index",
            ManipulatorError::Tessellation(_) => "manipulator.tessellation",
            ManipulatorError::RequiresBrep(_) => "manipulator.requires_brep",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            ManipulatorError::BadParameter { .. } => ErrorCategory::Config,
            ManipulatorError::BadIndex { .. } => ErrorCategory::Input,
            ManipulatorError::Tessellation(_) | ManipulatorError::RequiresBrep(_) => {
                ErrorCategory::Algorithm
            }
        }
    }
}
