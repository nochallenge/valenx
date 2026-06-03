//! Typed errors for `valenx-salome-bridge`.

use thiserror::Error;

/// Errors raised by the Salome-style platform bridge.
#[derive(Debug, Error)]
pub enum SalomeError {
    /// Bad parameter.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Node id not found in the study DAG.
    #[error("unknown node id {0}")]
    UnknownNode(u32),

    /// Dependency cycle detected during rebuild traversal.
    #[error("dependency cycle including node {0}")]
    Cycle(u32),

    /// Bridge module returned an error — wraps the original message.
    #[error("module `{module}`: {message}")]
    Module {
        /// Which module — `"geom"` / `"mesh"` / `"analysis"`.
        module: &'static str,
        /// Underlying message.
        message: String,
    },
}

/// Coarse category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// Study DAG / lookup.
    Study,
    /// Underlying bridged module.
    Module,
}

impl SalomeError {
    /// Stable kebab code.
    pub fn code(&self) -> &'static str {
        match self {
            SalomeError::BadParameter { .. } => "salome.bad_parameter",
            SalomeError::UnknownNode(_) => "salome.unknown_node",
            SalomeError::Cycle(_) => "salome.cycle",
            SalomeError::Module { .. } => "salome.module",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            SalomeError::BadParameter { .. } => ErrorCategory::Input,
            SalomeError::UnknownNode(_) | SalomeError::Cycle(_) => ErrorCategory::Study,
            SalomeError::Module { .. } => ErrorCategory::Module,
        }
    }
}
