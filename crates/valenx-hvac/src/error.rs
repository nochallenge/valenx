//! HVAC error taxonomy.

use thiserror::Error;

/// Errors raised by the HVAC crate.
#[derive(Debug, Error)]
pub enum HvacError {
    /// Bad parameter.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// CAD kernel wrapped.
    #[error("cad: {0}")]
    Cad(String),
}

/// Coarse error category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// Backend.
    Backend,
}

impl HvacError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            HvacError::BadParameter { .. } => "hvac.bad_parameter",
            HvacError::Cad(_) => "hvac.cad",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            HvacError::BadParameter { .. } => ErrorCategory::Input,
            HvacError::Cad(_) => ErrorCategory::Backend,
        }
    }
}
