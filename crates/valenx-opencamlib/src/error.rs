//! Typed error taxonomy for the OpenCamLib port.

use thiserror::Error;

/// Errors raised by `valenx-opencamlib`.
#[derive(Debug, Error)]
pub enum OpencamlibError {
    /// Bad parameter.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// The query point fell outside the spatial index extent — usually
    /// means the surface mesh is empty at that location.
    #[error("query out of extent at ({x}, {y})")]
    OutOfExtent {
        /// X coordinate of the query.
        x: f64,
        /// Y coordinate of the query.
        y: f64,
    },

    /// Algorithm produced a non-finite result.
    #[error("non-finite result from {algo}")]
    NonFinite {
        /// Algorithm label.
        algo: &'static str,
    },
}

/// Coarse category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// Algorithm / numerics.
    Algorithm,
}

impl OpencamlibError {
    /// Stable kebab code.
    pub fn code(&self) -> &'static str {
        match self {
            OpencamlibError::BadParameter { .. } => "opencamlib.bad_parameter",
            OpencamlibError::OutOfExtent { .. } => "opencamlib.out_of_extent",
            OpencamlibError::NonFinite { .. } => "opencamlib.non_finite",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            OpencamlibError::BadParameter { .. } => ErrorCategory::Input,
            OpencamlibError::OutOfExtent { .. } | OpencamlibError::NonFinite { .. } => {
                ErrorCategory::Algorithm
            }
        }
    }
}
