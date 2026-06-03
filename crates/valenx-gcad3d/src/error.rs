//! Typed errors for `valenx-gcad3d`.

use thiserror::Error;

/// Errors raised by the gCAD3D primitive constructors.
#[derive(Debug, Error)]
pub enum Gcad3dError {
    /// Bad parameter.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Degenerate input (e.g. colinear points for a 3-point arc).
    #[error("degenerate input: {0}")]
    Degenerate(String),

    /// Unsupported character in text extrude.
    #[error("unsupported text character `{0}`")]
    UnsupportedChar(char),
}

/// Coarse category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// Geometric degeneracy.
    Geometry,
}

impl Gcad3dError {
    /// Stable kebab code.
    pub fn code(&self) -> &'static str {
        match self {
            Gcad3dError::BadParameter { .. } => "gcad3d.bad_parameter",
            Gcad3dError::Degenerate(_) => "gcad3d.degenerate",
            Gcad3dError::UnsupportedChar(_) => "gcad3d.unsupported_char",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            Gcad3dError::BadParameter { .. } | Gcad3dError::UnsupportedChar(_) => {
                ErrorCategory::Input
            }
            Gcad3dError::Degenerate(_) => ErrorCategory::Geometry,
        }
    }
}
