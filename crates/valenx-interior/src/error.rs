//! Typed errors for `valenx-interior`.

use thiserror::Error;

/// Errors raised by the interior-design workbench.
#[derive(Debug, Error)]
pub enum InteriorError {
    /// Bad parameter.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Room id not found.
    #[error("room id `{0}` not found")]
    UnknownRoom(String),

    /// Furniture id not found.
    #[error("furniture id `{0}` not found")]
    UnknownFurniture(String),

    /// Geometry — empty polygon, degenerate, etc.
    #[error("geometry: {0}")]
    Geometry(String),
}

/// Coarse category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input / parameter.
    Input,
    /// Reference resolution.
    Lookup,
    /// Geometric problem.
    Geometry,
}

impl InteriorError {
    /// Stable kebab code.
    pub fn code(&self) -> &'static str {
        match self {
            InteriorError::BadParameter { .. } => "interior.bad_parameter",
            InteriorError::UnknownRoom(_) => "interior.unknown_room",
            InteriorError::UnknownFurniture(_) => "interior.unknown_furniture",
            InteriorError::Geometry(_) => "interior.geometry",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            InteriorError::BadParameter { .. } => ErrorCategory::Input,
            InteriorError::UnknownRoom(_) | InteriorError::UnknownFurniture(_) => {
                ErrorCategory::Lookup
            }
            InteriorError::Geometry(_) => ErrorCategory::Geometry,
        }
    }
}
