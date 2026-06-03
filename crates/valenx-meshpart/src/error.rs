//! MeshPart workbench error taxonomy.

use thiserror::Error;

/// Errors raised by mesh-part ops.
#[derive(Debug, Error)]
pub enum MeshPartError {
    /// Caller passed an empty mesh / empty list.
    #[error("empty input: {0}")]
    Empty(&'static str),

    /// Bad parameter (negative tolerance, etc).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Triangulation step rejected a non-simple polygon.
    #[error("polygon triangulation failed: {0}")]
    BadPolygon(String),

    /// Underlying CAD tessellation error.
    #[error("cad: {0}")]
    Cad(String),

    /// Sewing failed — the caller gets a mesh-backed solid instead.
    #[error("BRep sewing not available; returning mesh-backed solid: {0}")]
    SewingFallback(String),
}

/// Coarse error category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// Tunable knob.
    Config,
    /// Algorithm domain error.
    Algorithm,
    /// Not implemented in v1.
    NotImplemented,
}

impl MeshPartError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            MeshPartError::Empty(_) => "meshpart.empty",
            MeshPartError::BadParameter { .. } => "meshpart.bad_parameter",
            MeshPartError::BadPolygon(_) => "meshpart.bad_polygon",
            MeshPartError::Cad(_) => "meshpart.cad",
            MeshPartError::SewingFallback(_) => "meshpart.sewing_fallback",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            MeshPartError::Empty(_) => ErrorCategory::Input,
            MeshPartError::BadParameter { .. } => ErrorCategory::Config,
            MeshPartError::BadPolygon(_) | MeshPartError::Cad(_) => ErrorCategory::Algorithm,
            MeshPartError::SewingFallback(_) => ErrorCategory::NotImplemented,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_and_cats() {
        assert_eq!(MeshPartError::Empty("mesh").code(), "meshpart.empty");
        assert_eq!(
            MeshPartError::SewingFallback("x".into()).category(),
            ErrorCategory::NotImplemented
        );
    }
}
