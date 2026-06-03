//! Reverse-engineering workbench error taxonomy.

use thiserror::Error;

/// Errors raised by point-cloud loading, normal estimation, or
/// triangulation.
#[derive(Debug, Error)]
pub enum ReverseError {
    /// PLY parser couldn't make sense of the file.
    #[error("PLY parse error: {0}")]
    PlyParse(String),

    /// Caller passed an empty cloud / too-small k.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Name of the offending parameter (e.g. `"k"`).
        name: &'static str,
        /// Human-readable explanation.
        reason: String,
    },

    /// Triangulation produced no triangles (cloud too sparse, k too
    /// small).
    #[error("triangulation produced no triangles ({reason})")]
    EmptyTriangulation {
        /// Why the result was empty.
        reason: String,
    },

    /// Wrapped error from the Phase 23 mesh-to-brep stage.
    #[error("reconstruct: {0}")]
    Reconstruct(String),

    /// IO error.
    #[error("io: {0}")]
    Io(#[from] std::io::Error),

    /// RON parse / serialise error.
    #[error("ron: {0}")]
    Ron(String),
}

/// Coarse category an LLM / UI can branch on.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input — bad file, empty cloud.
    Input,
    /// Tunable knob (k, tolerance) out of range.
    Config,
    /// Transient / environmental — IO / parse.
    Runtime,
    /// Bug in reconstruct (file a report).
    Internal,
}

impl ReverseError {
    /// Stable kebab-cased identifier.
    pub fn code(&self) -> &'static str {
        match self {
            ReverseError::PlyParse(_) => "reverse.ply_parse",
            ReverseError::BadParameter { .. } => "reverse.bad_parameter",
            ReverseError::EmptyTriangulation { .. } => "reverse.empty_triangulation",
            ReverseError::Reconstruct(_) => "reverse.reconstruct",
            ReverseError::Io(_) => "reverse.io",
            ReverseError::Ron(_) => "reverse.ron",
        }
    }

    /// High-level classification.
    pub fn category(&self) -> ErrorCategory {
        match self {
            ReverseError::PlyParse(_) | ReverseError::EmptyTriangulation { .. } => {
                ErrorCategory::Input
            }
            ReverseError::BadParameter { .. } => ErrorCategory::Config,
            ReverseError::Reconstruct(_) => ErrorCategory::Internal,
            ReverseError::Io(_) | ReverseError::Ron(_) => ErrorCategory::Runtime,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_variant_has_stable_code_and_category() {
        let cases: Vec<(ReverseError, &'static str, ErrorCategory)> = vec![
            (
                ReverseError::PlyParse("missing header".into()),
                "reverse.ply_parse",
                ErrorCategory::Input,
            ),
            (
                ReverseError::BadParameter {
                    name: "k",
                    reason: "must be > 0".into(),
                },
                "reverse.bad_parameter",
                ErrorCategory::Config,
            ),
            (
                ReverseError::EmptyTriangulation {
                    reason: "too sparse".into(),
                },
                "reverse.empty_triangulation",
                ErrorCategory::Input,
            ),
            (
                ReverseError::Reconstruct("upstream".into()),
                "reverse.reconstruct",
                ErrorCategory::Internal,
            ),
            (
                ReverseError::Ron("bad token".into()),
                "reverse.ron",
                ErrorCategory::Runtime,
            ),
        ];
        for (err, code, cat) in cases {
            assert_eq!(err.code(), code);
            assert_eq!(err.category(), cat);
        }
    }
}
