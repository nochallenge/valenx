//! Draft workbench error taxonomy.
//!
//! Same shape as `valenx-fillet`'s `FilletError`: stable `code()`
//! strings an LLM can branch on plus a coarse [`ErrorCategory`] for
//! UI / escalation routing.

use thiserror::Error;

/// Errors raised by draft document construction, mutation, or
/// persistence.
#[derive(Debug, Error)]
pub enum DraftError {
    /// Tried to look up an entity that doesn't exist in this document.
    #[error("entity {0} not found in draft document")]
    UnknownEntity(usize),

    /// User-supplied parameter was out of range or nonsense (negative
    /// radius, zero-side polygon, empty text, …).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Name of the offending parameter (e.g. `"radius"`, `"sides"`).
        name: &'static str,
        /// Human-readable explanation.
        reason: String,
    },

    /// Tried to perform an operation that requires at least one entity
    /// on an empty document (e.g. compute a bounding box).
    #[error("draft document is empty")]
    EmptyDocument,

    /// IO error wrapping std::io.
    #[error("io: {0}")]
    Io(#[from] std::io::Error),

    /// RON parse / serialise error (string-wrapped because `ron` has
    /// distinct error types for ser / de).
    #[error("ron: {0}")]
    Ron(String),
}

/// Coarse category an LLM / UI can branch on to decide who to
/// escalate the error to.
///
/// Mirrors `valenx-fillet::ErrorCategory`.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User supplied bad data (fix the click coords / entity list).
    Input,
    /// User-tunable knob out of range (fix the radius / sides / size).
    Config,
    /// Transient or environmental failure (retry may help) — IO,
    /// parsing.
    Runtime,
    /// Bug in valenx-draft (file a report).
    Internal,
}

impl DraftError {
    /// Stable kebab-cased identifier; never changes across versions.
    pub fn code(&self) -> &'static str {
        match self {
            DraftError::UnknownEntity(_) => "draft.unknown_entity",
            DraftError::BadParameter { .. } => "draft.bad_parameter",
            DraftError::EmptyDocument => "draft.empty_document",
            DraftError::Io(_) => "draft.io",
            DraftError::Ron(_) => "draft.ron",
        }
    }

    /// High-level classification for LLM / UI routing.
    pub fn category(&self) -> ErrorCategory {
        match self {
            DraftError::UnknownEntity(_) | DraftError::EmptyDocument => ErrorCategory::Input,
            DraftError::BadParameter { .. } => ErrorCategory::Config,
            DraftError::Io(_) | DraftError::Ron(_) => ErrorCategory::Runtime,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_variant_has_stable_code_and_category() {
        let cases: Vec<(DraftError, &'static str, ErrorCategory)> = vec![
            (
                DraftError::UnknownEntity(7),
                "draft.unknown_entity",
                ErrorCategory::Input,
            ),
            (
                DraftError::BadParameter {
                    name: "radius",
                    reason: "negative".into(),
                },
                "draft.bad_parameter",
                ErrorCategory::Config,
            ),
            (
                DraftError::EmptyDocument,
                "draft.empty_document",
                ErrorCategory::Input,
            ),
            (
                DraftError::Ron("bad token".into()),
                "draft.ron",
                ErrorCategory::Runtime,
            ),
        ];
        for (err, expected_code, expected_cat) in cases {
            assert_eq!(err.code(), expected_code, "wrong code for {err:?}");
            assert_eq!(err.category(), expected_cat, "wrong category for {err:?}");
        }
    }
}
