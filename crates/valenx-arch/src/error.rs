//! Arch/BIM workbench error taxonomy.
//!
//! Same shape as `valenx-draft`'s `DraftError`: a stable
//! kebab-cased `code()` an LLM can branch on plus a coarse
//! [`ErrorCategory`] for UI / escalation routing.

use thiserror::Error;

/// Errors raised by arch entity construction, document mutation, IFC
/// or BCF export, or persistence.
#[derive(Debug, Error)]
pub enum ArchError {
    /// User-supplied dimension was out of range or nonsense (negative
    /// thickness, zero height, more end-points than start-points, …).
    #[error("bad dimension `{name}`: {reason}")]
    BadDimension {
        /// Name of the offending parameter (e.g. `"height"`,
        /// `"thickness"`).
        name: &'static str,
        /// Human-readable explanation.
        reason: String,
    },

    /// Tried to look up an entity that doesn't exist in this document.
    #[error("entity {0} not found in arch document")]
    UnknownEntity(usize),

    /// IFC writer failed to emit a well-formed file (I/O,
    /// schema-violation, or unsupported entity).
    #[error("ifc write failed: {0}")]
    IfcWriteFailed(String),

    /// BCF writer failed (I/O or invalid issue payload).
    #[error("bcf write failed: {0}")]
    BcfWriteFailed(String),

    /// IO error wrapping `std::io::Error`.
    #[error("io: {0}")]
    Io(#[from] std::io::Error),

    /// RON parse / serialise error (string-wrapped because `ron` has
    /// distinct error types for ser / de).
    #[error("ron: {0}")]
    Ron(String),

    /// A tessellation / CAD primitive call failed (negative size,
    /// degenerate profile, …). String-wrapped to break the dep
    /// cycle that an explicit `#[from] valenx_cad::CadError` would
    /// imply.
    #[error("cad: {0}")]
    Cad(String),
}

/// Coarse category an LLM / UI can branch on to decide who to
/// escalate the error to. Mirrors
/// `valenx_draft::ErrorCategory`.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User supplied bad data (fix the click coords / entity list).
    Input,
    /// User-tunable knob out of range (fix the dimension / count).
    Config,
    /// Transient or environmental failure (retry may help) — IO,
    /// parsing, CAD kernel state.
    Runtime,
    /// Bug in `valenx-arch` (file a report).
    Internal,
}

impl ArchError {
    /// Stable kebab-cased identifier; never changes across versions.
    pub fn code(&self) -> &'static str {
        match self {
            ArchError::BadDimension { .. } => "arch.bad_dimension",
            ArchError::UnknownEntity(_) => "arch.unknown_entity",
            ArchError::IfcWriteFailed(_) => "arch.ifc_write_failed",
            ArchError::BcfWriteFailed(_) => "arch.bcf_write_failed",
            ArchError::Io(_) => "arch.io",
            ArchError::Ron(_) => "arch.ron",
            ArchError::Cad(_) => "arch.cad",
        }
    }

    /// High-level classification for LLM / UI routing.
    pub fn category(&self) -> ErrorCategory {
        match self {
            ArchError::UnknownEntity(_) => ErrorCategory::Input,
            ArchError::BadDimension { .. } => ErrorCategory::Config,
            ArchError::IfcWriteFailed(_)
            | ArchError::BcfWriteFailed(_)
            | ArchError::Io(_)
            | ArchError::Ron(_)
            | ArchError::Cad(_) => ErrorCategory::Runtime,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_variant_has_stable_code_and_category() {
        let cases: Vec<(ArchError, &'static str, ErrorCategory)> = vec![
            (
                ArchError::BadDimension {
                    name: "height",
                    reason: "negative".into(),
                },
                "arch.bad_dimension",
                ErrorCategory::Config,
            ),
            (
                ArchError::UnknownEntity(7),
                "arch.unknown_entity",
                ErrorCategory::Input,
            ),
            (
                ArchError::IfcWriteFailed("schema mismatch".into()),
                "arch.ifc_write_failed",
                ErrorCategory::Runtime,
            ),
            (
                ArchError::BcfWriteFailed("empty title".into()),
                "arch.bcf_write_failed",
                ErrorCategory::Runtime,
            ),
            (
                ArchError::Ron("bad token".into()),
                "arch.ron",
                ErrorCategory::Runtime,
            ),
            (
                ArchError::Cad("degenerate".into()),
                "arch.cad",
                ErrorCategory::Runtime,
            ),
        ];
        for (err, expected_code, expected_cat) in cases {
            assert_eq!(err.code(), expected_code, "wrong code for {err:?}");
            assert_eq!(err.category(), expected_cat, "wrong category for {err:?}");
        }
    }

    #[test]
    fn display_includes_message() {
        let e = ArchError::BadDimension {
            name: "thickness",
            reason: "negative".into(),
        };
        let s = e.to_string();
        assert!(s.contains("thickness"));
        assert!(s.contains("negative"));
    }
}
