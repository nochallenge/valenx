//! TechDraw error taxonomy.
//!
//! Same shape as `valenx-draft`'s `DraftError` / `valenx-fillet`'s
//! `FilletError`: stable kebab-cased `code()` an LLM can branch on
//! plus a coarse [`ErrorCategory`] for UI / escalation routing.

use thiserror::Error;

/// Errors raised by drawing construction, projection, dimensioning,
/// export, or persistence.
#[derive(Debug, Error)]
pub enum TechDrawError {
    /// Tried to project from a solid with no faces / vertices / edges
    /// (BRep) or zero triangles (mesh). Nothing to draw.
    #[error("solid is empty — nothing to project")]
    EmptySolid,

    /// User-supplied view parameter is out of range or nonsense
    /// (negative scale, zero-length custom-camera axis, …).
    #[error("bad view parameter `{name}`: {reason}")]
    BadViewParameter {
        /// Name of the offending parameter (e.g. `"scale"`, `"eye"`).
        name: &'static str,
        /// Human-readable explanation.
        reason: String,
    },

    /// SVG / PDF / DXF writer hit an unexpected condition (IO,
    /// encoding, malformed input). The wrapped string is the formatter
    /// output of the underlying error so it survives across the
    /// `Send + Sync + 'static` boundary egui's error toasts assume.
    #[error("export failed: {0}")]
    ExportFailed(String),

    /// Tried to look up a view by index that doesn't exist in the
    /// drawing's `views` vector.
    #[error("view {0} not found in drawing")]
    UnknownView(usize),

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
/// Mirrors `valenx-draft::ErrorCategory`.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User supplied bad data (fix the input solid / view list).
    Input,
    /// User-tunable knob out of range (fix the scale / camera axis).
    Config,
    /// Transient or environmental failure (retry may help) — IO,
    /// parsing, export pipeline.
    Runtime,
    /// Bug in valenx-techdraw (file a report).
    Internal,
}

impl TechDrawError {
    /// Stable kebab-cased identifier; never changes across versions.
    pub fn code(&self) -> &'static str {
        match self {
            TechDrawError::EmptySolid => "techdraw.empty_solid",
            TechDrawError::BadViewParameter { .. } => "techdraw.bad_view_parameter",
            TechDrawError::ExportFailed(_) => "techdraw.export_failed",
            TechDrawError::UnknownView(_) => "techdraw.unknown_view",
            TechDrawError::Io(_) => "techdraw.io",
            TechDrawError::Ron(_) => "techdraw.ron",
        }
    }

    /// High-level classification for LLM / UI routing.
    pub fn category(&self) -> ErrorCategory {
        match self {
            TechDrawError::EmptySolid | TechDrawError::UnknownView(_) => ErrorCategory::Input,
            TechDrawError::BadViewParameter { .. } => ErrorCategory::Config,
            TechDrawError::Io(_) | TechDrawError::Ron(_) | TechDrawError::ExportFailed(_) => {
                ErrorCategory::Runtime
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_variant_has_stable_code_and_category() {
        let cases: Vec<(TechDrawError, &'static str, ErrorCategory)> = vec![
            (
                TechDrawError::EmptySolid,
                "techdraw.empty_solid",
                ErrorCategory::Input,
            ),
            (
                TechDrawError::BadViewParameter {
                    name: "scale",
                    reason: "negative".into(),
                },
                "techdraw.bad_view_parameter",
                ErrorCategory::Config,
            ),
            (
                TechDrawError::ExportFailed("disk full".into()),
                "techdraw.export_failed",
                ErrorCategory::Runtime,
            ),
            (
                TechDrawError::UnknownView(3),
                "techdraw.unknown_view",
                ErrorCategory::Input,
            ),
            (
                TechDrawError::Ron("bad token".into()),
                "techdraw.ron",
                ErrorCategory::Runtime,
            ),
        ];
        for (err, expected_code, expected_cat) in cases {
            assert_eq!(err.code(), expected_code, "wrong code for {err:?}");
            assert_eq!(err.category(), expected_cat, "wrong category for {err:?}");
        }
    }

    #[test]
    fn display_contains_useful_context() {
        let e = TechDrawError::BadViewParameter {
            name: "scale",
            reason: "must be > 0".into(),
        };
        let msg = e.to_string();
        assert!(msg.contains("scale"));
        assert!(msg.contains("must be > 0"));
    }
}
