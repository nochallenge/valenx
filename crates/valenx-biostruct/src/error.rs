//! Error taxonomy for `valenx-biostruct`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, BiostructError>`]. The variants are intentionally
//! coarse — a structure-analysis caller usually only cares about five
//! things:
//!
//! 1. Did a coordinate file (PDB / mmCIF) fail to parse
//!    ([`BiostructError::Parse`])?
//! 2. Is a constructed [`Structure`](crate::structure::Structure)
//!    itself malformed — an empty model, a residue with no atoms, a
//!    chain that does not exist ([`BiostructError::InvalidStructure`])?
//! 3. Did a selection string fail to parse or resolve
//!    ([`BiostructError::InvalidSelection`])?
//! 4. Did the caller pass nonsense arguments
//!    ([`BiostructError::Invalid`])?
//! 5. Is this a documented stub awaiting deeper work
//!    ([`BiostructError::NotYetImplemented`])?
//!
//! Use [`BiostructError::code`] for stable log / telemetry tagging and
//! [`BiostructError::category`] to bucket failures into Parse / Input /
//! Capability without matching every variant. The pattern mirrors
//! `valenx-bioseq`'s `BioseqError` and `valenx-cheminf`'s
//! `CheminfError`.

use std::fmt;

/// Errors produced by `valenx-biostruct`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BiostructError {
    /// A coordinate file failed to parse. `format` names the expected
    /// format (`"pdb"`, `"mmcif"`); `detail` is a human-readable
    /// reason surfaced verbatim in the UI.
    Parse {
        /// File format being parsed (`"pdb"`, `"mmcif"`).
        format: &'static str,
        /// 1-based line number, or `0` when not line-attributable.
        line: usize,
        /// Human-readable parse-failure reason.
        detail: String,
    },

    /// A [`Structure`](crate::structure::Structure) is structurally
    /// inconsistent: no models, a model with no chains, a residue with
    /// no atoms, a reference to a chain / residue / atom that does not
    /// exist. A property of the *hierarchy*, not of a parse or of a
    /// caller argument.
    InvalidStructure {
        /// Human-readable reason.
        reason: String,
    },

    /// A structure-selection string failed to parse, or it parsed but
    /// resolved to nothing / referenced an unknown token.
    InvalidSelection {
        /// The offending selection text (or a fragment of it).
        query: String,
        /// Human-readable reason.
        reason: String,
    },

    /// Caller passed an argument the algorithm cannot accept: an empty
    /// atom list, a non-positive cutoff, two atom sets of differing
    /// length for superposition, an out-of-range index, etc. A
    /// property of the *call*.
    Invalid {
        /// Logical parameter name (e.g. `"cutoff"`, `"atoms"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// Feature is part of this crate's public surface but not yet
    /// implemented as a real v1. The string identifies which algorithm
    /// the caller asked for.
    NotYetImplemented {
        /// Stable feature identifier (e.g. `"backbone_rebuild"`).
        feature: &'static str,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this
/// rather than on 5+ error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A coordinate file failed to parse.
    Parse,
    /// User-supplied input is wrong (bad structure, bad selection or
    /// bad argument).
    Input,
    /// Capability not available in v1 (documented stub).
    Capability,
}

impl BiostructError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"biostruct.<sub_id>"`. Codes never change
    /// across minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            BiostructError::Parse { .. } => "biostruct.parse",
            BiostructError::InvalidStructure { .. } => "biostruct.invalid_structure",
            BiostructError::InvalidSelection { .. } => "biostruct.invalid_selection",
            BiostructError::Invalid { .. } => "biostruct.invalid",
            BiostructError::NotYetImplemented { .. } => "biostruct.not_yet_implemented",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            BiostructError::Parse { .. } => "parse",
            BiostructError::InvalidStructure { .. }
            | BiostructError::InvalidSelection { .. }
            | BiostructError::Invalid { .. } => "input",
            BiostructError::NotYetImplemented { .. } => "capability",
        }
    }

    /// Typed category enum (for callers that want to `match` instead
    /// of comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            BiostructError::Parse { .. } => ErrorCategory::Parse,
            BiostructError::InvalidStructure { .. }
            | BiostructError::InvalidSelection { .. }
            | BiostructError::Invalid { .. } => ErrorCategory::Input,
            BiostructError::NotYetImplemented { .. } => ErrorCategory::Capability,
        }
    }

    /// Convenience constructor for [`BiostructError::Parse`].
    pub fn parse(format: &'static str, line: usize, detail: impl Into<String>) -> Self {
        BiostructError::Parse {
            format,
            line,
            detail: detail.into(),
        }
    }

    /// Convenience constructor for [`BiostructError::InvalidStructure`].
    pub fn invalid_structure(reason: impl Into<String>) -> Self {
        BiostructError::InvalidStructure {
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`BiostructError::InvalidSelection`].
    pub fn invalid_selection(query: impl Into<String>, reason: impl Into<String>) -> Self {
        BiostructError::InvalidSelection {
            query: query.into(),
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`BiostructError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        BiostructError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`BiostructError::NotYetImplemented`].
    pub fn not_yet(feature: &'static str) -> Self {
        BiostructError::NotYetImplemented { feature }
    }
}

impl fmt::Display for BiostructError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BiostructError::Parse {
                format,
                line,
                detail,
            } => {
                if *line == 0 {
                    write!(f, "{format} parse error: {detail}")
                } else {
                    write!(f, "{format} parse error (line {line}): {detail}")
                }
            }
            BiostructError::InvalidStructure { reason } => {
                write!(f, "invalid structure: {reason}")
            }
            BiostructError::InvalidSelection { query, reason } => {
                write!(f, "invalid selection `{query}`: {reason}")
            }
            BiostructError::Invalid { what, reason } => {
                write!(f, "invalid `{what}`: {reason}")
            }
            BiostructError::NotYetImplemented { feature } => {
                write!(
                    f,
                    "biostruct feature `{feature}` is not yet implemented (v1 scaffold)"
                )
            }
        }
    }
}

impl std::error::Error for BiostructError {}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, BiostructError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = BiostructError::parse("pdb", 12, "bad ATOM record");
        assert_eq!(err.code(), "biostruct.parse");
        assert_eq!(err.category(), "parse");
        assert_eq!(err.category_enum(), ErrorCategory::Parse);

        let err = BiostructError::invalid_structure("model has no chains");
        assert_eq!(err.code(), "biostruct.invalid_structure");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = BiostructError::invalid_selection("chain", "unknown token");
        assert_eq!(err.code(), "biostruct.invalid_selection");
        assert_eq!(err.category(), "input");

        let err = BiostructError::invalid("cutoff", "must be positive");
        assert_eq!(err.code(), "biostruct.invalid");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = BiostructError::not_yet("backbone_rebuild");
        assert_eq!(err.code(), "biostruct.not_yet_implemented");
        assert_eq!(err.category(), "capability");
        assert_eq!(err.category_enum(), ErrorCategory::Capability);
    }

    #[test]
    fn display_is_informative() {
        let msg = BiostructError::parse("mmcif", 4, "missing loop tag").to_string();
        assert!(msg.contains("mmcif"), "got: {msg}");
        assert!(msg.contains("line 4"), "got: {msg}");

        let msg = BiostructError::parse("pdb", 0, "empty file").to_string();
        assert!(!msg.contains("line"), "got: {msg}");

        let msg = BiostructError::not_yet("foo").to_string();
        assert!(msg.contains("foo"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(BiostructError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }
}
