//! Error taxonomy for `valenx-cheminf`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, CheminfError>`]. The variants are intentionally coarse —
//! a cheminformatics caller usually only cares about four things:
//!
//! 1. Did a notation string / file fail to parse
//!    ([`CheminfError::Parse`])?
//! 2. Is a constructed [`Molecule`](crate::molecule::Molecule) itself
//!    malformed — a dangling bond index, an unclosed ring, an
//!    impossible valence ([`CheminfError::InvalidMolecule`])?
//! 3. Did the caller pass nonsense arguments
//!    ([`CheminfError::Invalid`])?
//! 4. Is this a documented stub awaiting deeper work
//!    ([`CheminfError::NotYetImplemented`])?
//!
//! Use [`CheminfError::code`] for stable log/telemetry tagging and
//! [`CheminfError::category`] to bucket failures into Parse / Input /
//! Capability without matching every variant. The pattern mirrors
//! `valenx-bioseq`'s `BioseqError`.

use std::fmt;

/// Errors produced by `valenx-cheminf`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CheminfError {
    /// A notation string or file failed to parse. `format` names the
    /// expected notation (`"smiles"`, `"smarts"`, `"mol"`, `"sdf"`);
    /// `detail` is a human-readable reason surfaced verbatim in the UI.
    Parse {
        /// Notation being parsed (`"smiles"`, `"smarts"`, `"mol"`, …).
        format: &'static str,
        /// Human-readable parse-failure reason.
        detail: String,
    },

    /// A [`Molecule`](crate::molecule::Molecule) is structurally
    /// inconsistent: a bond referencing an out-of-range atom, an
    /// unclosed ring-bond pair, an atom whose declared bonding exceeds
    /// any plausible valence, etc. A property of the *graph*, not of a
    /// parse or of a caller argument.
    InvalidMolecule {
        /// Human-readable reason.
        reason: String,
    },

    /// Caller passed an argument the algorithm cannot accept: an empty
    /// input, a non-positive count, mismatched fingerprint lengths, an
    /// out-of-range atom index, etc. A property of the *call*.
    Invalid {
        /// Logical parameter name (e.g. `"radius"`, `"n_bits"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// Feature is part of this crate's public surface but not yet
    /// implemented as a real v1. The string identifies which algorithm
    /// the caller asked for.
    NotYetImplemented {
        /// Stable feature identifier (e.g. `"ringbond_stereo"`).
        feature: &'static str,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this
/// rather than on 4+ error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A notation string / file failed to parse.
    Parse,
    /// User-supplied input is wrong (bad molecule or bad argument).
    Input,
    /// Capability not available in v1 (documented stub).
    Capability,
}

impl CheminfError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"cheminf.<sub_id>"`. Codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            CheminfError::Parse { .. } => "cheminf.parse",
            CheminfError::InvalidMolecule { .. } => "cheminf.invalid_molecule",
            CheminfError::Invalid { .. } => "cheminf.invalid",
            CheminfError::NotYetImplemented { .. } => "cheminf.not_yet_implemented",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            CheminfError::Parse { .. } => "parse",
            CheminfError::InvalidMolecule { .. } | CheminfError::Invalid { .. } => "input",
            CheminfError::NotYetImplemented { .. } => "capability",
        }
    }

    /// Typed category enum (for callers that want to `match` instead of
    /// comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            CheminfError::Parse { .. } => ErrorCategory::Parse,
            CheminfError::InvalidMolecule { .. } | CheminfError::Invalid { .. } => {
                ErrorCategory::Input
            }
            CheminfError::NotYetImplemented { .. } => ErrorCategory::Capability,
        }
    }

    /// Convenience constructor for [`CheminfError::Parse`].
    pub fn parse(format: &'static str, detail: impl Into<String>) -> Self {
        CheminfError::Parse {
            format,
            detail: detail.into(),
        }
    }

    /// Convenience constructor for [`CheminfError::InvalidMolecule`].
    pub fn invalid_molecule(reason: impl Into<String>) -> Self {
        CheminfError::InvalidMolecule {
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`CheminfError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        CheminfError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`CheminfError::NotYetImplemented`].
    pub fn not_yet(feature: &'static str) -> Self {
        CheminfError::NotYetImplemented { feature }
    }
}

impl fmt::Display for CheminfError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CheminfError::Parse { format, detail } => {
                write!(f, "{format} parse error: {detail}")
            }
            CheminfError::InvalidMolecule { reason } => {
                write!(f, "invalid molecule: {reason}")
            }
            CheminfError::Invalid { what, reason } => {
                write!(f, "invalid `{what}`: {reason}")
            }
            CheminfError::NotYetImplemented { feature } => {
                write!(
                    f,
                    "cheminf feature `{feature}` is not yet implemented (v1 scaffold)"
                )
            }
        }
    }
}

impl std::error::Error for CheminfError {}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, CheminfError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = CheminfError::parse("smiles", "unclosed ring");
        assert_eq!(err.code(), "cheminf.parse");
        assert_eq!(err.category(), "parse");
        assert_eq!(err.category_enum(), ErrorCategory::Parse);

        let err = CheminfError::invalid_molecule("dangling bond");
        assert_eq!(err.code(), "cheminf.invalid_molecule");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = CheminfError::invalid("radius", "must be positive");
        assert_eq!(err.code(), "cheminf.invalid");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = CheminfError::not_yet("ringbond_stereo");
        assert_eq!(err.code(), "cheminf.not_yet_implemented");
        assert_eq!(err.category(), "capability");
        assert_eq!(err.category_enum(), ErrorCategory::Capability);
    }

    #[test]
    fn display_is_informative() {
        let msg = CheminfError::parse("smarts", "bad atom").to_string();
        assert!(msg.contains("smarts"), "got: {msg}");
        assert!(msg.contains("bad atom"), "got: {msg}");

        let msg = CheminfError::not_yet("foo").to_string();
        assert!(msg.contains("foo"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(CheminfError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }
}
