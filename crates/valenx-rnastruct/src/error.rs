//! Error taxonomy for `valenx-rnastruct`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, RnaStructError>`]. The variants are intentionally
//! coarse — an RNA-folding caller usually only cares about four things:
//!
//! 1. Did a structure / file string fail to parse
//!    ([`RnaStructError::Parse`])?
//! 2. Did the caller hand over an illegal sequence
//!    ([`RnaStructError::Sequence`]) or an inconsistent structure
//!    ([`RnaStructError::Structure`])?
//! 3. Did the caller pass nonsense arguments
//!    ([`RnaStructError::Invalid`])?
//! 4. Is this a documented stub awaiting deeper work
//!    ([`RnaStructError::NotYetImplemented`])?
//!
//! Use [`RnaStructError::code`] for stable log/telemetry tagging and
//! [`RnaStructError::category`] to bucket failures into Parse / Input /
//! Capability without matching every variant. The pattern mirrors
//! `valenx-bioseq`'s `BioseqError`.

use std::fmt;

/// Errors produced by `valenx-rnastruct`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RnaStructError {
    /// A structure string or a structure file (ct / bpseq /
    /// dot-bracket) failed to parse. `format` names the expected
    /// format; `detail` is a human-readable reason surfaced verbatim.
    Parse {
        /// Format being parsed (e.g. `"dot-bracket"`, `"ct"`, `"bpseq"`).
        format: &'static str,
        /// Human-readable parse-failure reason.
        detail: String,
    },

    /// The caller supplied a sequence the folder cannot accept — an
    /// empty sequence, a non-RNA alphabet, or a residue outside
    /// `A C G U`. A property of the *sequence*.
    Sequence {
        /// Human-readable reason.
        reason: String,
    },

    /// The caller supplied a base-pair set / dot-bracket structure that
    /// is internally inconsistent — an index out of range, a position
    /// paired twice, or unbalanced brackets.
    Structure {
        /// Human-readable reason.
        reason: String,
    },

    /// Caller passed an argument the algorithm cannot accept: a
    /// non-positive length, mismatched sequence / reactivity lengths,
    /// an out-of-range temperature, etc. A property of the *call*.
    Invalid {
        /// Logical parameter name (e.g. `"delta"`, `"temperature"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// Feature is part of this crate's public surface but not yet
    /// implemented as a real v1. The string identifies which algorithm
    /// the caller asked for.
    NotYetImplemented {
        /// Stable feature identifier.
        feature: &'static str,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this
/// rather than on 5 error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A structure file / string failed to parse.
    Parse,
    /// User-supplied input is wrong (bad sequence, bad structure, or
    /// bad argument).
    Input,
    /// Capability not available in v1 (documented stub).
    Capability,
}

impl RnaStructError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"rnastruct.<sub_id>"`. Codes never change
    /// across minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            RnaStructError::Parse { .. } => "rnastruct.parse",
            RnaStructError::Sequence { .. } => "rnastruct.sequence",
            RnaStructError::Structure { .. } => "rnastruct.structure",
            RnaStructError::Invalid { .. } => "rnastruct.invalid",
            RnaStructError::NotYetImplemented { .. } => "rnastruct.not_yet_implemented",
        }
    }

    /// Coarse category — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            RnaStructError::Parse { .. } => "parse",
            RnaStructError::Sequence { .. }
            | RnaStructError::Structure { .. }
            | RnaStructError::Invalid { .. } => "input",
            RnaStructError::NotYetImplemented { .. } => "capability",
        }
    }

    /// Typed category enum (for callers that want to `match` instead of
    /// comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            RnaStructError::Parse { .. } => ErrorCategory::Parse,
            RnaStructError::Sequence { .. }
            | RnaStructError::Structure { .. }
            | RnaStructError::Invalid { .. } => ErrorCategory::Input,
            RnaStructError::NotYetImplemented { .. } => ErrorCategory::Capability,
        }
    }

    /// Convenience constructor for [`RnaStructError::Parse`].
    pub fn parse(format: &'static str, detail: impl Into<String>) -> Self {
        RnaStructError::Parse {
            format,
            detail: detail.into(),
        }
    }

    /// Convenience constructor for [`RnaStructError::Sequence`].
    pub fn sequence(reason: impl Into<String>) -> Self {
        RnaStructError::Sequence {
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`RnaStructError::Structure`].
    pub fn structure(reason: impl Into<String>) -> Self {
        RnaStructError::Structure {
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`RnaStructError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        RnaStructError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`RnaStructError::NotYetImplemented`].
    pub fn not_yet(feature: &'static str) -> Self {
        RnaStructError::NotYetImplemented { feature }
    }
}

impl fmt::Display for RnaStructError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RnaStructError::Parse { format, detail } => {
                write!(f, "{format} parse error: {detail}")
            }
            RnaStructError::Sequence { reason } => {
                write!(f, "invalid RNA sequence: {reason}")
            }
            RnaStructError::Structure { reason } => {
                write!(f, "invalid structure: {reason}")
            }
            RnaStructError::Invalid { what, reason } => {
                write!(f, "invalid `{what}`: {reason}")
            }
            RnaStructError::NotYetImplemented { feature } => {
                write!(
                    f,
                    "rnastruct feature `{feature}` is not yet implemented (v1 scaffold)"
                )
            }
        }
    }
}

impl std::error::Error for RnaStructError {}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, RnaStructError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = RnaStructError::parse("dot-bracket", "unbalanced");
        assert_eq!(err.code(), "rnastruct.parse");
        assert_eq!(err.category(), "parse");
        assert_eq!(err.category_enum(), ErrorCategory::Parse);

        let err = RnaStructError::sequence("contains T");
        assert_eq!(err.code(), "rnastruct.sequence");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = RnaStructError::structure("index out of range");
        assert_eq!(err.code(), "rnastruct.structure");
        assert_eq!(err.category(), "input");

        let err = RnaStructError::invalid("delta", "must be >= 0");
        assert_eq!(err.code(), "rnastruct.invalid");
        assert_eq!(err.category(), "input");

        let err = RnaStructError::not_yet("tertiary");
        assert_eq!(err.code(), "rnastruct.not_yet_implemented");
        assert_eq!(err.category(), "capability");
        assert_eq!(err.category_enum(), ErrorCategory::Capability);
    }

    #[test]
    fn display_is_informative() {
        let msg = RnaStructError::sequence("contains a gap").to_string();
        assert!(msg.contains("gap"), "got: {msg}");

        let msg = RnaStructError::not_yet("foo").to_string();
        assert!(msg.contains("foo"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(RnaStructError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }
}
