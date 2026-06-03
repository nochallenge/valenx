//! Error taxonomy for `valenx-bioseq`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, BioseqError>`]. The variants are intentionally coarse —
//! a sequence-toolkit caller usually only cares about four things:
//!
//! 1. Did a file / format fail to parse ([`BioseqError::Parse`])?
//! 2. Did the caller mix up an alphabet or supply an illegal residue
//!    ([`BioseqError::Alphabet`])?
//! 3. Did the caller pass nonsense arguments ([`BioseqError::Invalid`])?
//! 4. Is this a documented stub awaiting deeper work
//!    ([`BioseqError::NotYetImplemented`])?
//!
//! Use [`BioseqError::code`] for stable log/telemetry tagging and
//! [`BioseqError::category`] to bucket failures into Parse / Input /
//! Capability without matching every variant. The pattern mirrors
//! `valenx-occt-surface`'s `OcctSurfaceError`.

use std::fmt;

/// Errors produced by `valenx-bioseq`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BioseqError {
    /// A file or in-memory string failed to parse. `format` names the
    /// expected format (`"fasta"`, `"genbank"`, …); `detail` is a
    /// human-readable reason surfaced verbatim in the UI.
    Parse {
        /// Format being parsed (e.g. `"fasta"`, `"fastq"`, `"genbank"`).
        format: &'static str,
        /// Human-readable parse-failure reason.
        detail: String,
    },

    /// An illegal residue or an alphabet mismatch — e.g. a `U` in a
    /// sequence declared DNA, or an unknown IUPAC code.
    Alphabet {
        /// The offending residue (uppercased), or a short descriptor.
        residue: char,
        /// Which alphabet rejected it (e.g. `"DNA"`, `"protein"`).
        alphabet: &'static str,
    },

    /// Caller passed an argument the algorithm cannot accept: an empty
    /// input, a non-positive length, mismatched lengths, an
    /// out-of-range coordinate, etc. A property of the *call*, not of
    /// the data being parsed.
    Invalid {
        /// Logical parameter name (e.g. `"k"`, `"min_len"`, `"primer"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// Feature is part of this crate's public surface but not yet
    /// implemented as a real v1. The string identifies which algorithm
    /// the caller asked for.
    NotYetImplemented {
        /// Stable feature identifier (e.g. `"rna_tertiary"`).
        feature: &'static str,
    },

    /// A GenBank / EMBL feature location references a different record
    /// (e.g. `accession:1..100`). Resolution requires the referenced
    /// record, which the single-record parsers cannot retrieve. The
    /// caller may catch this variant to fetch the cross-referenced
    /// record and re-parse the location in context, or skip the
    /// feature.
    CrossRecordLocation {
        /// The accession the location refers to (the text before the
        /// colon, e.g. `"J00194.1"`).
        accession: String,
        /// The raw location string, verbatim from the flat file
        /// (e.g. `"J00194.1:1..100"`).
        raw: String,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this
/// rather than on 4+ error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A file / string failed to parse.
    Parse,
    /// User-supplied input is wrong (bad alphabet or bad argument).
    Input,
    /// Capability not available in v1 (documented stub).
    Capability,
}

impl BioseqError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"bioseq.<sub_id>"`. Codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            BioseqError::Parse { .. } => "bioseq.parse",
            BioseqError::Alphabet { .. } => "bioseq.alphabet",
            BioseqError::Invalid { .. } => "bioseq.invalid",
            BioseqError::NotYetImplemented { .. } => "bioseq.not_yet_implemented",
            BioseqError::CrossRecordLocation { .. } => "bioseq.cross_record_location",
        }
    }

    /// Coarse category — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            BioseqError::Parse { .. } => "parse",
            BioseqError::Alphabet { .. } | BioseqError::Invalid { .. } => "input",
            BioseqError::NotYetImplemented { .. } | BioseqError::CrossRecordLocation { .. } => {
                "capability"
            }
        }
    }

    /// Typed category enum (for callers that want to `match` instead of
    /// comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            BioseqError::Parse { .. } => ErrorCategory::Parse,
            BioseqError::Alphabet { .. } | BioseqError::Invalid { .. } => ErrorCategory::Input,
            BioseqError::NotYetImplemented { .. } | BioseqError::CrossRecordLocation { .. } => {
                ErrorCategory::Capability
            }
        }
    }

    /// Convenience constructor for [`BioseqError::Parse`].
    pub fn parse(format: &'static str, detail: impl Into<String>) -> Self {
        BioseqError::Parse {
            format,
            detail: detail.into(),
        }
    }

    /// Convenience constructor for [`BioseqError::Alphabet`].
    pub fn alphabet(residue: char, alphabet: &'static str) -> Self {
        BioseqError::Alphabet { residue, alphabet }
    }

    /// Convenience constructor for [`BioseqError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        BioseqError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`BioseqError::NotYetImplemented`].
    pub fn not_yet(feature: &'static str) -> Self {
        BioseqError::NotYetImplemented { feature }
    }

    /// Convenience constructor for [`BioseqError::CrossRecordLocation`].
    pub fn cross_record_location(accession: impl Into<String>, raw: impl Into<String>) -> Self {
        BioseqError::CrossRecordLocation {
            accession: accession.into(),
            raw: raw.into(),
        }
    }
}

impl fmt::Display for BioseqError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BioseqError::Parse { format, detail } => {
                write!(f, "{format} parse error: {detail}")
            }
            BioseqError::Alphabet { residue, alphabet } => {
                write!(f, "invalid residue `{residue}` for {alphabet} alphabet")
            }
            BioseqError::Invalid { what, reason } => {
                write!(f, "invalid `{what}`: {reason}")
            }
            BioseqError::NotYetImplemented { feature } => {
                write!(
                    f,
                    "bioseq feature `{feature}` is not yet implemented (v1 scaffold)"
                )
            }
            BioseqError::CrossRecordLocation { accession, raw } => {
                write!(
                    f,
                    "cross-record location `{raw}` references accession `{accession}` — \
                     not resolvable without the referenced record"
                )
            }
        }
    }
}

impl std::error::Error for BioseqError {}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, BioseqError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = BioseqError::parse("fasta", "missing header");
        assert_eq!(err.code(), "bioseq.parse");
        assert_eq!(err.category(), "parse");
        assert_eq!(err.category_enum(), ErrorCategory::Parse);

        let err = BioseqError::alphabet('U', "DNA");
        assert_eq!(err.code(), "bioseq.alphabet");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = BioseqError::invalid("k", "must be positive");
        assert_eq!(err.code(), "bioseq.invalid");
        assert_eq!(err.category(), "input");

        let err = BioseqError::not_yet("rna_tertiary");
        assert_eq!(err.code(), "bioseq.not_yet_implemented");
        assert_eq!(err.category(), "capability");
        assert_eq!(err.category_enum(), ErrorCategory::Capability);
    }

    #[test]
    fn display_is_informative() {
        let msg = BioseqError::alphabet('Z', "protein").to_string();
        assert!(msg.contains('Z'), "got: {msg}");
        assert!(msg.contains("protein"), "got: {msg}");

        let msg = BioseqError::not_yet("foo").to_string();
        assert!(msg.contains("foo"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(BioseqError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }

    #[test]
    fn cross_record_location_error_is_typed_and_carries_accession() {
        let err = BioseqError::cross_record_location("J00194.1", "J00194.1:1..100");
        assert_eq!(err.code(), "bioseq.cross_record_location");
        assert_eq!(err.category_enum(), ErrorCategory::Capability);
        let msg = err.to_string();
        assert!(msg.contains("J00194.1"), "got: {msg}");
        assert!(msg.contains("1..100"), "got: {msg}");
        // The accession is recoverable from the typed variant.
        match err {
            BioseqError::CrossRecordLocation { accession, raw } => {
                assert_eq!(accession, "J00194.1");
                assert_eq!(raw, "J00194.1:1..100");
            }
            _ => panic!("wrong variant"),
        }
    }
}
