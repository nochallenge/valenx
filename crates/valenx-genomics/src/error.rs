//! Error taxonomy for `valenx-genomics`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, GenomicsError>`]. The variants are intentionally coarse
//! — an NGS-tooling caller usually only cares about four things:
//!
//! 1. Did a record / file fail to parse ([`GenomicsError::Parse`])?
//! 2. Is a constructed record itself structurally invalid — a CIGAR
//!    that does not consume the read, a VCF genotype outside the
//!    declared ploidy, a BED interval with `end < start`
//!    ([`GenomicsError::InvalidRecord`])?
//! 3. Did the caller pass nonsense arguments — an empty read set, a
//!    zero k-mer length, a negative error rate
//!    ([`GenomicsError::Invalid`])?
//! 4. Is this a documented stub awaiting deeper work
//!    ([`GenomicsError::NotYetImplemented`])?
//!
//! Use [`GenomicsError::code`] for stable log / telemetry tagging and
//! [`GenomicsError::category`] to bucket failures into Parse / Input /
//! Capability without matching every variant. The pattern mirrors
//! `valenx-bioseq`'s `BioseqError` and `valenx-cheminf`'s
//! `CheminfError`.

use std::fmt;

/// Errors produced by `valenx-genomics`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GenomicsError {
    /// A record string or file failed to parse. `format` names the
    /// expected format (`"sam"`, `"vcf"`, `"bed"`, `"gff"`, `"pileup"`);
    /// `detail` is a human-readable reason surfaced verbatim in the UI.
    Parse {
        /// Format being parsed (`"sam"`, `"vcf"`, `"bed"`, `"gff"`, …).
        format: &'static str,
        /// Human-readable parse-failure reason.
        detail: String,
    },

    /// A parsed record is structurally inconsistent: a SAM CIGAR whose
    /// query length disagrees with `SEQ`, a VCF genotype index past the
    /// allele list, a BED interval with `end < start`, etc. A property
    /// of the *record*, not of a parse position or a caller argument.
    InvalidRecord {
        /// Record kind (`"sam"`, `"vcf"`, `"bed"`, …).
        kind: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// Caller passed an argument the algorithm cannot accept: an empty
    /// input, a non-positive count, a probability outside `[0, 1]`, a
    /// k-mer length larger than the reads, etc. A property of the
    /// *call*.
    Invalid {
        /// Logical parameter name (e.g. `"k"`, `"error_rate"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// Feature is part of this crate's public surface but not yet
    /// implemented as a real v1. The string identifies which algorithm
    /// the caller asked for.
    NotYetImplemented {
        /// Stable feature identifier (e.g. `"bam_binary_io"`).
        feature: &'static str,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this
/// rather than on 4+ error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A record string / file failed to parse.
    Parse,
    /// User-supplied input is wrong (bad record or bad argument).
    Input,
    /// Capability not available in v1 (documented stub).
    Capability,
}

impl GenomicsError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"genomics.<sub_id>"`. Codes never change
    /// across minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            GenomicsError::Parse { .. } => "genomics.parse",
            GenomicsError::InvalidRecord { .. } => "genomics.invalid_record",
            GenomicsError::Invalid { .. } => "genomics.invalid",
            GenomicsError::NotYetImplemented { .. } => "genomics.not_yet_implemented",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            GenomicsError::Parse { .. } => "parse",
            GenomicsError::InvalidRecord { .. } | GenomicsError::Invalid { .. } => "input",
            GenomicsError::NotYetImplemented { .. } => "capability",
        }
    }

    /// Typed category enum (for callers that want to `match` instead of
    /// comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            GenomicsError::Parse { .. } => ErrorCategory::Parse,
            GenomicsError::InvalidRecord { .. } | GenomicsError::Invalid { .. } => {
                ErrorCategory::Input
            }
            GenomicsError::NotYetImplemented { .. } => ErrorCategory::Capability,
        }
    }

    /// Convenience constructor for [`GenomicsError::Parse`].
    pub fn parse(format: &'static str, detail: impl Into<String>) -> Self {
        GenomicsError::Parse {
            format,
            detail: detail.into(),
        }
    }

    /// Convenience constructor for [`GenomicsError::InvalidRecord`].
    pub fn invalid_record(kind: &'static str, reason: impl Into<String>) -> Self {
        GenomicsError::InvalidRecord {
            kind,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`GenomicsError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        GenomicsError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`GenomicsError::NotYetImplemented`].
    pub fn not_yet(feature: &'static str) -> Self {
        GenomicsError::NotYetImplemented { feature }
    }
}

impl fmt::Display for GenomicsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GenomicsError::Parse { format, detail } => {
                write!(f, "{format} parse error: {detail}")
            }
            GenomicsError::InvalidRecord { kind, reason } => {
                write!(f, "invalid {kind} record: {reason}")
            }
            GenomicsError::Invalid { what, reason } => {
                write!(f, "invalid `{what}`: {reason}")
            }
            GenomicsError::NotYetImplemented { feature } => {
                write!(
                    f,
                    "genomics feature `{feature}` is not yet implemented (v1 scaffold)"
                )
            }
        }
    }
}

impl std::error::Error for GenomicsError {}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, GenomicsError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = GenomicsError::parse("vcf", "bad header");
        assert_eq!(err.code(), "genomics.parse");
        assert_eq!(err.category(), "parse");
        assert_eq!(err.category_enum(), ErrorCategory::Parse);

        let err = GenomicsError::invalid_record("sam", "cigar mismatch");
        assert_eq!(err.code(), "genomics.invalid_record");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = GenomicsError::invalid("k", "must be positive");
        assert_eq!(err.code(), "genomics.invalid");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = GenomicsError::not_yet("bam_binary_io");
        assert_eq!(err.code(), "genomics.not_yet_implemented");
        assert_eq!(err.category(), "capability");
        assert_eq!(err.category_enum(), ErrorCategory::Capability);
    }

    #[test]
    fn display_is_informative() {
        let msg = GenomicsError::parse("bed", "too few columns").to_string();
        assert!(msg.contains("bed"), "got: {msg}");
        assert!(msg.contains("too few columns"), "got: {msg}");

        let msg = GenomicsError::not_yet("foo").to_string();
        assert!(msg.contains("foo"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(GenomicsError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }
}
