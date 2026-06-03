//! Error taxonomy for `valenx-popgen`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, PopgenError>`]. The variants are intentionally coarse — a
//! population-genetics caller usually only cares about five things:
//!
//! 1. Did a file / string fail to parse ([`PopgenError::Parse`])?
//! 2. Did the caller pass nonsense arguments — an empty sample, a
//!    negative rate, a probability outside `[0, 1]`
//!    ([`PopgenError::Invalid`])?
//! 3. Do two inputs disagree on a dimension — a genotype matrix whose
//!    rows differ in length, mismatched population labels
//!    ([`PopgenError::Dimension`])?
//! 4. Is a simulation / model in an inconsistent state — a genealogy
//!    with no MRCA, a tree-sequence with a dangling edge
//!    ([`PopgenError::Model`])?
//! 5. Is this a documented stub awaiting deeper work
//!    ([`PopgenError::NotYetImplemented`])?
//!
//! Use [`PopgenError::code`] for stable log / telemetry tagging and
//! [`PopgenError::category`] to bucket failures into Parse / Input /
//! Capability without matching every variant. The pattern mirrors
//! `valenx-bioseq`'s `BioseqError` and `valenx-phylo`'s `PhyloError`.

use std::fmt;

/// Errors produced by `valenx-popgen`.
#[derive(Debug, Clone, PartialEq)]
pub enum PopgenError {
    /// A file or in-memory string failed to parse. `format` names the
    /// expected format (`"vcf"`, `"ms"`, …); `detail` is a
    /// human-readable reason surfaced verbatim in the UI.
    Parse {
        /// Format being parsed (e.g. `"vcf"`, `"ms"`).
        format: &'static str,
        /// Human-readable parse-failure reason.
        detail: String,
    },

    /// Caller passed an argument the algorithm cannot accept: an empty
    /// sample, a non-positive population size, a rate or probability
    /// outside its domain, an out-of-range index, etc. A property of
    /// the *call*, not of data being parsed.
    Invalid {
        /// Logical parameter name (e.g. `"mutation_rate"`, `"n_demes"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// Two inputs disagree on a dimension — a genotype matrix with
    /// ragged rows, a per-deme size vector of the wrong length, a
    /// distance matrix that is not square, etc.
    Dimension {
        /// What was expected.
        expected: usize,
        /// What was actually supplied.
        actual: usize,
        /// Short context label (e.g. `"genotype matrix columns"`).
        context: &'static str,
    },

    /// A simulation, genealogy or tree-sequence is internally
    /// inconsistent: a coalescent that never reached a single MRCA, a
    /// tree-sequence edge referencing a missing node, a population with
    /// no individuals. A property of the *model state*, not of a
    /// call's arguments.
    Model {
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
/// rather than on 5+ error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A file / string failed to parse.
    Parse,
    /// User-supplied input is wrong (bad argument, bad shape, bad
    /// model state).
    Input,
    /// Capability not available in v1 (documented stub).
    Capability,
}

impl PopgenError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"popgen.<sub_id>"`. Codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            PopgenError::Parse { .. } => "popgen.parse",
            PopgenError::Invalid { .. } => "popgen.invalid",
            PopgenError::Dimension { .. } => "popgen.dimension",
            PopgenError::Model { .. } => "popgen.model",
            PopgenError::NotYetImplemented { .. } => "popgen.not_yet_implemented",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            PopgenError::Parse { .. } => "parse",
            PopgenError::Invalid { .. }
            | PopgenError::Dimension { .. }
            | PopgenError::Model { .. } => "input",
            PopgenError::NotYetImplemented { .. } => "capability",
        }
    }

    /// Typed category enum (for callers that want to `match` instead of
    /// comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            PopgenError::Parse { .. } => ErrorCategory::Parse,
            PopgenError::Invalid { .. }
            | PopgenError::Dimension { .. }
            | PopgenError::Model { .. } => ErrorCategory::Input,
            PopgenError::NotYetImplemented { .. } => ErrorCategory::Capability,
        }
    }

    /// Convenience constructor for [`PopgenError::Parse`].
    pub fn parse(format: &'static str, detail: impl Into<String>) -> Self {
        PopgenError::Parse {
            format,
            detail: detail.into(),
        }
    }

    /// Convenience constructor for [`PopgenError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        PopgenError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`PopgenError::Dimension`].
    pub fn dimension(expected: usize, actual: usize, context: &'static str) -> Self {
        PopgenError::Dimension {
            expected,
            actual,
            context,
        }
    }

    /// Convenience constructor for [`PopgenError::Model`].
    pub fn model(reason: impl Into<String>) -> Self {
        PopgenError::Model {
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`PopgenError::NotYetImplemented`].
    pub fn not_yet(feature: &'static str) -> Self {
        PopgenError::NotYetImplemented { feature }
    }
}

impl fmt::Display for PopgenError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PopgenError::Parse { format, detail } => {
                write!(f, "{format} parse error: {detail}")
            }
            PopgenError::Invalid { what, reason } => {
                write!(f, "invalid `{what}`: {reason}")
            }
            PopgenError::Dimension {
                expected,
                actual,
                context,
            } => {
                write!(
                    f,
                    "dimension mismatch for {context}: expected {expected}, got {actual}"
                )
            }
            PopgenError::Model { reason } => {
                write!(f, "inconsistent model: {reason}")
            }
            PopgenError::NotYetImplemented { feature } => {
                write!(
                    f,
                    "popgen feature `{feature}` is not yet implemented (v1 scaffold)"
                )
            }
        }
    }
}

impl std::error::Error for PopgenError {}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, PopgenError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = PopgenError::parse("vcf", "missing #CHROM header");
        assert_eq!(err.code(), "popgen.parse");
        assert_eq!(err.category(), "parse");
        assert_eq!(err.category_enum(), ErrorCategory::Parse);

        let err = PopgenError::invalid("mutation_rate", "must be non-negative");
        assert_eq!(err.code(), "popgen.invalid");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = PopgenError::dimension(4, 3, "genotype columns");
        assert_eq!(err.code(), "popgen.dimension");
        assert_eq!(err.category(), "input");

        let err = PopgenError::model("genealogy has no MRCA");
        assert_eq!(err.code(), "popgen.model");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = PopgenError::not_yet("spatial_continuous");
        assert_eq!(err.code(), "popgen.not_yet_implemented");
        assert_eq!(err.category(), "capability");
        assert_eq!(err.category_enum(), ErrorCategory::Capability);
    }

    #[test]
    fn display_is_informative() {
        let msg = PopgenError::dimension(4, 9, "genotype columns").to_string();
        assert!(msg.contains('4') && msg.contains('9'), "got: {msg}");
        assert!(msg.contains("genotype columns"), "got: {msg}");

        let msg = PopgenError::not_yet("foo").to_string();
        assert!(msg.contains("foo"), "got: {msg}");

        let msg = PopgenError::model("dangling edge").to_string();
        assert!(msg.contains("dangling edge"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> =
            Box::new(PopgenError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }
}
