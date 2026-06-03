//! Error taxonomy for `valenx-phylo`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, PhyloError>`]. The variants are intentionally coarse — a
//! phylogenetics caller usually only cares about five things:
//!
//! 1. Did a tree / alignment file fail to parse ([`PhyloError::Parse`])?
//! 2. Is a tree structurally invalid — a cycle, a dangling edge, a
//!    duplicate label ([`PhyloError::InvalidTree`])?
//! 3. Did the caller pass nonsense arguments ([`PhyloError::Invalid`])?
//! 4. Do two matrices / sequences disagree on a dimension
//!    ([`PhyloError::Dimension`])?
//! 5. Is this a documented stub awaiting deeper work
//!    ([`PhyloError::NotYetImplemented`])?
//!
//! Use [`PhyloError::code`] for stable log/telemetry tagging and
//! [`PhyloError::category`] to bucket failures into Parse / Input /
//! Capability without matching every variant. The pattern mirrors
//! `valenx-bioseq`'s `BioseqError`.

use std::fmt;

/// Errors produced by `valenx-phylo`.
#[derive(Debug, Clone, PartialEq)]
pub enum PhyloError {
    /// A file or in-memory string failed to parse. `format` names the
    /// expected format (`"newick"`, `"nexus"`, …); `detail` is a
    /// human-readable reason surfaced verbatim in the UI.
    Parse {
        /// Format being parsed (e.g. `"newick"`, `"nexus"`).
        format: &'static str,
        /// Human-readable parse-failure reason.
        detail: String,
    },

    /// A tree is structurally invalid: a cycle, a node with no parent
    /// that is not the root, a duplicate leaf label, a missing child,
    /// etc. A property of the *tree*, not of a call's arguments.
    InvalidTree {
        /// Human-readable reason.
        reason: String,
    },

    /// Caller passed an argument the algorithm cannot accept: an empty
    /// input, a non-positive count, an out-of-range index, a model
    /// parameter outside its domain, etc.
    Invalid {
        /// Logical parameter name (e.g. `"kappa"`, `"n_taxa"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// Two inputs disagree on a dimension — e.g. a 4×3 rate matrix, an
    /// alignment whose rows differ in length, or a distance matrix that
    /// is not square.
    Dimension {
        /// What was expected.
        expected: usize,
        /// What was actually supplied.
        actual: usize,
        /// Short context label (e.g. `"alignment rows"`).
        context: &'static str,
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
    /// User-supplied input is wrong (bad tree, bad argument, bad shape).
    Input,
    /// Capability not available in v1 (documented stub).
    Capability,
}

impl PhyloError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"phylo.<sub_id>"`. Codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            PhyloError::Parse { .. } => "phylo.parse",
            PhyloError::InvalidTree { .. } => "phylo.invalid_tree",
            PhyloError::Invalid { .. } => "phylo.invalid",
            PhyloError::Dimension { .. } => "phylo.dimension",
            PhyloError::NotYetImplemented { .. } => "phylo.not_yet_implemented",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            PhyloError::Parse { .. } => "parse",
            PhyloError::InvalidTree { .. }
            | PhyloError::Invalid { .. }
            | PhyloError::Dimension { .. } => "input",
            PhyloError::NotYetImplemented { .. } => "capability",
        }
    }

    /// Typed category enum (for callers that want to `match` instead of
    /// comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            PhyloError::Parse { .. } => ErrorCategory::Parse,
            PhyloError::InvalidTree { .. }
            | PhyloError::Invalid { .. }
            | PhyloError::Dimension { .. } => ErrorCategory::Input,
            PhyloError::NotYetImplemented { .. } => ErrorCategory::Capability,
        }
    }

    /// Convenience constructor for [`PhyloError::Parse`].
    pub fn parse(format: &'static str, detail: impl Into<String>) -> Self {
        PhyloError::Parse {
            format,
            detail: detail.into(),
        }
    }

    /// Convenience constructor for [`PhyloError::InvalidTree`].
    pub fn invalid_tree(reason: impl Into<String>) -> Self {
        PhyloError::InvalidTree {
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`PhyloError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        PhyloError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`PhyloError::Dimension`].
    pub fn dimension(expected: usize, actual: usize, context: &'static str) -> Self {
        PhyloError::Dimension {
            expected,
            actual,
            context,
        }
    }

    /// Convenience constructor for [`PhyloError::NotYetImplemented`].
    pub fn not_yet(feature: &'static str) -> Self {
        PhyloError::NotYetImplemented { feature }
    }
}

impl fmt::Display for PhyloError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PhyloError::Parse { format, detail } => {
                write!(f, "{format} parse error: {detail}")
            }
            PhyloError::InvalidTree { reason } => {
                write!(f, "invalid tree: {reason}")
            }
            PhyloError::Invalid { what, reason } => {
                write!(f, "invalid `{what}`: {reason}")
            }
            PhyloError::Dimension {
                expected,
                actual,
                context,
            } => {
                write!(
                    f,
                    "dimension mismatch for {context}: expected {expected}, got {actual}"
                )
            }
            PhyloError::NotYetImplemented { feature } => {
                write!(
                    f,
                    "phylo feature `{feature}` is not yet implemented (v1 scaffold)"
                )
            }
        }
    }
}

impl std::error::Error for PhyloError {}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, PhyloError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = PhyloError::parse("newick", "unbalanced parentheses");
        assert_eq!(err.code(), "phylo.parse");
        assert_eq!(err.category(), "parse");
        assert_eq!(err.category_enum(), ErrorCategory::Parse);

        let err = PhyloError::invalid_tree("cycle detected");
        assert_eq!(err.code(), "phylo.invalid_tree");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = PhyloError::invalid("kappa", "must be positive");
        assert_eq!(err.code(), "phylo.invalid");
        assert_eq!(err.category(), "input");

        let err = PhyloError::dimension(4, 3, "rate matrix");
        assert_eq!(err.code(), "phylo.dimension");
        assert_eq!(err.category(), "input");

        let err = PhyloError::not_yet("bayesian_mcmc");
        assert_eq!(err.code(), "phylo.not_yet_implemented");
        assert_eq!(err.category(), "capability");
        assert_eq!(err.category_enum(), ErrorCategory::Capability);
    }

    #[test]
    fn display_is_informative() {
        let msg = PhyloError::dimension(4, 9, "alignment rows").to_string();
        assert!(msg.contains('4') && msg.contains('9'), "got: {msg}");
        assert!(msg.contains("alignment rows"), "got: {msg}");

        let msg = PhyloError::not_yet("foo").to_string();
        assert!(msg.contains("foo"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(PhyloError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }
}
