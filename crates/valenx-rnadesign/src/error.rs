//! Error taxonomy for `valenx-rnadesign`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, RnaDesignError>`]. The variants are intentionally coarse
//! — a synthetic-RNA-design caller usually only cares about four things:
//!
//! 1. Did the caller hand over an invalid design *goal* or *constraint*
//!    set — an empty / pseudoknotted target, a non-protein "encode this
//!    protein" input, a GC range that excludes itself
//!    ([`RnaDesignError::Goal`])?
//! 2. Did the caller pass nonsense arguments to a workflow step — a
//!    zero iteration budget, a stage transition out of order, a batch
//!    size of zero ([`RnaDesignError::Invalid`])?
//! 3. Did the design / optimisation search run to completion without a
//!    usable candidate — no inverse-fold solution reaches the target,
//!    no synonymous variant clears every constraint
//!    ([`RnaDesignError::NoDesign`])?
//! 4. Is this a documented stub awaiting deeper work
//!    ([`RnaDesignError::NotYetImplemented`])?
//!
//! A fifth variant, [`RnaDesignError::Upstream`], wraps a failure that
//! bubbled up from one of the three building-block crates
//! (`valenx-rnastruct`, `valenx-genediting`, `valenx-bioseq`) — this
//! crate is an orchestration layer, so a folding or codon-optimisation
//! error is surfaced verbatim with its originating crate named.
//!
//! Use [`RnaDesignError::code`] for stable log / telemetry tagging and
//! [`RnaDesignError::category`] to bucket failures into Input / NoDesign
//! / Upstream / Capability without matching every variant. The pattern
//! mirrors `valenx-genediting`'s `GeneditingError` and `valenx-bioseq`'s
//! `BioseqError`.

use std::fmt;

/// Errors produced by `valenx-rnadesign`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RnaDesignError {
    /// The design *goal* or the *constraint* set is itself invalid: an
    /// empty or pseudoknotted target dot-bracket, a protein-encoding
    /// goal whose input is not a protein sequence, a GC range whose
    /// lower bound exceeds its upper bound, a length window that
    /// excludes every length. A property of the supplied design intent,
    /// not of a generic argument.
    Goal {
        /// What kind of goal / constraint was rejected (`"target"`,
        /// `"protein"`, `"gc_range"`, `"length"`, `"constraints"`, …).
        what: &'static str,
        /// Human-readable reason surfaced verbatim in the UI.
        reason: String,
    },

    /// Caller passed an argument a workflow step cannot accept: a
    /// zero-iteration optimisation budget, a batch size of zero, a
    /// stage transition attempted out of order, a weight outside
    /// `[0, 1]`. A property of the *call*.
    Invalid {
        /// Logical parameter name (e.g. `"iterations"`, `"stage"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// The design or optimisation search ran to completion but found
    /// nothing usable — no inverse-fold candidate reaches the target
    /// structure, no synonymous CDS variant clears every constraint,
    /// the two-state design cannot satisfy both target structures.
    /// Distinct from [`Invalid`] (the request was well-formed) — the
    /// *design space* simply offers no solution under the constraints.
    ///
    /// [`Invalid`]: RnaDesignError::Invalid
    NoDesign {
        /// What was being designed (`"structural"`, `"coding"`,
        /// `"riboswitch"`, `"optimization"`, …).
        what: &'static str,
        /// Human-readable reason — why the search came up empty.
        reason: String,
    },

    /// A building-block crate (`valenx-rnastruct`, `valenx-genediting`
    /// or `valenx-bioseq`) returned an error while this crate was
    /// orchestrating it. The originating crate is named and its
    /// human-readable message is surfaced verbatim — this crate adds
    /// the workflow context, it does not swallow the cause.
    Upstream {
        /// The building-block crate the failure came from
        /// (`"valenx-rnastruct"`, `"valenx-genediting"`,
        /// `"valenx-bioseq"`).
        crate_name: &'static str,
        /// The originating error's message.
        message: String,
    },

    /// Feature is part of this crate's public surface but not yet
    /// implemented as a real v1. The string identifies which workflow
    /// the caller asked for.
    NotYetImplemented {
        /// Stable feature identifier.
        feature: &'static str,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this
/// rather than on the five error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// User-supplied input is wrong (a bad goal / constraint or a bad
    /// argument).
    Input,
    /// The request was well-formed but the design search found no
    /// solution.
    NoDesign,
    /// A building-block crate failed underneath this orchestration
    /// layer.
    Upstream,
    /// Capability not available in v1 (a documented stub).
    Capability,
}

impl RnaDesignError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"rnadesign.<sub_id>"`. Codes never change
    /// across minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            RnaDesignError::Goal { .. } => "rnadesign.goal",
            RnaDesignError::Invalid { .. } => "rnadesign.invalid",
            RnaDesignError::NoDesign { .. } => "rnadesign.no_design",
            RnaDesignError::Upstream { .. } => "rnadesign.upstream",
            RnaDesignError::NotYetImplemented { .. } => "rnadesign.not_yet_implemented",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            RnaDesignError::Goal { .. } | RnaDesignError::Invalid { .. } => "input",
            RnaDesignError::NoDesign { .. } => "no_design",
            RnaDesignError::Upstream { .. } => "upstream",
            RnaDesignError::NotYetImplemented { .. } => "capability",
        }
    }

    /// Typed category enum (for callers that want to `match` instead of
    /// comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            RnaDesignError::Goal { .. } | RnaDesignError::Invalid { .. } => ErrorCategory::Input,
            RnaDesignError::NoDesign { .. } => ErrorCategory::NoDesign,
            RnaDesignError::Upstream { .. } => ErrorCategory::Upstream,
            RnaDesignError::NotYetImplemented { .. } => ErrorCategory::Capability,
        }
    }

    /// Convenience constructor for [`RnaDesignError::Goal`].
    pub fn goal(what: &'static str, reason: impl Into<String>) -> Self {
        RnaDesignError::Goal {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`RnaDesignError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        RnaDesignError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`RnaDesignError::NoDesign`].
    pub fn no_design(what: &'static str, reason: impl Into<String>) -> Self {
        RnaDesignError::NoDesign {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`RnaDesignError::Upstream`].
    pub fn upstream(crate_name: &'static str, message: impl Into<String>) -> Self {
        RnaDesignError::Upstream {
            crate_name,
            message: message.into(),
        }
    }

    /// Convenience constructor for [`RnaDesignError::NotYetImplemented`].
    pub fn not_yet(feature: &'static str) -> Self {
        RnaDesignError::NotYetImplemented { feature }
    }
}

// --- `From` conversions for the three building-block crates -----------
//
// These let `?` lift an upstream error straight into `RnaDesignError`,
// always tagging the originating crate so the failure stays traceable.

impl From<valenx_rnastruct::RnaStructError> for RnaDesignError {
    fn from(e: valenx_rnastruct::RnaStructError) -> Self {
        RnaDesignError::upstream("valenx-rnastruct", e.to_string())
    }
}

impl From<valenx_genediting::GeneditingError> for RnaDesignError {
    fn from(e: valenx_genediting::GeneditingError) -> Self {
        RnaDesignError::upstream("valenx-genediting", e.to_string())
    }
}

impl From<valenx_bioseq::BioseqError> for RnaDesignError {
    fn from(e: valenx_bioseq::BioseqError) -> Self {
        RnaDesignError::upstream("valenx-bioseq", e.to_string())
    }
}

impl fmt::Display for RnaDesignError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RnaDesignError::Goal { what, reason } => {
                write!(f, "invalid design goal ({what}): {reason}")
            }
            RnaDesignError::Invalid { what, reason } => {
                write!(f, "invalid `{what}`: {reason}")
            }
            RnaDesignError::NoDesign { what, reason } => {
                write!(f, "no valid {what} design found: {reason}")
            }
            RnaDesignError::Upstream {
                crate_name,
                message,
            } => {
                write!(f, "{crate_name} error: {message}")
            }
            RnaDesignError::NotYetImplemented { feature } => {
                write!(
                    f,
                    "rnadesign feature `{feature}` is not yet implemented (v1 scaffold)"
                )
            }
        }
    }
}

impl std::error::Error for RnaDesignError {}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, RnaDesignError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let e = RnaDesignError::goal("target", "pseudoknotted");
        assert_eq!(e.code(), "rnadesign.goal");
        assert_eq!(e.category(), "input");
        assert_eq!(e.category_enum(), ErrorCategory::Input);

        let e = RnaDesignError::invalid("iterations", "must be positive");
        assert_eq!(e.code(), "rnadesign.invalid");
        assert_eq!(e.category(), "input");
        assert_eq!(e.category_enum(), ErrorCategory::Input);

        let e = RnaDesignError::no_design("structural", "target unreachable");
        assert_eq!(e.code(), "rnadesign.no_design");
        assert_eq!(e.category(), "no_design");
        assert_eq!(e.category_enum(), ErrorCategory::NoDesign);

        let e = RnaDesignError::upstream("valenx-rnastruct", "fold failed");
        assert_eq!(e.code(), "rnadesign.upstream");
        assert_eq!(e.category(), "upstream");
        assert_eq!(e.category_enum(), ErrorCategory::Upstream);

        let e = RnaDesignError::not_yet("tertiary_design");
        assert_eq!(e.code(), "rnadesign.not_yet_implemented");
        assert_eq!(e.category(), "capability");
        assert_eq!(e.category_enum(), ErrorCategory::Capability);
    }

    #[test]
    fn display_is_informative() {
        let msg = RnaDesignError::goal("gc_range", "min exceeds max").to_string();
        assert!(msg.contains("min exceeds max"), "got: {msg}");

        let msg = RnaDesignError::not_yet("foo").to_string();
        assert!(msg.contains("foo"), "got: {msg}");
    }

    #[test]
    fn upstream_errors_convert() {
        // Each building-block error lifts into an Upstream variant that
        // names the originating crate.
        let rs: RnaDesignError = valenx_rnastruct::RnaStructError::sequence("bad").into();
        assert!(matches!(
            rs,
            RnaDesignError::Upstream {
                crate_name: "valenx-rnastruct",
                ..
            }
        ));
        let ge: RnaDesignError = valenx_genediting::GeneditingError::invalid("x", "y").into();
        assert!(matches!(
            ge,
            RnaDesignError::Upstream {
                crate_name: "valenx-genediting",
                ..
            }
        ));
        let bs: RnaDesignError = valenx_bioseq::BioseqError::invalid("x", "y").into();
        assert!(matches!(
            bs,
            RnaDesignError::Upstream {
                crate_name: "valenx-bioseq",
                ..
            }
        ));
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(RnaDesignError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }
}
