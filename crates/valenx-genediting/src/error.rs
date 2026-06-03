//! Error taxonomy for `valenx-genediting`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, GeneditingError>`]. The variants are intentionally
//! coarse — a gene-editing / mRNA-design caller usually only cares
//! about four things:
//!
//! 1. Is the *editing target* itself invalid — a target window
//!    outside the supplied sequence, a desired edit that does not
//!    match the reference base, a CDS whose length is not a multiple
//!    of three ([`GeneditingError::InvalidTarget`])?
//! 2. Did the caller pass nonsense arguments — an empty sequence, a
//!    zero homology-arm length, a PBS length outside the scanned
//!    range ([`GeneditingError::Invalid`])?
//! 3. Was the design search exhausted without a usable result — no
//!    PAM-adjacent guide reaches the locus, no base editor can install
//!    the requested transition, no pegRNA passes the constraints
//!    ([`GeneditingError::NoValidDesign`])?
//! 4. Is this a documented stub awaiting deeper work
//!    ([`GeneditingError::NotYetImplemented`])?
//!
//! Use [`GeneditingError::code`] for stable log / telemetry tagging
//! and [`GeneditingError::category`] to bucket failures into
//! Input / NoDesign / Capability without matching every variant. The
//! pattern mirrors `valenx-genomics`'s `GenomicsError` and
//! `valenx-bioseq`'s `BioseqError`.

use std::fmt;

/// Errors produced by `valenx-genediting`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GeneditingError {
    /// The editing *target* (or an input sequence interpreted as a
    /// target) is structurally invalid: a target span past the end of
    /// the sequence, a non-ACGT base in a region that must be precise,
    /// a CDS length not divisible by three, a desired edit whose
    /// "from" base disagrees with the reference. A property of the
    /// supplied biology, not of a generic argument.
    InvalidTarget {
        /// What kind of target was rejected (`"region"`, `"cds"`,
        /// `"variant"`, `"locus"`, …).
        kind: &'static str,
        /// Human-readable reason surfaced verbatim in the UI.
        reason: String,
    },

    /// Caller passed an argument the algorithm cannot accept: an empty
    /// input, a non-positive length, a probability outside `[0, 1]`, a
    /// scan range whose lower bound exceeds its upper bound, etc. A
    /// property of the *call*.
    Invalid {
        /// Logical parameter name (e.g. `"homology_arm_len"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// The design search ran to completion but found nothing usable —
    /// no guide reaches the locus within the PAM constraints, no base
    /// editor installs the requested transition, no pegRNA satisfies
    /// the length / structure filters. Distinct from [`Invalid`] (the
    /// request was well-formed) — the *biology* simply offers no
    /// solution.
    ///
    /// [`Invalid`]: GeneditingError::Invalid
    NoValidDesign {
        /// What was being designed (`"guide"`, `"base_edit"`,
        /// `"pegrna"`, `"donor"`, …).
        what: &'static str,
        /// Human-readable reason — why the search came up empty.
        reason: String,
    },

    /// Feature is part of this crate's public surface but not yet
    /// implemented as a real v1. The string identifies which workflow
    /// the caller asked for.
    NotYetImplemented {
        /// Stable feature identifier (e.g. `"crispr_a_to_t_edit"`).
        feature: &'static str,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this
/// rather than on 4+ error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// User-supplied input is wrong (a bad target or a bad argument).
    Input,
    /// The request was well-formed but the design search found no
    /// solution.
    NoDesign,
    /// Capability not available in v1 (a documented stub).
    Capability,
}

impl GeneditingError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"genediting.<sub_id>"`. Codes never change
    /// across minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            GeneditingError::InvalidTarget { .. } => "genediting.invalid_target",
            GeneditingError::Invalid { .. } => "genediting.invalid",
            GeneditingError::NoValidDesign { .. } => "genediting.no_valid_design",
            GeneditingError::NotYetImplemented { .. } => "genediting.not_yet_implemented",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            GeneditingError::InvalidTarget { .. } | GeneditingError::Invalid { .. } => "input",
            GeneditingError::NoValidDesign { .. } => "no_design",
            GeneditingError::NotYetImplemented { .. } => "capability",
        }
    }

    /// Typed category enum (for callers that want to `match` instead
    /// of comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            GeneditingError::InvalidTarget { .. } | GeneditingError::Invalid { .. } => {
                ErrorCategory::Input
            }
            GeneditingError::NoValidDesign { .. } => ErrorCategory::NoDesign,
            GeneditingError::NotYetImplemented { .. } => ErrorCategory::Capability,
        }
    }

    /// Convenience constructor for [`GeneditingError::InvalidTarget`].
    pub fn invalid_target(kind: &'static str, reason: impl Into<String>) -> Self {
        GeneditingError::InvalidTarget {
            kind,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`GeneditingError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        GeneditingError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`GeneditingError::NoValidDesign`].
    pub fn no_valid_design(what: &'static str, reason: impl Into<String>) -> Self {
        GeneditingError::NoValidDesign {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`GeneditingError::NotYetImplemented`].
    pub fn not_yet(feature: &'static str) -> Self {
        GeneditingError::NotYetImplemented { feature }
    }
}

impl fmt::Display for GeneditingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GeneditingError::InvalidTarget { kind, reason } => {
                write!(f, "invalid {kind} target: {reason}")
            }
            GeneditingError::Invalid { what, reason } => {
                write!(f, "invalid `{what}`: {reason}")
            }
            GeneditingError::NoValidDesign { what, reason } => {
                write!(f, "no valid {what} design found: {reason}")
            }
            GeneditingError::NotYetImplemented { feature } => {
                write!(
                    f,
                    "gene-editing feature `{feature}` is not yet implemented (v1 scaffold)"
                )
            }
        }
    }
}

impl std::error::Error for GeneditingError {}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, GeneditingError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let e = GeneditingError::invalid_target("region", "span past end");
        assert_eq!(e.code(), "genediting.invalid_target");
        assert_eq!(e.category(), "input");
        assert_eq!(e.category_enum(), ErrorCategory::Input);

        let e = GeneditingError::invalid("homology_arm_len", "must be positive");
        assert_eq!(e.code(), "genediting.invalid");
        assert_eq!(e.category(), "input");
        assert_eq!(e.category_enum(), ErrorCategory::Input);

        let e = GeneditingError::no_valid_design("guide", "no PAM reaches the locus");
        assert_eq!(e.code(), "genediting.no_valid_design");
        assert_eq!(e.category(), "no_design");
        assert_eq!(e.category_enum(), ErrorCategory::NoDesign);

        let e = GeneditingError::not_yet("foo");
        assert_eq!(e.code(), "genediting.not_yet_implemented");
        assert_eq!(e.category(), "capability");
        assert_eq!(e.category_enum(), ErrorCategory::Capability);
    }

    #[test]
    fn display_is_informative() {
        let msg = GeneditingError::no_valid_design("pegrna", "PBS too short").to_string();
        assert!(msg.contains("pegrna"), "got: {msg}");
        assert!(msg.contains("PBS too short"), "got: {msg}");

        let msg = GeneditingError::not_yet("bar").to_string();
        assert!(msg.contains("bar"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> =
            Box::new(GeneditingError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }
}
