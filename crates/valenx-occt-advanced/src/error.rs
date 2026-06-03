//! Error taxonomy for the OCCT advanced-operations feature set.
//!
//! Every public function in this crate returns
//! [`Result<_, OcctAdvancedError>`]. The variants are intentionally
//! coarse: this crate stitches 30 OCCT-equivalent advanced builders /
//! analyzers / upgraders / geometric-query APIs onto Rust-native
//! backends (truck-modeling for topology, valenx-surface for NURBS
//! curve/surface evaluation). Most callers only care about three
//! things:
//!
//! 1. Did the caller pass nonsense ([`OcctAdvancedError::BadInput`])?
//! 2. Is this feature a documented stub awaiting deep work
//!    ([`OcctAdvancedError::NotYetImplemented`])?
//! 3. Did the analysis report a defect that prevents the op
//!    ([`OcctAdvancedError::Defect`])?
//!
//! Use [`OcctAdvancedError::code`] for log/telemetry tagging and
//! [`OcctAdvancedError::category`] to classify failures into Input /
//! Capability / Io buckets without matching every variant.
//!
//! ## Why `Defect` is its own variant
//!
//! Shape-analysis APIs ([`crate::shape_analysis_freebounds()`],
//! [`crate::shape_analysis_orientedclosedsolid()`], etc.) commonly hit
//! "input shape itself is broken" cases that aren't I/O or
//! NotYetImplemented — they're real findings about the geometry.
//! Callers usually want to either repair the shape (via the matching
//! `shape_upgrade_*` op) or surface the defect to the user. Carrying
//! a structured defect payload keeps that branch typed.

use std::io;

use thiserror::Error;

/// Errors produced by `valenx-occt-advanced`.
#[derive(Debug, Error)]
pub enum OcctAdvancedError {
    /// Feature is documented in this crate's public API surface but
    /// not yet implemented. The string identifies which OCCT-equivalent
    /// API the caller asked for so the UI / telemetry can suggest the
    /// concrete follow-up phase that will deliver it (typically Phase
    /// `N.5` where `N` is the originating phase index).
    #[error("occt-advanced feature `{feature}` is not yet implemented (v1 scaffold; deep impl tracked in Phase 131.5+)")]
    NotYetImplemented {
        /// Stable feature identifier (e.g.
        /// `"offset_api_thru_sections_with_guides"`,
        /// `"shape_upgrade_unifysamedomain"`).
        feature: &'static str,
    },

    /// Caller passed a parameter the underlying kernel cannot accept.
    /// Use this for shape-of-input violations: empty input lists,
    /// non-finite dimensions, mismatched array lengths, out-of-range
    /// parameters. Anything that's a property of the *call* rather
    /// than the kernel state.
    #[error("bad input: `{field}` — {reason}")]
    BadInput {
        /// Logical parameter name (e.g. `"u_param"`, `"profile_xy"`,
        /// `"face_index"`).
        field: &'static str,
        /// Human-readable reason, surfaced verbatim in the UI.
        reason: String,
    },

    /// truck (or the downstream `valenx-cad` / `valenx-surface`
    /// wrapper) refused the requested op. This is the "known kernel
    /// limitation" channel — surface it verbatim so users know they
    /// hit a truck issue, not a Valenx bug. Real OCCT would typically
    /// succeed where truck returns this.
    #[error("backend limitation: {0}")]
    Backend(String),

    /// Shape-analysis op completed and reported a geometric defect in
    /// the *input* shape. Distinct from `BadInput` (caller's shape-of-
    /// argument bug) — the input is well-formed at the type level but
    /// the geometry itself violates an invariant (open shell, reversed
    /// face, degenerate edge, etc.). Pair the matching
    /// `shape_upgrade_*` op to attempt a repair.
    #[error("shape defect at {locus}: {kind}")]
    Defect {
        /// Where in the shape the defect was found (face index, edge
        /// index, wire id — caller's choice). Stable for round-tripping
        /// back into a repair call.
        locus: String,
        /// What went wrong, surfaced verbatim in the UI.
        kind: String,
    },

    /// I/O failure during persist/restore round-trips. Most modules
    /// don't touch disk, but the diagnostic dumps in `shape_analysis_*`
    /// optionally write a defect-report file for the UI to render.
    #[error("io: {0}")]
    Io(#[from] io::Error),
}

/// Coarse category for routing / display purposes.
///
/// Use this to switch a single `match` against three buckets rather
/// than 5+ variants. Stable across crate versions.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// User-supplied input is wrong.
    Input,
    /// Feature/capability not available in v1 (either stub or
    /// backend-limitation) — or the input shape itself is defective
    /// in a way the analyzer flagged but the op cannot recover from.
    Capability,
    /// I/O subsystem failure.
    Io,
}

impl OcctAdvancedError {
    /// Stable kebab-cased error code suitable for log/telemetry
    /// tagging. Format: `"occt_advanced.<sub_id>"`. Codes never change
    /// across minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            OcctAdvancedError::NotYetImplemented { .. } => "occt_advanced.not_yet_implemented",
            OcctAdvancedError::BadInput { .. } => "occt_advanced.bad_input",
            OcctAdvancedError::Backend(_) => "occt_advanced.backend",
            OcctAdvancedError::Defect { .. } => "occt_advanced.defect",
            OcctAdvancedError::Io(_) => "occt_advanced.io",
        }
    }

    /// Coarse category — see [`ErrorCategory`] for the meaning of
    /// each bucket. `Defect` lands in `Capability` (the op cannot
    /// proceed against the given geometry), `BadInput` lands in
    /// `Input`.
    pub fn category(&self) -> ErrorCategory {
        match self {
            OcctAdvancedError::BadInput { .. } => ErrorCategory::Input,
            OcctAdvancedError::NotYetImplemented { .. }
            | OcctAdvancedError::Backend(_)
            | OcctAdvancedError::Defect { .. } => ErrorCategory::Capability,
            OcctAdvancedError::Io(_) => ErrorCategory::Io,
        }
    }

    /// Convenience constructor — most modules build this once at the
    /// top of their stub function body.
    pub fn not_yet(feature: &'static str) -> Self {
        OcctAdvancedError::NotYetImplemented { feature }
    }

    /// Convenience constructor for the BadInput variant.
    pub fn bad_input(field: &'static str, reason: impl Into<String>) -> Self {
        OcctAdvancedError::BadInput {
            field,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for the Defect variant.
    pub fn defect(locus: impl Into<String>, kind: impl Into<String>) -> Self {
        OcctAdvancedError::Defect {
            locus: locus.into(),
            kind: kind.into(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = OcctAdvancedError::not_yet("offset_api_thru_sections_with_guides");
        assert_eq!(err.code(), "occt_advanced.not_yet_implemented");
        assert_eq!(err.category(), ErrorCategory::Capability);

        let err = OcctAdvancedError::bad_input("u_param", "must be in [0, 1]");
        assert_eq!(err.code(), "occt_advanced.bad_input");
        assert_eq!(err.category(), ErrorCategory::Input);

        let err = OcctAdvancedError::Backend("truck refused split".into());
        assert_eq!(err.code(), "occt_advanced.backend");
        assert_eq!(err.category(), ErrorCategory::Capability);

        let err = OcctAdvancedError::defect("edge[3]", "degenerate (length < 1e-9)");
        assert_eq!(err.code(), "occt_advanced.defect");
        assert_eq!(err.category(), ErrorCategory::Capability);

        let err: OcctAdvancedError = io::Error::other("disk full").into();
        assert_eq!(err.code(), "occt_advanced.io");
        assert_eq!(err.category(), ErrorCategory::Io);
    }

    #[test]
    fn display_includes_feature_name() {
        let err = OcctAdvancedError::not_yet("shape_upgrade_unifysamedomain");
        let msg = err.to_string();
        assert!(msg.contains("shape_upgrade_unifysamedomain"), "got: {msg}");
    }

    #[test]
    fn defect_display_includes_locus_and_kind() {
        let err = OcctAdvancedError::defect("face[2]", "reversed orientation");
        let msg = err.to_string();
        assert!(msg.contains("face[2]"), "got: {msg}");
        assert!(msg.contains("reversed orientation"), "got: {msg}");
    }
}
