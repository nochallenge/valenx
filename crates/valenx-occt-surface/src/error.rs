//! Error taxonomy for the OCCT surface-modeling feature set.
//!
//! Every public function in this crate returns
//! [`Result<_, OcctSurfaceError>`]. The variants are intentionally
//! coarse: this crate stitches 31 OCCT-equivalent features onto a
//! Rust-native kernel, so most callers only care about three things:
//!
//! 1. Did the caller pass nonsense ([`OcctSurfaceError::BadInput`])?
//! 2. Did the underlying truck/valenx-cad layer refuse the request
//!    ([`OcctSurfaceError::TruckLimit`])?
//! 3. Is this feature a documented stub awaiting deep work
//!    ([`OcctSurfaceError::NotYetImplemented`])?
//!
//! Use [`OcctSurfaceError::code`] for log/telemetry tagging and
//! [`OcctSurfaceError::category`] to classify failures into Input /
//! Capability / Io buckets without matching every variant.

use std::io;

use thiserror::Error;

/// Errors produced by `valenx-occt-surface`.
#[derive(Debug, Error)]
pub enum OcctSurfaceError {
    /// Feature is documented in this crate's public API surface but
    /// not yet implemented. The string identifies which OCCT-equivalent
    /// API the caller asked for so the UI / telemetry can suggest the
    /// concrete follow-up phase that will deliver it (typically Phase
    /// `N.5` where `N` is the originating phase index).
    #[error("occt-surface feature `{feature}` is not yet implemented (v1 scaffold; deep impl tracked in Phase 70.5+)")]
    NotYetImplemented {
        /// Stable feature identifier (e.g. `"pipe_shell"`, `"feat_make_revol"`).
        feature: &'static str,
    },

    /// Caller passed a parameter the underlying kernel cannot accept.
    /// Use this for shape-of-input violations: empty input lists,
    /// non-finite dimensions, mismatched array lengths, etc. Anything
    /// that's a property of the *call* rather than the kernel state.
    #[error("bad input: `{field}` — {reason}")]
    BadInput {
        /// Logical parameter name (e.g. `"radius"`, `"profiles"`).
        field: &'static str,
        /// Human-readable reason, surfaced verbatim in the UI.
        reason: String,
    },

    /// truck (or the downstream `valenx-cad` wrapper) refused the
    /// requested op. This is the "known kernel limitation" channel —
    /// surface it verbatim so users know they hit a truck issue, not
    /// a Valenx bug. Real OCCT would typically succeed where truck
    /// returns this.
    #[error("truck/kernel limitation: {0}")]
    TruckLimit(String),

    /// I/O failure during persist/restore round-trips. Most modules
    /// don't touch disk, but the ones that do (e.g. fit-curve LSQ
    /// dump for debugging) need a typed I/O variant for consistency.
    #[error("io: {0}")]
    Io(#[from] io::Error),
}

/// Coarse category for routing / display purposes.
///
/// Use this to switch a single `match` against three buckets rather
/// than 4+ variants. Stable across crate versions.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// User-supplied input is wrong.
    Input,
    /// Feature/capability not available in v1 (either stub or
    /// kernel-limitation).
    Capability,
    /// I/O subsystem failure.
    Io,
}

impl OcctSurfaceError {
    /// Stable kebab-cased error code suitable for log/telemetry
    /// tagging. Format: `"occt_surface.<sub_id>"`. Codes never change
    /// across minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            OcctSurfaceError::NotYetImplemented { .. } => "occt_surface.not_yet_implemented",
            OcctSurfaceError::BadInput { .. } => "occt_surface.bad_input",
            OcctSurfaceError::TruckLimit(_) => "occt_surface.truck_limit",
            OcctSurfaceError::Io(_) => "occt_surface.io",
        }
    }

    /// Coarse category — see [`ErrorCategory`] for the meaning of
    /// each bucket.
    pub fn category(&self) -> ErrorCategory {
        match self {
            OcctSurfaceError::BadInput { .. } => ErrorCategory::Input,
            OcctSurfaceError::NotYetImplemented { .. }
            | OcctSurfaceError::TruckLimit(_) => ErrorCategory::Capability,
            OcctSurfaceError::Io(_) => ErrorCategory::Io,
        }
    }

    /// Convenience constructor — most modules build this once at the
    /// top of their stub function body.
    pub fn not_yet(feature: &'static str) -> Self {
        OcctSurfaceError::NotYetImplemented { feature }
    }

    /// Convenience constructor for the BadInput variant.
    pub fn bad_input(field: &'static str, reason: impl Into<String>) -> Self {
        OcctSurfaceError::BadInput {
            field,
            reason: reason.into(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = OcctSurfaceError::not_yet("pipe_shell");
        assert_eq!(err.code(), "occt_surface.not_yet_implemented");
        assert_eq!(err.category(), ErrorCategory::Capability);

        let err = OcctSurfaceError::bad_input("radius", "must be positive");
        assert_eq!(err.code(), "occt_surface.bad_input");
        assert_eq!(err.category(), ErrorCategory::Input);

        let err = OcctSurfaceError::TruckLimit("no fillet algo".into());
        assert_eq!(err.code(), "occt_surface.truck_limit");
        assert_eq!(err.category(), ErrorCategory::Capability);

        let err: OcctSurfaceError = io::Error::other("disk full").into();
        assert_eq!(err.code(), "occt_surface.io");
        assert_eq!(err.category(), ErrorCategory::Io);
    }

    #[test]
    fn display_includes_feature_name() {
        let err = OcctSurfaceError::not_yet("approx_curve_fit");
        let msg = err.to_string();
        assert!(msg.contains("approx_curve_fit"), "got: {msg}");
    }
}
