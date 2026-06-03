//! Error taxonomy for the OCCT data-exchange feature set.
//!
//! Every public function in this crate returns
//! [`Result<_, OcctExchangeError>`]. The variants are intentionally
//! coarse: this crate stitches 30 OCCT-equivalent importers /
//! exporters onto Rust-native backends (truck-stepio, the hand-rolled
//! IGES writer, valenx-mesh's format modules). Most callers only care
//! about three things:
//!
//! 1. Did the caller pass nonsense ([`OcctExchangeError::BadInput`])?
//! 2. Is this format/feature a documented stub awaiting deep work
//!    ([`OcctExchangeError::NotYetImplemented`])?
//! 3. Did disk / parsing fail ([`OcctExchangeError::Io`] /
//!    [`OcctExchangeError::Parse`])?
//!
//! Use [`OcctExchangeError::code`] for log/telemetry tagging and
//! [`OcctExchangeError::category`] to classify failures into Input /
//! Capability / Io buckets without matching every variant.
//!
//! ## Why `Parse` is its own variant
//!
//! Mesh interchange formats (OBJ, PLY, STL, glTF) commonly hit
//! malformed-input cases that aren't I/O failures and aren't the
//! caller's fault — they're the *file*'s fault. We want to distinguish
//! "the disk is broken" (`Io`) from "the file is broken" (`Parse`) so
//! the UI can blame the right thing.

use std::io;

use thiserror::Error;

/// Errors produced by `valenx-occt-exchange`.
#[derive(Debug, Error)]
pub enum OcctExchangeError {
    /// Feature is documented in this crate's public API surface but
    /// not yet implemented. The string identifies which OCCT-equivalent
    /// importer / exporter the caller asked for so the UI / telemetry
    /// can suggest the concrete follow-up phase that will deliver it
    /// (typically Phase `N.5` where `N` is the phase index from
    /// `docs/GOALS.md`).
    #[error("occt-exchange feature `{feature}` is not yet implemented (v1 scaffold; deep impl tracked in Phase 101.5+)")]
    NotYetImplemented {
        /// Stable feature identifier (e.g. `"step_ap242_full_writer"`,
        /// `"gltf2_writer"`).
        feature: &'static str,
    },

    /// Caller passed a parameter the underlying kernel cannot accept.
    /// Use this for shape-of-input violations: empty input lists,
    /// non-finite vertex coordinates, mismatched array lengths, etc.
    /// Anything that's a property of the *call* rather than the kernel
    /// state or the file contents.
    #[error("bad input: `{field}` — {reason}")]
    BadInput {
        /// Logical parameter name (e.g. `"path"`, `"solids"`,
        /// `"vertex_colors"`).
        field: &'static str,
        /// Human-readable reason, surfaced verbatim in the UI.
        reason: String,
    },

    /// The downstream `valenx-step-iges` or `valenx-mesh` backend
    /// refused the requested op. Surface verbatim so users know they
    /// hit a known back-end limitation, not a Valenx bug.
    #[error("backend limitation: {0}")]
    Backend(String),

    /// File parsing failed: malformed line, unexpected token, wrong
    /// magic bytes, etc. Use this when *the file itself* is the
    /// problem, distinct from disk I/O failures.
    #[error("parse error at {context}: {reason}")]
    Parse {
        /// Where in the file the error happened (line number,
        /// byte-offset, entity index — caller's choice).
        context: String,
        /// What went wrong, surfaced verbatim in the UI.
        reason: String,
    },

    /// I/O failure during persist/restore round-trips: file not
    /// found, permission denied, disk full, etc.
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
    /// backend-limitation) — or the file we tried to read is itself
    /// malformed in a way the parser cannot recover from.
    Capability,
    /// I/O subsystem failure.
    Io,
}

impl OcctExchangeError {
    /// Stable kebab-cased error code suitable for log/telemetry
    /// tagging. Format: `"occt_exchange.<sub_id>"`. Codes never change
    /// across minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            OcctExchangeError::NotYetImplemented { .. } => "occt_exchange.not_yet_implemented",
            OcctExchangeError::BadInput { .. } => "occt_exchange.bad_input",
            OcctExchangeError::Backend(_) => "occt_exchange.backend",
            OcctExchangeError::Parse { .. } => "occt_exchange.parse",
            OcctExchangeError::Io(_) => "occt_exchange.io",
        }
    }

    /// Coarse category — see [`ErrorCategory`] for the meaning of
    /// each bucket. `Parse` lands in `Capability` (it's a file-format
    /// limitation surfaced as a runtime failure), `BadInput` lands in
    /// `Input`.
    pub fn category(&self) -> ErrorCategory {
        match self {
            OcctExchangeError::BadInput { .. } => ErrorCategory::Input,
            OcctExchangeError::NotYetImplemented { .. }
            | OcctExchangeError::Backend(_)
            | OcctExchangeError::Parse { .. } => ErrorCategory::Capability,
            OcctExchangeError::Io(_) => ErrorCategory::Io,
        }
    }

    /// Convenience constructor — most modules build this once at the
    /// top of their stub function body.
    pub fn not_yet(feature: &'static str) -> Self {
        OcctExchangeError::NotYetImplemented { feature }
    }

    /// Convenience constructor for the BadInput variant.
    pub fn bad_input(field: &'static str, reason: impl Into<String>) -> Self {
        OcctExchangeError::BadInput {
            field,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for the Parse variant.
    pub fn parse(context: impl Into<String>, reason: impl Into<String>) -> Self {
        OcctExchangeError::Parse {
            context: context.into(),
            reason: reason.into(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = OcctExchangeError::not_yet("step_ap214_writer");
        assert_eq!(err.code(), "occt_exchange.not_yet_implemented");
        assert_eq!(err.category(), ErrorCategory::Capability);

        let err = OcctExchangeError::bad_input("path", "extension must be .step");
        assert_eq!(err.code(), "occt_exchange.bad_input");
        assert_eq!(err.category(), ErrorCategory::Input);

        let err = OcctExchangeError::Backend("ruststep refused".into());
        assert_eq!(err.code(), "occt_exchange.backend");
        assert_eq!(err.category(), ErrorCategory::Capability);

        let err = OcctExchangeError::parse("line 4", "expected `magic`");
        assert_eq!(err.code(), "occt_exchange.parse");
        assert_eq!(err.category(), ErrorCategory::Capability);

        let err: OcctExchangeError = io::Error::other("disk full").into();
        assert_eq!(err.code(), "occt_exchange.io");
        assert_eq!(err.category(), ErrorCategory::Io);
    }

    #[test]
    fn display_includes_feature_name() {
        let err = OcctExchangeError::not_yet("jt_writer");
        let msg = err.to_string();
        assert!(msg.contains("jt_writer"), "got: {msg}");
    }
}
