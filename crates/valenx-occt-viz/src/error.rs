//! Error taxonomy for the OCCT visualization-pattern feature set.
//!
//! Every public function in this crate returns
//! [`Result<_, OcctVizError>`]. The variants are intentionally coarse:
//! this crate adapts 40 OCCT-equivalent V3d / AIS / Prs3d / animation
//! APIs onto Valenx's existing egui + wgpu pipeline rather than direct-
//! porting the C++ Aspect/OpenGl layer. Most callers only care about
//! three things:
//!
//! 1. Did the caller pass nonsense ([`OcctVizError::BadInput`])?
//! 2. Is this feature a documented stub awaiting deep work
//!    ([`OcctVizError::NotYetImplemented`])?
//! 3. Did the renderer / windowing backend refuse the op
//!    ([`OcctVizError::Render`])?
//!
//! Use [`OcctVizError::code`] for log/telemetry tagging and
//! [`OcctVizError::category`] to classify failures into Input /
//! Capability / Io buckets without matching every variant.
//!
//! ## Why `Render` is its own variant
//!
//! Visualization APIs commonly hit "the GPU / surface refused" cases
//! (lost device, out-of-memory, unsupported swap chain format,
//! incompatible texture sampling) that aren't I/O *and* aren't user-
//! input bugs. Carrying a typed payload keeps that branch routable so
//! the UI can surface a "GPU device lost — please restart" toast
//! distinct from "you passed a bad rectangle" or "this feature isn't
//! implemented yet".

use std::io;

use thiserror::Error;

/// Errors produced by `valenx-occt-viz`.
#[derive(Debug, Error)]
pub enum OcctVizError {
    /// Feature is documented in this crate's public API surface but
    /// not yet implemented. The string identifies which OCCT-equivalent
    /// API the caller asked for so the UI / telemetry can suggest the
    /// concrete follow-up phase that will deliver it (typically Phase
    /// `N.5` where `N` is the originating phase index).
    #[error("occt-viz feature `{feature}` is not yet implemented (v1 scaffold; deep impl tracked in Phase 161.5+)")]
    NotYetImplemented {
        /// Stable feature identifier (e.g.
        /// `"ais_select_polygon"`,
        /// `"transformation_rotation_widget"`).
        feature: &'static str,
    },

    /// Caller passed a parameter the windowing / rendering layer
    /// cannot accept. Use this for shape-of-input violations: empty
    /// selection rectangles, non-finite camera angles, out-of-range
    /// transparency, mismatched array lengths. Anything that's a
    /// property of the *call* rather than backend state.
    #[error("bad input: `{field}` — {reason}")]
    BadInput {
        /// Logical parameter name (e.g. `"distance"`, `"transparency"`,
        /// `"clipping_planes"`).
        field: &'static str,
        /// Human-readable reason, surfaced verbatim in the UI.
        reason: String,
    },

    /// The egui + wgpu backend refused the requested op. This is the
    /// "known renderer limitation" channel — surface it verbatim so
    /// users know they hit a backend issue (lost device, unsupported
    /// surface format, etc.), not a Valenx bug.
    #[error("render backend: {0}")]
    Render(String),

    /// I/O failure during screenshot / video-frame persistence. Most
    /// modules don't touch disk, but
    /// [`crate::view_screenshot::view_screenshot()`] and
    /// [`crate::view_video_export::view_video_export()`] write image
    /// data so callers see this variant when the file system refuses.
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
    /// renderer-backend limitation).
    Capability,
    /// I/O subsystem failure.
    Io,
}

impl OcctVizError {
    /// Stable kebab-cased error code suitable for log/telemetry
    /// tagging. Format: `"occt_viz.<sub_id>"`. Codes never change
    /// across minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            OcctVizError::NotYetImplemented { .. } => "occt_viz.not_yet_implemented",
            OcctVizError::BadInput { .. } => "occt_viz.bad_input",
            OcctVizError::Render(_) => "occt_viz.render",
            OcctVizError::Io(_) => "occt_viz.io",
        }
    }

    /// Coarse category — see [`ErrorCategory`] for the meaning of
    /// each bucket. `Render` lands in `Capability` (the renderer
    /// cannot proceed against the given state), `BadInput` lands in
    /// `Input`.
    pub fn category(&self) -> ErrorCategory {
        match self {
            OcctVizError::BadInput { .. } => ErrorCategory::Input,
            OcctVizError::NotYetImplemented { .. } | OcctVizError::Render(_) => {
                ErrorCategory::Capability
            }
            OcctVizError::Io(_) => ErrorCategory::Io,
        }
    }

    /// Convenience constructor — most modules build this once at the
    /// top of their stub function body.
    pub fn not_yet(feature: &'static str) -> Self {
        OcctVizError::NotYetImplemented { feature }
    }

    /// Convenience constructor for the BadInput variant.
    pub fn bad_input(field: &'static str, reason: impl Into<String>) -> Self {
        OcctVizError::BadInput {
            field,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for the Render variant.
    pub fn render(msg: impl Into<String>) -> Self {
        OcctVizError::Render(msg.into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = OcctVizError::not_yet("ais_select_polygon");
        assert_eq!(err.code(), "occt_viz.not_yet_implemented");
        assert_eq!(err.category(), ErrorCategory::Capability);

        let err = OcctVizError::bad_input("transparency", "must be in [0, 1]");
        assert_eq!(err.code(), "occt_viz.bad_input");
        assert_eq!(err.category(), ErrorCategory::Input);

        let err = OcctVizError::render("wgpu device lost");
        assert_eq!(err.code(), "occt_viz.render");
        assert_eq!(err.category(), ErrorCategory::Capability);

        let err: OcctVizError = io::Error::other("disk full").into();
        assert_eq!(err.code(), "occt_viz.io");
        assert_eq!(err.category(), ErrorCategory::Io);
    }

    #[test]
    fn display_includes_feature_name() {
        let err = OcctVizError::not_yet("transformation_rotation_widget");
        let msg = err.to_string();
        assert!(msg.contains("transformation_rotation_widget"), "got: {msg}");
    }
}
