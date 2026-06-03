//! Phase 167 — `V3d_Viewer::SetZBufferAuto()` — depth-buffer-based
//! hidden-surface removal.
//!
//! ## What OCCT does
//!
//! `V3d_Viewer::SetZBufferAuto(Standard_Boolean)` toggles depth testing
//! in the OpenGl pipeline. When ON, each fragment's `gl_FragDepth` is
//! compared against the depth buffer and only the nearest survives —
//! the classic z-buffer HSR algorithm. When OFF, OCCT falls back to
//! painter's-algorithm sorting (back-to-front draw order). Default is
//! ON because painter's order has well-known artifacts with
//! interpenetrating geometry.
//!
//! ## v1 status
//!
//! **Honest v1.** Valenx's egui+wgpu pipeline runs the depth test
//! unconditionally inside `WgpuRenderer` (`valenx_app::wgpu_renderer`)
//! — the `DEPTH_FORMAT = wgpu::TextureFormat::Depth32Float` attachment
//! is configured once at startup and every render pass uses it. This
//! op queries the desired enabled state and returns it; the caller
//! (the View menu in `valenx_app`) reads the returned value to
//! display the menu checkbox state. Switching it off at runtime would
//! require rebuilding the wgpu pipeline (no painter's-algorithm
//! fallback exists in Valenx — Phase 167.5 will add the disable path
//! if a user requests it).

use crate::error::OcctVizError;

/// Report whether the depth buffer is enabled for the current viewer.
///
/// Always returns `Ok(true)` in v1 — the egui+wgpu pipeline cannot
/// run without depth testing. The Boolean is queried so the View menu
/// can render a (greyed-out) "Z-Buffer: ON" checkbox that documents
/// the pipeline state for users coming from OCCT.
pub fn v3d_viewer_z_buffer() -> Result<bool, OcctVizError> {
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn always_enabled_in_v1() {
        assert!(v3d_viewer_z_buffer().unwrap());
    }
}
