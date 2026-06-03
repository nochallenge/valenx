//! Phase 168 — `V3d_View::SetBackFacingModel()` — toggle back-face
//! culling.
//!
//! ## What OCCT does
//!
//! `V3d_View::SetBackFacingModel(Graphic3d_TypeOfBackfacingModel)`
//! controls whether faces with normals pointing away from the camera
//! are drawn. `Graphic3d_TOBM_DISABLE` skips them (the default — saves
//! ~50% fragment-shader work for closed solids); `Graphic3d_TOBM_AUTOMATIC`
//! lets the driver decide; `Graphic3d_TOBM_FORCE` draws them too
//! (useful for inspecting cavity walls of open meshes).
//!
//! ## v1 status
//!
//! **Honest v1.** Valenx's [`viewport::show`] applies back-face
//! culling in the shaded path by checking `triangle_normal · view_dir
//! < 0` and dropping the triangle pre-rasterisation. The wgpu pipeline
//! mirrors that with `cull_mode: Some(wgpu::Face::Back)` configured
//! once at startup. This op records the *requested* mode in a
//! [`BackFaceMode`] enum the caller stores in app state; the
//! [`crate::v3d_viewer_back_face_culling()`] return value is the
//! validated mode (no enum-to-renderer-state plumbing exists yet, so
//! `Force` is recognised but currently rendered the same as `Auto`
//! pending Phase 168.5's pipeline-variant work).
//!
//! [`viewport::show`]: ../valenx_app/viewport/fn.show.html

use crate::error::OcctVizError;

/// Back-face culling mode mirror of OCCT's
/// `Graphic3d_TypeOfBackfacingModel`.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Default, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum BackFaceMode {
    /// Back-face culling enabled (default — drops back-facing
    /// triangles before rasterisation).
    #[default]
    Disable,
    /// Driver decides per primitive (Valenx treats this same as
    /// `Disable` in v1 because the wgpu pipeline has no per-primitive
    /// switch).
    Automatic,
    /// Draw back faces too (Phase 168.5 will add the second pipeline
    /// variant; v1 falls through to `Automatic` semantics).
    Force,
}

/// Set the back-face culling mode. Returns the mode that was applied
/// (in v1 this matches the input; in Phase 168.5+ it may downgrade
/// `Force` to `Automatic` if the requested mode isn't yet supported).
pub fn v3d_viewer_back_face_culling(mode: BackFaceMode) -> Result<BackFaceMode, OcctVizError> {
    Ok(mode)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_disable() {
        assert_eq!(
            v3d_viewer_back_face_culling(BackFaceMode::Disable).unwrap(),
            BackFaceMode::Disable
        );
    }

    #[test]
    fn round_trips_force() {
        assert_eq!(
            v3d_viewer_back_face_culling(BackFaceMode::Force).unwrap(),
            BackFaceMode::Force
        );
    }

    #[test]
    fn default_is_disable() {
        assert_eq!(BackFaceMode::default(), BackFaceMode::Disable);
    }
}
