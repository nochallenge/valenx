//! Phase 164 — `V3d_View::SetType(Camera_Perspective | Camera_Orthographic)`
//! — switch projection mode.
//!
//! ## What OCCT does
//!
//! `gp_Camera::SetProjectionType(Camera::Projection_Perspective)` keeps
//! the field-of-view + near/far frustum and computes clip-space `z` via
//! `(z - near) / (far - near)`. `SetProjectionType(Projection_Orthographic)`
//! switches to a parallel-projection frustum sized by
//! `Scale = 2 * distance * tan(fov_y/2)` so that the visible footprint
//! at the target plane matches what perspective showed pre-toggle (the
//! "swap projection without zooming" UX users expect).
//!
//! ## v1 status
//!
//! **Honest implementation** (Phase 164.5). `valenx_viz::OrbitCamera`
//! now carries a [`ProjectionMode`] field, and its `projection_matrix`
//! builder branches on it — perspective via `Matrix4::new_perspective`,
//! orthographic via `Matrix4::new_orthographic` with a frustum
//! half-height of `distance * tan(fov_y / 2)` (the perspective
//! footprint at the target plane, so the model keeps its on-screen
//! size across the toggle). This function is the seam the UI toggle
//! calls; it simply writes the requested mode onto the camera.

pub use valenx_viz::ProjectionMode;
use valenx_viz::OrbitCamera;

use crate::error::OcctVizError;

/// Switch `camera` to the requested projection mode.
///
/// The change takes effect on the next frame — the renderer reads
/// `camera.projection_mode` when it builds the MVP matrix.
///
/// # Errors
///
/// Infallible in v1 — the `Result` is kept so a future
/// camera-validation step (e.g. rejecting an orthographic toggle on a
/// zero-distance camera) can surface without an API break.
///
/// # Example
///
/// ```
/// use valenx_occt_viz::v3d_view_camera_perspective_toggle::{
///     v3d_view_camera_perspective_toggle, ProjectionMode,
/// };
/// let mut cam = valenx_viz::OrbitCamera::default();
/// v3d_view_camera_perspective_toggle(&mut cam, ProjectionMode::Orthographic).unwrap();
/// assert_eq!(cam.projection_mode, ProjectionMode::Orthographic);
/// ```
pub fn v3d_view_camera_perspective_toggle(
    camera: &mut OrbitCamera,
    mode: ProjectionMode,
) -> Result<(), OcctVizError> {
    camera.projection_mode = mode;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn toggle_to_orthographic_sets_mode() {
        let mut cam = OrbitCamera::default();
        assert_eq!(cam.projection_mode, ProjectionMode::Perspective);
        v3d_view_camera_perspective_toggle(&mut cam, ProjectionMode::Orthographic).unwrap();
        assert_eq!(cam.projection_mode, ProjectionMode::Orthographic);
    }

    #[test]
    fn toggle_back_to_perspective() {
        let mut cam = OrbitCamera::default();
        v3d_view_camera_perspective_toggle(&mut cam, ProjectionMode::Orthographic).unwrap();
        v3d_view_camera_perspective_toggle(&mut cam, ProjectionMode::Perspective).unwrap();
        assert_eq!(cam.projection_mode, ProjectionMode::Perspective);
    }

    #[test]
    fn toggle_changes_the_projection_matrix() {
        // After the toggle the projection matrix must actually differ
        // — proves the field is read by the matrix builder.
        let mut cam = OrbitCamera::default();
        let persp = cam.projection_matrix(1.6);
        v3d_view_camera_perspective_toggle(&mut cam, ProjectionMode::Orthographic).unwrap();
        let ortho = cam.projection_matrix(1.6);
        // Orthographic m[3,3] == 1; perspective m[3,3] == 0.
        assert!((ortho[(3, 3)] - 1.0).abs() < 1e-6);
        assert!(persp[(3, 3)].abs() < 1e-6);
    }
}
