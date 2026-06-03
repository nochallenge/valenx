//! Phase 162 — `V3d_View::Pan()` — translate camera + target in the
//! screen-aligned plane.
//!
//! ## What OCCT does
//!
//! `V3d_View::Pan(dx_px, dy_px)` shifts both the eye and the look-at
//! target by the same vector, expressed in screen-pixel deltas and
//! projected into world space using the camera's current view matrix.
//! The effect is "drag the scene around without rotating it" — eye and
//! target move together so the orbit centre stays put relative to the
//! geometry. The screen-pixel → world scale factor depends on the
//! camera distance and the viewport size; OCCT computes it from
//! `Convert(dx, dy)` (a viewport projection helper).
//!
//! ## v1 status
//!
//! **Honest v1.** Computes the camera-local screen-right and
//! screen-up axes from the [`valenx_viz::OrbitCamera::view_matrix`]
//! (transpose of the rotation block) and applies the world-space
//! translation to `target`. Because `OrbitCamera::eye` is *derived*
//! from `target + spherical(distance, azimuth, elevation)`, moving the
//! target shifts the eye by the same vector — exactly the OCCT
//! semantics. The screen→world scale uses `distance * fov_y * dy_px /
//! viewport_h` so panning is distance-aware (further-away scenes need
//! bigger world deltas per pixel).

use nalgebra::Vector3;
use valenx_viz::OrbitCamera;

use crate::error::OcctVizError;

/// Pan `camera` by `(dx_px, dy_px)` screen-pixel deltas. `viewport_h`
/// is the viewport height in pixels (must be > 0).
///
/// Positive `dx_px` drags the scene to the right (target moves left);
/// positive `dy_px` drags the scene up (target moves down). Sign
/// matches OCCT and the standard egui drag-pan convention.
///
/// # Errors
///
/// - [`OcctVizError::BadInput`] if any input is non-finite, or
///   `viewport_h <= 0`.
pub fn v3d_view_camera_pan(
    camera: &mut OrbitCamera,
    dx_px: f32,
    dy_px: f32,
    viewport_h: f32,
) -> Result<(), OcctVizError> {
    if !dx_px.is_finite() {
        return Err(OcctVizError::bad_input("dx_px", "must be finite"));
    }
    if !dy_px.is_finite() {
        return Err(OcctVizError::bad_input("dy_px", "must be finite"));
    }
    if !viewport_h.is_finite() || viewport_h <= 0.0 {
        return Err(OcctVizError::bad_input("viewport_h", "must be > 0"));
    }

    // World-units-per-pixel at the target plane: half-height of the
    // frustum at `distance` is `distance * tan(fov_y/2)`, full height
    // is twice that. Per-pixel scale = full / viewport_h.
    let half = (camera.fov_y_deg * 0.5).to_radians().tan();
    let world_per_pixel = 2.0 * camera.distance * half / viewport_h;

    // Derive screen-right + screen-up in world coords from the view
    // matrix's rotation block (rows are world axes mapped into view;
    // transpose gives view axes in world).
    let view = camera.view_matrix();
    let right = Vector3::new(view[(0, 0)], view[(0, 1)], view[(0, 2)]);
    let up = Vector3::new(view[(1, 0)], view[(1, 1)], view[(1, 2)]);

    // Drag-the-scene convention: positive dx moves target leftward
    // (so scene appears to move right). Same for dy.
    let translation = -right * (dx_px * world_per_pixel) + up * (dy_px * world_per_pixel);
    camera.target += translation;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_zero_viewport_height() {
        let mut cam = OrbitCamera::default();
        let err = v3d_view_camera_pan(&mut cam, 1.0, 0.0, 0.0).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_nan_dx() {
        let mut cam = OrbitCamera::default();
        let err = v3d_view_camera_pan(&mut cam, f32::NAN, 0.0, 600.0).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn moves_target() {
        let mut cam = OrbitCamera::default();
        let before = cam.target;
        v3d_view_camera_pan(&mut cam, 10.0, 0.0, 600.0).unwrap();
        assert!(
            (cam.target - before).norm() > 1e-6,
            "pan should move target"
        );
    }

    #[test]
    fn pan_distance_scales_with_camera_distance() {
        let mut near = OrbitCamera {
            distance: 1.0,
            ..Default::default()
        };
        let mut far = OrbitCamera {
            distance: 10.0,
            ..Default::default()
        };
        let before_near = near.target;
        let before_far = far.target;
        v3d_view_camera_pan(&mut near, 100.0, 0.0, 600.0).unwrap();
        v3d_view_camera_pan(&mut far, 100.0, 0.0, 600.0).unwrap();
        let near_delta = (near.target - before_near).norm();
        let far_delta = (far.target - before_far).norm();
        assert!(far_delta > near_delta * 9.0, "far camera should pan ~10x");
    }
}
