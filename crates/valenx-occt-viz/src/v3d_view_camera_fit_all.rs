//! Phase 165 — `V3d_View::FitAll()` — auto-frame the camera around
//! the bounded scene volume so every visible object fits.
//!
//! ## What OCCT does
//!
//! `V3d_View::FitAll(margin)` walks every interactive object the
//! viewer knows about, computes a union AABB, then calls
//! `Camera::SetCenter(bbox_centre)` and `SetDistance(half_diag /
//! tan(fov_y/2) * (1 + margin))`. The `margin` parameter (default
//! 1% in OCCT, 15% in Valenx for visual breathing room) extends the
//! framing distance so the geometry doesn't crowd the edge of the
//! viewport.
//!
//! ## v1 status
//!
//! **Honest v1.** Delegates to
//! [`valenx_viz::OrbitCamera::frame_bounds`] which implements the
//! identical "centre target + back the camera off by `diag/2 /
//! tan(fov_y/2)` with a 15% margin" formula. The caller passes the
//! union AABB pre-computed; this crate doesn't enumerate scene
//! objects (that's the responsibility of `valenx-app`'s
//! `ViewportState`).

use valenx_viz::OrbitCamera;

use crate::error::OcctVizError;

/// Fit `camera` to the world-space AABB `[min, max]`.
///
/// # Errors
///
/// - [`OcctVizError::BadInput`] if any component is non-finite, or
///   `min[i] > max[i]` on any axis.
pub fn v3d_view_camera_fit_all(
    camera: &mut OrbitCamera,
    min: [f32; 3],
    max: [f32; 3],
) -> Result<(), OcctVizError> {
    for (i, v) in min.iter().chain(max.iter()).enumerate() {
        if !v.is_finite() {
            return Err(OcctVizError::bad_input(
                "aabb",
                format!("component {i} is non-finite"),
            ));
        }
    }
    for i in 0..3 {
        if min[i] > max[i] {
            return Err(OcctVizError::bad_input(
                "aabb",
                format!("min[{i}]={} > max[{i}]={}", min[i], max[i]),
            ));
        }
    }
    camera.frame_bounds(min, max);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_inverted_aabb() {
        let mut cam = OrbitCamera::default();
        let err = v3d_view_camera_fit_all(&mut cam, [1.0; 3], [0.0; 3]).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_non_finite() {
        let mut cam = OrbitCamera::default();
        let err = v3d_view_camera_fit_all(&mut cam, [0.0; 3], [f32::INFINITY, 1.0, 1.0])
            .unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn centers_target_on_aabb_centre() {
        let mut cam = OrbitCamera::default();
        v3d_view_camera_fit_all(&mut cam, [0.0, 0.0, 0.0], [10.0, 10.0, 10.0]).unwrap();
        assert!((cam.target.x - 5.0).abs() < 1e-5);
        assert!((cam.target.y - 5.0).abs() < 1e-5);
        assert!((cam.target.z - 5.0).abs() < 1e-5);
    }
}
