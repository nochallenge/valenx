//! Phase 161 — `V3d_View::Rotate()` family — orbit the camera around
//! its target point.
//!
//! ## What OCCT does
//!
//! `V3d_View::Rotate(dAz, dEl, dRoll)` increments the view's azimuth,
//! elevation, and roll angles (radians) about the camera's pivot point
//! (defaults to the bounded scene centre). `V3d_View::Turn` is the
//! same op in screen-pixel deltas pre-scaled by a sensitivity factor.
//! Both wrap `gp_Camera::SetEye` internally to recompute the eye
//! position from the spherical coordinates.
//!
//! ## v1 status
//!
//! **Honest v1.** Delegates to [`valenx_viz::OrbitCamera::orbit`] which
//! implements the same spherical accumulation with elevation clamped
//! to ±89.9° (OCCT clamps the same way — looking straight up/down
//! through the pole causes the azimuth math to gimbal-lock). Roll is
//! not implemented (Valenx's CAD-style orbit camera has no roll
//! degree-of-freedom; OCCT exposes it but the Valenx viewport
//! conventionally locks the world-up vector to +Y).

use valenx_viz::OrbitCamera;

use crate::error::OcctVizError;

/// Orbit `camera` by `(dx_deg, dy_deg)` screen-pixel-scaled angle
/// increments. Positive `dx_deg` rotates clockwise looking down; positive
/// `dy_deg` raises the elevation (camera goes up). Both are in *degrees*
/// — convert from screen pixels at the call site using whatever
/// sensitivity factor the UI prefers (the egui viewport uses 0.5°/px).
///
/// # Errors
///
/// - [`OcctVizError::BadInput`] if either delta is non-finite.
pub fn v3d_view_camera_orbit(
    camera: &mut OrbitCamera,
    dx_deg: f32,
    dy_deg: f32,
) -> Result<(), OcctVizError> {
    if !dx_deg.is_finite() {
        return Err(OcctVizError::bad_input("dx_deg", "must be finite"));
    }
    if !dy_deg.is_finite() {
        return Err(OcctVizError::bad_input("dy_deg", "must be finite"));
    }
    camera.orbit(dx_deg, dy_deg);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_nan_dx() {
        let mut cam = OrbitCamera::default();
        let err = v3d_view_camera_orbit(&mut cam, f32::NAN, 0.0).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn orbits_azimuth_by_dx() {
        let mut cam = OrbitCamera::default();
        let before = cam.azimuth_deg;
        v3d_view_camera_orbit(&mut cam, 10.0, 0.0).unwrap();
        let after = cam.azimuth_deg;
        // OrbitCamera::orbit applies the raw degree increment.
        assert!((after - before - 10.0).abs() < 1e-4);
    }

    #[test]
    fn elevation_clamps_at_pole() {
        let mut cam = OrbitCamera::default();
        v3d_view_camera_orbit(&mut cam, 0.0, 1000.0).unwrap();
        assert!(cam.elevation_deg <= 89.9);
    }
}
