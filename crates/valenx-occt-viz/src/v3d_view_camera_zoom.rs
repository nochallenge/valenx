//! Phase 163 — `V3d_View::SetZoom()` / `Zoom()` — dolly and framing
//! zoom.
//!
//! ## What OCCT does
//!
//! OCCT exposes three zoom flavours: `SetZoom(coef)` scales the
//! orthographic scale by `coef` (or, in perspective mode, divides the
//! camera distance by `coef`); `Zoom(dx, dy)` interprets a pixel delta
//! using a viewer-configured sensitivity; `FitAll` reframes around the
//! whole scene's bounding box (delegated to Phase 165). All three end
//! up calling `gp_Camera::SetDistance` for perspective cameras.
//!
//! ## v1 status
//!
//! **Honest v1.** Delegates to [`valenx_viz::OrbitCamera::zoom`]
//! which treats `frac` as a fractional dolly: 0.1 = "zoom in 10%",
//! −0.1 = "zoom out 10%". The clamp at distance ≥ 1e-4 mirrors OCCT's
//! anti-inversion guard (`SetDistance` floors at machine epsilon).

use valenx_viz::OrbitCamera;

use crate::error::OcctVizError;

/// Zoom `camera` by a fractional dolly amount.
///
/// Positive `frac` zooms in (camera approaches target); negative
/// zooms out. `0.0` is a no-op. Clamped internally so the camera
/// can't invert through the target.
///
/// # Errors
///
/// - [`OcctVizError::BadInput`] if `frac` is non-finite or ≥ 1.0
///   (which would attempt to drive `distance` to zero or negative).
pub fn v3d_view_camera_zoom(camera: &mut OrbitCamera, frac: f32) -> Result<(), OcctVizError> {
    if !frac.is_finite() {
        return Err(OcctVizError::bad_input("frac", "must be finite"));
    }
    if frac >= 1.0 {
        return Err(OcctVizError::bad_input(
            "frac",
            format!("must be < 1.0 (got {frac}); would invert through target"),
        ));
    }
    camera.zoom(frac);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_nan() {
        let mut cam = OrbitCamera::default();
        let err = v3d_view_camera_zoom(&mut cam, f32::NAN).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_at_one() {
        let mut cam = OrbitCamera::default();
        let err = v3d_view_camera_zoom(&mut cam, 1.0).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn zoom_in_reduces_distance() {
        let mut cam = OrbitCamera::default();
        let before = cam.distance;
        v3d_view_camera_zoom(&mut cam, 0.5).unwrap();
        assert!(cam.distance < before, "zoom in should reduce distance");
    }

    #[test]
    fn zoom_out_increases_distance() {
        let mut cam = OrbitCamera::default();
        let before = cam.distance;
        v3d_view_camera_zoom(&mut cam, -0.5).unwrap();
        assert!(cam.distance > before, "zoom out should increase distance");
    }
}
