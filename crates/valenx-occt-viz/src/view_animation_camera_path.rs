//! Phase 189 — `AIS_Animation_Camera::SetView()` — animate the camera
//! along a parametric path.
//!
//! ## What OCCT does
//!
//! `AIS_Animation_Camera::SetView(Start, End)` registers a tween
//! between two `gp_Camera` poses. On each `Update(time)` call, OCCT
//! linearly interpolates the eye position, the target position, and
//! the up vector — but interpolates the field-of-view + projection
//! type with step-discontinuity at `t=1` (so the destination FoV is
//! used the whole time, avoiding mid-tween FoV warble). The view
//! matrix is recomputed every frame from the interpolated camera.
//!
//! ## v1 status
//!
//! **Honest v1.** Linearly interpolates azimuth (along the shorter
//! arc — wraps around the 0° / 360° boundary correctly), elevation,
//! distance, and the target position. FoV uses step-at-1 to mirror
//! OCCT. Caller passes `t ∈ [0, 1]` (clamped) and gets a new
//! [`valenx_viz::OrbitCamera`] that's safe to swap into the viewport
//! immediately.

use nalgebra::Point3;
use valenx_viz::OrbitCamera;

use crate::error::OcctVizError;

/// Compute the interpolated camera at parameter `t ∈ [0, 1]` along
/// the path `start → end`.
///
/// `t = 0.0` returns a camera equal to `start`. `t = 1.0` returns
/// a camera equal to `end`. Intermediate values lerp in a way that
/// avoids the gimbal flip at the 0°/360° azimuth boundary.
///
/// # Errors
///
/// - [`OcctVizError::BadInput`] if `t` is not finite.
pub fn view_animation_camera_path(
    start: &OrbitCamera,
    end: &OrbitCamera,
    t: f32,
) -> Result<OrbitCamera, OcctVizError> {
    if !t.is_finite() {
        return Err(OcctVizError::bad_input("t", "must be finite"));
    }
    let t = t.clamp(0.0, 1.0);

    // Shortest-arc azimuth lerp: wrap delta to [-180, 180].
    let mut d_az = end.azimuth_deg - start.azimuth_deg;
    while d_az > 180.0 {
        d_az -= 360.0;
    }
    while d_az < -180.0 {
        d_az += 360.0;
    }
    let az = start.azimuth_deg + d_az * t;

    let el = start.elevation_deg + (end.elevation_deg - start.elevation_deg) * t;
    let distance = start.distance + (end.distance - start.distance) * t;
    let target = Point3::new(
        start.target.x + (end.target.x - start.target.x) * t,
        start.target.y + (end.target.y - start.target.y) * t,
        start.target.z + (end.target.z - start.target.z) * t,
    );
    // FoV: step at t=1 (OCCT semantics).
    let fov = if t >= 1.0 {
        end.fov_y_deg
    } else {
        start.fov_y_deg
    };

    // Projection mode is discrete — switch at t=1 (same step
    // semantics as FoV).
    let projection_mode = if t >= 1.0 {
        end.projection_mode
    } else {
        start.projection_mode
    };

    Ok(OrbitCamera {
        target,
        distance,
        azimuth_deg: az,
        elevation_deg: el,
        fov_y_deg: fov,
        near: start.near.min(end.near),
        far: start.far.max(end.far),
        projection_mode,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cam(az: f32, el: f32, dist: f32) -> OrbitCamera {
        OrbitCamera {
            azimuth_deg: az,
            elevation_deg: el,
            distance: dist,
            ..Default::default()
        }
    }

    #[test]
    fn rejects_nan_t() {
        let a = cam(0.0, 0.0, 5.0);
        let b = cam(90.0, 0.0, 5.0);
        let err = view_animation_camera_path(&a, &b, f32::NAN).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn t_zero_returns_start() {
        let a = cam(0.0, 30.0, 5.0);
        let b = cam(90.0, 60.0, 10.0);
        let c = view_animation_camera_path(&a, &b, 0.0).unwrap();
        assert!((c.azimuth_deg - a.azimuth_deg).abs() < 1e-4);
        assert!((c.elevation_deg - a.elevation_deg).abs() < 1e-4);
        assert!((c.distance - a.distance).abs() < 1e-4);
    }

    #[test]
    fn t_one_returns_end() {
        let a = cam(0.0, 30.0, 5.0);
        let b = cam(90.0, 60.0, 10.0);
        let c = view_animation_camera_path(&a, &b, 1.0).unwrap();
        assert!((c.azimuth_deg - b.azimuth_deg).abs() < 1e-4);
        assert!((c.elevation_deg - b.elevation_deg).abs() < 1e-4);
        assert!((c.distance - b.distance).abs() < 1e-4);
    }

    #[test]
    fn midpoint_is_midway() {
        let a = cam(0.0, 30.0, 5.0);
        let b = cam(90.0, 60.0, 10.0);
        let c = view_animation_camera_path(&a, &b, 0.5).unwrap();
        assert!((c.azimuth_deg - 45.0).abs() < 1e-4);
        assert!((c.elevation_deg - 45.0).abs() < 1e-4);
        assert!((c.distance - 7.5).abs() < 1e-4);
    }

    #[test]
    fn shortest_arc_wraps_through_zero() {
        // 350° → 10° should lerp through 0°, not the long way around.
        let a = cam(350.0, 0.0, 5.0);
        let b = cam(10.0, 0.0, 5.0);
        let c = view_animation_camera_path(&a, &b, 0.5).unwrap();
        // Midpoint should be ~360° (≡ 0°), not 180°.
        let az_norm = ((c.azimuth_deg % 360.0) + 360.0) % 360.0;
        assert!(
            !(5.0..=355.0).contains(&az_norm),
            "midpoint az={az_norm} should be near 0/360"
        );
    }

    #[test]
    fn clamps_t_out_of_range() {
        let a = cam(0.0, 0.0, 5.0);
        let b = cam(90.0, 0.0, 10.0);
        let c = view_animation_camera_path(&a, &b, 5.0).unwrap();
        // Clamped to 1.0 → matches end.
        assert!((c.azimuth_deg - b.azimuth_deg).abs() < 1e-4);
    }
}
