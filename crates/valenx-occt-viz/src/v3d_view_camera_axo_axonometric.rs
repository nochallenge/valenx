//! Phase 166 — `V3d_View::SetProj()` for standard axonometric views.
//!
//! ## What OCCT does
//!
//! `V3d_View::SetProj(V3d_TypeOfOrientation)` snaps the camera direction
//! vector to one of the canonical orientations (`V3d_Xpos`, `V3d_Ypos`,
//! `V3d_Zpos`, `V3d_XposYposZpos` for isometric). OCCT then resets the
//! roll to zero and re-runs `Camera::SetUp` to align world-up with the
//! viewport's vertical axis (CAD-standard with no tilt).
//!
//! ## v1 status
//!
//! **Honest v1.** Delegates to [`valenx_viz::OrbitCamera::set_view`]
//! which handles all six axis-aligned views plus `Iso`. The mapping
//! from OCCT names → Valenx [`ViewDirection`] is direct: Xpos=Right,
//! Ypos=Top, Zpos=Front (OCCT uses Z-up by default; Valenx uses Y-up
//! to match Fusion 360 convention — the axis swap is baked into the
//! [`AxoView`] mapping below).
//!
//! [`ViewDirection`]: ../valenx_viz/camera/enum.ViewDirection.html

use valenx_viz::{OrbitCamera, ViewDirection};

use crate::error::OcctVizError;

/// One of the standard axonometric / orthographic snaps. Mirrors
/// OCCT's `V3d_TypeOfOrientation` enum with Y-up axis convention to
/// match Valenx's [`OrbitCamera`] world.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum AxoView {
    /// Look down +Z axis at the XY plane (top-down).
    Top,
    /// Look up −Z axis at the XY plane (bottom-up).
    Bottom,
    /// Look down +X axis at the YZ plane (right side).
    Right,
    /// Look up −X axis at the YZ plane (left side).
    Left,
    /// Look down +Y axis at the XZ plane (front).
    Front,
    /// Look up −Y axis at the XZ plane (back).
    Back,
    /// Standard isometric (45° azimuth, ~35.264° elevation —
    /// arctan(1/√2), the cube-corner viewing angle).
    Iso,
}

impl AxoView {
    fn to_view_direction(self) -> ViewDirection {
        match self {
            AxoView::Top => ViewDirection::Top,
            AxoView::Bottom => ViewDirection::Bottom,
            AxoView::Right => ViewDirection::Right,
            AxoView::Left => ViewDirection::Left,
            AxoView::Front => ViewDirection::Front,
            AxoView::Back => ViewDirection::Back,
            AxoView::Iso => ViewDirection::Iso,
        }
    }
}

/// Snap `camera` to the given axonometric / orthographic view.
///
/// This op cannot fail under valid input — but takes the standard
/// `Result` for API consistency so the caller's match-arm pattern
/// matches the rest of this crate.
pub fn v3d_view_camera_axo_axonometric(
    camera: &mut OrbitCamera,
    view: AxoView,
) -> Result<(), OcctVizError> {
    camera.set_view(view.to_view_direction());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn top_view_pins_elevation_to_90() {
        let mut cam = OrbitCamera::default();
        v3d_view_camera_axo_axonometric(&mut cam, AxoView::Top).unwrap();
        assert!((cam.elevation_deg - 90.0).abs() < 1e-4);
    }

    #[test]
    fn front_view_pins_azimuth_to_0() {
        let mut cam = OrbitCamera::default();
        v3d_view_camera_axo_axonometric(&mut cam, AxoView::Front).unwrap();
        assert!(cam.azimuth_deg.abs() < 1e-4);
        assert!(cam.elevation_deg.abs() < 1e-4);
    }

    #[test]
    fn iso_view_uses_arctan_half_sqrt2() {
        let mut cam = OrbitCamera::default();
        v3d_view_camera_axo_axonometric(&mut cam, AxoView::Iso).unwrap();
        assert!((cam.azimuth_deg - 45.0).abs() < 1e-3);
        // arctan(1/sqrt(2)) ≈ 35.264°
        assert!((cam.elevation_deg - 35.264).abs() < 0.1);
    }
}
