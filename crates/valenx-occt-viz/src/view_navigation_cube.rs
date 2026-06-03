//! Phase 197 — `AIS_ViewCube` — 3D nav-cube widget in the corner
//! (click face → align view).
//!
//! ## What OCCT does
//!
//! `AIS_ViewCube` draws a small 3D cube in a viewport corner (typically
//! upper-right) with face labels (`Front`/`Back`/`Top`/`Bottom`/
//! `Left`/`Right`) and edge / corner hot-zones. Clicking a face snaps
//! the main view to that orthographic direction; clicking an edge
//! gives a two-axis isometric; clicking a corner gives the full
//! `Iso` view. Animated transitions via [`crate::view_animation_camera_path()`].
//!
//! ## v1 status
//!
//! **Honest v1.** The cube's geometry + rendering is the caller's
//! responsibility (egui's `Painter` can draw the small 2D projection
//! of the cube directly, no wgpu pass needed for the widget itself).
//! This op maps a clicked face name to the corresponding
//! [`crate::v3d_view_camera_axo_axonometric::AxoView`] so the caller's
//! click handler can call the existing axonometric snap to drive the
//! view change. The 7-region mapping (6 faces + iso for any corner)
//! is baked into the [`NavCubeRegion::to_axo_view()`] helper.

use crate::error::OcctVizError;
use crate::v3d_view_camera_axo_axonometric::AxoView;

/// One of the regions on the navigation cube the user can click.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum NavCubeRegion {
    /// Front face (looks down +Y).
    Front,
    /// Back face (looks down −Y).
    Back,
    /// Top face (looks down +Z).
    Top,
    /// Bottom face (looks down −Z).
    Bottom,
    /// Right face (looks down +X).
    Right,
    /// Left face (looks down −X).
    Left,
    /// Any of the 8 corners → isometric view.
    Corner,
    /// Any of the 12 edges → snap to nearest axonometric
    /// (Phase 197.5 will distinguish edges per-axis; v1 collapses
    /// to Iso).
    Edge,
}

impl NavCubeRegion {
    /// Map this region to the corresponding [`AxoView`] for
    /// [`crate::v3d_view_camera_axo_axonometric()`].
    pub fn to_axo_view(self) -> AxoView {
        match self {
            NavCubeRegion::Front => AxoView::Front,
            NavCubeRegion::Back => AxoView::Back,
            NavCubeRegion::Top => AxoView::Top,
            NavCubeRegion::Bottom => AxoView::Bottom,
            NavCubeRegion::Right => AxoView::Right,
            NavCubeRegion::Left => AxoView::Left,
            NavCubeRegion::Corner | NavCubeRegion::Edge => AxoView::Iso,
        }
    }
}

/// Resolve a clicked nav-cube region to the [`AxoView`] the main view
/// should snap to.
///
/// This op cannot fail under valid input — wraps the enum mapping
/// for API consistency.
pub fn view_navigation_cube(region: NavCubeRegion) -> Result<AxoView, OcctVizError> {
    Ok(region.to_axo_view())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn front_face_maps_to_front_view() {
        assert_eq!(
            view_navigation_cube(NavCubeRegion::Front).unwrap(),
            AxoView::Front
        );
    }

    #[test]
    fn top_face_maps_to_top_view() {
        assert_eq!(
            view_navigation_cube(NavCubeRegion::Top).unwrap(),
            AxoView::Top
        );
    }

    #[test]
    fn corner_and_edge_both_map_to_iso() {
        assert_eq!(
            view_navigation_cube(NavCubeRegion::Corner).unwrap(),
            AxoView::Iso
        );
        assert_eq!(
            view_navigation_cube(NavCubeRegion::Edge).unwrap(),
            AxoView::Iso
        );
    }

    #[test]
    fn all_six_faces_unique() {
        let faces = [
            NavCubeRegion::Front,
            NavCubeRegion::Back,
            NavCubeRegion::Top,
            NavCubeRegion::Bottom,
            NavCubeRegion::Right,
            NavCubeRegion::Left,
        ];
        let mapped: Vec<_> = faces.iter().map(|r| r.to_axo_view()).collect();
        // 6 unique outputs.
        for (i, a) in mapped.iter().enumerate() {
            for (j, b) in mapped.iter().enumerate() {
                if i != j {
                    assert_ne!(a, b, "regions {i} and {j} map to the same view");
                }
            }
        }
    }
}
