//! Phase 173 — `AIS_InteractiveContext::SelectRectangle()` —
//! rubber-band box selection of every object whose screen-projected
//! bounding box falls inside the rectangle.
//!
//! ## What OCCT does
//!
//! After the user releases a drag-rectangle (the visual painted by
//! [`crate::v3d_viewer_xor_drag()`]), OCCT projects every selectable
//! object's AABB into screen space and tests for inclusion in the
//! drag rectangle. With the "inclusive" mode, partial overlap counts;
//! with "exclusive" mode, the AABB must lie *entirely* inside.
//! Selection promotion follows the same path as
//! [`crate::ais_select_single()`].
//!
//! ## v1 status — real screen-AABB box selection
//!
//! With the picking substrate from Phase 171 in place, this is a real
//! selector: every object's registered geometry vertices are
//! projected to screen pixels with the installed camera, their screen
//! AABB is compared to the drag rectangle per `mode`
//! ([`BoxSelectionMode::Inclusive`] = overlap, `Exclusive` = fully
//! enclosed), and the matching objects are promoted to `Selected`.
//!
//! Objects with no registered geometry, or behind the camera, are
//! skipped. With no camera installed nothing is selected.

use crate::ais_interactive_context::{InteractiveContext, ObjectState};
use crate::error::OcctVizError;
use crate::v3d_viewer_xor_drag::v3d_viewer_xor_drag;

/// Selection mode mirror of OCCT's `Aspect_TypeOfHighlightMethod`.
#[derive(
    Copy, Clone, Debug, Eq, PartialEq, Hash, Default, serde::Serialize, serde::Deserialize,
)]
#[serde(rename_all = "lowercase")]
pub enum BoxSelectionMode {
    /// Object's screen-AABB must overlap the rectangle by ≥ 1 pixel
    /// (default — least surprising for users coming from CAD apps).
    #[default]
    Inclusive,
    /// Object's screen-AABB must lie entirely inside the rectangle
    /// (matches FreeCAD's "fully enclosed" toggle).
    Exclusive,
}

/// Select every object in `ctx` whose screen-projected AABB satisfies
/// `mode` against the rectangle defined by corners `p1` and `p2`.
/// Returns the list of newly-selected IDs.
///
/// # Errors
///
/// [`OcctVizError::BadInput`] if either corner contains non-finite
/// coordinates (delegated to [`crate::v3d_viewer_xor_drag()`]).
pub fn ais_select_box(
    ctx: &mut InteractiveContext,
    p1: [f32; 2],
    p2: [f32; 2],
    mode: BoxSelectionMode,
) -> Result<Vec<usize>, OcctVizError> {
    // Normalise the rectangle (handles any drag direction).
    let rect = v3d_viewer_xor_drag(p1, p2)?;
    let (rx0, ry0) = (rect.min[0], rect.min[1]);
    let (rx1, ry1) = (rect.max[0], rect.max[1]);

    let mut selected: Vec<usize> = Vec::new();
    for id in ctx.geometry_ids() {
        if ctx.state(id) == Some(ObjectState::Hidden) {
            continue;
        }
        let Some(geom) = ctx.geometry(id) else {
            continue;
        };
        // Project every geometry vertex to screen and accumulate the
        // screen-space AABB.
        let mut sx0 = f32::INFINITY;
        let mut sy0 = f32::INFINITY;
        let mut sx1 = f32::NEG_INFINITY;
        let mut sy1 = f32::NEG_INFINITY;
        let mut any = false;
        for v in geometry_vertices(geom) {
            if let Some((sx, sy)) = ctx.project(v) {
                sx0 = sx0.min(sx);
                sy0 = sy0.min(sy);
                sx1 = sx1.max(sx);
                sy1 = sy1.max(sy);
                any = true;
            }
        }
        if !any {
            continue;
        }
        let inside = match mode {
            BoxSelectionMode::Inclusive => {
                // Screen AABBs overlap.
                sx0 <= rx1 && sx1 >= rx0 && sy0 <= ry1 && sy1 >= ry0
            }
            BoxSelectionMode::Exclusive => {
                // Object AABB fully inside the rectangle.
                sx0 >= rx0 && sx1 <= rx1 && sy0 >= ry0 && sy1 <= ry1
            }
        };
        if inside {
            ctx.set_state(id, ObjectState::Selected);
            selected.push(id);
        }
    }
    selected.sort();
    Ok(selected)
}

/// All world-space vertices of a [`Pickable`] — used to compute its
/// screen AABB.
fn geometry_vertices(
    geom: &crate::ais_interactive_context::Pickable,
) -> Vec<nalgebra::Vector3<f64>> {
    use crate::ais_interactive_context::Pickable;
    let v3 = |p: &[f64; 3]| nalgebra::Vector3::new(p[0], p[1], p[2]);
    match geom {
        Pickable::Mesh { triangles } | Pickable::TaggedMesh { triangles, .. } => {
            triangles.iter().map(v3).collect()
        }
        Pickable::Polyline { points } => points.iter().map(v3).collect(),
        Pickable::Point { position } => vec![v3(position)],
        Pickable::TaggedEdges { edges } => edges
            .iter()
            .flat_map(|(_, pts)| pts.iter().map(v3))
            .collect(),
        Pickable::TaggedVertices { vertices } => vertices.iter().map(|(_, p)| v3(p)).collect(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ais_interactive_context::{ais_interactive_context, Pickable};
    use nalgebra::{Matrix4, Point3, Vector3};

    fn test_view() -> Matrix4<f64> {
        let view = Matrix4::look_at_rh(
            &Point3::new(0.0, 0.0, 10.0),
            &Point3::origin(),
            &Vector3::y(),
        );
        let proj = Matrix4::new_perspective(1.0, 60f64.to_radians(), 0.1, 100.0);
        proj * view
    }

    fn axis_triangle() -> Pickable {
        Pickable::Mesh {
            triangles: vec![[-0.5, -0.5, 0.0], [0.5, -0.5, 0.0], [0.0, 0.5, 0.0]],
        }
    }

    #[test]
    fn rejects_nan_corner() {
        let mut ctx = ais_interactive_context().unwrap();
        let err = ais_select_box(
            &mut ctx,
            [0.0, 0.0],
            [f32::NAN, 10.0],
            BoxSelectionMode::Inclusive,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn default_mode_is_inclusive() {
        assert_eq!(BoxSelectionMode::default(), BoxSelectionMode::Inclusive);
    }

    #[test]
    fn no_camera_selects_nothing() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.display_geometry(axis_triangle());
        let sel = ais_select_box(
            &mut ctx,
            [0.0, 0.0],
            [100.0, 100.0],
            BoxSelectionMode::Inclusive,
        )
        .unwrap();
        assert!(sel.is_empty());
    }

    #[test]
    fn full_screen_rectangle_selects_the_centred_object() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        let id = ctx.display_geometry(axis_triangle());
        let sel = ais_select_box(
            &mut ctx,
            [0.0, 0.0],
            [100.0, 100.0],
            BoxSelectionMode::Inclusive,
        )
        .unwrap();
        assert_eq!(sel, vec![id]);
        assert_eq!(ctx.state(id), Some(ObjectState::Selected));
    }

    #[test]
    fn corner_rectangle_misses_the_centred_object() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        ctx.display_geometry(axis_triangle());
        // A tiny rectangle in the top-left corner — the centred
        // triangle projects to the middle of the screen.
        let sel = ais_select_box(
            &mut ctx,
            [0.0, 0.0],
            [5.0, 5.0],
            BoxSelectionMode::Inclusive,
        )
        .unwrap();
        assert!(sel.is_empty());
    }

    #[test]
    fn exclusive_mode_needs_full_enclosure() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        ctx.display_geometry(axis_triangle());
        // Inclusive: a half-screen box that overlaps the centre.
        let inc = ais_select_box(
            &mut ctx,
            [40.0, 40.0],
            [100.0, 100.0],
            BoxSelectionMode::Inclusive,
        )
        .unwrap();
        // Reset for the exclusive run. The centred triangle projects
        // to a screen AABB of roughly x∈[45.7,54.3], y∈[45.7,54.3]; an
        // exclusive box starting at (48,48) clips off its lower-left
        // corner, so it is not fully enclosed.
        let mut ctx2 = ais_interactive_context().unwrap();
        ctx2.set_view(test_view(), 100.0, 100.0);
        ctx2.display_geometry(axis_triangle());
        let exc = ais_select_box(
            &mut ctx2,
            [48.0, 48.0],
            [60.0, 60.0],
            BoxSelectionMode::Exclusive,
        )
        .unwrap();
        // The object straddles the exclusive box edge → not fully
        // enclosed → not selected; but the larger inclusive box
        // overlaps it → selected.
        assert_eq!(inc.len(), 1, "inclusive overlap should select");
        assert!(exc.is_empty(), "exclusive needs full enclosure");
    }
}
