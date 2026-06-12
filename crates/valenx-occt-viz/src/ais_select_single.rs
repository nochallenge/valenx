//! Phase 172 — `AIS_InteractiveContext::Select()` / `MoveTo` then
//! single-pick — pick the topmost object under the cursor.
//!
//! ## What OCCT does
//!
//! On click, OCCT performs a per-pixel ray-cast (with the camera's
//! current view matrix) through every selectable object's `SelectMgr`-
//! registered primitive list (triangles for face owners, polyline
//! segments for edge owners, points for vertex owners). The first
//! intersection along the ray wins; the owner is then promoted to
//! `Selected`. ID-based hit tests (the GPU read-back picking pass) are
//! an alternative path enabled with `SelectionMode > 0`.
//!
//! ## v1 status — real CPU ray-cast pick
//!
//! This implements OCCT's default `SelectionMode 0` path: a CPU
//! ray-cast. [`InteractiveContext::ray_at`] unprojects the cursor
//! pixel into a world ray (the renderer keeps the camera current via
//! [`InteractiveContext::set_view`]); the ray is tested against every
//! object's registered [`crate::ais_interactive_context::Pickable`]
//! geometry and the nearest hit's owner is promoted to `Selected`.
//!
//! Objects registered geometry-less (via plain `display()`) are not
//! pickable — only those registered with `display_geometry` are. If
//! no camera has been installed, or the ray misses everything, the
//! pick returns `None` and the selection set is unchanged.

use crate::ais_interactive_context::{InteractiveContext, ObjectState};
use crate::error::OcctVizError;

/// World-space proximity (in model units) for edge / vertex picks. A
/// triangle pick is exact; a polyline or point pick counts as a hit
/// when the ray passes within this distance.
pub const PICK_TOLERANCE: f64 = 0.05;

/// Pick the topmost object at screen-pixel `(x, y)` in `ctx` and
/// promote it to the `Selected` state. On hit, the returned
/// `Option<usize>` carries the picked object's ID.
///
/// The previous selection is **not** cleared — call
/// [`InteractiveContext::set_state`] to reset states for a
/// single-select-replaces-all UX, or keep it for additive selection.
///
/// # Errors
///
/// [`OcctVizError::BadInput`] if either coordinate is non-finite or
/// negative.
pub fn ais_select_single(
    ctx: &mut InteractiveContext,
    x: f32,
    y: f32,
) -> Result<Option<usize>, OcctVizError> {
    if !x.is_finite() || x < 0.0 {
        return Err(OcctVizError::bad_input(
            "x",
            format!("must be finite and >= 0 (got {x})"),
        ));
    }
    if !y.is_finite() || y < 0.0 {
        return Err(OcctVizError::bad_input(
            "y",
            format!("must be finite and >= 0 (got {y})"),
        ));
    }
    // Unproject + ray-cast. No camera installed ⇒ nothing to pick.
    let Some(ray) = ctx.ray_at(x, y) else {
        return Ok(None);
    };
    let Some(hit) = ctx.pick_nearest(&ray, PICK_TOLERANCE) else {
        return Ok(None);
    };
    ctx.set_state(hit.id, ObjectState::Selected);
    Ok(Some(hit.id))
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
            triangles: vec![[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 1.0, 0.0]],
        }
    }

    #[test]
    fn rejects_negative_x() {
        let mut ctx = ais_interactive_context().unwrap();
        let err = ais_select_single(&mut ctx, -1.0, 10.0).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_nan_y() {
        let mut ctx = ais_interactive_context().unwrap();
        let err = ais_select_single(&mut ctx, 10.0, f32::NAN).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn no_camera_picks_nothing() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.display_geometry(axis_triangle());
        // No set_view → ray_at returns None → no pick.
        assert_eq!(ais_select_single(&mut ctx, 50.0, 50.0).unwrap(), None);
    }

    #[test]
    fn picks_and_selects_the_object_under_the_cursor() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        let id = ctx.display_geometry(axis_triangle());
        let picked = ais_select_single(&mut ctx, 50.0, 50.0).unwrap();
        assert_eq!(picked, Some(id));
        assert_eq!(ctx.state(id), Some(ObjectState::Selected));
    }

    #[test]
    fn miss_returns_none_and_changes_nothing() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        let id = ctx.display_geometry(axis_triangle());
        // A corner pixel — the centred triangle isn't there.
        let picked = ais_select_single(&mut ctx, 1.0, 1.0).unwrap();
        assert_eq!(picked, None);
        assert_eq!(ctx.state(id), Some(ObjectState::Visible));
    }
}
