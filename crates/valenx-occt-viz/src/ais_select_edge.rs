//! Phase 178 — edge-level selection inside a solid.
//!
//! ## What OCCT does
//!
//! Activation mode 2 (`StdSelect_BRepOwner_for_Edges`) makes picks
//! return individual `TopoDS_Edge` owners. The picking pass tests each
//! edge polyline against a screen-space cylinder (radius = the pick
//! tolerance, typically 4 px) so users don't need to click *exactly*
//! on the 1-pixel-wide edge line.
//!
//! ## v1 status — real CPU subshape ray-cast
//!
//! With the picking substrate from Phase 171 in place, this is a real
//! edge picker. The parent solid is registered with a
//! [`Pickable::TaggedEdges`](crate::ais_interactive_context::Pickable::TaggedEdges)
//! — each edge a labelled polyline. The cursor pixel is unprojected to
//! a world ray and the ray's closest approach to each edge polyline is
//! compared against [`crate::ais_select_single::PICK_TOLERANCE`] — the
//! screen-space cylinder OCCT uses so users don't have to click the
//! 1-pixel edge line exactly. The nearest hit's edge index is returned.
//!
//! With no camera installed the pick returns `None`.

use crate::ais_interactive_context::InteractiveContext;
use crate::ais_select_single::PICK_TOLERANCE;
use crate::error::OcctVizError;

/// Pick the edge at screen-pixel `(x, y)` inside the solid identified
/// by `parent_id`. Returns `Option<usize>` — the edge index within
/// the parent's edge list, or `None` if no edge was hit within the
/// pick tolerance.
///
/// # Errors
///
/// [`OcctVizError::BadInput`] for non-finite or negative coordinates,
/// or unknown `parent_id`.
pub fn ais_select_edge(
    ctx: &mut InteractiveContext,
    parent_id: usize,
    x: f32,
    y: f32,
) -> Result<Option<usize>, OcctVizError> {
    if !x.is_finite() || x < 0.0 || !y.is_finite() || y < 0.0 {
        return Err(OcctVizError::bad_input(
            "x_or_y",
            format!("must be finite and >= 0 (got {x}, {y})"),
        ));
    }
    if ctx.state(parent_id).is_none() {
        return Err(OcctVizError::bad_input(
            "parent_id",
            format!("not registered: {parent_id}"),
        ));
    }
    let Some(ray) = ctx.ray_at(x, y) else {
        return Ok(None);
    };
    // Edges are picked with a generous tolerance — they are 1-px-wide
    // lines, so the screen-space pick cylinder needs slack.
    Ok(ctx.pick_subshape(parent_id, &ray, PICK_TOLERANCE * 4.0))
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

    /// A solid whose edge 3 passes through the world origin and edge 5
    /// sits far away.
    fn two_edge_solid() -> Pickable {
        Pickable::TaggedEdges {
            edges: vec![
                (3, vec![[-2.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
                (5, vec![[40.0, 40.0, 0.0], [44.0, 40.0, 0.0]]),
            ],
        }
    }

    #[test]
    fn rejects_nan_y() {
        let mut ctx = ais_interactive_context().unwrap();
        let id = ctx.display();
        let err = ais_select_edge(&mut ctx, id, 10.0, f32::NAN).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn no_camera_returns_none() {
        let mut ctx = ais_interactive_context().unwrap();
        let id = ctx.display_geometry(two_edge_solid());
        assert_eq!(ais_select_edge(&mut ctx, id, 50.0, 50.0).unwrap(), None);
    }

    #[test]
    fn picks_the_edge_under_the_cursor() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        let id = ctx.display_geometry(two_edge_solid());
        // The centre ray passes through the origin → edge 3.
        let edge = ais_select_edge(&mut ctx, id, 50.0, 50.0).unwrap();
        assert_eq!(edge, Some(3));
    }

    #[test]
    fn miss_returns_none() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        let id = ctx.display_geometry(two_edge_solid());
        // The top-left corner hits neither edge.
        assert_eq!(ais_select_edge(&mut ctx, id, 2.0, 2.0).unwrap(), None);
    }
}
