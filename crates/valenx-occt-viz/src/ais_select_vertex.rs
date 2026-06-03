//! Phase 179 — vertex-level selection inside a solid.
//!
//! ## What OCCT does
//!
//! Activation mode 1 (`StdSelect_BRepOwner_for_Vertices`) makes picks
//! return individual `TopoDS_Vertex` owners. The picking pass tests
//! each vertex's screen-space position against a small square (typical
//! tolerance = 6 px) — vertices are zero-extent geometry so the
//! tolerance is generous compared to edges or faces.
//!
//! ## v1 status — real CPU subshape ray-cast
//!
//! With the picking substrate from Phase 171 in place, this is a real
//! vertex picker. The parent solid is registered with a
//! [`Pickable::TaggedVertices`](crate::ais_interactive_context::Pickable::TaggedVertices)
//! — each labelled with its vertex index. The cursor pixel is
//! unprojected to a world ray; the ray's closest approach to each
//! vertex is compared against a generous tolerance (vertices are
//! zero-extent geometry, so the pick square is wider than for edges).
//! The nearest hit's vertex index is returned.
//!
//! With no camera installed the pick returns `None`.

use crate::ais_interactive_context::InteractiveContext;
use crate::ais_select_single::PICK_TOLERANCE;
use crate::error::OcctVizError;

/// Pick the vertex at screen-pixel `(x, y)` inside the solid identified
/// by `parent_id`. Returns `Option<usize>` — the vertex index within
/// the parent's vertex list, or `None` if no vertex was hit within the
/// pick tolerance.
///
/// # Errors
///
/// [`OcctVizError::BadInput`] for non-finite or negative coordinates,
/// or unknown `parent_id`.
pub fn ais_select_vertex(
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
    // Vertices get the most generous tolerance — they have zero extent.
    Ok(ctx.pick_subshape(parent_id, &ray, PICK_TOLERANCE * 6.0))
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

    /// A solid whose vertex 2 sits at the world origin and vertex 8
    /// sits far away.
    fn two_vertex_solid() -> Pickable {
        Pickable::TaggedVertices {
            vertices: vec![
                (2, [0.0, 0.0, 0.0]),
                (8, [40.0, 40.0, 0.0]),
            ],
        }
    }

    #[test]
    fn rejects_unknown_parent() {
        let mut ctx = ais_interactive_context().unwrap();
        let err = ais_select_vertex(&mut ctx, 99, 10.0, 20.0).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn no_camera_returns_none() {
        let mut ctx = ais_interactive_context().unwrap();
        let id = ctx.display_geometry(two_vertex_solid());
        assert_eq!(ais_select_vertex(&mut ctx, id, 50.0, 50.0).unwrap(), None);
    }

    #[test]
    fn picks_the_vertex_under_the_cursor() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        let id = ctx.display_geometry(two_vertex_solid());
        // The centre ray passes through the origin → vertex 2.
        let vtx = ais_select_vertex(&mut ctx, id, 50.0, 50.0).unwrap();
        assert_eq!(vtx, Some(2));
    }

    #[test]
    fn miss_returns_none() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        let id = ctx.display_geometry(two_vertex_solid());
        assert_eq!(ais_select_vertex(&mut ctx, id, 2.0, 2.0).unwrap(), None);
    }
}
