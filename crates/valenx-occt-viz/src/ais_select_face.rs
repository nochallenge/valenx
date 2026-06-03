//! Phase 177 — face-level selection inside a solid.
//!
//! ## What OCCT does
//!
//! `AIS_InteractiveContext::Activate(StdSelect_BRepOwner_for_Faces)`
//! switches the selection mode so picks return individual `TopoDS_Face`
//! owners rather than the parent solid. The picking pass then tests
//! each face's triangulation against the ray and returns the first hit
//! along with the parent shape ID.
//!
//! ## v1 status — real CPU subshape ray-cast
//!
//! With the picking substrate from Phase 171 in place, this is a real
//! face picker. The parent solid is registered with a
//! [`Pickable::TaggedMesh`](crate::ais_interactive_context::Pickable::TaggedMesh)
//! — each triangle labelled with its source face index, exactly the
//! "per-face tri-range tagging" OCCT does. The cursor pixel is
//! unprojected to a world ray, the ray is tested against the parent's
//! triangles, and the nearest hit's face tag is returned.
//!
//! A parent registered without geometry (or with an untagged
//! [`Pickable::Mesh`](crate::ais_interactive_context::Pickable::Mesh))
//! has no per-face structure — the pick still works but returns face
//! index `0`. With no camera installed the pick returns `None`.

use crate::ais_interactive_context::InteractiveContext;
use crate::ais_select_single::PICK_TOLERANCE;
use crate::error::OcctVizError;

/// Pick the face at screen-pixel `(x, y)` inside the solid identified
/// by `parent_id`. Returns `Option<usize>` — the face index within
/// the parent's face list, or `None` if no face was hit.
///
/// # Errors
///
/// [`OcctVizError::BadInput`] for non-finite or negative coordinates,
/// or unknown `parent_id`.
pub fn ais_select_face(
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
    Ok(ctx.pick_subshape(parent_id, &ray, PICK_TOLERANCE))
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

    /// A solid with two faces: face 7 near (z=0), face 9 off to the
    /// side. Triangles tagged with their face index.
    fn two_face_solid() -> Pickable {
        Pickable::TaggedMesh {
            triangles: vec![
                // Face 7 — centred on the axis.
                [-1.0, -1.0, 0.0],
                [1.0, -1.0, 0.0],
                [0.0, 1.0, 0.0],
                // Face 9 — far off to the side.
                [50.0, 50.0, 0.0],
                [52.0, 50.0, 0.0],
                [51.0, 52.0, 0.0],
            ],
            face_tags: vec![7, 9],
        }
    }

    #[test]
    fn rejects_negative_coords() {
        let mut ctx = ais_interactive_context().unwrap();
        let id = ctx.display();
        let err = ais_select_face(&mut ctx, id, -1.0, 5.0).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_unknown_parent() {
        let mut ctx = ais_interactive_context().unwrap();
        let err = ais_select_face(&mut ctx, 99, 10.0, 20.0).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn no_camera_returns_none() {
        let mut ctx = ais_interactive_context().unwrap();
        let id = ctx.display_geometry(two_face_solid());
        assert_eq!(ais_select_face(&mut ctx, id, 50.0, 50.0).unwrap(), None);
    }

    #[test]
    fn picks_the_face_under_the_cursor() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        let id = ctx.display_geometry(two_face_solid());
        // The centre ray hits face 7.
        let face = ais_select_face(&mut ctx, id, 50.0, 50.0).unwrap();
        assert_eq!(face, Some(7));
    }

    #[test]
    fn miss_returns_none() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        let id = ctx.display_geometry(two_face_solid());
        // A corner pixel hits neither tagged face.
        assert_eq!(ais_select_face(&mut ctx, id, 1.0, 1.0).unwrap(), None);
    }
}
