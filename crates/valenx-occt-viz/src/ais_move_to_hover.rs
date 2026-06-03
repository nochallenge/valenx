//! Sub-face hover preview — `AIS_InteractiveContext::MoveTo` with
//! sub-element resolution (Phase 171-180 polish).
//!
//! ## What OCCT does
//!
//! Every time the cursor moves, OCCT calls
//! `AIS_InteractiveContext::MoveTo(x, y, view)`. It ray-casts the
//! pixel against the registered selection primitives and *dynamically
//! highlights* whatever is under the cursor **before** the user
//! clicks — and it highlights at the granularity of the active
//! selection mode: the whole object in mode 0, an individual face in
//! mode 4, an edge in mode 2, a vertex in mode 1. That live "this is
//! what you would pick" feedback is what makes sub-face selection
//! usable.
//!
//! The bare [`ais_highlight_dynamic`](fn@crate::ais_highlight_dynamic)
//! op already flips an object's [`ObjectState`] to `Hovered`, but it
//! requires the caller to *already know which object* — it does no
//! picking. And it has no
//! notion of *which sub-element* is under the cursor.
//!
//! ## v1 status — real sub-element hover state machine
//!
//! This module ships the missing piece: a [`HoverPreview`] that holds
//! the cursor-resolved hover target, and [`move_to`], the real
//! `MoveTo` op:
//!
//! 1. Unproject the cursor pixel to a world ray
//!    ([`InteractiveContext::ray_at`]).
//! 2. Ray-cast every pickable object; the nearest hit wins.
//! 3. If the hit object carries a `Tagged*` [`Pickable`], also
//!    resolve **which sub-element** (face / edge / vertex index) the
//!    ray struck — exactly OCCT's subshape owner.
//! 4. Update the object's [`ObjectState`] to `Hovered`, clearing the
//!    hover on whatever was hovered last frame. A `Selected` object
//!    is left `Selected` (hover never demotes a selection).
//!
//! [`HoverPreview::confirm`] then promotes the current hover to
//! `Selected` — that is the click, and it returns the resolved
//! [`HoverTarget`] so the caller knows whether a whole object, a
//! face, an edge, or a vertex was picked.
//!
//! Together [`move_to`] + [`HoverPreview::confirm`] are the
//! hover-then-click loop OCCT's selection UX is built on, now with
//! genuine sub-face granularity.

use crate::ais_interactive_context::{InteractiveContext, ObjectState, Pickable};
use crate::error::OcctVizError;

/// World-space proximity (model units) for edge / vertex hover picks —
/// matches [`crate::ais_select_single::PICK_TOLERANCE`]. A triangle
/// hover is exact; a polyline / point hover counts when the ray passes
/// within this distance.
pub const HOVER_TOLERANCE: f64 = 0.05;

/// Which kind of sub-element the cursor resolved to.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, serde::Serialize, serde::Deserialize)]
pub enum SubKind {
    /// The whole object (the picked [`Pickable`] carries no
    /// sub-element tags — `Mesh` / `Polyline` / `Point`).
    Whole,
    /// An individual face of a tagged solid.
    Face,
    /// An individual edge of a tagged solid.
    Edge,
    /// An individual vertex of a tagged solid.
    Vertex,
}

/// The fully-resolved hover / selection target: which object, which
/// sub-element kind, and (for a sub-element) its index.
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct HoverTarget {
    /// ID of the object under the cursor.
    pub object: usize,
    /// The granularity that resolved.
    pub kind: SubKind,
    /// Sub-element index — the face / edge / vertex index for a
    /// `Tagged*` pick, `0` for a `Whole` pick.
    pub sub_index: usize,
    /// Ray parameter `t` at the hit (distance along the cursor ray).
    pub distance: f64,
}

/// The live hover-preview state machine — OCCT's "what would I pick"
/// feedback.
///
/// Construct one per viewport, feed it [`move_to`] every time the
/// cursor moves, and call [`HoverPreview::confirm`] on click.
#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct HoverPreview {
    /// The current hover target, or `None` when the cursor is over
    /// empty space.
    current: Option<HoverTarget>,
}

impl HoverPreview {
    /// A fresh preview with nothing hovered.
    pub fn new() -> Self {
        Self::default()
    }

    /// The currently-hovered target, if any.
    pub fn current(&self) -> Option<HoverTarget> {
        self.current
    }

    /// Whether anything is under the cursor right now.
    pub fn is_hovering(&self) -> bool {
        self.current.is_some()
    }

    /// Promote the current hover to a `Selected` object in `ctx` and
    /// return the resolved [`HoverTarget`].
    ///
    /// This is the click. Returns `None` (and changes nothing) when
    /// the cursor is over empty space. The previous selection is *not*
    /// cleared — the caller decides single-select-replaces-all vs.
    /// additive selection (call [`InteractiveContext::set_state`] to
    /// reset other objects first for the former).
    pub fn confirm(&self, ctx: &mut InteractiveContext) -> Option<HoverTarget> {
        let target = self.current?;
        ctx.set_state(target.object, ObjectState::Selected);
        Some(target)
    }
}

/// Resolve what the cursor at screen-pixel `(x, y)` is over and update
/// the hover preview + object states accordingly.
///
/// This is OCCT's `AIS_InteractiveContext::MoveTo`. It:
///
/// - ray-casts the cursor pixel against every pickable object,
/// - resolves the nearest hit down to its sub-element (face / edge /
///   vertex) when the object carries `Tagged*` geometry,
/// - flips the newly-hovered object to [`ObjectState::Hovered`] and
///   clears the hover on whatever was hovered before (a `Selected`
///   object keeps its `Selected` state — hover never demotes it),
/// - records the resolved target in `preview`.
///
/// Returns the new [`HoverTarget`] (`None` when the cursor is over
/// empty space, or when no camera has been installed via
/// [`InteractiveContext::set_view`]).
///
/// # Errors
///
/// [`OcctVizError::BadInput`] if either coordinate is non-finite or
/// negative.
pub fn move_to(
    ctx: &mut InteractiveContext,
    preview: &mut HoverPreview,
    x: f32,
    y: f32,
) -> Result<Option<HoverTarget>, OcctVizError> {
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

    // Resolve the hit. No camera ⇒ no ray ⇒ nothing hovered.
    let new_target = ctx.ray_at(x, y).and_then(|ray| {
        let hit = ctx.pick_nearest(&ray, HOVER_TOLERANCE)?;
        // Resolve the sub-element of the hit object.
        let (kind, sub_index) = match ctx.geometry(hit.id) {
            Some(geom) => {
                let kind = sub_kind_of(geom);
                let sub = if kind == SubKind::Whole {
                    0
                } else {
                    // pick_subshape re-casts against this object only.
                    ctx.pick_subshape(hit.id, &ray, HOVER_TOLERANCE)
                        .unwrap_or(0)
                };
                (kind, sub)
            }
            None => (SubKind::Whole, 0),
        };
        Some(HoverTarget {
            object: hit.id,
            kind,
            sub_index,
            distance: hit.distance,
        })
    });

    // Clear the previous hover (unless it became Selected meanwhile).
    if let Some(prev) = preview.current {
        let still_hovering_same = new_target.map(|t| t.object) == Some(prev.object);
        if !still_hovering_same && ctx.state(prev.object) == Some(ObjectState::Hovered) {
            ctx.set_state(prev.object, ObjectState::Visible);
        }
    }

    // Apply the new hover (don't demote a Selected object).
    if let Some(target) = new_target {
        if ctx.state(target.object) != Some(ObjectState::Selected) {
            ctx.set_state(target.object, ObjectState::Hovered);
        }
    }

    preview.current = new_target;
    Ok(new_target)
}

/// The sub-element granularity a [`Pickable`] supports.
fn sub_kind_of(geom: &Pickable) -> SubKind {
    match geom {
        Pickable::Mesh { .. } | Pickable::Polyline { .. } | Pickable::Point { .. } => {
            SubKind::Whole
        }
        Pickable::TaggedMesh { .. } => SubKind::Face,
        Pickable::TaggedEdges { .. } => SubKind::Edge,
        Pickable::TaggedVertices { .. } => SubKind::Vertex,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ais_interactive_context::ais_interactive_context;
    use nalgebra::{Matrix4, Point3, Vector3};

    /// A camera looking down -Z at the origin from z = +10.
    fn test_view() -> Matrix4<f64> {
        let view = Matrix4::look_at_rh(
            &Point3::new(0.0, 0.0, 10.0),
            &Point3::origin(),
            &Vector3::y(),
        );
        let proj = Matrix4::new_perspective(1.0, 60f64.to_radians(), 0.1, 100.0);
        proj * view
    }

    /// A triangle facing the camera, centred on the axis.
    fn axis_triangle() -> Pickable {
        Pickable::Mesh {
            triangles: vec![
                [-1.0, -1.0, 0.0],
                [1.0, -1.0, 0.0],
                [0.0, 1.0, 0.0],
            ],
        }
    }

    /// A tagged two-face mesh: triangle 0 is face 7, triangle 1 is
    /// face 9. Both lie at z = 0 in front of the camera; face 7 fills
    /// the left half, face 9 the right.
    fn two_face_solid() -> Pickable {
        Pickable::TaggedMesh {
            triangles: vec![
                // face 7 — a triangle on the -X side.
                [-2.0, -1.0, 0.0],
                [0.0, -1.0, 0.0],
                [-1.0, 1.0, 0.0],
                // face 9 — a triangle on the +X side.
                [0.0, -1.0, 0.0],
                [2.0, -1.0, 0.0],
                [1.0, 1.0, 0.0],
            ],
            face_tags: vec![7, 9],
        }
    }

    #[test]
    fn rejects_negative_coordinates() {
        let mut ctx = ais_interactive_context().unwrap();
        let mut hp = HoverPreview::new();
        assert!(move_to(&mut ctx, &mut hp, -1.0, 5.0).is_err());
        assert!(move_to(&mut ctx, &mut hp, 5.0, f32::NAN).is_err());
    }

    #[test]
    fn fresh_preview_hovers_nothing() {
        let hp = HoverPreview::new();
        assert!(!hp.is_hovering());
        assert!(hp.current().is_none());
    }

    #[test]
    fn no_camera_hovers_nothing() {
        let mut ctx = ais_interactive_context().unwrap();
        let mut hp = HoverPreview::new();
        ctx.display_geometry(axis_triangle());
        // No set_view → ray_at returns None.
        assert_eq!(move_to(&mut ctx, &mut hp, 50.0, 50.0).unwrap(), None);
        assert!(!hp.is_hovering());
    }

    #[test]
    fn move_to_hovers_the_object_under_the_cursor() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        let id = ctx.display_geometry(axis_triangle());
        let mut hp = HoverPreview::new();
        let target = move_to(&mut ctx, &mut hp, 50.0, 50.0).unwrap();
        let target = target.expect("centre cursor should hover the triangle");
        assert_eq!(target.object, id);
        assert_eq!(target.kind, SubKind::Whole);
        // The object's state is now Hovered.
        assert_eq!(ctx.state(id), Some(ObjectState::Hovered));
    }

    #[test]
    fn moving_off_the_object_clears_the_hover() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        let id = ctx.display_geometry(axis_triangle());
        let mut hp = HoverPreview::new();
        // Hover the triangle.
        move_to(&mut ctx, &mut hp, 50.0, 50.0).unwrap();
        assert_eq!(ctx.state(id), Some(ObjectState::Hovered));
        // Move to a corner — empty space.
        let target = move_to(&mut ctx, &mut hp, 1.0, 1.0).unwrap();
        assert!(target.is_none());
        // The hover was released back to Visible.
        assert_eq!(ctx.state(id), Some(ObjectState::Visible));
        assert!(!hp.is_hovering());
    }

    #[test]
    fn move_to_resolves_the_face_under_the_cursor() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        let _id = ctx.display_geometry(two_face_solid());
        let mut hp = HoverPreview::new();
        // Cursor left-of-centre → over face 7 (the -X triangle). Face
        // 7's two-triangle solid projects to screen-x ≈ [32.7, 50] —
        // (40, 55) lands well inside that projected triangle.
        let left = move_to(&mut ctx, &mut hp, 40.0, 55.0).unwrap();
        let left = left.expect("left cursor should hover a face");
        assert_eq!(left.kind, SubKind::Face);
        assert_eq!(left.sub_index, 7, "left cursor should resolve face 7");
        // Cursor right-of-centre → over face 9 (the +X triangle),
        // projected to screen-x ≈ [50, 67.3].
        let right = move_to(&mut ctx, &mut hp, 60.0, 55.0).unwrap();
        let right = right.expect("right cursor should hover a face");
        assert_eq!(right.kind, SubKind::Face);
        assert_eq!(right.sub_index, 9, "right cursor should resolve face 9");
    }

    #[test]
    fn confirm_promotes_the_hover_to_selected() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        let id = ctx.display_geometry(two_face_solid());
        let mut hp = HoverPreview::new();
        // Hover face 9 (its projected triangle is screen-x ≈ [50, 67.3]).
        move_to(&mut ctx, &mut hp, 60.0, 55.0).unwrap();
        // Click.
        let picked = hp.confirm(&mut ctx).expect("confirm should pick");
        assert_eq!(picked.object, id);
        assert_eq!(picked.kind, SubKind::Face);
        assert_eq!(picked.sub_index, 9);
        assert_eq!(ctx.state(id), Some(ObjectState::Selected));
    }

    #[test]
    fn confirm_over_empty_space_picks_nothing() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        ctx.display_geometry(axis_triangle());
        let mut hp = HoverPreview::new();
        // Cursor over a corner — empty.
        move_to(&mut ctx, &mut hp, 1.0, 1.0).unwrap();
        assert!(hp.confirm(&mut ctx).is_none());
    }

    #[test]
    fn hover_does_not_demote_a_selected_object() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        let id = ctx.display_geometry(axis_triangle());
        let mut hp = HoverPreview::new();
        // Select the object outright.
        ctx.set_state(id, ObjectState::Selected);
        // Hover it — it must stay Selected, not drop to Hovered.
        move_to(&mut ctx, &mut hp, 50.0, 50.0).unwrap();
        assert_eq!(ctx.state(id), Some(ObjectState::Selected));
        // Moving away must not demote the still-Selected object.
        move_to(&mut ctx, &mut hp, 1.0, 1.0).unwrap();
        assert_eq!(ctx.state(id), Some(ObjectState::Selected));
    }

    #[test]
    fn re_hovering_the_same_object_keeps_it_hovered() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        let id = ctx.display_geometry(axis_triangle());
        let mut hp = HoverPreview::new();
        // Two consecutive move_to calls over the same triangle.
        move_to(&mut ctx, &mut hp, 48.0, 50.0).unwrap();
        move_to(&mut ctx, &mut hp, 52.0, 50.0).unwrap();
        assert_eq!(ctx.state(id), Some(ObjectState::Hovered));
        assert!(hp.is_hovering());
    }
}
