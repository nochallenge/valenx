//! Orthographic projection groups + isometric layout.
//!
//! A [`ProjectionGroup`] auto-arranges Front / Top / Right (and an
//! optional isometric) views from a single [`valenx_cad::Solid`] at a
//! base position on the sheet with a chosen scale. The relative
//! placement follows the chosen [`Projection`] convention:
//!
//! - **Third-angle** (ASME default, common in North America) — Top is
//!   placed **above** the Front view; Right is placed to the **right**
//!   of Front. The mental model: the view is the projection seen
//!   through a glass cube *in front of* the object.
//! - **First-angle** (ISO default, common in Europe) — Top is placed
//!   **below** Front; Right is placed to the **left**. The mental
//!   model: the view is what you'd see by *unfolding* the box around
//!   the object outward.
//!
//! The isometric view, when included, goes in the upper-right corner
//! of the group, offset from the Front view by Front's width plus a
//! gap so it doesn't overlap with Right.
//!
//! Layout math operates on the source solid's world-space bounding
//! box: width / height in each view are derived from the BBox extents
//! of the *projected* edges, scaled by [`ProjectionGroup::scale`], then
//! padded by [`ProjectionGroup::gap_mm`].
//!
//! Auto-labels (`"Front"`, `"Top"`, `"Right"`, `"Isometric"`) are
//! emitted as `(x_mm, y_mm, label_text)` tuples — the renderer places
//! them above each view's bbox top.
//!
//! When `feature_id` is `Some`, [`crate::Drawing::regenerate_all`]
//! treats every view in the group as parametric: changing the source
//! feature reruns the whole group through [`ProjectionGroup::rebuild`].
//!
//! # Example
//!
//! ```
//! use valenx_techdraw::{Drawing, Sheet};
//! use valenx_techdraw::projection_group::{ProjectionGroup, Projection};
//! use valenx_cad::primitives::box_solid;
//!
//! let mut d = Drawing::new(Sheet::a3_landscape("Bracket", "A. Engineer", "A"));
//! let cube = box_solid(50.0, 30.0, 20.0).unwrap();
//! let mut group = ProjectionGroup::new(
//!     [80.0, 80.0],        // base position (Front's lower-left)
//!     1.0,                 // scale
//!     Projection::ThirdAngle,
//! );
//! group.gap_mm = 30.0;
//! group.include_isometric = true;
//! let view_ids = group.build_into(&mut d, &cube).unwrap();
//! assert_eq!(view_ids.len(), 4); // Front, Top, Right, Iso
//! ```

use serde::{Deserialize, Serialize};

use crate::error::TechDrawError;
use crate::view::{View, ViewKind};
use crate::Drawing;
use valenx_feature_tree::feature::FeatureId;

/// First- vs third-angle orthographic projection convention.
///
/// Affects the *placement* of Top and Right relative to Front. The
/// view content (the projected edges) is identical between the two —
/// only the layout differs.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Projection {
    /// ASME / US default: Top above Front, Right to the right.
    ThirdAngle,
    /// ISO / European default: Top below Front, Right to the left.
    FirstAngle,
}

impl Projection {
    /// Short label for UI.
    pub fn label(self) -> &'static str {
        match self {
            Projection::ThirdAngle => "Third-angle",
            Projection::FirstAngle => "First-angle",
        }
    }
}

/// A group of auto-arranged projection views (Front / Top / Right /
/// optional Isometric) for a single source solid.
///
/// Carries layout knobs (base position, scale, gap, isometric toggle,
/// projection convention) and an optional parametric link
/// (`feature_id`). The actual views and their on-sheet positions are
/// computed on demand by [`Self::build_into`] (creates the views) or
/// [`Self::rebuild`] (regenerates pre-existing views).
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct ProjectionGroup {
    /// Position of the Front view's lower-left corner on the sheet
    /// (mm).
    pub base_position: [f64; 2],
    /// Common scale applied to every view in the group.
    pub scale: f64,
    /// Projection convention (controls Top + Right placement).
    pub projection: Projection,
    /// Gap between adjacent view bounding boxes (mm). Default 20 mm
    /// matches the standard "1 inch breathing room" most drawings use.
    pub gap_mm: f64,
    /// Whether to add an isometric view in the upper-right corner.
    pub include_isometric: bool,
    /// Optional [`FeatureId`] backing this group. When `Some`,
    /// [`crate::Drawing::regenerate_all`] will rebuild every view in
    /// this group when the feature tree replays.
    pub feature_id: Option<FeatureId>,
    /// Indices of the views that belong to this group in the drawing's
    /// `views` vector. Populated by [`Self::build_into`]; the same
    /// slots are reused by [`Self::rebuild`].
    pub view_indices: Vec<usize>,
}

impl ProjectionGroup {
    /// Construct a new group with sensible defaults: 20 mm gap, no
    /// isometric, no parametric link.
    pub fn new(base_position: [f64; 2], scale: f64, projection: Projection) -> Self {
        Self {
            base_position,
            scale,
            projection,
            gap_mm: 20.0,
            include_isometric: false,
            feature_id: None,
            view_indices: Vec::new(),
        }
    }

    /// Build the four (or three) views, generate their edges from
    /// `solid`, append them to `drawing`, store their indices, and
    /// return the index list. The returned vector is ordered:
    /// `[front, top, right]` or `[front, top, right, iso]`.
    ///
    /// Errors: forwarded from [`View::generate`] (scale ≤ 0,
    /// empty solid, etc.). On failure the drawing is left in whatever
    /// state succeeded so far — callers can `remove_view` to clean up.
    pub fn build_into(
        &mut self,
        drawing: &mut Drawing,
        solid: &valenx_cad::Solid,
    ) -> Result<Vec<usize>, TechDrawError> {
        if !self.scale.is_finite() || self.scale <= 0.0 {
            return Err(TechDrawError::BadViewParameter {
                name: "scale",
                reason: format!("must be > 0, got {}", self.scale),
            });
        }
        // First pass: generate Front edges so we know its bbox; the
        // bbox drives Top + Right placement.
        let mut front = View::new(ViewKind::Front, self.scale, self.base_position);
        front.generate(solid)?;
        let (front_min, front_max) = front.bbox().ok_or(TechDrawError::EmptySolid)?;
        let front_w = (front_max[0] - front_min[0]) * self.scale;
        let front_h = (front_max[1] - front_min[1]) * self.scale;

        let positions = self.layout_positions(front_w, front_h);
        front.position = positions.front;
        let front_idx = drawing.add_view(front);

        let mut top = View::new(ViewKind::Top, self.scale, positions.top);
        top.generate(solid)?;
        let top_idx = drawing.add_view(top);

        let mut right = View::new(ViewKind::Right, self.scale, positions.right);
        right.generate(solid)?;
        let right_idx = drawing.add_view(right);

        let mut out = vec![front_idx, top_idx, right_idx];
        if self.include_isometric {
            let mut iso = View::new(ViewKind::Isometric, self.scale, positions.iso);
            iso.generate(solid)?;
            out.push(drawing.add_view(iso));
        }
        self.view_indices = out.clone();
        Ok(out)
    }

    /// Re-extract edges for every view in the group from `solid`. Used
    /// by [`crate::Drawing::regenerate_all`] when the source feature
    /// changes. The view positions are left as-is so user-edited
    /// layouts survive a regenerate.
    ///
    /// Returns a list of per-view errors (view index + error) for any
    /// view that fails to regenerate. An empty vec means everything
    /// succeeded.
    pub fn rebuild(
        &self,
        drawing: &mut Drawing,
        solid: &valenx_cad::Solid,
    ) -> Vec<(usize, TechDrawError)> {
        let mut errs = Vec::new();
        for &idx in &self.view_indices {
            if let Some(v) = drawing.views.get_mut(idx) {
                if let Err(e) = v.generate(solid) {
                    errs.push((idx, e));
                }
            }
        }
        errs
    }

    /// Compute per-view target positions on the sheet given the Front
    /// view's already-known projected width + height (in mm at this
    /// group's scale). Public so tests can verify the layout without
    /// instantiating a [`Drawing`].
    pub fn layout_positions(&self, front_w: f64, front_h: f64) -> GroupLayout {
        let [bx, by] = self.base_position;
        let gap = self.gap_mm;
        match self.projection {
            Projection::ThirdAngle => {
                let front = [bx, by];
                // Top placed above Front.
                let top = [bx, by + front_h + gap];
                // Right placed to the right of Front.
                let right = [bx + front_w + gap, by];
                // Iso in the upper-right corner (above Right).
                let iso = [bx + front_w + gap, by + front_h + gap];
                GroupLayout {
                    front,
                    top,
                    right,
                    iso,
                }
            }
            Projection::FirstAngle => {
                let front = [bx, by];
                // Top placed below Front (so Front goes above Top).
                let top = [bx, by - (front_h + gap)];
                // Right placed to the left of Front.
                let right = [bx - (front_w + gap), by];
                // Iso in the (relative) upper-right — at +x, +y from Front.
                let iso = [bx + front_w + gap, by + front_h + gap];
                GroupLayout {
                    front,
                    top,
                    right,
                    iso,
                }
            }
        }
    }

    /// Return suggested `(x_mm, y_mm, label_text)` label positions for
    /// the views in this group. Placed below each view's bbox (the
    /// standard "view caption" position in engineering drawings).
    ///
    /// Front_h is required so the caption can be placed below Front
    /// regardless of how the per-view local frame is positioned.
    pub fn label_positions(&self, front_w: f64, front_h: f64) -> Vec<(f64, f64, String)> {
        let l = self.layout_positions(front_w, front_h);
        let cap_dy = -8.0; // label 8 mm below view bottom
        let mut out = vec![
            (
                l.front[0] + front_w * 0.5,
                l.front[1] + cap_dy,
                "Front".into(),
            ),
            (l.top[0] + front_w * 0.5, l.top[1] + cap_dy, "Top".into()),
            (
                l.right[0] + front_w * 0.5, // approximated — right view is roughly the same width
                l.right[1] + cap_dy,
                "Right".into(),
            ),
        ];
        if self.include_isometric {
            out.push((
                l.iso[0] + front_w * 0.5,
                l.iso[1] + cap_dy,
                "Isometric".into(),
            ));
        }
        out
    }
}

/// Computed positions for every view in a [`ProjectionGroup`].
///
/// `iso` is populated regardless of [`ProjectionGroup::include_isometric`]
/// — callers consult that flag to know whether to use it.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct GroupLayout {
    /// Front view position (mm).
    pub front: [f64; 2],
    /// Top view position (mm).
    pub top: [f64; 2],
    /// Right view position (mm).
    pub right: [f64; 2],
    /// Isometric view position (mm).
    pub iso: [f64; 2],
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sheet::Sheet;
    use valenx_cad::primitives::box_solid;

    /// Layout math sanity: Top is *above* Front for third-angle.
    #[test]
    fn third_angle_top_is_above_front() {
        let g = ProjectionGroup::new([100.0, 50.0], 1.0, Projection::ThirdAngle);
        let l = g.layout_positions(40.0, 30.0);
        assert_eq!(l.front, [100.0, 50.0]);
        assert!(l.top[1] > l.front[1], "Top should be above Front");
        assert!(
            l.right[0] > l.front[0],
            "Right should be to the right of Front"
        );
        // Iso in the upper-right corner.
        assert!(l.iso[0] > l.front[0] && l.iso[1] > l.front[1]);
    }

    /// Layout math sanity: Top is *below* Front for first-angle.
    #[test]
    fn first_angle_top_is_below_front() {
        let g = ProjectionGroup::new([100.0, 100.0], 1.0, Projection::FirstAngle);
        let l = g.layout_positions(40.0, 30.0);
        assert_eq!(l.front, [100.0, 100.0]);
        assert!(l.top[1] < l.front[1], "Top should be below Front");
        assert!(
            l.right[0] < l.front[0],
            "Right should be to the left of Front"
        );
    }

    /// Build a third-angle group from a 50×30×20 box; verify the three
    /// view kinds are right and Top ends up above Front.
    #[test]
    fn build_third_angle_group_places_views_correctly() {
        let cube = box_solid(50.0, 30.0, 20.0).unwrap();
        let mut d = Drawing::new(Sheet::a3_landscape("X", "Y", "Z"));
        let mut g = ProjectionGroup::new([100.0, 100.0], 1.0, Projection::ThirdAngle);
        let ids = g.build_into(&mut d, &cube).unwrap();
        assert_eq!(ids.len(), 3);
        assert_eq!(d.views[ids[0]].kind, ViewKind::Front);
        assert_eq!(d.views[ids[1]].kind, ViewKind::Top);
        assert_eq!(d.views[ids[2]].kind, ViewKind::Right);
        let front_pos = d.views[ids[0]].position;
        let top_pos = d.views[ids[1]].position;
        let right_pos = d.views[ids[2]].position;
        assert!(top_pos[1] > front_pos[1], "Top above Front");
        assert!(right_pos[0] > front_pos[0], "Right right of Front");
    }

    /// Including the iso adds a fourth view kind.
    #[test]
    fn build_with_isometric_adds_fourth_view() {
        let cube = box_solid(20.0, 20.0, 20.0).unwrap();
        let mut d = Drawing::new(Sheet::a3_landscape("X", "Y", "Z"));
        let mut g = ProjectionGroup::new([100.0, 100.0], 1.0, Projection::ThirdAngle);
        g.include_isometric = true;
        let ids = g.build_into(&mut d, &cube).unwrap();
        assert_eq!(ids.len(), 4);
        assert_eq!(d.views[ids[3]].kind, ViewKind::Isometric);
        // Iso should be in the upper-right (offset from Front in both
        // axes).
        let front_pos = d.views[ids[0]].position;
        let iso_pos = d.views[ids[3]].position;
        assert!(iso_pos[0] > front_pos[0]);
        assert!(iso_pos[1] > front_pos[1]);
    }

    /// Sanity: every projected view actually has edges (we projected a
    /// real cube; the HLR pass should give us at least one visible
    /// edge per view).
    #[test]
    fn all_views_in_group_have_edges() {
        let cube = box_solid(40.0, 20.0, 10.0).unwrap();
        let mut d = Drawing::new(Sheet::a3_landscape("X", "Y", "Z"));
        let mut g = ProjectionGroup::new([100.0, 100.0], 1.0, Projection::ThirdAngle);
        g.include_isometric = true;
        g.build_into(&mut d, &cube).unwrap();
        for v in &d.views {
            assert!(
                !v.visible_edges.is_empty(),
                "{} view should have visible edges",
                v.kind.label(),
            );
        }
    }

    /// rebuild reruns edge extraction in place — change the solid,
    /// rebuild, the new edge bbox should match the new solid.
    #[test]
    fn rebuild_updates_views_in_place() {
        let small = box_solid(10.0, 10.0, 10.0).unwrap();
        let big = box_solid(40.0, 40.0, 40.0).unwrap();
        let mut d = Drawing::new(Sheet::a3_landscape("X", "Y", "Z"));
        let mut g = ProjectionGroup::new([150.0, 150.0], 1.0, Projection::ThirdAngle);
        g.build_into(&mut d, &small).unwrap();
        let front_bbox_small = d.views[g.view_indices[0]].bbox().unwrap();
        let w_small = front_bbox_small.1[0] - front_bbox_small.0[0];

        let errs = g.rebuild(&mut d, &big);
        assert!(errs.is_empty(), "rebuild errors: {errs:?}");
        let front_bbox_big = d.views[g.view_indices[0]].bbox().unwrap();
        let w_big = front_bbox_big.1[0] - front_bbox_big.0[0];
        assert!(
            w_big > w_small * 2.0,
            "front view should grow with new solid ({w_small} → {w_big})"
        );
    }

    /// Right edges of Top, Front, Right should be self-consistent for
    /// a known cube — verify Front's x-extent ≈ Top's x-extent (both
    /// project the cube's x-axis onto the drawing x), and Front's
    /// y-extent ≈ Right's y-extent (both project z onto drawing y).
    /// Catches bogus camera orientations.
    #[test]
    fn projection_group_view_extents_are_consistent_for_a_known_cube() {
        let cube = box_solid(50.0, 30.0, 20.0).unwrap(); // x=50, y=30, z=20
        let mut d = Drawing::new(Sheet::a3_landscape("X", "Y", "Z"));
        let mut g = ProjectionGroup::new([100.0, 100.0], 1.0, Projection::ThirdAngle);
        g.build_into(&mut d, &cube).unwrap();
        let front_bbox = d.views[g.view_indices[0]].bbox().unwrap();
        let top_bbox = d.views[g.view_indices[1]].bbox().unwrap();
        let right_bbox = d.views[g.view_indices[2]].bbox().unwrap();
        let fw = front_bbox.1[0] - front_bbox.0[0];
        let fh = front_bbox.1[1] - front_bbox.0[1];
        let tw = top_bbox.1[0] - top_bbox.0[0];
        let th = top_bbox.1[1] - top_bbox.0[1];
        let rw = right_bbox.1[0] - right_bbox.0[0];
        let rh = right_bbox.1[1] - right_bbox.0[1];
        // Front: world x→draw x (50), world z→draw y (20).
        assert!((fw - 50.0).abs() < 1e-3, "front width {fw}");
        assert!((fh - 20.0).abs() < 1e-3, "front height {fh}");
        // Top: world x→draw x (50), world y→draw y (30).
        assert!((tw - 50.0).abs() < 1e-3, "top width {tw}");
        assert!((th - 30.0).abs() < 1e-3, "top height {th}");
        // Right: world y→draw x (30), world z→draw y (20).
        assert!((rw - 30.0).abs() < 1e-3, "right width {rw}");
        assert!((rh - 20.0).abs() < 1e-3, "right height {rh}");
    }

    /// Bad scale → BadViewParameter error before any view is added.
    #[test]
    fn build_with_bad_scale_returns_bad_view_parameter() {
        let cube = box_solid(10.0, 10.0, 10.0).unwrap();
        let mut d = Drawing::new(Sheet::a3_landscape("X", "Y", "Z"));
        let mut g = ProjectionGroup::new([100.0, 100.0], 0.0, Projection::ThirdAngle);
        let e = g.build_into(&mut d, &cube).unwrap_err();
        match e {
            TechDrawError::BadViewParameter { name, .. } => assert_eq!(name, "scale"),
            other => panic!("wrong error: {other:?}"),
        }
        assert!(d.views.is_empty(), "no views added when scale is bad");
    }

    /// `label_positions` returns one entry per view in the group.
    #[test]
    fn label_positions_count_matches_view_count() {
        let mut g = ProjectionGroup::new([100.0, 100.0], 1.0, Projection::ThirdAngle);
        assert_eq!(g.label_positions(40.0, 30.0).len(), 3);
        g.include_isometric = true;
        assert_eq!(g.label_positions(40.0, 30.0).len(), 4);
    }

    /// Projection labels are distinct strings.
    #[test]
    fn projection_labels_are_distinct() {
        assert_ne!(
            Projection::FirstAngle.label(),
            Projection::ThirdAngle.label()
        );
    }
}
