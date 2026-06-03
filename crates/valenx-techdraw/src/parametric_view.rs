//! Parametric views — drawing views that re-generate themselves when
//! the source feature tree changes.
//!
//! A [`ParametricView`] is a small wrapper around a regular [`View`]
//! that remembers which [`FeatureId`] in the feature tree produced
//! the source solid, plus an `auto_update` flag. When the feature
//! tree replays and the resulting solid changes, calling
//! [`regenerate`] re-extracts visible / hidden edges in place.
//!
//! Why a wrapper instead of a field on [`View`]? Keeping it separate
//! avoids forcing every existing call site through a parametric link;
//! a hand-placed view (Phase 5) still works exactly as before.

use serde::{Deserialize, Serialize};

use crate::error::TechDrawError;
use crate::view::View;
use valenx_feature_tree::feature::FeatureId;

/// Marker that ties a [`View`] back to a feature in a feature tree.
///
/// `feature_id` is the producing feature in the source
/// [`valenx_feature_tree::FeatureTree`]. `auto_update` toggles whether
/// `Drawing::regenerate_all` re-runs this view's edge extraction on
/// each replay — users can pin a view to a specific snapshot by
/// flipping this off.
///
/// The view itself lives in the parent drawing; this struct only
/// carries the link metadata. The drawing layer joins
/// `ParametricView` to `View` by index (the `view_idx` field).
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ParametricView {
    /// Index of the [`View`] in the drawing's `views` vector.
    pub view_idx: usize,
    /// Feature in the source feature tree whose output solid we draw.
    /// Persisted so the link survives a serialise / deserialise round
    /// trip.
    pub feature_id: FeatureId,
    /// When `true`, [`crate::Drawing::regenerate_all`] re-runs
    /// [`regenerate`] for this view; when `false`, the view is frozen
    /// at its last manual generation.
    pub auto_update: bool,
}

impl ParametricView {
    /// Construct a new parametric link. `auto_update` defaults to
    /// `true` — callers that want a frozen snapshot flip it off after
    /// construction.
    pub fn new(view_idx: usize, feature_id: FeatureId) -> Self {
        Self {
            view_idx,
            feature_id,
            auto_update: true,
        }
    }
}

/// Regenerate a single view from `solid`. Thin wrapper around
/// [`View::generate`] — pulled out as a free function so callers that
/// already have a `&mut View` + `&Solid` don't need to construct a
/// [`ParametricView`] just to re-extract edges.
///
/// Returns [`TechDrawError::BadViewParameter`] when `view.scale` is
/// non-positive (same contract as [`View::generate`]).
pub fn regenerate(view: &mut View, solid: &valenx_cad::Solid) -> Result<(), TechDrawError> {
    view.generate(solid)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::view::ViewKind;
    use valenx_cad::primitives::box_solid;
    use valenx_feature_tree::feature::FeatureId;

    #[test]
    fn parametric_view_defaults_to_auto_update() {
        let pv = ParametricView::new(0, FeatureId(3));
        assert_eq!(pv.view_idx, 0);
        assert_eq!(pv.feature_id, FeatureId(3));
        assert!(pv.auto_update);
    }

    #[test]
    fn regenerate_populates_edges_for_cube() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let mut v = View::new(ViewKind::Front, 1.0, [0.0, 0.0]);
        regenerate(&mut v, &cube).unwrap();
        // Even for a coarse cube the front view should produce >= 4
        // visible edge segments (the four sides of the visible face).
        assert!(
            !v.visible_edges.is_empty(),
            "front view of cube should yield visible edges"
        );
    }

    #[test]
    fn regenerate_with_bad_scale_returns_bad_view_parameter() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let mut v = View::new(ViewKind::Front, -1.0, [0.0, 0.0]);
        let e = regenerate(&mut v, &cube).unwrap_err();
        match e {
            TechDrawError::BadViewParameter { name, .. } => assert_eq!(name, "scale"),
            other => panic!("wrong error: {other:?}"),
        }
    }
}
