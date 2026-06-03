//! [`FeatureTree`] — the ordered list of features + the sketches they
//! reference.
//!
//! The tree owns both the sketches (so a save round-trip is
//! self-contained) and the feature entries (each with a suppression
//! flag and a display name).

use serde::{Deserialize, Serialize};
use valenx_cad::Solid;
use valenx_spreadsheet::Spreadsheet;

use crate::feature::{Feature, FeatureId, SketchRef};
use crate::FeatureError;

/// Entry in the tree — a feature plus its suppression flag and display
/// name. Keeping these in a wrapper means [`Feature`] stays free of UI
/// concerns; the tree decides what to render and what to skip.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TreeEntry {
    /// The underlying operation.
    pub feature: Feature,
    /// `true` = skip during replay (greyed out in UI, no solid produced).
    pub suppressed: bool,
    /// User-editable display name (e.g. `"Bottom plate Pad"`).
    pub name: String,
}

/// Ordered list of features + the sketches they reference.
///
/// `sketches` and `features` are independent: a sketch can exist without
/// any feature consuming it (the sketcher created it but the user hasn't
/// pad'd it yet), and a feature can reference an earlier feature
/// without touching any sketch (e.g. Mirror, patterns).
///
/// `PartialEq` is derived structurally so the host app can wrap a
/// `FeatureTree` in an `undo::History<FeatureTree>` for the Part
/// Design panel's `↶ ↷` controls. Sub-types carrying `f64` (lengths,
/// angles, vectors) use IEEE 754 semantics; the snapshot dedupe
/// silently fails on NaN-valued params, which the panel's drag-value
/// widgets prevent from happening in practice.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct FeatureTree {
    /// All sketches owned by this part. Indexed by [`SketchRef`].
    pub sketches: Vec<valenx_sketch::Sketch>,
    /// All features in document order. Indexed by [`FeatureId`].
    pub features: Vec<TreeEntry>,
}

impl FeatureTree {
    /// Create an empty tree (no sketches, no features).
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a sketch and return its new [`SketchRef`].
    pub fn add_sketch(&mut self, sketch: valenx_sketch::Sketch) -> SketchRef {
        let id = SketchRef(self.sketches.len());
        self.sketches.push(sketch);
        id
    }

    /// Append a feature and return its new [`FeatureId`]. Defaults to
    /// non-suppressed with the given display name.
    pub fn add_feature(&mut self, feature: Feature, name: impl Into<String>) -> FeatureId {
        let id = FeatureId(self.features.len());
        self.features.push(TreeEntry {
            feature,
            suppressed: false,
            name: name.into(),
        });
        id
    }

    /// Toggle the suppressed flag of a feature.
    ///
    /// # Errors
    ///
    /// Returns [`FeatureError::UnknownFeature`] when `id` does not index
    /// into [`Self::features`].
    pub fn set_suppressed(&mut self, id: FeatureId, suppressed: bool) -> Result<(), FeatureError> {
        self.features
            .get_mut(id.0)
            .ok_or(FeatureError::UnknownFeature(id.0))?
            .suppressed = suppressed;
        Ok(())
    }

    /// Remove a feature by id and return the [`Feature`] payload.
    ///
    /// Note: this re-indexes every later feature in the tree. Callers
    /// holding live [`FeatureId`]s for those entries will need to refresh
    /// them.
    ///
    /// # Errors
    ///
    /// Returns [`FeatureError::UnknownFeature`] when `id` does not index
    /// into [`Self::features`].
    pub fn delete_feature(&mut self, id: FeatureId) -> Result<Feature, FeatureError> {
        if id.0 >= self.features.len() {
            return Err(FeatureError::UnknownFeature(id.0));
        }
        Ok(self.features.remove(id.0).feature)
    }

    /// Borrow a sketch by reference.
    ///
    /// # Errors
    ///
    /// Returns [`FeatureError::UnknownSketch`] when `r` does not index
    /// into [`Self::sketches`].
    pub fn get_sketch(&self, r: SketchRef) -> Result<&valenx_sketch::Sketch, FeatureError> {
        self.sketches
            .get(r.0)
            .ok_or(FeatureError::UnknownSketch(r.0))
    }

    /// Validate every sketch owned by this tree (R33 H1).
    ///
    /// Each [`valenx_sketch::Sketch`] is run through
    /// [`valenx_sketch::Sketch::validate`], which rejects any entity
    /// carrying a variable handle past the end of its `vars`. A
    /// hand-edited / corrupt / version-skewed `.valenx` could otherwise
    /// deserialize a sketch with an out-of-range handle that panics
    /// "index out of bounds" the first time a feature consumes it
    /// (e.g. `Hole` indexing `sketch.vars[pt.x_var]`, or the solver's
    /// Jacobian assembly) during replay on the UI thread — which has no
    /// `catch_unwind`, so the whole app crashes.
    ///
    /// Called by the project-load path
    /// (`crate::persist::ValenxProject::from_ron`) so a bad document is
    /// rejected at load instead of crashing during replay.
    ///
    /// # Errors
    ///
    /// Returns [`FeatureError::SketchError`] wrapping the first
    /// [`valenx_sketch::SketchError::CorruptHandle`] encountered.
    pub fn validate(&self) -> Result<(), FeatureError> {
        for sketch in &self.sketches {
            sketch.validate()?;
        }
        Ok(())
    }

    /// Replay the tree against a [`Spreadsheet`] context — resolves
    /// any [`crate::feature::Value::Expression`] params via `ss`
    /// before dispatching to the per-op evaluators.
    ///
    /// New in Phase 16. Convenience wrapper around the free function
    /// [`crate::replay::replay_with_spreadsheet`].
    ///
    /// # Errors
    ///
    /// See [`crate::replay::replay`].
    pub fn replay_with_spreadsheet(&self, ss: &Spreadsheet) -> Result<Option<Solid>, FeatureError> {
        crate::replay::replay_with_spreadsheet(self, ss)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::{PadParams, PocketParams};

    fn pad(sketch: SketchRef) -> Feature {
        Feature::Pad(PadParams {
            sketch,
            depth: 1.0.into(),
            direction_positive: true,
        })
    }

    #[test]
    fn add_sketch_and_get_sketch_round_trip() {
        let mut tree = FeatureTree::new();
        let s0 = tree.add_sketch(valenx_sketch::Sketch::new());
        let s1 = tree.add_sketch(valenx_sketch::Sketch::new());
        assert_eq!(s0, SketchRef(0));
        assert_eq!(s1, SketchRef(1));
        assert!(tree.get_sketch(s0).is_ok());
        assert!(tree.get_sketch(s1).is_ok());
    }

    #[test]
    fn get_sketch_out_of_bounds_returns_unknown_sketch() {
        let tree = FeatureTree::new();
        let err = tree.get_sketch(SketchRef(42)).unwrap_err();
        assert_eq!(err.code(), "feature_tree.unknown_sketch");
    }

    #[test]
    fn add_feature_assigns_sequential_ids() {
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(valenx_sketch::Sketch::new());
        let f0 = tree.add_feature(pad(s), "Pad 1");
        let f1 = tree.add_feature(pad(s), "Pad 2");
        assert_eq!(f0, FeatureId(0));
        assert_eq!(f1, FeatureId(1));
        assert_eq!(tree.features.len(), 2);
        assert_eq!(tree.features[0].name, "Pad 1");
        assert!(!tree.features[0].suppressed);
    }

    #[test]
    fn set_suppressed_toggles_flag() {
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(valenx_sketch::Sketch::new());
        let f = tree.add_feature(pad(s), "Pad");
        tree.set_suppressed(f, true).unwrap();
        assert!(tree.features[0].suppressed);
        tree.set_suppressed(f, false).unwrap();
        assert!(!tree.features[0].suppressed);
    }

    #[test]
    fn set_suppressed_on_unknown_feature_errors() {
        let mut tree = FeatureTree::new();
        let err = tree.set_suppressed(FeatureId(99), true).unwrap_err();
        assert_eq!(err.code(), "feature_tree.unknown_feature");
    }

    #[test]
    fn delete_feature_removes_and_returns_payload() {
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(valenx_sketch::Sketch::new());
        let _f0 = tree.add_feature(pad(s), "First");
        let _f1 = tree.add_feature(
            Feature::Pocket(PocketParams {
                sketch: s,
                depth: 0.5.into(),
                direction_positive: false,
            }),
            "Second",
        );
        assert_eq!(tree.features.len(), 2);
        let removed = tree.delete_feature(FeatureId(0)).unwrap();
        assert!(matches!(removed, Feature::Pad(_)));
        assert_eq!(tree.features.len(), 1);
        // What was at index 1 is now at index 0.
        assert!(matches!(tree.features[0].feature, Feature::Pocket(_)));
    }

    #[test]
    fn delete_feature_out_of_bounds_errors() {
        let mut tree = FeatureTree::new();
        let err = tree.delete_feature(FeatureId(0)).unwrap_err();
        assert_eq!(err.code(), "feature_tree.unknown_feature");
    }
}
