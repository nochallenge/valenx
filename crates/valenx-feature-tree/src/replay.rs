//! Replay: walk the feature tree top-to-bottom and produce a final
//! [`valenx_cad::Solid`].
//!
//! Each feature is dispatched to its `ops::*::evaluate` function with
//! enough context to honour cross-feature references:
//!
//! - `Pad` and `Revolve` consume just a sketch.
//! - `Pocket` also needs the *preceding* solid (the boolean base).
//! - `Mirror` and the patterns reference an *earlier* feature by id
//!   and need access to its evaluated result, so they get a borrow
//!   of the `prior` results map.
//!
//! Suppressed features short-circuit: their slot in the results map is
//! `FeatureResult::Suppressed` (a crate-private enum variant) and they
//! do not advance `last_solid`.

use std::collections::HashMap;

use valenx_cad::Solid;
use valenx_spreadsheet::Spreadsheet;

use crate::feature::{Feature, FeatureId};
use crate::tree::FeatureTree;
use crate::FeatureError;

/// Per-feature evaluated result, built incrementally during replay.
///
/// Crate-private because callers outside the crate should reach for
/// [`replay`] which returns just the final solid; the per-feature
/// intermediates are an implementation detail of cross-feature
/// references (`Mirror`, patterns).
///
/// The `Solid` variant is ~200 bytes (truck `Solid` carries shells
/// with edge / face / vertex tables). The Suppressed variant is
/// zero-sized. Clippy's `large_enum_variant` lint flags that, but
/// boxing every result would balloon the heap traffic during replay
/// for a marginal stack-size win — replay holds a `HashMap` of these,
/// so it's already on the heap anyway.
#[allow(clippy::large_enum_variant)]
#[derive(Clone)]
pub(crate) enum FeatureResult {
    /// Feature produced a solid (most cases).
    #[allow(dead_code)]
    Solid(Solid),
    /// Feature was suppressed; pass-through for the next feature that
    /// needs to combine with "the previous solid".
    Suppressed,
}

/// Walk the tree top-to-bottom, evaluate each non-suppressed feature
/// in order, and return the final solid (the last solid-producing
/// feature). Returns `Ok(None)` if the tree is empty or every feature
/// is suppressed.
///
/// Calls through to [`replay_with_spreadsheet`] with an empty
/// [`Spreadsheet`] — Phase 16 features that hold
/// [`crate::feature::Value::Expression`] will error during evaluation
/// because the empty workbook has no cells to resolve against. Tree
/// authors who want to plug expressions in must call
/// [`FeatureTree::replay_with_spreadsheet`] directly with a populated
/// [`Spreadsheet`].
///
/// # Errors
///
/// Propagates any per-op evaluation error verbatim. The first failing
/// feature aborts the rest of the tree — partial-replay semantics are
/// deferred to a future task.
pub fn replay(tree: &FeatureTree) -> Result<Option<Solid>, FeatureError> {
    replay_with_spreadsheet(tree, &Spreadsheet::new())
}

/// Like [`replay`] but resolves
/// [`crate::feature::Value::Expression`] params against `ss`.
///
/// New in Phase 16. Callers that don't use parametric expressions
/// (the common case) can keep calling the original [`replay`].
pub fn replay_with_spreadsheet(
    tree: &FeatureTree,
    ss: &Spreadsheet,
) -> Result<Option<Solid>, FeatureError> {
    let mut results: HashMap<FeatureId, FeatureResult> = HashMap::new();
    let mut last_solid: Option<Solid> = None;
    for (idx, entry) in tree.features.iter().enumerate() {
        let id = FeatureId(idx);
        if entry.suppressed {
            results.insert(id, FeatureResult::Suppressed);
            continue;
        }
        let solid = dispatch(&entry.feature, tree, &results, last_solid.as_ref(), ss)?;
        results.insert(id, FeatureResult::Solid(solid.clone()));
        last_solid = Some(solid);
    }
    Ok(last_solid)
}

fn dispatch(
    feature: &Feature,
    tree: &FeatureTree,
    prior: &HashMap<FeatureId, FeatureResult>,
    last_solid: Option<&Solid>,
    ss: &Spreadsheet,
) -> Result<Solid, FeatureError> {
    match feature {
        Feature::Pad(p) => crate::ops::pad::evaluate(tree, p, ss),
        Feature::Pocket(p) => crate::ops::pocket::evaluate(tree, p, last_solid, ss),
        Feature::Revolve(p) => crate::ops::revolve::evaluate(tree, p, ss),
        Feature::Mirror(p) => crate::ops::mirror::evaluate(tree, p, prior),
        Feature::LinearPattern(p) => crate::ops::pattern_linear::evaluate(tree, p, prior),
        Feature::CircularPattern(p) => crate::ops::pattern_circular::evaluate(tree, p, prior),
        Feature::Fillet(p) => crate::ops::fillet::evaluate(tree, p, prior),
        Feature::Chamfer(p) => crate::ops::chamfer::evaluate(tree, p, prior),
        Feature::ImportedSolid(p) => crate::ops::imported_solid::evaluate(p),
        Feature::Hole(p) => crate::ops::hole::evaluate(tree, p, last_solid),
        Feature::Loft(p) => crate::ops::loft::evaluate(tree, p),
        Feature::Sweep(p) => crate::ops::sweep::evaluate(tree, p),
        Feature::Pipe(p) => crate::ops::pipe::evaluate(tree, p),
        Feature::Helix(p) => crate::ops::helix::evaluate(tree, p),
        Feature::MultiTransform(p) => crate::ops::multi_transform::evaluate(tree, p, prior),
        Feature::DraftAngle(p) => crate::ops::draft_angle::evaluate(tree, p, prior),
        Feature::Shell(p) => crate::ops::shell::evaluate(tree, p, prior),
        Feature::Thickness(p) => crate::ops::thickness::evaluate(tree, p, prior),
        Feature::BooleanHistory(p) => crate::ops::boolean_history::evaluate(tree, p, prior),
        Feature::ImportedAdvanced(p) => crate::ops::imported_advanced::evaluate(p),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::PadParams;

    #[test]
    fn empty_tree_returns_none() {
        let tree = FeatureTree::new();
        let result = replay(&tree).unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn suppressed_only_tree_returns_none() {
        // A single suppressed Pad — would otherwise hit the stub error
        // for Pad evaluation. Replay must short-circuit before
        // dispatching it.
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(valenx_sketch::Sketch::new());
        let f = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Pad",
        );
        tree.set_suppressed(f, true).unwrap();
        let result = replay(&tree).unwrap();
        assert!(result.is_none());
    }
}
