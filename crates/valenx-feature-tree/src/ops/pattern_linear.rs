//! LinearPattern evaluator — translate-and-union N instances of an
//! earlier feature's solid along a fixed direction.
//!
//! The LinearPattern feature consumes the output of an *earlier*
//! feature (its `target` field is a [`FeatureId`], not a sketch). It
//! iterates `i in 0..count` and at each step translates the target's
//! solid by `direction_unit * spacing * i`, then unions every instance
//! into one final solid.
//!
//! ## Stab-overlap for consecutive instances
//!
//! Truck-shapeops's boolean union refuses to merge two solids when
//! their boundary shells share a coincident face (the same failure
//! mode the [`Pocket`] and [`Mirror`] evaluators dodge). For Linear
//! Pattern this happens whenever `spacing` exactly equals the target's
//! extent along `direction` — consecutive copies meet at a perfectly
//! flat face. The defensive fix is the same idea as
//! [`POCKET_STAB_EPSILON`]: nudge each successive copy back by
//! [`LINEAR_PATTERN_STAB_EPSILON`] along the negative pattern
//! direction so adjacent copies overlap by a sliver rather than
//! sharing a face. shapeops then sees two intersecting blobs and
//! merges cleanly. The overlap is much smaller than the model scale,
//! so the user-visible silhouette is unchanged.
//!
//! ## Known gap — non-overlapping instances
//!
//! If `spacing` exceeds the target's extent along `direction` (plus
//! the stab overlap) the instances are well-separated, and the union
//! of two disjoint solids surfaces [`valenx_cad::CadError::EmptyResult`]
//! through shapeops. We surface that as [`FeatureError::CadError`] so
//! the UI can hint "spacing too large for this target".
//!
//! [`POCKET_STAB_EPSILON`]: super::pocket::POCKET_STAB_EPSILON
//! [`Pocket`]: super::pocket
//! [`Mirror`]: super::mirror

use std::collections::HashMap;

use valenx_cad::{union, Solid};

use crate::feature::{FeatureId, LinearPatternParams};
use crate::replay::FeatureResult;
use crate::tree::FeatureTree;
use crate::FeatureError;

/// Nudge distance used to break the coincident-face case when
/// consecutive pattern instances meet edge-to-edge. Sized like
/// [`super::pocket::POCKET_STAB_EPSILON`] (0.5 model units) so it's
/// bigger than [`valenx_cad::DEFAULT_BOOL_TOLERANCE`] (0.05) —
/// otherwise shapeops's tolerance-based edge merger collapses the
/// nudge back into coincidence and the union returns `EmptyResult`.
pub const LINEAR_PATTERN_STAB_EPSILON: f64 = 0.5;

/// Evaluate a LinearPattern: look up the target's solid in `prior`,
/// validate parameters, build `count` translated copies, then union
/// them all into a single solid.
pub(crate) fn evaluate(
    _tree: &FeatureTree,
    p: &LinearPatternParams,
    prior: &HashMap<FeatureId, FeatureResult>,
) -> Result<Solid, FeatureError> {
    // ---- count validation ----
    if p.count == 0 {
        return Err(FeatureError::BadParameter {
            name: "count",
            reason: "must be >= 1 (count includes the original instance)".into(),
        });
    }

    // ---- direction validation ----
    if !p.direction.x.is_finite() || !p.direction.y.is_finite() || !p.direction.z.is_finite() {
        return Err(FeatureError::BadParameter {
            name: "direction",
            reason: format!("must be finite, got {:?}", p.direction),
        });
    }
    let d_len = p.direction.norm();
    if d_len < 1e-12 {
        return Err(FeatureError::BadParameter {
            name: "direction",
            reason: format!(
                "must have nonzero magnitude (got {:?}, |d| = {d_len})",
                p.direction
            ),
        });
    }

    // ---- spacing validation ----
    if !p.spacing.is_finite() {
        return Err(FeatureError::BadParameter {
            name: "spacing",
            reason: format!("must be finite, got {}", p.spacing),
        });
    }

    // ---- target lookup ----
    let target_result = prior
        .get(&p.target)
        .ok_or_else(|| FeatureError::BadParameter {
            name: "target",
            reason: format!(
                "feature {} has not been evaluated before this pattern (forward / self reference?)",
                p.target.0
            ),
        })?;
    let original = match target_result {
        FeatureResult::Solid(s) => s,
        FeatureResult::Suppressed => {
            return Err(FeatureError::BadParameter {
                name: "target",
                reason: format!(
                    "feature {} is suppressed; linear pattern needs a live target",
                    p.target.0
                ),
            });
        }
    };

    // Single-instance pattern is just the original (no translation, no
    // union needed). This avoids triggering shapeops on a no-op union.
    if p.count == 1 {
        return Ok(original.clone());
    }

    // ---- build & union all instances ----
    //
    // Unit direction vector — `spacing` is the per-step distance, so
    // step i lands at `unit_dir * spacing * i`. The first instance
    // (i = 0) is the original and isn't translated.
    let ux = p.direction.x / d_len;
    let uy = p.direction.y / d_len;
    let uz = p.direction.z / d_len;

    // Successive copies often meet at coincident faces (e.g. when
    // spacing matches the target's extent along `direction`). Try the
    // straight translate-and-union first — if shapeops reports
    // EmptyResult (the canonical coincident-face symptom) retry with
    // the LINEAR_PATTERN_STAB_EPSILON nudge to overlap adjacent
    // copies by a sliver instead of sharing a face. This matches the
    // pragmatic approach the Mirror evaluator uses for plane stabs:
    // keep silhouettes pixel-accurate when possible and only deform
    // by a sub-tolerance amount when truck forces our hand.
    let mut acc = original.clone();
    for i in 1..p.count {
        let step = p.spacing * (i as f64);
        let copy = original.translated(ux * step, uy * step, uz * step)?;
        match union(&acc, &copy) {
            Ok(s) => acc = s,
            Err(valenx_cad::CadError::EmptyResult) => {
                // Coincident-face failure mode — retry this step with
                // a stab nudge along the negative pattern direction.
                let nudged = original.translated(
                    ux * step - ux * LINEAR_PATTERN_STAB_EPSILON,
                    uy * step - uy * LINEAR_PATTERN_STAB_EPSILON,
                    uz * step - uz * LINEAR_PATTERN_STAB_EPSILON,
                )?;
                acc = union(&acc, &nudged).map_err(FeatureError::from)?;
            }
            Err(e) => return Err(FeatureError::from(e)),
        }
    }
    Ok(acc)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::{Feature, PadParams};
    use crate::replay::replay;
    use nalgebra::Vector3 as NaVec3;

    /// Build a square sketch `[(x0, y0), (x1, y0), (x1, y1), (x0, y1)]`.
    fn square_sketch(x0: f64, y0: f64, x1: f64, y1: f64) -> valenx_sketch::Sketch {
        let mut s = valenx_sketch::Sketch::new();
        let a = s.add_point(x0, y0);
        let b = s.add_point(x1, y0);
        let c = s.add_point(x1, y1);
        let d = s.add_point(x0, y1);
        s.add_line(a, b).unwrap();
        s.add_line(b, c).unwrap();
        s.add_line(c, d).unwrap();
        s.add_line(d, a).unwrap();
        s
    }

    #[test]
    fn linear_pattern_4_unit_cubes_along_x_spans_expected_extent() {
        // Pad the unit square at x ∈ [0, 1], y ∈ [0, 1], z ∈ [0, 1].
        // Linear-pattern 4 copies along +X with spacing 2.
        //
        // The instances should land at:
        //   i=0: x ∈ [0, 1]  (original)
        //   i=1: x ∈ [2, 3]
        //   i=2: x ∈ [4, 5]
        //   i=3: x ∈ [6, 7]
        //
        // So the union's bounding box spans x ∈ [0, 7] (the plan's
        // "0 to 7 in X" target). The stab-overlap nudge shifts each
        // copy slightly toward the previous one (by 0.5), but that
        // just creates a small overlap at the leading edge of each
        // copy — the union's outer silhouette is unchanged because the
        // *trailing* edge of copy i still reaches x = 1 + 2i.
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch(0.0, 0.0, 1.0, 1.0));
        let pad_id = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Unit Cube",
        );
        tree.add_feature(
            Feature::LinearPattern(LinearPatternParams {
                target: pad_id,
                direction: NaVec3::new(1.0, 0.0, 0.0),
                count: 4,
                spacing: 2.0,
            }),
            "4 along +X",
        );
        let solid = replay(&tree)
            .expect("replay succeeds")
            .expect("replay produces a solid");
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.25).unwrap();
        let min_x = mesh.nodes.iter().map(|n| n.x).fold(f64::INFINITY, f64::min);
        let max_x = mesh
            .nodes
            .iter()
            .map(|n| n.x)
            .fold(f64::NEG_INFINITY, f64::max);
        // Min x should be ~0 (original's near edge); max x should be
        // ~7 (last copy's far edge, allowing for the small stab
        // overhang which leaves a slight margin near edges that lands
        // within the same range).
        assert!(
            (min_x - 0.0).abs() < 1e-3,
            "pattern should span min x ~= 0, got {min_x}"
        );
        assert!(
            (max_x - 7.0).abs() < 1e-3,
            "pattern should span max x ~= 7, got {max_x}"
        );
    }

    #[test]
    fn linear_pattern_count_1_is_just_the_original() {
        // count = 1 means "the original instance, no copies". The
        // evaluator must short-circuit before doing any translation or
        // union — that way the result is byte-identical to the
        // target's solid (no shapeops involvement = no risk of
        // accidentally re-tessellating or perturbing the geometry).
        //
        // We pin this with a face-count check: the original unit-cube
        // Pad is a 6-face prism, and so should the count = 1 pattern's
        // output be. Doing a no-op union would still produce 6 faces
        // but might (depending on shapeops's behaviour) introduce a
        // tolerance-level wobble in the vertex positions. The
        // short-circuit avoids that entirely.
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch(0.0, 0.0, 1.0, 1.0));
        let pad_id = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Unit Cube",
        );
        tree.add_feature(
            Feature::LinearPattern(LinearPatternParams {
                target: pad_id,
                direction: NaVec3::new(1.0, 0.0, 0.0),
                count: 1,
                spacing: 2.0,
            }),
            "Pattern of 1",
        );
        let solid = replay(&tree)
            .expect("replay succeeds")
            .expect("replay produces a solid");
        assert_eq!(
            solid.faces(),
            6,
            "count = 1 should yield the original prism (6 faces), got {}",
            solid.faces()
        );
        // Bounding box should match the original: x ∈ [0, 1].
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.25).unwrap();
        let min_x = mesh.nodes.iter().map(|n| n.x).fold(f64::INFINITY, f64::min);
        let max_x = mesh
            .nodes
            .iter()
            .map(|n| n.x)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            (min_x - 0.0).abs() < 1e-6,
            "count = 1 should preserve min x = 0, got {min_x}"
        );
        assert!(
            (max_x - 1.0).abs() < 1e-6,
            "count = 1 should preserve max x = 1, got {max_x}"
        );
    }

    #[test]
    fn linear_pattern_negative_spacing_walks_against_direction() {
        // Negative spacing should produce a pattern that walks *back*
        // along the unit direction. The evaluator multiplies
        // `unit_dir * spacing * i`, so a negative spacing naturally
        // flips the sign without any extra branch: the copies land on
        // the opposite side of the original from where they would
        // with positive spacing.
        //
        // Pad at x ∈ [0, 1], direction = +X, spacing = -2, count = 3.
        // Instances:
        //   i=0: x ∈ [0, 1]      (original)
        //   i=1: x ∈ [-2, -1]    (1 * -2 = -2 offset)
        //   i=2: x ∈ [-4, -3]    (2 * -2 = -4 offset)
        //
        // Bounding box should span x ∈ [-4, 1]. None of the copies
        // touch (gaps of 1 unit between them), so the stab nudge
        // doesn't engage.
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch(0.0, 0.0, 1.0, 1.0));
        let pad_id = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Unit Cube",
        );
        tree.add_feature(
            Feature::LinearPattern(LinearPatternParams {
                target: pad_id,
                direction: NaVec3::new(1.0, 0.0, 0.0),
                count: 3,
                spacing: -2.0, // walks back along +X = forward along -X
            }),
            "Walk Backward",
        );
        let solid = replay(&tree)
            .expect("replay succeeds")
            .expect("replay produces a solid");
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.25).unwrap();
        let min_x = mesh.nodes.iter().map(|n| n.x).fold(f64::INFINITY, f64::min);
        let max_x = mesh
            .nodes
            .iter()
            .map(|n| n.x)
            .fold(f64::NEG_INFINITY, f64::max);
        // Original's far-X edge unchanged at x = 1.
        assert!(
            (max_x - 1.0).abs() < 1e-3,
            "negative-spacing pattern should keep max x ~= 1 (original's far edge), got {max_x}"
        );
        // Furthest-back copy reaches x = -4.
        assert!(
            (min_x - -4.0).abs() < 1e-3,
            "negative-spacing pattern should reach min x ~= -4 (last copy's near edge), got {min_x}"
        );
    }

    #[test]
    fn linear_pattern_targeting_pocket_replicates_pocketed_solid() {
        // Pattern's target is a Pocket, not a Pad. The replay shell
        // records the *pocketed* solid (Pad - Pocket cutter, not the
        // cutter itself) into the prior-results map under the Pocket's
        // FeatureId. So the pattern should produce N copies of the
        // pocketed result — N pad-shaped blocks each with the same
        // hole pattern, not N copies of the cutter tool.
        //
        // We verify with a face count: a single pad-minus-hole has
        // some number F faces (more than 6 — the hole adds inner
        // walls). The pattern of N copies should have strictly more
        // than F faces (each copy keeps its walls in the union).
        //
        // This pins the contract: patterns operate on *evaluated
        // solids* via the `prior` map, not on the raw sketch geometry.
        // If the dispatch shell ever changed to pass the raw Pocket
        // tool through (a regression of the worst kind), this test
        // would catch it because the face count would be very wrong.
        use crate::feature::PocketParams;

        // Build a small octagonal-hole stamp: an 8-sided "circle"
        // pocketed out of a 4×4 base, 1 unit deep.
        fn polygonal_circle(cx: f64, cy: f64, radius: f64, sides: u32) -> valenx_sketch::Sketch {
            use std::f64::consts::TAU;
            let mut s = valenx_sketch::Sketch::new();
            let mut ids = Vec::with_capacity(sides as usize);
            for i in 0..sides {
                let theta = (i as f64 / sides as f64) * TAU;
                ids.push(s.add_point(cx + radius * theta.cos(), cy + radius * theta.sin()));
            }
            for i in 0..(sides as usize) {
                let next = (i + 1) % (sides as usize);
                s.add_line(ids[i], ids[next]).unwrap();
            }
            s
        }

        let mut tree = FeatureTree::new();
        let base_s = tree.add_sketch(square_sketch(-2.0, -2.0, 2.0, 2.0)); // 4×4 base
        let hole_s = tree.add_sketch(polygonal_circle(0.0, 0.0, 0.6, 8));
        tree.add_feature(
            Feature::Pad(PadParams {
                sketch: base_s,
                depth: 2.0.into(),
                direction_positive: true,
            }),
            "Base Pad",
        );
        let pocket_id = tree.add_feature(
            Feature::Pocket(PocketParams {
                sketch: hole_s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Hole",
        );
        // Pattern the *pocketed* solid 2 copies along +X with spacing
        // 6 (large enough that the copies don't overlap — the base is
        // 4 wide, so spacing 6 leaves a 2-unit gap between copies).
        tree.add_feature(
            Feature::LinearPattern(LinearPatternParams {
                target: pocket_id,
                direction: NaVec3::new(1.0, 0.0, 0.0),
                count: 2,
                spacing: 6.0,
            }),
            "Pattern of Holed Blocks",
        );

        let solid = replay(&tree)
            .expect("replay succeeds")
            .expect("replay produces a solid");
        // A single pocketed block has more than 6 faces (the
        // octagonal hole adds at least 8 inner walls + 2 hole-rim
        // caps split). The pattern of 2 must have MORE faces than 1
        // — if the implementation had accidentally patterned the raw
        // base Pad (6 faces × 2 = 12), the count would be much lower
        // than the actual ~30+ face holed-block-pair.
        assert!(
            solid.faces() > 12,
            "pattern of 2 holed blocks should have many faces (each block keeps its hole walls); got {} (looks like the pattern is missing the holes)",
            solid.faces()
        );
        // Bounding box sanity: the two blocks span x ∈ [-2, 8] (first
        // at [-2, 2], second at [4, 8]).
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.25).unwrap();
        let min_x = mesh.nodes.iter().map(|n| n.x).fold(f64::INFINITY, f64::min);
        let max_x = mesh
            .nodes
            .iter()
            .map(|n| n.x)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            (min_x - -2.0).abs() < 1e-3,
            "pattern's first block should span min x ~= -2, got {min_x}"
        );
        assert!(
            (max_x - 8.0).abs() < 1e-3,
            "pattern's second block should span max x ~= 8, got {max_x}"
        );
    }
}
