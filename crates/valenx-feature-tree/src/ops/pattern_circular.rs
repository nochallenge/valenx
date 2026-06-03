//! CircularPattern evaluator — rotate-and-union N instances of an
//! earlier feature's solid around an axis.
//!
//! The CircularPattern feature consumes the output of an *earlier*
//! feature (its `target` field is a [`FeatureId`], not a sketch). It
//! iterates `i in 0..count` and at each step rotates the target's
//! solid about `(axis_origin, axis_direction)` by `angle_step * i`,
//! then unions every instance into one final solid.
//!
//! ## Step angle: full vs partial sweep
//!
//! There are two natural interpretations of "N copies around a
//! `total_angle` sweep":
//!
//! - **Full circle (`total_angle ≈ 2π`):** the last copy lands on top
//!   of the original, so we want `count` evenly spaced copies and the
//!   step is `total_angle / count`. The wrap-around makes the
//!   "endpoint" implicit.
//! - **Partial sweep (`total_angle < 2π`):** the user typically wants
//!   the last copy to land *at* the sweep endpoint, not just before
//!   it. With `count` copies the step is `total_angle / (count - 1)`.
//!
//! The evaluator picks between these two formulas based on whether
//! `total_angle` is within 1e-6 of `2π` (TAU). This matches the
//! plan's spec and is the convention used by most parametric CAD
//! systems (FreeCAD, SolidWorks).
//!
//! ## Stab-overlap for touching instances
//!
//! As with [`LinearPattern`], consecutive rotated copies can meet at
//! a coincident face when the target straddles the rotation axis.
//! We use the same EmptyResult-fallback strategy: try the straight
//! union first, and on EmptyResult retry with a small rotation nudge
//! ([`CIRCULAR_PATTERN_STAB_EPSILON_RAD`], 0.01 radians ≈ 0.57°)
//! that breaks the coincident-face case. The overlap is much smaller
//! than the model scale, so the user-visible silhouette is unchanged.
//!
//! [`LinearPattern`]: super::pattern_linear

use std::collections::HashMap;
use std::f64::consts::TAU;

use valenx_cad::{union, Solid};

use crate::feature::{CircularPatternParams, FeatureId};
use crate::replay::FeatureResult;
use crate::tree::FeatureTree;
use crate::FeatureError;

/// Rotational nudge applied when shapeops returns EmptyResult on a
/// circular-pattern union step — typically the symptom of two
/// rotated copies sharing a coincident face. 0.01 rad ≈ 0.57° is
/// much smaller than any sweep we'd resolve at sketch tolerance and
/// well above [`valenx_cad::DEFAULT_BOOL_TOLERANCE`]'s angular
/// equivalent.
pub const CIRCULAR_PATTERN_STAB_EPSILON_RAD: f64 = 0.01;

/// Tolerance for detecting "full circle" sweep. If `|total_angle - TAU| < this`,
/// the step formula is `total_angle / count` (endpoint wraps around);
/// otherwise it's `total_angle / (count - 1)` (last copy lands at the
/// sweep endpoint).
const FULL_CIRCLE_EPS: f64 = 1e-6;

/// Evaluate a CircularPattern: look up the target's solid in `prior`,
/// validate parameters, build `count` rotated copies, then union
/// them all into a single solid.
pub(crate) fn evaluate(
    _tree: &FeatureTree,
    p: &CircularPatternParams,
    prior: &HashMap<FeatureId, FeatureResult>,
) -> Result<Solid, FeatureError> {
    // ---- count validation ----
    if p.count == 0 {
        return Err(FeatureError::BadParameter {
            name: "count",
            reason: "must be >= 1 (count includes the original instance)".into(),
        });
    }

    // ---- axis validation ----
    if !p.axis_direction.x.is_finite()
        || !p.axis_direction.y.is_finite()
        || !p.axis_direction.z.is_finite()
    {
        return Err(FeatureError::BadParameter {
            name: "axis_direction",
            reason: format!("must be finite, got {:?}", p.axis_direction),
        });
    }
    let ax_len = p.axis_direction.norm();
    if ax_len < 1e-12 {
        return Err(FeatureError::BadParameter {
            name: "axis_direction",
            reason: format!(
                "must have nonzero magnitude (got {:?}, |axis| = {ax_len})",
                p.axis_direction
            ),
        });
    }
    if !p.axis_origin.x.is_finite() || !p.axis_origin.y.is_finite() || !p.axis_origin.z.is_finite()
    {
        return Err(FeatureError::BadParameter {
            name: "axis_origin",
            reason: format!("must be finite, got {:?}", p.axis_origin),
        });
    }

    // ---- total_angle validation ----
    if !p.total_angle.is_finite() {
        return Err(FeatureError::BadParameter {
            name: "total_angle",
            reason: format!("must be finite, got {}", p.total_angle),
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
                    "feature {} is suppressed; circular pattern needs a live target",
                    p.target.0
                ),
            });
        }
    };

    // Single-instance pattern is just the original (no rotation, no
    // union needed) — avoids no-op shapeops work.
    if p.count == 1 {
        return Ok(original.clone());
    }

    // ---- step angle: full circle vs partial sweep ----
    //
    // Full circle: `count` evenly spaced copies, step = total_angle /
    // count (so the wrap-around closes the ring without doubling-up
    // on the start position).
    //
    // Partial sweep: last copy lands AT the sweep endpoint, step =
    // total_angle / (count - 1).
    let is_full_circle = (p.total_angle.abs() - TAU).abs() < FULL_CIRCLE_EPS;
    let step_angle = if is_full_circle {
        p.total_angle / (p.count as f64)
    } else {
        p.total_angle / ((p.count - 1) as f64)
    };

    // ---- build & union all instances ----
    //
    // Pass non-normalised axis to `Solid::rotated` — truck's
    // `builder::rotated` handles normalisation internally.
    let axis_origin = (p.axis_origin.x, p.axis_origin.y, p.axis_origin.z);
    let axis_dir = (p.axis_direction.x, p.axis_direction.y, p.axis_direction.z);

    let mut acc = original.clone();
    for i in 1..p.count {
        let angle = step_angle * (i as f64);
        let copy = original.rotated(axis_origin, axis_dir, angle)?;
        match union(&acc, &copy) {
            Ok(s) => acc = s,
            Err(valenx_cad::CadError::EmptyResult) => {
                // Coincident-face failure mode — retry this step with
                // a small angular nudge so the two rotated copies
                // overlap by a sliver instead of sharing a face. The
                // sub-degree nudge is invisible in user-visible output.
                let nudged_angle = angle - CIRCULAR_PATTERN_STAB_EPSILON_RAD;
                let nudged = original.rotated(axis_origin, axis_dir, nudged_angle)?;
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
    fn circular_pattern_6_cubes_around_z_at_radius_5_spans_full_circle() {
        // Pad a unit-sized cube at x ∈ [4.5, 5.5], y ∈ [-0.5, 0.5],
        // z ∈ [0, 1] — i.e. its centre is at (5, 0, _), radius 5
        // from the Z-axis. Circular-pattern 6 copies around +Z with
        // total_angle = 2π (full circle).
        //
        // The 6 copies sit at angles 0°, 60°, 120°, 180°, 240°, 300°.
        // The bounding box should span x ∈ [-5.5, 5.5] and y ∈
        // [-5.5, 5.5] — i.e. roughly ±(5 + 0.5) on each axis.
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch(4.5, -0.5, 5.5, 0.5));
        let pad_id = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Cube at +X radius 5",
        );
        tree.add_feature(
            Feature::CircularPattern(CircularPatternParams {
                target: pad_id,
                axis_origin: NaVec3::new(0.0, 0.0, 0.0),
                axis_direction: NaVec3::new(0.0, 0.0, 1.0),
                count: 6,
                total_angle: TAU,
            }),
            "6 around Z",
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
        let min_y = mesh.nodes.iter().map(|n| n.y).fold(f64::INFINITY, f64::min);
        let max_y = mesh
            .nodes
            .iter()
            .map(|n| n.y)
            .fold(f64::NEG_INFINITY, f64::max);
        // 6 unit cubes at radius 5 around the Z-axis. The cubes sit
        // at angles 0°, 60°, 120°, 180°, 240°, 300° — no cube lies
        // exactly on +Y or -Y. The X-extreme cubes (at 0° and 180°)
        // have their outer-radial corners at x = ±5.5. The Y extremes
        // come from the cubes at 60°/120° (and 240°/300°): for the
        // original cube's outer-radial corner (5.5, 0.5) rotated 60°
        // → y = 5.5*sin60° + 0.5*cos60° = 5.013. So max y ≈ 5.013,
        // min y ≈ -5.013. The plan called for "X and Y both span -6
        // to +6" as a rough sanity check; the actual bounds are
        // tighter and depend on count.
        assert!(
            (5.3..=5.6).contains(&max_x),
            "max x should be ~5.5 (cube at angle 0°'s far edge), got {max_x}"
        );
        assert!(
            (-5.6..=-5.3).contains(&min_x),
            "min x should be ~-5.5 (cube at angle 180°'s far edge), got {min_x}"
        );
        // No cube lies exactly on Y — the Y extreme is the far corner
        // of the 60°/120° cubes, ≈ 5.013.
        assert!(
            (max_y - 5.013).abs() < 0.1,
            "max y should be ~5.013 (cube at angle 60°'s far corner), got {max_y}"
        );
        assert!(
            (min_y - -5.013).abs() < 0.1,
            "min y should be ~-5.013 (cube at angle 240°'s far corner), got {min_y}"
        );
    }

    #[test]
    fn circular_pattern_partial_90_with_4_instances_endpoint_lands_at_90deg() {
        // Partial sweep test: 4 copies over 90° total_angle. The
        // partial-sweep formula divides by (count - 1) = 3 so the
        // step is 30°, putting copies at 0°, 30°, 60°, 90° — the
        // last copy lands exactly at the sweep endpoint.
        //
        // Pad a thin block at x ∈ [4.5, 5.5], y ∈ [-0.1, 0.1] —
        // narrow in Y so its rotated extents are dominated by the
        // radial position. The copy at 90° should sit centred at
        // (0, 5, _) with extents x ∈ [-0.1, 0.1], y ∈ [4.5, 5.5].
        //
        // Two outermost extents to verify:
        //   - max x ≈ 5.5: cube at angle 0° has its far-radial edge
        //     at x = 5.5.
        //   - max y ≈ 5.5: cube at angle 90° has its far-radial edge
        //     at y = 5.5 (this is the partial-sweep endpoint and is
        //     the headline contract — full-circle formula would put
        //     it at 22.5° instead and miss y = 5.5 entirely).
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch(4.5, -0.1, 5.5, 0.1));
        let pad_id = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Thin Block at +X",
        );
        tree.add_feature(
            Feature::CircularPattern(CircularPatternParams {
                target: pad_id,
                axis_origin: NaVec3::new(0.0, 0.0, 0.0),
                axis_direction: NaVec3::new(0.0, 0.0, 1.0),
                count: 4,
                total_angle: std::f64::consts::FRAC_PI_2, // 90°
            }),
            "Quarter Fan",
        );
        let solid = replay(&tree)
            .expect("replay succeeds")
            .expect("replay produces a solid");
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.25).unwrap();
        let max_x = mesh
            .nodes
            .iter()
            .map(|n| n.x)
            .fold(f64::NEG_INFINITY, f64::max);
        let max_y = mesh
            .nodes
            .iter()
            .map(|n| n.y)
            .fold(f64::NEG_INFINITY, f64::max);
        // Original cube at angle 0° contributes max x = 5.5.
        assert!(
            (max_x - 5.5).abs() < 0.05,
            "max x should be ~5.5 (cube at angle 0°'s far edge), got {max_x}"
        );
        // Cube at angle 90° (the sweep endpoint) contributes max y =
        // 5.5. If the evaluator wrongly used the full-circle formula
        // (step = 90 / 4 = 22.5°), the last copy would land at 67.5°
        // and max y would be < 5.5 — closer to 5.5 * sin(67.5) ≈
        // 5.08. The 0.05 tolerance rules that out cleanly.
        assert!(
            (max_y - 5.5).abs() < 0.05,
            "max y should be ~5.5 (cube at endpoint angle 90°'s far edge); \
             got {max_y} — partial-sweep step formula may be wrong"
        );
    }

    #[test]
    fn circular_pattern_off_axis_origin_translates_ring_centre() {
        // Rotation axis offset from the world origin. The original
        // cube must be placed relative to the *offset* axis, and the
        // resulting ring's centre must be at the offset axis position
        // — not at the world origin.
        //
        // Setup:
        //   - Axis through (10, 10, 0) along +Z.
        //   - Cube at x ∈ [12.5, 13.5], y ∈ [9.5, 10.5] — i.e. the
        //     cube's centre is at (13, 10, _), which is 3 units along
        //     +X from the axis position (10, 10).
        //   - 4 copies, total_angle = 2π (full circle, step = 90°).
        //
        // Copies sit at world positions:
        //   i=0:   centre = (10+3, 10  , _) = (13, 10, _)
        //   i=90°: centre = (10  , 10+3, _) = (10, 13, _)
        //   i=180: centre = (10-3, 10  , _) = ( 7, 10, _)
        //   i=270: centre = (10  , 10-3, _) = (10,  7, _)
        //
        // Bounding box should be x ∈ [7-0.5, 13+0.5] = [6.5, 13.5]
        // and y ∈ [7-0.5, 13+0.5] = [6.5, 13.5]. If the evaluator
        // wrongly used the world origin instead of axis_origin the
        // result would be centred at (0, 0) and span x ∈ [-13.5,
        // 13.5] — a wildly different bounding box.
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch(12.5, 9.5, 13.5, 10.5));
        let pad_id = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Cube near (13, 10)",
        );
        tree.add_feature(
            Feature::CircularPattern(CircularPatternParams {
                target: pad_id,
                axis_origin: NaVec3::new(10.0, 10.0, 0.0),
                axis_direction: NaVec3::new(0.0, 0.0, 1.0),
                count: 4,
                total_angle: TAU,
            }),
            "4 around offset axis",
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
        let min_y = mesh.nodes.iter().map(|n| n.y).fold(f64::INFINITY, f64::min);
        let max_y = mesh
            .nodes
            .iter()
            .map(|n| n.y)
            .fold(f64::NEG_INFINITY, f64::max);
        // The ring's bounding box: the cube at angle 0° contributes
        // max x = 13.5 (its far-radial edge); the cube at 180° gives
        // min x = 6.5; symmetric for Y. Allow 0.05 slop.
        assert!(
            (max_x - 13.5).abs() < 0.05,
            "off-axis ring's max x should be ~13.5 (cube at 0° around axis 10,10); got {max_x} — axis_origin ignored?"
        );
        assert!(
            (min_x - 6.5).abs() < 0.05,
            "off-axis ring's min x should be ~6.5 (cube at 180° around axis 10,10); got {min_x}"
        );
        assert!(
            (max_y - 13.5).abs() < 0.05,
            "off-axis ring's max y should be ~13.5 (cube at 90° around axis 10,10); got {max_y}"
        );
        assert!(
            (min_y - 6.5).abs() < 0.05,
            "off-axis ring's min y should be ~6.5 (cube at 270° around axis 10,10); got {min_y}"
        );
    }

    #[test]
    fn circular_pattern_count_1_is_just_the_original() {
        // count = 1 short-circuits the rotation loop entirely — no
        // shapeops union, no rotation matrix applied, the result is
        // byte-identical to the target's solid. Same contract as
        // LinearPattern's count = 1 fast path.
        //
        // We pin this with a face-count check (a Pad of a unit
        // square is a 6-face prism; count = 1 pattern thereof should
        // also be 6) and a bounding-box check (cube at radius 5
        // stays at x ∈ [4.5, 5.5] — no rotation applied).
        //
        // The face-count check is the more rigorous half: a no-op
        // union (acc = original; for _ in 1..1 {}) would also give
        // 6 faces, but might tolerance-wobble the vertices. The
        // short-circuit in the evaluator avoids both shapeops and
        // any matrix construction, so the output is the cloned
        // input.
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch(4.5, -0.5, 5.5, 0.5));
        let pad_id = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Unit Cube at +X radius 5",
        );
        tree.add_feature(
            Feature::CircularPattern(CircularPatternParams {
                target: pad_id,
                axis_origin: NaVec3::new(0.0, 0.0, 0.0),
                axis_direction: NaVec3::new(0.0, 0.0, 1.0),
                count: 1,
                total_angle: TAU,
            }),
            "Singleton Pattern",
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
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.25).unwrap();
        let min_x = mesh.nodes.iter().map(|n| n.x).fold(f64::INFINITY, f64::min);
        let max_x = mesh
            .nodes
            .iter()
            .map(|n| n.x)
            .fold(f64::NEG_INFINITY, f64::max);
        // Cube should be at exactly x ∈ [4.5, 5.5] — no rotation,
        // no shapeops, no tolerance wobble.
        assert!(
            (min_x - 4.5).abs() < 1e-6,
            "count = 1 should preserve min x = 4.5, got {min_x}"
        );
        assert!(
            (max_x - 5.5).abs() < 1e-6,
            "count = 1 should preserve max x = 5.5, got {max_x}"
        );
    }
}
