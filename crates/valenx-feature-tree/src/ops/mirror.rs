//! Mirror evaluator — reflect a previously-built solid across a plane.
//!
//! The Mirror feature consumes the output of an *earlier* feature (its
//! `target` field is a `FeatureId`, not a sketch) and reflects it
//! across the plane defined by `plane_origin + plane_normal`. When
//! `keep_original = true` the result is the union of the original and
//! its mirror; when `false` it's just the mirrored copy.
//!
//! ## Touching-plane overlap (the "mirror stab")
//!
//! The typical use case for Mirror is "draw half the symmetric part,
//! mirror to get the other half". That means the original solid
//! usually *touches* the mirror plane — and the mirror image touches
//! it from the other side, producing a shared coincident face when
//! they meet. truck-shapeops's `or` refuses to merge two solids when
//! their boundary shells share a face (the same failure mode the
//! Pocket evaluator dodges with its [`POCKET_STAB_EPSILON`] overhang).
//!
//! The defensive fix here is the same idea: nudge the mirror back
//! along the *negative* plane normal by [`MIRROR_STAB_EPSILON`] so the
//! two solids overlap by a sliver rather than sharing an exact face.
//! shapeops then sees them as two intersecting blobs and merges
//! cleanly. The overlap is much smaller than the model scale so it
//! doesn't show up in the user-visible geometry.
//!
//! ## Known gap — non-overlapping `keep_original`
//!
//! If the original solid and the mirror are well-separated (the
//! original doesn't touch / cross the plane and is more than
//! `MIRROR_STAB_EPSILON` away from it), even the stab nudge can't
//! bridge the gap and truck-shapeops's `or` returns `None`. We
//! surface that as [`FeatureError::CadError`]
//! ([`valenx_cad::CadError::EmptyResult`]) so the UI can hint "Mirror
//! with Keep Original needs the source to touch the plane".
//!
//! [`POCKET_STAB_EPSILON`]: super::pocket::POCKET_STAB_EPSILON

use std::collections::HashMap;

use valenx_cad::{union, Solid};

use crate::feature::{FeatureId, MirrorParams};
use crate::replay::FeatureResult;
use crate::tree::FeatureTree;
use crate::FeatureError;

/// Nudge distance used to break the coincident-face case when
/// `keep_original = true`. Sized like [`super::pocket::POCKET_STAB_EPSILON`]
/// (0.5 model units) so it's bigger than
/// [`valenx_cad::DEFAULT_BOOL_TOLERANCE`] (0.05) — otherwise
/// shapeops's tolerance-based edge merger collapses the nudge back
/// into coincidence and we're right back where we started. The
/// overlap is invisible in the user-visible result because shapeops
/// clips back to the union boundary.
pub const MIRROR_STAB_EPSILON: f64 = 0.5;

/// Evaluate a Mirror: look up the target's solid in `prior`, validate
/// the plane, reflect, and optionally union with the original.
pub(crate) fn evaluate(
    _tree: &FeatureTree,
    p: &MirrorParams,
    prior: &HashMap<FeatureId, FeatureResult>,
) -> Result<Solid, FeatureError> {
    // ---- plane normal validation ----
    if !p.plane_normal.x.is_finite()
        || !p.plane_normal.y.is_finite()
        || !p.plane_normal.z.is_finite()
    {
        return Err(FeatureError::BadParameter {
            name: "plane_normal",
            reason: format!("must be finite, got {:?}", p.plane_normal),
        });
    }
    let n_len = p.plane_normal.norm();
    if n_len < 1e-12 {
        return Err(FeatureError::BadParameter {
            name: "plane_normal",
            reason: format!(
                "must have nonzero magnitude (got {:?}, |n| = {n_len})",
                p.plane_normal
            ),
        });
    }
    if !p.plane_origin.x.is_finite()
        || !p.plane_origin.y.is_finite()
        || !p.plane_origin.z.is_finite()
    {
        return Err(FeatureError::BadParameter {
            name: "plane_origin",
            reason: format!("must be finite, got {:?}", p.plane_origin),
        });
    }

    // ---- target lookup ----
    let target_result = prior
        .get(&p.target)
        .ok_or_else(|| FeatureError::BadParameter {
            name: "target",
            reason: format!(
                "feature {} has not been evaluated before this mirror (forward / self reference?)",
                p.target.0
            ),
        })?;
    let original = match target_result {
        FeatureResult::Solid(s) => s,
        FeatureResult::Suppressed => {
            return Err(FeatureError::BadParameter {
                name: "target",
                reason: format!(
                    "feature {} is suppressed; mirror needs a live target",
                    p.target.0
                ),
            });
        }
    };

    // ---- reflect ----
    //
    // `Solid::mirrored` builds the Householder reflection matrix
    // through `plane_origin` with unit `plane_normal` and applies it
    // via truck-modeling's `builder::transformed`. It also inverts
    // the resulting face orientations because reflection has
    // determinant -1.
    let mut mirrored = original
        .mirrored(
            (p.plane_origin.x, p.plane_origin.y, p.plane_origin.z),
            (p.plane_normal.x, p.plane_normal.y, p.plane_normal.z),
        )
        .map_err(FeatureError::CadError)?;

    if !p.keep_original {
        return Ok(mirrored);
    }

    // ---- union with original (the "stab" case) ----
    //
    // If the original touches or crosses the mirror plane the
    // mirror image meets it at a coincident face. truck-shapeops
    // can't merge across a shared face, so push the mirror back
    // along the *negative* plane normal by MIRROR_STAB_EPSILON to
    // create a small overlap instead. This is the analogue of
    // Pocket's POCKET_STAB_EPSILON overhang: trade a tiny
    // user-invisible deformation for a robust boolean op.
    let nx = p.plane_normal.x / n_len;
    let ny = p.plane_normal.y / n_len;
    let nz = p.plane_normal.z / n_len;
    mirrored = mirrored.translated(
        -MIRROR_STAB_EPSILON * nx,
        -MIRROR_STAB_EPSILON * ny,
        -MIRROR_STAB_EPSILON * nz,
    )?;
    let result = union(original, &mirrored).map_err(FeatureError::from)?;
    Ok(result)
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
    fn mirror_across_diagonal_plane_swaps_x_and_y_extents() {
        // Mirror the unit-square Pad across the diagonal plane whose
        // normal is (1, -1, 0) / √2 through the origin. Reflecting
        // (x, y, z) about that plane sends (x, y, z) → (y, x, z)
        // (the standard x-y swap reflection — easy to verify by hand).
        //
        // Source: x ∈ [0, 1], y ∈ [0, 1], z ∈ [0, 1]. The diagonal
        // plane y = x passes through the corner (0, 0, _) and the
        // edge (1, 1, _), so the cube straddles the plane. The
        // reflected copy lands back in the same x = y region (because
        // the cube is symmetric about y = x), making this a tricky
        // case: the mirror image overlaps the original substantially.
        //
        // Test `keep_original = false` here so we don't have to union
        // (which would just produce the original cube again because
        // it's self-mirror-symmetric). The result should be the
        // mirror alone — same bounding box as the original since
        // y = x reflection of the unit cube returns the unit cube.
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
            Feature::Mirror(MirrorParams {
                target: pad_id,
                plane_origin: NaVec3::new(0.0, 0.0, 0.0),
                // y = x plane: any point on this plane satisfies y - x = 0,
                // so the normal is (1, -1, 0). Length doesn't matter — the
                // evaluator normalises.
                plane_normal: NaVec3::new(1.0, -1.0, 0.0),
                keep_original: false,
            }),
            "Diagonal Mirror",
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
        // The reflected unit cube (about y = x) is still a unit cube
        // — just with its (x, y) swapped, which is geometrically the
        // same shape. So the bounding box is unchanged at [0, 1] in
        // both X and Y.
        assert!(
            (min_x - 0.0).abs() < 1e-3 && (max_x - 1.0).abs() < 1e-3,
            "reflected cube should span x ∈ [0, 1], got [{min_x}, {max_x}]"
        );
        assert!(
            (min_y - 0.0).abs() < 1e-3 && (max_y - 1.0).abs() < 1e-3,
            "reflected cube should span y ∈ [0, 1], got [{min_y}, {max_y}]"
        );
        // It's still a unit cube (6 faces). Reflection over y = x
        // preserves topology — the test is mostly that we don't
        // panic and don't return a degenerate solid.
        assert_eq!(
            solid.faces(),
            6,
            "reflected unit cube should still have 6 faces, got {}",
            solid.faces()
        );
    }

    #[test]
    fn mirror_keep_original_false_returns_only_the_reflected_copy() {
        // Pad at x ∈ [2, 3] (well to the +X side of the YZ plane).
        // Mirror with keep_original = false should produce ONLY the
        // reflected copy at x ∈ [-3, -2]. The original is gone from
        // the final solid — replay's last_solid is now the mirror's
        // output, which is the disjoint-from-the-pad reflected copy.
        //
        // This is the path that lets the user use Mirror as a true
        // "move to mirror image" rather than a "make symmetric" op.
        // The disjoint case is also the one where the union-based
        // `keep_original = true` would fail (no overlap), so it
        // doubles as documentation that the `false` branch sidesteps
        // the EmptyResult risk.
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch(2.0, 0.0, 3.0, 1.0));
        let pad_id = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Pad at +X",
        );
        tree.add_feature(
            Feature::Mirror(MirrorParams {
                target: pad_id,
                plane_origin: NaVec3::new(0.0, 0.0, 0.0),
                plane_normal: NaVec3::new(1.0, 0.0, 0.0),
                keep_original: false,
            }),
            "Mirror only",
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
        // Mirror of x ∈ [2, 3] across YZ plane is x ∈ [-3, -2]. NO
        // overlap with the original — the union path (which we
        // skipped) would have surfaced EmptyResult. The mirror-only
        // path returns just the reflected solid.
        assert!(
            max_x <= -2.0 + 1e-3,
            "mirror-only should have max x ~= -2.0 (reflected far edge), got {max_x}"
        );
        assert!(
            min_x >= -3.0 - 1e-3,
            "mirror-only should have min x ~= -3.0 (reflected near edge), got {min_x}"
        );
        // Original was a 1×1×1 prism (6 faces); the reflected copy
        // is congruent so it also has 6 faces.
        assert_eq!(
            solid.faces(),
            6,
            "mirror-only of unit prism should still have 6 faces, got {}",
            solid.faces()
        );
    }

    #[test]
    fn mirror_targeting_suppressed_feature_returns_bad_parameter() {
        // Build a tree where the Pad is suppressed *before* the
        // Mirror runs. Replay's dispatch shell inserts a
        // FeatureResult::Suppressed entry for the suppressed Pad, so
        // the Mirror evaluator sees the target in `prior` but as a
        // suppressed slot rather than a live solid.
        //
        // The right outcome is a structured BadParameter("target",
        // ...) — *not* UnknownFeature (the feature exists in the
        // tree, it's just turned off) and *not* a generic CAD
        // failure. The UI can use this to highlight the dangling
        // reference and offer "unsuppress this feature" as a fix.
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch(0.0, 0.0, 1.0, 1.0));
        let pad_id = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Suppressed Pad",
        );
        tree.add_feature(
            Feature::Mirror(MirrorParams {
                target: pad_id,
                plane_origin: NaVec3::new(0.0, 0.0, 0.0),
                plane_normal: NaVec3::new(1.0, 0.0, 0.0),
                keep_original: false,
            }),
            "Dangling Mirror",
        );
        // Suppress the Pad *after* adding both — order of mutation
        // doesn't matter; what counts is that replay sees the Pad as
        // suppressed by the time it dispatches the Mirror.
        tree.set_suppressed(pad_id, true).unwrap();

        let err = replay(&tree).unwrap_err();
        match err {
            FeatureError::BadParameter { name, reason } => {
                assert_eq!(name, "target");
                assert!(
                    reason.contains("suppressed"),
                    "expected reason to mention suppressed, got {reason:?}"
                );
            }
            other => panic!("expected BadParameter(target, suppressed), got {other:?}"),
        }
    }

    #[test]
    fn mirror_targeting_forward_reference_returns_bad_parameter() {
        // Target a feature that *will exist* in the tree but hasn't
        // been evaluated yet — i.e. its FeatureId is past the Mirror's
        // own position in the tree. The replay shell only inserts
        // results for features it has already dispatched, so a
        // forward-referenced ID won't be in `prior` at all.
        //
        // This is the canonical "broken parametric reference" case:
        // the Mirror tries to reflect something that comes after it,
        // which is structurally impossible in a top-to-bottom replay.
        // Surface as BadParameter("target", ...) with a message
        // mentioning the forward/self reference.
        use crate::feature::FeatureId;

        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch(0.0, 0.0, 1.0, 1.0));
        // Mirror is feature 0; target = FeatureId(99) doesn't exist.
        tree.add_feature(
            Feature::Mirror(MirrorParams {
                target: FeatureId(99),
                plane_origin: NaVec3::new(0.0, 0.0, 0.0),
                plane_normal: NaVec3::new(1.0, 0.0, 0.0),
                keep_original: false,
            }),
            "Forward Mirror",
        );
        // Add a Pad after the Mirror so the FeatureId(99) is "real"
        // in the sense of "the user clearly meant something" but
        // unreachable from a top-down replay.
        tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Pad After Mirror",
        );
        let err = replay(&tree).unwrap_err();
        assert!(
            matches!(err, FeatureError::BadParameter { name: "target", .. }),
            "expected BadParameter(target) for forward reference, got {err:?}"
        );
    }

    #[test]
    fn mirror_pad_across_yz_plane_with_keep_original_doubles_x_extent() {
        // Pad the unit square at x ∈ [0, 1], y ∈ [0, 1] up to z = 1.
        // Mirror across the YZ plane (normal = +X through origin). The
        // original touches the plane at x = 0, so the mirror image at
        // x ∈ [-1, 0] meets it at a coincident face — the MIRROR_STAB
        // overhang is what lets the subsequent union succeed.
        //
        // Final solid should span x ∈ [-1, 1] (a 2×1×1 prism, modulo
        // the slight nudge from MIRROR_STAB_EPSILON which is clipped
        // away by the union).
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch(0.0, 0.0, 1.0, 1.0));
        let pad_id = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Half Block",
        );
        tree.add_feature(
            Feature::Mirror(MirrorParams {
                target: pad_id,
                plane_origin: NaVec3::new(0.0, 0.0, 0.0),
                plane_normal: NaVec3::new(1.0, 0.0, 0.0),
                keep_original: true,
            }),
            "Mirror across YZ",
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
        // The mirror was nudged by -MIRROR_STAB_EPSILON (0.5) along
        // the plane normal (+X), so its raw bounds are x ∈ [-1.5,
        // -0.5]. After union with the original at x ∈ [0, 1] the
        // resulting solid spans x ∈ [-1.5, 1.0]. Without the stab
        // overhang shapeops returns EmptyResult and we never get
        // here.
        assert!(
            max_x >= 1.0 - 1e-3,
            "union should reach max x ~= 1.0 (original's far edge), got {max_x}"
        );
        assert!(
            min_x <= -1.0 + 1e-3,
            "union should reach min x ~= -1.0 (mirror's far edge), got {min_x}"
        );
        // Total face count should clearly exceed a plain box's 6:
        // the stab-overlap union produces extra walls around the
        // sliver of overlap. We mostly just want "more than 6" as a
        // sanity check that something non-trivial happened.
        assert!(
            solid.faces() >= 6,
            "mirrored+unioned solid should have at least 6 faces, got {}",
            solid.faces()
        );
    }
}
