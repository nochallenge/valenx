//! Pocket evaluator — boolean-subtracts an extruded profile from the
//! preceding solid.
//!
//! A Pocket is the inverse of a Pad: same "sketch + depth → tool" recipe,
//! but instead of returning the tool we feed it to
//! [`valenx_cad::difference`] against the *base* solid that comes
//! through from the prior step of the feature tree. The base is the
//! `last_solid` argument forwarded by [`crate::replay::replay`]; if it's
//! `None` the pocket has nothing to cut into and we bail with a
//! `BadParameter` instead of fabricating a result.
//!
//! ## Stabbing overhang
//!
//! truck-shapeops's boolean operator refuses to subtract two solids when
//! they share a coincident face (e.g. both extruded from the same XY
//! plane). Since Pad and Pocket both extrude from z = 0, a naive
//! `difference` would always share at least one cap face and surface
//! [`valenx_cad::CadError::EmptyResult`].
//!
//! The standard CAD workaround is "stabbing": grow the cutting tool a
//! small epsilon **past the open face of the cut** and translate it so
//! the cut still starts at the working plane. We use
//! [`POCKET_STAB_EPSILON`] = `0.5` model units; that's larger than
//! [`valenx_cad::DEFAULT_BOOL_TOLERANCE`] (0.05) so shapeops's
//! tolerance-based edge merger doesn't collapse the cap planes back
//! into coincidence.
//!
//! The overhang is applied to the **open end only**. The blind end of
//! the cut (the pocket bottom) is a fresh interior face with no
//! coincidence to clear, so it must sit exactly at the requested
//! depth — overhanging it too would cut `epsilon` deeper than asked.
//! A *blind* pocket of depth `d` therefore removes exactly `d` of
//! material; a *through* pocket is specified by setting `d` beyond the
//! part (the conventional "through all" idiom), so its far end clears
//! the back face on its own.

use valenx_cad::{difference, Solid};
use valenx_spreadsheet::Spreadsheet;

use crate::feature::PocketParams;
use crate::tree::FeatureTree;
use crate::FeatureError;

/// Overshoot added to the pocket cutter past the **open** cap plane so
/// the tool's open-end cap never coincides with the base's face. Must
/// be larger than [`valenx_cad::DEFAULT_BOOL_TOLERANCE`] (0.05) or
/// shapeops's edge-merge step collapses the overhang back into
/// coincidence and the difference returns `EmptyResult`. The blind end
/// is *not* over-extended — see the module docs.
pub const POCKET_STAB_EPSILON: f64 = 0.5;

/// Evaluate a Pocket: build the cutting tool from the sketch + depth,
/// then subtract it from `base`.
///
/// `direction_positive == true` extrudes the tool along +Z; `false`
/// flips the sweep vector so the pocket cuts downwards (-Z) from the
/// working plane.
pub(crate) fn evaluate(
    tree: &FeatureTree,
    p: &PocketParams,
    base: Option<&Solid>,
    ss: &Spreadsheet,
) -> Result<Solid, FeatureError> {
    let base = base.ok_or(FeatureError::BadParameter {
        name: "pocket",
        reason: "pocket requires a base solid (must be preceded by a solid-producing feature)"
            .into(),
    })?;
    let depth = p
        .depth
        .resolve(ss)
        .map_err(|e| FeatureError::BadParameter {
            name: "depth",
            reason: format!("could not resolve depth expression: {e}"),
        })?;
    if !depth.is_finite() || depth.abs() < 1e-12 {
        return Err(FeatureError::BadParameter {
            name: "depth",
            reason: format!("must be nonzero and finite, got {depth}"),
        });
    }
    let sketch = tree.get_sketch(p.sketch)?;
    let signed_depth = if p.direction_positive { depth } else { -depth };

    // Stabbing: overhang the tool past the *open* end of the cut only.
    //
    // The pocket profile lives on the working plane (z = 0); a Pad
    // typically places the base solid's face flush with that plane.
    // The tool's *open*-end cap therefore coincides with a base face,
    // and `truck_shapeops` returns no solid for a flush-face
    // difference — so the tool is grown by one [`POCKET_STAB_EPSILON`]
    // *past the open face* and translated back, clearing the
    // coincidence.
    //
    // Crucially the overhang is applied to the **open end only**. An
    // earlier revision over-extended *both* ends, which made a blind
    // pocket of depth `d` actually cut `d + epsilon` deep — the far
    // overhang ate `epsilon` of material past the requested bottom.
    // The blind end of the cut is a fresh interior face (no
    // coincidence to clear), so it must sit exactly at `d`. For a
    // *through* pocket the user specifies `d` beyond the part (the
    // standard "through all" idiom) so the far end clears the back
    // face naturally — no far overhang needed there either.
    //
    // signed_depth ≥ 0 → tool extrudes +Z; open end at z = 0, so the
    // tool spans [-epsilon, depth]: extrude `depth + epsilon`, then
    // translate by -epsilon. signed_depth < 0 mirrors it.
    let (stabbed_depth, stab_shift) = if signed_depth >= 0.0 {
        (signed_depth + POCKET_STAB_EPSILON, -POCKET_STAB_EPSILON)
    } else {
        (signed_depth - POCKET_STAB_EPSILON, POCKET_STAB_EPSILON)
    };
    let tool = sketch.extrude(stabbed_depth)?;
    let tool = tool.translated(0.0, 0.0, stab_shift)?;

    // `difference` uses `valenx_cad::DEFAULT_BOOL_TOLERANCE` (0.05)
    // internally — the same default the standalone Part toolbox uses,
    // so user-visible behaviour stays consistent across surfaces.
    let result = difference(base, &tool).map_err(FeatureError::from)?;
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::{Feature, PadParams};
    use crate::replay::replay;

    /// 4 × 4 square base profile centred at the origin (corners at ±2).
    fn square_sketch(half: f64) -> valenx_sketch::Sketch {
        let mut s = valenx_sketch::Sketch::new();
        let a = s.add_point(-half, -half);
        let b = s.add_point(half, -half);
        let c = s.add_point(half, half);
        let d = s.add_point(-half, half);
        s.add_line(a, b).unwrap();
        s.add_line(b, c).unwrap();
        s.add_line(c, d).unwrap();
        s.add_line(d, a).unwrap();
        s
    }

    /// Build an 8-sided regular polygon centred at `(cx, cy)` with
    /// `radius`. Used in place of a circle for Pad / Pocket profiles
    /// because Phase 1's `extract_profile_lines` only handles
    /// straight-line entities (arc support is a Phase 2.5 follow-up).
    fn polygonal_circle(cx: f64, cy: f64, radius: f64, sides: u32) -> valenx_sketch::Sketch {
        use std::f64::consts::TAU;
        assert!(sides >= 3);
        let mut s = valenx_sketch::Sketch::new();
        let mut ids = Vec::with_capacity(sides as usize);
        for i in 0..sides {
            let theta = (i as f64 / sides as f64) * TAU;
            let x = cx + radius * theta.cos();
            let y = cy + radius * theta.sin();
            ids.push(s.add_point(x, y));
        }
        for i in 0..(sides as usize) {
            let next = (i + 1) % (sides as usize);
            s.add_line(ids[i], ids[next]).unwrap();
        }
        s
    }

    #[test]
    fn pocket_subtracts_polygonal_hole_from_padded_box() {
        // Box minus circular hole: Pad a 4x4 square 2 units tall, then
        // Pocket an 8-sided "circle" (radius 1) part-way into it (depth
        // 1 < box height 2 so the cut produces a blind pocket rather
        // than punching all the way through).
        //
        // The polygonal stand-in approximates a circle because Phase 1's
        // `extract_profile_lines` doesn't yet handle arcs — once arc
        // support lands (Phase 2.5) we can swap in `Sketch::add_circle`.
        //
        // The Pocket evaluator adds POCKET_STAB_EPSILON overhang on
        // both ends of the cutting tool, so the depth-1 pocket actually
        // extrudes a depth-(1 + 1.0) cutter translated by -0.5 — that
        // overhang is what lets truck-shapeops produce a clean result
        // instead of EmptyResult.
        let mut tree = FeatureTree::new();
        let base_s = tree.add_sketch(square_sketch(2.0));
        let hole_s = tree.add_sketch(polygonal_circle(0.0, 0.0, 1.0, 8));
        tree.add_feature(
            Feature::Pad(PadParams {
                sketch: base_s,
                depth: 2.0.into(),
                direction_positive: true,
            }),
            "Base Pad",
        );
        tree.add_feature(
            Feature::Pocket(PocketParams {
                sketch: hole_s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Blind Hole",
        );
        let solid = replay(&tree)
            .expect("replay succeeds")
            .expect("replay produces a solid");
        // The plain Pad of a square gives a 6-face prism. After
        // pocketing an 8-sided cylinder through it, the solid must
        // have strictly more faces (top + bottom keep their hole rim,
        // plus 8 new inner walls). Don't pin the exact count — truck
        // may simplify or split — but it MUST be more than 6.
        assert!(
            solid.faces() > 6,
            "pocketed box should have more than 6 faces, got {}",
            solid.faces()
        );
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.25).unwrap();
        // A plain box tessellates to ~12 triangles; the pocketed box
        // gains 16+ from the inner walls and the punched-through caps.
        assert!(
            mesh.total_elements() > 12,
            "pocketed tessellation should have more triangles than a plain box, got {}",
            mesh.total_elements()
        );
    }

    #[test]
    fn pocket_without_base_solid_errors() {
        // Pocket called with `base = None` — the standard situation
        // when Pocket is the first feature in a tree with no Pad ahead
        // of it. Must surface as BadParameter("pocket", ...) so the UI
        // can offer "add a Pad above this Pocket" instead of a generic
        // CAD failure.
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(polygonal_circle(0.0, 0.0, 1.0, 8));
        let params = PocketParams {
            sketch: s,
            depth: 1.0.into(),
            direction_positive: true,
        };
        let err = super::evaluate(&tree, &params, None, &Spreadsheet::new()).unwrap_err();
        match err {
            FeatureError::BadParameter { name, reason } => {
                assert_eq!(name, "pocket");
                assert!(
                    reason.contains("base solid"),
                    "expected reason to mention 'base solid', got {reason:?}"
                );
            }
            other => panic!("expected BadParameter, got {other:?}"),
        }
    }

    #[test]
    fn pocket_with_negative_direction_cuts_downward() {
        // Pocket with direction_positive = false flips the sweep so the
        // cutting tool extrudes along -Z from the working plane. To
        // overlap with the base in this configuration the base must
        // extend below z = 0, which means Pad-down: extrude depth 2
        // with direction_positive = false → base at z=-2..0. The
        // downward Pocket tool then sits at z=-1..0 (overlapping the
        // upper half of the base). The Pocket's stab-overhang carries
        // the right sign automatically so shapeops still gets a clean
        // disjoint setup.
        let mut tree = FeatureTree::new();
        let base_s = tree.add_sketch(square_sketch(2.0));
        let hole_s = tree.add_sketch(polygonal_circle(0.0, 0.0, 1.0, 8));
        tree.add_feature(
            Feature::Pad(PadParams {
                sketch: base_s,
                depth: 2.0.into(),
                direction_positive: false, // base z = -2..0
            }),
            "Base Pad Down",
        );
        tree.add_feature(
            Feature::Pocket(PocketParams {
                sketch: hole_s,
                depth: 1.0.into(),
                direction_positive: false, // tool z = -1..0
            }),
            "Down Pocket",
        );
        let solid = replay(&tree)
            .expect("replay succeeds")
            .expect("replay produces a solid");
        assert!(
            solid.faces() > 6,
            "down-pocketed box should have > 6 faces, got {}",
            solid.faces()
        );
        // Verify the result sits below z = 0 + epsilon (matches the
        // base's z range, not flipped to +Z).
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.25).unwrap();
        let max_z = mesh
            .nodes
            .iter()
            .map(|n| n.z)
            .fold(f64::NEG_INFINITY, f64::max);
        // Allow tessellation slop. Even with stab overhang the result
        // is clipped to the base, so max z should not exceed ~0.
        assert!(
            max_z <= 1e-3,
            "down pocket should keep the solid below z = 0, got max z = {max_z}"
        );
    }

    #[test]
    fn pocket_with_oversized_profile_is_handled_gracefully() {
        // Pocket profile is larger than the base's footprint (the
        // "hole" extends beyond the base in X and Y). Conceptually this
        // would erase the entire base — truck-shapeops may handle that
        // cleanly, may degenerate, or may surface EmptyResult; any of
        // those is acceptable, what we're checking is that the op
        // doesn't panic and that any error comes back through
        // FeatureError rather than tearing down the calling thread.
        //
        // Per Task 17 in the plan: "Pocket with sketch profile larger
        // than base footprint — result handled gracefully (truck may or
        // may not flag it; document)."
        let mut tree = FeatureTree::new();
        // Base is 1x1 at half = 0.5.
        let base_s = tree.add_sketch(square_sketch(0.5));
        // Tool is a 6x6 octagon centred at the origin — completely
        // engulfs the base in X and Y.
        let huge_s = tree.add_sketch(polygonal_circle(0.0, 0.0, 3.0, 8));
        tree.add_feature(
            Feature::Pad(PadParams {
                sketch: base_s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Small Base",
        );
        tree.add_feature(
            Feature::Pocket(PocketParams {
                sketch: huge_s,
                depth: 0.5.into(),
                direction_positive: true,
            }),
            "Oversized Pocket",
        );
        // Whatever happens, we must NOT panic. Either:
        //   (a) replay succeeds → either Some(solid) (truck preserved
        //       something) or None (impossible here since base is
        //       non-empty, but plausible after a future shapeops bump)
        //   (b) replay returns FeatureError::CadError wrapping shapeops's
        //       EmptyResult / similar
        // Both are documented graceful outcomes; the regression we are
        // guarding against is a process-killing assertion inside truck.
        match replay(&tree) {
            Ok(Some(_solid)) => {
                // Successful subtraction — fine, truck handled it.
            }
            Ok(None) => {
                // Replay returned no solid — unexpected for a non-empty
                // tree, but document & accept as graceful.
            }
            Err(FeatureError::CadError(_)) => {
                // shapeops surfaced an EmptyResult or similar — accept.
            }
            Err(other) => panic!(
                "oversized pocket should either succeed or surface a CadError, got {other:?}"
            ),
        }
    }

    #[test]
    fn pocket_inside_pad_pocket_pad_chain_advances_last_solid() {
        // Pad → Pocket → Pad chain. The replay shell forwards
        // `last_solid` to Pocket (so the pocket lands ON the first
        // pad) but Pad ignores `last_solid` entirely — it produces a
        // fresh extrusion from its own sketch. So a third Pad after
        // Pocket *replaces* the pocketed result in `last_solid` rather
        // than welding onto it. Union semantics live behind the Mirror
        // / Pattern ops; Pad-as-add is a future Phase 2.5 enhancement.
        //
        // What this test really pins down: the chain runs to
        // completion, the final Pad becomes the result of replay, and
        // its topology matches a plain extruded prism (not a leftover
        // shape from the middle of the chain).
        let mut tree = FeatureTree::new();
        let big_s = tree.add_sketch(square_sketch(2.0)); // 4x4 base
        let hole_s = tree.add_sketch(polygonal_circle(0.0, 0.0, 0.6, 8));
        let final_s = tree.add_sketch(square_sketch(1.0)); // 2x2 capping pad
        tree.add_feature(
            Feature::Pad(PadParams {
                sketch: big_s,
                depth: 2.0.into(),
                direction_positive: true,
            }),
            "Base Pad",
        );
        tree.add_feature(
            Feature::Pocket(PocketParams {
                sketch: hole_s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Middle Pocket",
        );
        tree.add_feature(
            Feature::Pad(PadParams {
                sketch: final_s,
                depth: 1.5.into(),
                direction_positive: true,
            }),
            "Final Pad",
        );
        let solid = replay(&tree)
            .expect("3-step chain replays")
            .expect("produces a final solid");
        // The final pad is a plain 2x2x1.5 prism — six faces, the same
        // as a `valenx_cad::box_solid` of equal dimensions. If the
        // chain accidentally fed the pocketed solid back in, the face
        // count would be much higher (8 walls + 2 caps + bottom hole
        // ring from the polygonal cutter).
        assert_eq!(
            solid.faces(),
            6,
            "final pad should be a plain 6-face prism, got {} faces (chain ordering issue?)",
            solid.faces()
        );
    }
}
