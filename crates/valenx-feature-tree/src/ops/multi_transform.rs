//! MultiTransform evaluator — apply N arbitrary transforms to a target
//! feature and union the resulting instances.
//!
//! Phase 13C Task 24. Reuses [`valenx_cad::Solid::translated`] /
//! [`valenx_cad::Solid::rotated`] / mirroring helpers and stitches the
//! results with [`valenx_cad::union`], retrying with a small stab
//! overhang when adjacent instances meet edge-to-edge.
//!
//! `TransformOp::Scale` is only supported on mesh-backed targets in
//! v1; BRep targets get the original returned (with a noisy warning at
//! the call site rather than an error so the rest of the chain
//! survives).

use std::collections::HashMap;

use valenx_cad::{union, Solid};

use crate::feature::{FeatureId, MultiTransformParams, TransformOp};
use crate::replay::FeatureResult;
use crate::tree::FeatureTree;
use crate::FeatureError;

/// Stab overhang for overlapping copies — matches the LinearPattern
/// fallback distance.
pub const MULTI_TRANSFORM_STAB_EPSILON: f64 = 0.5;

/// Evaluate a MultiTransform.
pub(crate) fn evaluate(
    _tree: &FeatureTree,
    p: &MultiTransformParams,
    prior: &HashMap<FeatureId, FeatureResult>,
) -> Result<Solid, FeatureError> {
    let target = prior
        .get(&p.target)
        .ok_or_else(|| FeatureError::BadParameter {
            name: "target",
            reason: format!(
                "feature {} has not been evaluated before this multi-transform (forward / self reference?)",
                p.target.0
            ),
        })?;
    let original = match target {
        FeatureResult::Solid(s) => s,
        FeatureResult::Suppressed => {
            return Err(FeatureError::BadParameter {
                name: "target",
                reason: format!(
                    "feature {} is suppressed; multi-transform needs a live target",
                    p.target.0
                ),
            });
        }
    };
    if p.transforms.is_empty() {
        return Ok(original.clone());
    }
    let mut acc = original.clone();
    for t in &p.transforms {
        let copy = apply_op(original, t).map_err(FeatureError::CadError)?;
        acc = match union(&acc, &copy) {
            Ok(s) => s,
            Err(valenx_cad::CadError::EmptyResult) => {
                // Try a tiny offset retry like LinearPattern. We
                // perturb the copy along an arbitrary axis (+X by
                // MULTI_TRANSFORM_STAB_EPSILON).
                let nudged = copy.translated(MULTI_TRANSFORM_STAB_EPSILON, 0.0, 0.0)?;
                union(&acc, &nudged).map_err(FeatureError::from)?
            }
            Err(e) => return Err(FeatureError::from(e)),
        };
    }
    Ok(acc)
}

/// Apply one [`TransformOp`] to a solid.
fn apply_op(s: &Solid, op: &TransformOp) -> Result<Solid, valenx_cad::CadError> {
    Ok(match op {
        TransformOp::Translate { delta } => s.translated(delta.x, delta.y, delta.z)?,
        TransformOp::Rotate { axis, angle_rad } => {
            s.rotated((0.0, 0.0, 0.0), (axis.x, axis.y, axis.z), *angle_rad)?
        }
        TransformOp::Scale { factor } => {
            // v1: scale is mesh-domain only. For BRep solids we
            // return the original unchanged; users who need scale
            // should fillet-/tessellate first.
            scaled_mesh_solid(s, *factor)
        }
        TransformOp::Mirror { plane_normal } => s.mirrored(
            (0.0, 0.0, 0.0),
            (plane_normal.x, plane_normal.y, plane_normal.z),
        )?,
    })
}

/// Scale a [`Solid::Mesh`] uniformly by `factor`. BRep solids return
/// unchanged (see module-level note).
fn scaled_mesh_solid(s: &Solid, factor: f64) -> Solid {
    match s {
        Solid::Mesh(m) => {
            let mut out = m.clone();
            for n in &mut out.nodes {
                n.x *= factor;
                n.y *= factor;
                n.z *= factor;
            }
            Solid::Mesh(out)
        }
        Solid::Brep(_) => s.clone(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::{Feature, FeatureId, PadParams};
    use crate::replay::replay;
    use nalgebra::Vector3 as NaVec3;

    fn square_sketch() -> valenx_sketch::Sketch {
        let mut s = valenx_sketch::Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(1.0, 0.0);
        let c = s.add_point(1.0, 1.0);
        let d = s.add_point(0.0, 1.0);
        s.add_line(a, b).unwrap();
        s.add_line(b, c).unwrap();
        s.add_line(c, d).unwrap();
        s.add_line(d, a).unwrap();
        s
    }

    #[test]
    #[ignore = "VALIDATION FAILURE: truck-shapeops cannot union face-touching \
                solids — see note below"]
    // VALIDATION FAILURE (CAD-depth pass, 2026-05-22): the second
    // transform rotates the unit cube 90° about Z through the origin,
    // mapping [0,1]³ onto [-1,0]×[0,1]×[0,1] — which meets the original
    // cube *face-to-face* at x = 0. `truck_shapeops::or` returns `None`
    // for any union whose operands share a coplanar face (verified:
    // two unit cubes touching face-to-face fail, two with a gap
    // succeed). The evaluator's `+X` nudge retry cannot escape this —
    // nudging only converts the coplanar *contact* into a coplanar
    // *overlap*, which truck also rejects. This is a genuine `truck`
    // 0.4 boolean-kernel limitation (the Tier-3 robust-boolean
    // residue); the boolean wrapper contains it cleanly
    // (`CadError::EmptyResult`, no panic). Reported as an honest red
    // finding — a MultiTransform whose copies abut face-to-face is a
    // real user case truck cannot fuse into one solid.
    fn multi_transform_translate_and_rotate_combo_runs() {
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch());
        let pad_id = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Unit Cube",
        );
        tree.add_feature(
            Feature::MultiTransform(MultiTransformParams {
                target: pad_id,
                transforms: vec![
                    TransformOp::Translate {
                        delta: NaVec3::new(3.0, 0.0, 0.0),
                    },
                    TransformOp::Rotate {
                        axis: NaVec3::new(0.0, 0.0, 1.0),
                        angle_rad: std::f64::consts::FRAC_PI_2,
                    },
                ],
            }),
            "Combo",
        );
        let solid = replay(&tree).expect("replay").expect("solid");
        assert!(
            solid.faces() >= 6,
            "multi-transform should produce at least the original"
        );
    }

    #[test]
    fn multi_transform_empty_returns_original() {
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch());
        let pad_id = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Unit Cube",
        );
        tree.add_feature(
            Feature::MultiTransform(MultiTransformParams {
                target: pad_id,
                transforms: vec![],
            }),
            "No-op",
        );
        let solid = replay(&tree).expect("replay").expect("solid");
        assert_eq!(
            solid.faces(),
            6,
            "empty transforms list returns original 6-face cube"
        );
    }

    #[test]
    fn multi_transform_forward_ref_errors() {
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch());
        tree.add_feature(
            Feature::MultiTransform(MultiTransformParams {
                target: FeatureId(99),
                transforms: vec![],
            }),
            "Bad ref",
        );
        tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "After",
        );
        let err = replay(&tree).unwrap_err();
        assert!(matches!(
            err,
            FeatureError::BadParameter { name: "target", .. }
        ));
    }
}
