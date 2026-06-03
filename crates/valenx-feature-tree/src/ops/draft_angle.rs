//! DraftAngle evaluator — tilt selected faces of a target solid.
//!
//! Phase 13D Task 28. v1: tessellation-based approximation. We
//! tessellate the target into a mesh, identify the vertices belonging
//! to the listed triangle indices, and rotate each one about the
//! neutral plane axis by `draft_angle_deg`.
//!
//! Mesh-backed output — see Loft.

use std::collections::{HashMap, HashSet};

use nalgebra::Vector3;
use valenx_cad::{solid_to_mesh, Solid};
use valenx_mesh::element::ElementType;

use crate::feature::{DraftAngleParams, FeatureId};
use crate::replay::FeatureResult;
use crate::tree::FeatureTree;
use crate::FeatureError;

/// Evaluate a DraftAngle.
pub(crate) fn evaluate(
    _tree: &FeatureTree,
    p: &DraftAngleParams,
    prior: &HashMap<FeatureId, FeatureResult>,
) -> Result<Solid, FeatureError> {
    let target = lookup_target(prior, p.target)?;
    let raw = solid_to_mesh(target, valenx_cad::DEFAULT_TESS_TOLERANCE)?;

    let normal_len = p.neutral_plane_normal.norm();
    if normal_len < 1e-12 {
        return Err(FeatureError::BadParameter {
            name: "neutral_plane_normal",
            reason: "must have nonzero magnitude".into(),
        });
    }
    let normal = p.neutral_plane_normal / normal_len;
    let angle = p.draft_angle_deg.to_radians();
    let (cos_a, sin_a) = (angle.cos(), angle.sin());

    let mut mesh = raw;
    // Build a set of vertex indices touched by the selected triangles.
    let mut affected: HashSet<u32> = HashSet::new();
    let mut tri_cursor = 0usize;
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            tri_cursor += block.count();
            continue;
        }
        let tris = block.count();
        for local in 0..tris {
            let global = tri_cursor + local;
            if p.face_indices.contains(&global) {
                let base = local * 3;
                for offset in 0..3 {
                    affected.insert(block.connectivity[base + offset]);
                }
            }
        }
        tri_cursor += tris;
    }

    // Rotate each affected vertex about the world origin around `normal`.
    // For v1 the neutral plane goes through the origin; future work
    // could add a `neutral_origin` field.
    for &idx in &affected {
        let pos = mesh.nodes[idx as usize];
        let rotated = rotate_about_axis(&pos, &normal, cos_a, sin_a);
        mesh.nodes[idx as usize] = rotated;
    }
    mesh.recompute_stats();
    Ok(Solid::from_mesh(mesh))
}

fn lookup_target(
    prior: &HashMap<FeatureId, FeatureResult>,
    id: FeatureId,
) -> Result<&Solid, FeatureError> {
    let res = prior.get(&id).ok_or_else(|| FeatureError::BadParameter {
        name: "target",
        reason: format!(
            "feature {} has not been evaluated before this draft (forward / self reference?)",
            id.0
        ),
    })?;
    match res {
        FeatureResult::Solid(s) => Ok(s),
        FeatureResult::Suppressed => Err(FeatureError::BadParameter {
            name: "target",
            reason: format!("feature {} is suppressed; draft needs a live target", id.0),
        }),
    }
}

fn rotate_about_axis(
    p: &Vector3<f64>,
    axis: &Vector3<f64>,
    cos_a: f64,
    sin_a: f64,
) -> Vector3<f64> {
    // Rodrigues' rotation formula.
    let term1 = p * cos_a;
    let term2 = axis.cross(p) * sin_a;
    let term3 = axis * (axis.dot(p)) * (1.0 - cos_a);
    term1 + term2 + term3
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::{Feature, PadParams};
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
    fn draft_angle_zero_angle_preserves_geometry_node_count() {
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch());
        let pad_id = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Cube",
        );
        tree.add_feature(
            Feature::DraftAngle(DraftAngleParams {
                target: pad_id,
                face_indices: vec![0, 1],
                neutral_plane_normal: NaVec3::new(0.0, 0.0, 1.0),
                draft_angle_deg: 0.0,
            }),
            "No tilt",
        );
        let solid = replay(&tree).expect("replay").expect("solid");
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.5).unwrap();
        assert!(mesh.total_elements() > 0);
    }

    #[test]
    fn draft_rejects_zero_normal() {
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch());
        let pad_id = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 1.0.into(),
                direction_positive: true,
            }),
            "Cube",
        );
        tree.add_feature(
            Feature::DraftAngle(DraftAngleParams {
                target: pad_id,
                face_indices: vec![],
                neutral_plane_normal: NaVec3::zeros(),
                draft_angle_deg: 5.0,
            }),
            "Bad",
        );
        let err = replay(&tree).unwrap_err();
        assert!(matches!(
            err,
            FeatureError::BadParameter {
                name: "neutral_plane_normal",
                ..
            }
        ));
    }
}
