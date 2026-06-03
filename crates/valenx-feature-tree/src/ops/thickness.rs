//! Thickness evaluator — extrude a single face of a target solid along
//! its normal to produce a thin slab.
//!
//! Phase 13D Task 32. v1 implementation: tessellate the target,
//! identify the triangle at `face_index`, build a prism by translating
//! the triangle's three vertices along the triangle normal by
//! `thickness`, emit the resulting wedge as a mesh-backed Solid.

use std::collections::HashMap;

use nalgebra::Vector3;
use valenx_cad::{solid_to_mesh, Solid};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::feature::{FeatureId, ThicknessParams};
use crate::replay::FeatureResult;
use crate::tree::FeatureTree;
use crate::FeatureError;

/// Evaluate a Thickness.
pub(crate) fn evaluate(
    _tree: &FeatureTree,
    p: &ThicknessParams,
    prior: &HashMap<FeatureId, FeatureResult>,
) -> Result<Solid, FeatureError> {
    if !p.thickness.is_finite() || p.thickness <= 0.0 {
        return Err(FeatureError::BadParameter {
            name: "thickness",
            reason: format!("must be > 0 and finite, got {}", p.thickness),
        });
    }
    let target = lookup_target(prior, p.target)?;
    let raw = solid_to_mesh(target, valenx_cad::DEFAULT_TESS_TOLERANCE)?;
    let mut tri_cursor = 0usize;
    // (v0, v1, v2, normal) for the matched triangle.
    type FoundTri = (Vector3<f64>, Vector3<f64>, Vector3<f64>, Vector3<f64>);
    let mut found: Option<FoundTri> = None;
    for block in &raw.element_blocks {
        if block.element_type != ElementType::Tri3 {
            tri_cursor += block.count();
            continue;
        }
        let tris = block.count();
        if p.face_index >= tri_cursor && p.face_index < tri_cursor + tris {
            let local = p.face_index - tri_cursor;
            let base = local * 3;
            let i0 = block.connectivity[base] as usize;
            let i1 = block.connectivity[base + 1] as usize;
            let i2 = block.connectivity[base + 2] as usize;
            let v0 = raw.nodes[i0];
            let v1 = raw.nodes[i1];
            let v2 = raw.nodes[i2];
            let n = (v1 - v0).cross(&(v2 - v0));
            let normal = if n.norm() > 1e-12 {
                n.normalize()
            } else {
                Vector3::z()
            };
            found = Some((v0, v1, v2, normal));
            break;
        }
        tri_cursor += tris;
    }
    let (v0, v1, v2, normal) = found.ok_or(FeatureError::BadParameter {
        name: "face_index",
        reason: format!(
            "no triangle at index {} (total {})",
            p.face_index,
            tri_cursor
                + raw
                    .element_blocks
                    .iter()
                    .filter(|b| b.element_type == ElementType::Tri3)
                    .map(|b| b.count())
                    .sum::<usize>(),
        ),
    })?;

    // Build wedge: 6 vertices (3 base + 3 offset), 8 triangles (2 caps
    // + 6 side triangles for the 3 quads).
    let offset = normal * p.thickness;
    let nodes: Vec<Vector3<f64>> = vec![v0, v1, v2, v0 + offset, v1 + offset, v2 + offset];
    // Front cap (CCW from outside)
    let mut conn: Vec<u32> = vec![0, 1, 2];
    // Back cap (reversed)
    conn.extend_from_slice(&[5, 4, 3]);
    // Side quads as triangle pairs (0-1-4, 0-4-3) etc.
    conn.extend_from_slice(&[0, 1, 4, 0, 4, 3]);
    conn.extend_from_slice(&[1, 2, 5, 1, 5, 4]);
    conn.extend_from_slice(&[2, 0, 3, 2, 3, 5]);

    let mut mesh = Mesh::new("thickness");
    mesh.nodes = nodes;
    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = conn;
    mesh.element_blocks.push(block);
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
            "feature {} has not been evaluated before this thickness (forward / self reference?)",
            id.0
        ),
    })?;
    match res {
        FeatureResult::Solid(s) => Ok(s),
        FeatureResult::Suppressed => Err(FeatureError::BadParameter {
            name: "target",
            reason: format!(
                "feature {} is suppressed; thickness needs a live target",
                id.0
            ),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::{Feature, PadParams};
    use crate::replay::replay;

    fn square_sketch() -> valenx_sketch::Sketch {
        let mut s = valenx_sketch::Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(2.0, 0.0);
        let c = s.add_point(2.0, 2.0);
        let d = s.add_point(0.0, 2.0);
        s.add_line(a, b).unwrap();
        s.add_line(b, c).unwrap();
        s.add_line(c, d).unwrap();
        s.add_line(d, a).unwrap();
        s
    }

    #[test]
    fn thickness_on_first_triangle_emits_wedge() {
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
            Feature::Thickness(ThicknessParams {
                target: pad_id,
                face_index: 0,
                thickness: 0.5,
            }),
            "Wedge",
        );
        let solid = replay(&tree).expect("replay").expect("solid");
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.5).unwrap();
        // 8 tris in our wedge.
        assert!(mesh.total_elements() >= 8);
    }

    #[test]
    fn thickness_rejects_zero_thickness() {
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
            Feature::Thickness(ThicknessParams {
                target: pad_id,
                face_index: 0,
                thickness: 0.0,
            }),
            "Bad",
        );
        let err = replay(&tree).unwrap_err();
        assert!(matches!(
            err,
            FeatureError::BadParameter {
                name: "thickness",
                ..
            }
        ));
    }
}
