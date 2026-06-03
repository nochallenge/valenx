//! Shell evaluator — hollow out a solid into a thin-walled shell.
//!
//! Phase 13D Task 30. v1: tessellation-based — tessellate the target,
//! offset every triangle's vertices inward (or outward) by `thickness`
//! along the triangle's normal, and emit both the original outer mesh
//! and the offset inner mesh as a combined shell. Triangles in
//! `face_indices_to_remove` are skipped on both surfaces, leaving an
//! opening.
//!
//! The result is a mesh-backed [`valenx_cad::Solid::Mesh`] containing
//! both surfaces plus a strip of "wall" triangles connecting the
//! boundaries of removed faces.

use std::collections::{HashMap, HashSet};

use nalgebra::Vector3;
use valenx_cad::{solid_to_mesh, Solid};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::feature::{FeatureId, ShellParams, ShellSide};
use crate::replay::FeatureResult;
use crate::tree::FeatureTree;
use crate::FeatureError;

/// Evaluate a Shell.
pub(crate) fn evaluate(
    _tree: &FeatureTree,
    p: &ShellParams,
    prior: &HashMap<FeatureId, FeatureResult>,
) -> Result<Solid, FeatureError> {
    if !p.thickness.is_finite() || p.thickness <= 0.0 {
        return Err(FeatureError::BadParameter {
            name: "thickness",
            reason: format!("must be > 0 and finite, got {}", p.thickness),
        });
    }
    let target = lookup_target(prior, p.target)?;
    let mut raw = solid_to_mesh(target, valenx_cad::DEFAULT_TESS_TOLERANCE)?;
    raw = valenx_mesh::boolean::merge_coincident_nodes(&raw, 1e-4);

    let remove: HashSet<usize> = p.face_indices_to_remove.iter().copied().collect();
    let sign: f64 = match p.inward_or_outward {
        ShellSide::Inward => -1.0,
        ShellSide::Outward => 1.0,
    };

    // Build per-vertex normal by averaging triangle normals.
    let n_nodes = raw.nodes.len();
    let mut vertex_normal: Vec<Vector3<f64>> = vec![Vector3::zeros(); n_nodes];
    for block in &raw.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        let tris = block.count();
        for tri in 0..tris {
            let base = tri * 3;
            let i0 = block.connectivity[base] as usize;
            let i1 = block.connectivity[base + 1] as usize;
            let i2 = block.connectivity[base + 2] as usize;
            let p0 = raw.nodes[i0];
            let p1 = raw.nodes[i1];
            let p2 = raw.nodes[i2];
            let n = (p1 - p0).cross(&(p2 - p0));
            if n.norm() > 1e-12 {
                let n = n.normalize();
                vertex_normal[i0] += n;
                vertex_normal[i1] += n;
                vertex_normal[i2] += n;
            }
        }
    }
    for v in &mut vertex_normal {
        if v.norm() > 1e-12 {
            *v = v.normalize();
        }
    }

    // Build the offset inner mesh: new vertex positions = original +
    // sign * thickness * vertex_normal.
    let mut nodes_out: Vec<Vector3<f64>> = Vec::with_capacity(n_nodes * 2);
    nodes_out.extend_from_slice(&raw.nodes);
    let offset_base = nodes_out.len();
    for (i, normal_i) in vertex_normal.iter().enumerate().take(n_nodes) {
        nodes_out.push(raw.nodes[i] + *normal_i * sign * p.thickness);
    }

    // Emit triangles. For each triangle:
    //   - If its global index is in `remove`, skip both faces (leaves
    //     the opening).
    //   - Otherwise emit the outer triangle as-is and the inner
    //     triangle with reversed winding (opposite normal).
    let mut conn: Vec<u32> = Vec::new();
    let mut tri_cursor = 0usize;
    let mut removed_edges: Vec<[u32; 2]> = Vec::new();
    for block in &raw.element_blocks {
        if block.element_type != ElementType::Tri3 {
            tri_cursor += block.count();
            continue;
        }
        let tris = block.count();
        for tri in 0..tris {
            let global = tri_cursor + tri;
            let base = tri * 3;
            let i0 = block.connectivity[base];
            let i1 = block.connectivity[base + 1];
            let i2 = block.connectivity[base + 2];
            if remove.contains(&global) {
                // Track its three edges as boundary edges for wall stitching.
                removed_edges.push([i0, i1]);
                removed_edges.push([i1, i2]);
                removed_edges.push([i2, i0]);
                continue;
            }
            conn.push(i0);
            conn.push(i1);
            conn.push(i2);
            // Inner copy with reversed winding.
            conn.push(offset_base as u32 + i2);
            conn.push(offset_base as u32 + i1);
            conn.push(offset_base as u32 + i0);
        }
        tri_cursor += tris;
    }

    // Wall stitching: for each removed-face boundary edge (i, j) emit
    // a quad (two tris) bridging outer (i, j) and inner (i, j).
    for [i, j] in removed_edges {
        let i_out = i;
        let j_out = j;
        let i_in = offset_base as u32 + i;
        let j_in = offset_base as u32 + j;
        // Quad: i_out → j_out → j_in → i_in
        conn.push(i_out);
        conn.push(j_out);
        conn.push(j_in);
        conn.push(i_out);
        conn.push(j_in);
        conn.push(i_in);
    }

    let mut mesh = Mesh::new("shell");
    mesh.nodes = nodes_out;
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
            "feature {} has not been evaluated before this shell (forward / self reference?)",
            id.0
        ),
    })?;
    match res {
        FeatureResult::Solid(s) => Ok(s),
        FeatureResult::Suppressed => Err(FeatureError::BadParameter {
            name: "target",
            reason: format!("feature {} is suppressed; shell needs a live target", id.0),
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
    fn shell_of_cube_produces_more_triangles_than_original() {
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch());
        let pad_id = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 2.0.into(),
                direction_positive: true,
            }),
            "Cube",
        );
        tree.add_feature(
            Feature::Shell(ShellParams {
                target: pad_id,
                face_indices_to_remove: vec![],
                thickness: 0.1,
                inward_or_outward: ShellSide::Inward,
            }),
            "Hollow",
        );
        let solid = replay(&tree).expect("replay").expect("solid");
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.5).unwrap();
        assert!(
            mesh.total_elements() > 6,
            "shell should have many triangles"
        );
    }

    #[test]
    fn shell_rejects_zero_thickness() {
        let mut tree = FeatureTree::new();
        let s = tree.add_sketch(square_sketch());
        let pad_id = tree.add_feature(
            Feature::Pad(PadParams {
                sketch: s,
                depth: 2.0.into(),
                direction_positive: true,
            }),
            "Cube",
        );
        tree.add_feature(
            Feature::Shell(ShellParams {
                target: pad_id,
                face_indices_to_remove: vec![],
                thickness: 0.0,
                inward_or_outward: ShellSide::Inward,
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
