//! Apply a [`ManipulateOp`] to a [`valenx_cad::Solid`].
//!
//! v1 pipeline: tessellate the BRep (chord error 0.5 mm — caller
//! controls precision via the input mesh resolution if it pre-meshed),
//! mutate the resulting triangle mesh, return [`Solid::Mesh`].
//! Mesh-backed solids skip the tessellation and mutate directly.

use nalgebra::Vector3;
use valenx_cad::{solid_to_mesh, Solid};
use valenx_mesh::element::ElementType;
use valenx_mesh::Mesh;

use crate::error::ManipulatorError;
use crate::op::ManipulateOp;

/// Default tessellation tolerance (mm) when the input is a BRep.
pub const DEFAULT_TOLERANCE_MM: f64 = 0.5;

/// Apply a single op to `solid` and return the resulting mesh-
/// backed solid.
pub fn apply(solid: &Solid, op: &ManipulateOp) -> Result<Solid, ManipulatorError> {
    let mut mesh = into_mesh(solid)?;
    apply_to_mesh(&mut mesh, op)?;
    mesh.recompute_stats();
    Ok(Solid::from_mesh(mesh))
}

/// Apply a sequence of ops in order.
pub fn apply_sequence(solid: &Solid, ops: &[ManipulateOp]) -> Result<Solid, ManipulatorError> {
    let mut mesh = into_mesh(solid)?;
    for op in ops {
        apply_to_mesh(&mut mesh, op)?;
    }
    mesh.recompute_stats();
    Ok(Solid::from_mesh(mesh))
}

fn into_mesh(solid: &Solid) -> Result<Mesh, ManipulatorError> {
    match solid {
        Solid::Mesh(m) => Ok(m.clone()),
        Solid::Brep(_) => solid_to_mesh(solid, DEFAULT_TOLERANCE_MM)
            .map_err(|e| ManipulatorError::Tessellation(e.to_string())),
    }
}

fn apply_to_mesh(mesh: &mut Mesh, op: &ManipulateOp) -> Result<(), ManipulatorError> {
    match op {
        ManipulateOp::MoveVertex { vertex_idx, delta } => {
            let n = mesh.nodes.len();
            if *vertex_idx >= n {
                return Err(ManipulatorError::BadIndex {
                    got: *vertex_idx,
                    n,
                });
            }
            mesh.nodes[*vertex_idx] += delta;
            Ok(())
        }
        ManipulateOp::MoveFace { face_idx, delta } => {
            let (block_idx, tri_idx) = locate_face(mesh, *face_idx)?;
            let block = &mesh.element_blocks[block_idx];
            let a = block.connectivity[tri_idx * 3] as usize;
            let b = block.connectivity[tri_idx * 3 + 1] as usize;
            let c = block.connectivity[tri_idx * 3 + 2] as usize;
            // Duplicate vertices so neighbouring faces don't drag
            // along — this is the canonical push/pull behaviour.
            let new_a = duplicate_vertex(mesh, a, *delta);
            let new_b = duplicate_vertex(mesh, b, *delta);
            let new_c = duplicate_vertex(mesh, c, *delta);
            let block = &mut mesh.element_blocks[block_idx];
            block.connectivity[tri_idx * 3] = new_a as u32;
            block.connectivity[tri_idx * 3 + 1] = new_b as u32;
            block.connectivity[tri_idx * 3 + 2] = new_c as u32;
            Ok(())
        }
        ManipulateOp::RotateFace {
            face_idx,
            axis,
            angle_deg,
        } => {
            let (block_idx, tri_idx) = locate_face(mesh, *face_idx)?;
            let block = &mesh.element_blocks[block_idx];
            let a = block.connectivity[tri_idx * 3] as usize;
            let b = block.connectivity[tri_idx * 3 + 1] as usize;
            let c = block.connectivity[tri_idx * 3 + 2] as usize;
            let centroid = (mesh.nodes[a] + mesh.nodes[b] + mesh.nodes[c]) / 3.0;
            let ax = axis
                .try_normalize(1e-12)
                .ok_or_else(|| ManipulatorError::BadParameter {
                    name: "axis",
                    reason: "zero-length axis vector".into(),
                })?;
            let rad = angle_deg.to_radians();
            let cs = rad.cos();
            let sn = rad.sin();
            for &v in &[a, b, c] {
                let rel = mesh.nodes[v] - centroid;
                let rot = rotate_rodrigues(rel, ax, cs, sn);
                mesh.nodes[v] = centroid + rot;
            }
            Ok(())
        }
        ManipulateOp::MoveEdge { edge_idx, delta } => {
            let (block_idx, tri_idx, edge_in_tri) = locate_edge(mesh, *edge_idx)?;
            let block = &mesh.element_blocks[block_idx];
            let i0 = block.connectivity[tri_idx * 3 + edge_in_tri] as usize;
            let i1 = block.connectivity[tri_idx * 3 + (edge_in_tri + 1) % 3] as usize;
            mesh.nodes[i0] += delta;
            mesh.nodes[i1] += delta;
            Ok(())
        }
        ManipulateOp::ExtrudeFace { face_idx, distance } => {
            let (block_idx, tri_idx) = locate_face(mesh, *face_idx)?;
            let block = &mesh.element_blocks[block_idx];
            let a = block.connectivity[tri_idx * 3] as usize;
            let b = block.connectivity[tri_idx * 3 + 1] as usize;
            let c = block.connectivity[tri_idx * 3 + 2] as usize;
            let normal = face_normal(mesh, a, b, c);
            let delta = normal * (*distance);
            let new_a = duplicate_vertex(mesh, a, delta);
            let new_b = duplicate_vertex(mesh, b, delta);
            let new_c = duplicate_vertex(mesh, c, delta);
            // Move the face to the extruded position.
            let block = &mut mesh.element_blocks[block_idx];
            block.connectivity[tri_idx * 3] = new_a as u32;
            block.connectivity[tri_idx * 3 + 1] = new_b as u32;
            block.connectivity[tri_idx * 3 + 2] = new_c as u32;
            // Add 3 side-wall quads (2 tris each = 6 new tris).
            let aa = a as u32;
            let bb = b as u32;
            let cc = c as u32;
            let na = new_a as u32;
            let nb = new_b as u32;
            let nc = new_c as u32;
            block.connectivity.extend_from_slice(&[aa, bb, nb]);
            block.connectivity.extend_from_slice(&[aa, nb, na]);
            block.connectivity.extend_from_slice(&[bb, cc, nc]);
            block.connectivity.extend_from_slice(&[bb, nc, nb]);
            block.connectivity.extend_from_slice(&[cc, aa, na]);
            block.connectivity.extend_from_slice(&[cc, na, nc]);
            Ok(())
        }
        ManipulateOp::OffsetFace { face_idx, distance } => {
            let (block_idx, tri_idx) = locate_face(mesh, *face_idx)?;
            let block = &mesh.element_blocks[block_idx];
            let a = block.connectivity[tri_idx * 3] as usize;
            let b = block.connectivity[tri_idx * 3 + 1] as usize;
            let c = block.connectivity[tri_idx * 3 + 2] as usize;
            let normal = face_normal(mesh, a, b, c);
            let delta = normal * (*distance);
            mesh.nodes[a] += delta;
            mesh.nodes[b] += delta;
            mesh.nodes[c] += delta;
            Ok(())
        }
    }
}

fn duplicate_vertex(mesh: &mut Mesh, idx: usize, delta: Vector3<f64>) -> usize {
    let new = mesh.nodes[idx] + delta;
    mesh.nodes.push(new);
    mesh.nodes.len() - 1
}

fn face_normal(mesh: &Mesh, a: usize, b: usize, c: usize) -> Vector3<f64> {
    let pa = mesh.nodes[a];
    let pb = mesh.nodes[b];
    let pc = mesh.nodes[c];
    let n = (pb - pa).cross(&(pc - pa));
    n.try_normalize(1e-12).unwrap_or_else(Vector3::z)
}

fn rotate_rodrigues(v: Vector3<f64>, k: Vector3<f64>, cs: f64, sn: f64) -> Vector3<f64> {
    v * cs + k.cross(&v) * sn + k * k.dot(&v) * (1.0 - cs)
}

fn locate_face(mesh: &Mesh, face_idx: usize) -> Result<(usize, usize), ManipulatorError> {
    let mut running = 0usize;
    for (bi, block) in mesh.element_blocks.iter().enumerate() {
        if !matches!(block.element_type, ElementType::Tri3) {
            continue;
        }
        let tris_in_block = block.connectivity.len() / 3;
        if face_idx < running + tris_in_block {
            return Ok((bi, face_idx - running));
        }
        running += tris_in_block;
    }
    Err(ManipulatorError::BadIndex {
        got: face_idx,
        n: running,
    })
}

fn locate_edge(mesh: &Mesh, edge_idx: usize) -> Result<(usize, usize, usize), ManipulatorError> {
    // Unique-edge id = sum over preceding triangles × 3.
    // Each triangle contributes 3 edges; edge_idx is `tri × 3 + e`.
    let mut running = 0usize;
    for (bi, block) in mesh.element_blocks.iter().enumerate() {
        if !matches!(block.element_type, ElementType::Tri3) {
            continue;
        }
        let tris = block.connectivity.len() / 3;
        let edges_in_block = tris * 3;
        if edge_idx < running + edges_in_block {
            let local = edge_idx - running;
            return Ok((bi, local / 3, local % 3));
        }
        running += edges_in_block;
    }
    Err(ManipulatorError::BadIndex {
        got: edge_idx,
        n: running,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_mesh::element::ElementBlock;

    fn unit_tri_mesh() -> Mesh {
        let mut mesh = Mesh::new("unit_tri");
        mesh.nodes.push(Vector3::new(0.0, 0.0, 0.0));
        mesh.nodes.push(Vector3::new(1.0, 0.0, 0.0));
        mesh.nodes.push(Vector3::new(0.0, 1.0, 0.0));
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity.extend_from_slice(&[0, 1, 2]);
        mesh.element_blocks.push(block);
        mesh.recompute_stats();
        mesh
    }

    #[test]
    fn move_vertex_translates_single_node() {
        let solid = Solid::from_mesh(unit_tri_mesh());
        let op = ManipulateOp::MoveVertex {
            vertex_idx: 0,
            delta: Vector3::new(0.0, 0.0, 5.0),
        };
        let s2 = apply(&solid, &op).unwrap();
        match s2 {
            Solid::Mesh(m) => assert!((m.nodes[0].z - 5.0).abs() < 1e-9),
            _ => panic!(),
        }
    }

    #[test]
    fn move_face_duplicates_vertices() {
        let solid = Solid::from_mesh(unit_tri_mesh());
        let op = ManipulateOp::MoveFace {
            face_idx: 0,
            delta: Vector3::new(0.0, 0.0, 1.0),
        };
        let s2 = apply(&solid, &op).unwrap();
        match s2 {
            Solid::Mesh(m) => {
                // 3 original + 3 duplicates.
                assert_eq!(m.nodes.len(), 6);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn extrude_face_adds_six_side_tris() {
        let solid = Solid::from_mesh(unit_tri_mesh());
        let op = ManipulateOp::ExtrudeFace {
            face_idx: 0,
            distance: 1.0,
        };
        let s2 = apply(&solid, &op).unwrap();
        match s2 {
            Solid::Mesh(m) => {
                // 1 face tri + 6 side wall tris.
                assert_eq!(m.total_elements(), 7);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn offset_face_translates_in_place() {
        let solid = Solid::from_mesh(unit_tri_mesh());
        let op = ManipulateOp::OffsetFace {
            face_idx: 0,
            distance: 1.0,
        };
        let s2 = apply(&solid, &op).unwrap();
        match s2 {
            Solid::Mesh(m) => {
                // No new vertices, just translated up by 1 (normal = +z).
                assert_eq!(m.nodes.len(), 3);
                assert!((m.nodes[0].z - 1.0).abs() < 1e-9);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn rotate_face_keeps_centroid_fixed() {
        let solid = Solid::from_mesh(unit_tri_mesh());
        let op = ManipulateOp::RotateFace {
            face_idx: 0,
            axis: Vector3::new(0.0, 0.0, 1.0),
            angle_deg: 90.0,
        };
        let s2 = apply(&solid, &op).unwrap();
        match s2 {
            Solid::Mesh(m) => {
                let centroid_before = Vector3::new(1.0 / 3.0, 1.0 / 3.0, 0.0);
                let centroid_after = (m.nodes[0] + m.nodes[1] + m.nodes[2]) / 3.0;
                assert!((centroid_after - centroid_before).norm() < 1e-9);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn move_edge_translates_both_endpoints() {
        let solid = Solid::from_mesh(unit_tri_mesh());
        let op = ManipulateOp::MoveEdge {
            edge_idx: 0,
            delta: Vector3::new(0.0, 0.0, 2.0),
        };
        let s2 = apply(&solid, &op).unwrap();
        match s2 {
            Solid::Mesh(m) => {
                assert!((m.nodes[0].z - 2.0).abs() < 1e-9);
                assert!((m.nodes[1].z - 2.0).abs() < 1e-9);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn out_of_range_index_errors() {
        let solid = Solid::from_mesh(unit_tri_mesh());
        let op = ManipulateOp::MoveVertex {
            vertex_idx: 99,
            delta: Vector3::zeros(),
        };
        assert!(matches!(
            apply(&solid, &op),
            Err(ManipulatorError::BadIndex { .. })
        ));
    }

    #[test]
    fn apply_sequence_chains_ops() {
        let solid = Solid::from_mesh(unit_tri_mesh());
        let ops = vec![
            ManipulateOp::OffsetFace {
                face_idx: 0,
                distance: 1.0,
            },
            ManipulateOp::OffsetFace {
                face_idx: 0,
                distance: 1.0,
            },
        ];
        let s2 = apply_sequence(&solid, &ops).unwrap();
        match s2 {
            Solid::Mesh(m) => assert!((m.nodes[0].z - 2.0).abs() < 1e-9),
            _ => panic!(),
        }
    }
}
