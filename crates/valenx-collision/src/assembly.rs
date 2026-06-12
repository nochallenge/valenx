//! Pairwise assembly collision check.

use valenx_assembly::Assembly;
use valenx_cad::{solid_to_mesh, DEFAULT_TESS_TOLERANCE};
use valenx_mesh::Mesh;

use crate::error::CollisionError;
use crate::mesh_pair::{collide, CollisionInfo};

/// Run a pairwise collision check across every distinct pair of
/// parts in `assembly`. Returns `(part_a_id, part_b_id,
/// CollisionInfo)` tuples — one entry per detected pair.
pub fn check_collisions(
    assembly: &Assembly,
) -> Result<Vec<(usize, usize, CollisionInfo)>, CollisionError> {
    // Tessellate each part once.
    let mut tess: Vec<(usize, Mesh)> = Vec::with_capacity(assembly.parts.len());
    for p in &assembly.parts {
        let mesh = solid_to_mesh(&p.solid, DEFAULT_TESS_TOLERANCE)
            .map_err(|e| CollisionError::Tessellation(e.to_string()))?;
        tess.push((p.id, mesh));
    }

    let mut out = Vec::new();
    for i in 0..tess.len() {
        for j in (i + 1)..tess.len() {
            if let Some(info) = collide(&tess[i].1, &tess[j].1) {
                out.push((tess[i].0, tess[j].0, info));
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;
    use valenx_assembly::part::Part;
    use valenx_cad::Solid;
    use valenx_mesh::element::{ElementBlock, ElementType};

    fn unit_tri_solid(translate: Vector3<f64>) -> Solid {
        let mut m = Mesh::new("t");
        m.nodes.push(translate);
        m.nodes.push(translate + Vector3::new(1.0, 0.0, 0.0));
        m.nodes.push(translate + Vector3::new(0.0, 1.0, 0.0));
        let mut blk = ElementBlock::new(ElementType::Tri3);
        blk.connectivity.extend_from_slice(&[0, 1, 2]);
        m.element_blocks.push(blk);
        m.recompute_stats();
        Solid::from_mesh(m)
    }

    #[test]
    fn check_collisions_finds_overlapping_pair() {
        let mut asm = Assembly::new();
        asm.add_part(Part::new(0, "A", unit_tri_solid(Vector3::zeros())));
        asm.add_part(Part::new(1, "B", unit_tri_solid(Vector3::zeros())));
        let result = check_collisions(&asm).unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].0, 0);
        assert_eq!(result[0].1, 1);
    }

    #[test]
    fn check_collisions_skips_disjoint_pair() {
        let mut asm = Assembly::new();
        asm.add_part(Part::new(0, "A", unit_tri_solid(Vector3::zeros())));
        asm.add_part(Part::new(
            1,
            "B",
            unit_tri_solid(Vector3::new(100.0, 0.0, 0.0)),
        ));
        let result = check_collisions(&asm).unwrap();
        assert!(result.is_empty());
    }
}
