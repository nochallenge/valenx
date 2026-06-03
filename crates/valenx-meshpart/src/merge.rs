//! Concatenate meshes into one.

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::error::MeshPartError;

/// Concatenate the nodes + element blocks of `meshes` into one new
/// mesh. Connectivity indices in each input block are rebased by the
/// running node offset so the merged connectivity points at the
/// correct global node indices.
///
/// Output mesh `id` is `merged_<n>` where `n` is the input count.
/// Output `stats` is recomputed (counts only — no quality metrics).
///
/// Errors:
/// - [`MeshPartError::Empty`] when `meshes` is empty.
pub fn merge_meshes(meshes: &[Mesh]) -> Result<Mesh, MeshPartError> {
    if meshes.is_empty() {
        return Err(MeshPartError::Empty("meshes"));
    }
    let mut out = Mesh::new(format!("merged_{}", meshes.len()));
    // Aggregate blocks keyed by ElementType so the result stays in the
    // canonical "one block per type" shape rather than producing N
    // single-mesh Tri3 blocks.
    let mut tri3 = ElementBlock::new(ElementType::Tri3);
    let mut tet4 = ElementBlock::new(ElementType::Tet4);
    let mut quad4 = ElementBlock::new(ElementType::Quad4);
    let mut hex8 = ElementBlock::new(ElementType::Hex8);

    for m in meshes {
        let offset = out.nodes.len() as u32;
        out.nodes.extend_from_slice(&m.nodes);
        for block in &m.element_blocks {
            let bumped: Vec<u32> = block.connectivity.iter().map(|i| i + offset).collect();
            match block.element_type {
                ElementType::Tri3 => tri3.connectivity.extend(bumped),
                ElementType::Tet4 => tet4.connectivity.extend(bumped),
                ElementType::Quad4 => quad4.connectivity.extend(bumped),
                ElementType::Hex8 => hex8.connectivity.extend(bumped),
                _ => {
                    // Other element types stay as their own block — we
                    // preserve them rather than dropping.
                    let mut b = ElementBlock::new(block.element_type);
                    b.connectivity = bumped;
                    out.element_blocks.push(b);
                }
            }
        }
    }
    if !tri3.connectivity.is_empty() {
        out.element_blocks.push(tri3);
    }
    if !tet4.connectivity.is_empty() {
        out.element_blocks.push(tet4);
    }
    if !quad4.connectivity.is_empty() {
        out.element_blocks.push(quad4);
    }
    if !hex8.connectivity.is_empty() {
        out.element_blocks.push(hex8);
    }
    out.recompute_stats();
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    fn make_one_tri(start: f64) -> Mesh {
        let mut m = Mesh::new("a");
        m.nodes = vec![
            Vector3::new(start, 0.0, 0.0),
            Vector3::new(start + 1.0, 0.0, 0.0),
            Vector3::new(start, 1.0, 0.0),
        ];
        let mut b = ElementBlock::new(ElementType::Tri3);
        b.connectivity = vec![0, 1, 2];
        m.element_blocks.push(b);
        m
    }

    #[test]
    fn empty_errors() {
        assert!(matches!(merge_meshes(&[]), Err(MeshPartError::Empty(_))));
    }

    #[test]
    fn two_tris_merged() {
        let m1 = make_one_tri(0.0);
        let m2 = make_one_tri(5.0);
        let m = merge_meshes(&[m1, m2]).unwrap();
        assert_eq!(m.nodes.len(), 6);
        // One aggregated Tri3 block, connectivity rebased.
        let blk = m
            .element_blocks
            .iter()
            .find(|b| b.element_type == ElementType::Tri3)
            .unwrap();
        assert_eq!(blk.connectivity, vec![0, 1, 2, 3, 4, 5]);
    }
}
