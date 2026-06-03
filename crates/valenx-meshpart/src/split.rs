//! Mesh splitting by a plane list.

use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::error::MeshPartError;

/// One cutting plane: `n · (x - p0) = 0`.
#[derive(Copy, Clone, Debug)]
pub struct Plane {
    /// A point on the plane.
    pub point: Vector3<f64>,
    /// Plane normal (need not be unit length — normalised internally).
    pub normal: Vector3<f64>,
}

/// Split `mesh` into regions by the supplied plane list. Each
/// triangle is classified by the signs of its three vertices relative
/// to every plane → that yields a bitmask region id. Each region with
/// ≥1 triangle becomes its own [`Mesh`] in the output Vec; the
/// triangles whose vertices straddle a plane are kept whole (no
/// edge-cut splitting in v1).
///
/// `planes`'s order is significant — region id bit `i` is set iff the
/// triangle lies on the positive side of `planes[i]`.
///
/// Errors:
/// - [`MeshPartError::Empty`] when `planes` is empty or the mesh has
///   no Tri3 block.
/// - [`MeshPartError::BadParameter`] when a plane has a zero-length
///   normal.
pub fn split_mesh_by_planes(mesh: &Mesh, planes: &[Plane]) -> Result<Vec<Mesh>, MeshPartError> {
    if planes.is_empty() {
        return Err(MeshPartError::Empty("planes"));
    }
    let mut norms = Vec::with_capacity(planes.len());
    for (i, p) in planes.iter().enumerate() {
        let n = p.normal.try_normalize(1e-12).ok_or_else(|| {
            MeshPartError::BadParameter {
                name: "planes[i].normal",
                reason: format!("zero-length at index {i}"),
            }
        })?;
        norms.push((p.point, n));
    }
    let block = mesh
        .element_blocks
        .iter()
        .find(|b| b.element_type == ElementType::Tri3)
        .ok_or(MeshPartError::Empty("Tri3 block"))?;
    if block.count() == 0 {
        return Err(MeshPartError::Empty("triangles"));
    }

    let mut groups: std::collections::BTreeMap<u32, ElementBlock> = std::collections::BTreeMap::new();
    for tri in block.connectivity.chunks_exact(3) {
        let centroid = (mesh.nodes[tri[0] as usize]
            + mesh.nodes[tri[1] as usize]
            + mesh.nodes[tri[2] as usize])
            / 3.0;
        let mut bits: u32 = 0;
        for (i, (p0, n)) in norms.iter().enumerate() {
            if i >= 32 {
                break;
            }
            if (centroid - p0).dot(n) > 0.0 {
                bits |= 1 << i;
            }
        }
        let entry = groups
            .entry(bits)
            .or_insert_with(|| ElementBlock::new(ElementType::Tri3));
        entry.connectivity.extend_from_slice(tri);
    }

    let mut out = Vec::with_capacity(groups.len());
    for (region_id, blk) in groups {
        let mut m = Mesh::new(format!("split_region_{region_id}"));
        m.nodes = mesh.nodes.clone();
        m.element_blocks.push(blk);
        m.recompute_stats();
        out.push(m);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn two_tri_strip() -> Mesh {
        let mut m = Mesh::new("strip");
        m.nodes = vec![
            Vector3::new(-1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, -1.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
        ];
        let mut b = ElementBlock::new(ElementType::Tri3);
        b.connectivity = vec![0, 1, 2, 1, 3, 2];
        m.element_blocks.push(b);
        m
    }

    #[test]
    fn split_by_yz_plane_makes_two_regions() {
        let m = two_tri_strip();
        let p = Plane {
            point: Vector3::zeros(),
            normal: Vector3::x(),
        };
        let regions = split_mesh_by_planes(&m, &[p]).unwrap();
        assert_eq!(regions.len(), 2);
    }

    #[test]
    fn empty_plane_list_errors() {
        let m = two_tri_strip();
        assert!(matches!(
            split_mesh_by_planes(&m, &[]),
            Err(MeshPartError::Empty(_))
        ));
    }
}
