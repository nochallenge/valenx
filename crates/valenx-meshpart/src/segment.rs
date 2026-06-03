//! Triangle-group segmentation by face normal.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use valenx_mesh::element::ElementType;
use valenx_mesh::Mesh;

use crate::error::MeshPartError;

/// A connected group of triangle indices sharing similar face
/// normals.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct TriangleGroup {
    /// Indices into the first Tri3 element block's triangle array
    /// (zero-based triangle index, not vertex index).
    pub triangle_indices: Vec<u32>,
    /// The representative (first triangle's) normal — useful for
    /// downstream tools (e.g. extract a sketch on the plane it spans).
    pub representative_normal: Vector3<f64>,
}

/// Group every triangle in `mesh`'s first Tri3 block by face normal,
/// where two triangles join the same group iff their face normals
/// differ by at most `angle_threshold_deg`. The grouping is purely
/// normal-based — triangles in non-adjacent regions can land in the
/// same group if they happen to face the same direction.
///
/// Returns one [`TriangleGroup`] per cluster (order: discovery order).
///
/// Errors:
/// - [`MeshPartError::BadParameter`] when `angle_threshold_deg` is
///   negative or non-finite.
/// - [`MeshPartError::Empty`] when the mesh has no Tri3 elements.
pub fn segment_by_normal(
    mesh: &Mesh,
    angle_threshold_deg: f64,
) -> Result<Vec<TriangleGroup>, MeshPartError> {
    if !angle_threshold_deg.is_finite() || angle_threshold_deg < 0.0 {
        return Err(MeshPartError::BadParameter {
            name: "angle_threshold_deg",
            reason: format!("must be finite >= 0, got {angle_threshold_deg}"),
        });
    }
    let block = mesh
        .element_blocks
        .iter()
        .find(|b| b.element_type == ElementType::Tri3)
        .ok_or(MeshPartError::Empty("Tri3 block"))?;
    let tri_count = block.count();
    if tri_count == 0 {
        return Err(MeshPartError::Empty("triangles"));
    }
    // Pre-compute per-triangle face normals.
    let mut normals: Vec<Vector3<f64>> = Vec::with_capacity(tri_count);
    for tri in block.connectivity.chunks_exact(3) {
        let (a, b, c) = (tri[0] as usize, tri[1] as usize, tri[2] as usize);
        let pa = mesh.nodes[a];
        let pb = mesh.nodes[b];
        let pc = mesh.nodes[c];
        let n = (pb - pa).cross(&(pc - pa));
        normals.push(n.try_normalize(1e-12).unwrap_or_else(Vector3::z));
    }
    let cos_threshold = angle_threshold_deg.to_radians().cos();

    let mut groups: Vec<TriangleGroup> = Vec::new();
    'outer: for (i, n) in normals.iter().enumerate() {
        for g in groups.iter_mut() {
            // Compare against representative.
            if n.dot(&g.representative_normal) >= cos_threshold {
                g.triangle_indices.push(i as u32);
                continue 'outer;
            }
        }
        groups.push(TriangleGroup {
            triangle_indices: vec![i as u32],
            representative_normal: *n,
        });
    }
    Ok(groups)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_mesh::element::ElementBlock;

    fn cube_top_and_bottom() -> Mesh {
        // 4 tris: 2 facing +Z (top), 2 facing -Z (bottom).
        let mut m = Mesh::new("c");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(0.0, 1.0, 1.0),
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut b = ElementBlock::new(ElementType::Tri3);
        b.connectivity = vec![
            // Top.
            0, 1, 2, 0, 2, 3, // +Z
            // Bottom.
            4, 5, 6, 4, 7, 5, // -Z
        ];
        m.element_blocks.push(b);
        m
    }

    #[test]
    fn segments_into_two_groups() {
        let m = cube_top_and_bottom();
        let groups = segment_by_normal(&m, 5.0).unwrap();
        assert_eq!(groups.len(), 2);
    }

    #[test]
    fn bad_threshold_errors() {
        let m = cube_top_and_bottom();
        assert!(matches!(
            segment_by_normal(&m, -1.0),
            Err(MeshPartError::BadParameter { name: "angle_threshold_deg", .. })
        ));
    }
}
