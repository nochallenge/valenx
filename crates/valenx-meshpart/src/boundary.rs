//! Boundary loop extraction + planar flattening.

use std::collections::HashMap;

use nalgebra::Vector3;

use valenx_mesh::element::ElementType;
use valenx_mesh::Mesh;

use crate::error::MeshPartError;
use crate::segment::TriangleGroup;

/// Extract the boundary polyline of the triangle group `group` within
/// `mesh`. The boundary is the chain of edges that appear in
/// exactly one of the group's triangles.
///
/// Returns the loop as an ordered list of vertex coordinates (closed
/// loop — first and last vertex are the same point only when the
/// underlying triangle ring forms a closed polygon).
///
/// Errors:
/// - [`MeshPartError::Empty`] when the group has no triangles or the
///   mesh has no Tri3 block.
/// - [`MeshPartError::BadPolygon`] when the boundary doesn't form a
///   single closed loop (e.g. a torus group has two boundary loops —
///   only the first one is returned; non-manifold groups can produce
///   chains that don't close).
pub fn extract_boundary_loop(
    mesh: &Mesh,
    group: &TriangleGroup,
) -> Result<Vec<Vector3<f64>>, MeshPartError> {
    if group.triangle_indices.is_empty() {
        return Err(MeshPartError::Empty("triangle group"));
    }
    let block = mesh
        .element_blocks
        .iter()
        .find(|b| b.element_type == ElementType::Tri3)
        .ok_or(MeshPartError::Empty("Tri3 block"))?;
    crate::check_connectivity(block, mesh.nodes.len())?;
    let tri_count = block.count();
    for &ti in &group.triangle_indices {
        if ti as usize >= tri_count {
            return Err(MeshPartError::BadParameter {
                name: "triangle_indices",
                reason: format!("triangle index {ti} >= triangle count {tri_count}"),
            });
        }
    }

    // Count each (a,b) directed edge across the group.
    let mut edge_count: HashMap<(u32, u32), i32> = HashMap::new();
    for &ti in &group.triangle_indices {
        let base = ti as usize * 3;
        let tri = [
            block.connectivity[base],
            block.connectivity[base + 1],
            block.connectivity[base + 2],
        ];
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let und = if a < b { (a, b) } else { (b, a) };
            *edge_count.entry(und).or_insert(0) += 1;
        }
    }
    // Boundary edges have count 1.
    let mut next: HashMap<u32, u32> = HashMap::new();
    for &ti in &group.triangle_indices {
        let base = ti as usize * 3;
        let tri = [
            block.connectivity[base],
            block.connectivity[base + 1],
            block.connectivity[base + 2],
        ];
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let und = if a < b { (a, b) } else { (b, a) };
            if edge_count[&und] == 1 {
                next.insert(a, b);
            }
        }
    }
    if next.is_empty() {
        return Err(MeshPartError::Empty("boundary edges"));
    }
    // Walk one loop starting at any node.
    let start = *next.keys().next().unwrap();
    let mut loop_ids = vec![start];
    let mut cur = start;
    loop {
        let Some(&n) = next.get(&cur) else {
            return Err(MeshPartError::BadPolygon(
                "boundary chain broke before closing".into(),
            ));
        };
        if n == start {
            break;
        }
        if loop_ids.len() > next.len() + 2 {
            return Err(MeshPartError::BadPolygon(
                "boundary walk did not terminate".into(),
            ));
        }
        loop_ids.push(n);
        cur = n;
    }
    Ok(loop_ids.into_iter().map(|i| mesh.nodes[i as usize]).collect())
}

/// Project a 3D polyline `loop_3d` onto the 2D plane normal to
/// `plane_normal`. The basis is built by orthogonalising against the
/// closest world axis to `plane_normal`. Returns one `[u, v]` per
/// input point.
///
/// Errors:
/// - [`MeshPartError::Empty`] when `loop_3d` is empty.
/// - [`MeshPartError::BadParameter`] when `plane_normal` is zero-
///   length.
pub fn flatten_boundary(
    loop_3d: &[Vector3<f64>],
    plane_normal: Vector3<f64>,
) -> Result<Vec<[f64; 2]>, MeshPartError> {
    if loop_3d.is_empty() {
        return Err(MeshPartError::Empty("loop_3d"));
    }
    let n = plane_normal
        .try_normalize(1e-12)
        .ok_or(MeshPartError::BadParameter {
            name: "plane_normal",
            reason: "zero-length normal".into(),
        })?;
    // Pick the world axis least-aligned with n.
    let helper = if n.x.abs() < n.y.abs() && n.x.abs() < n.z.abs() {
        Vector3::x()
    } else if n.y.abs() < n.z.abs() {
        Vector3::y()
    } else {
        Vector3::z()
    };
    let u = n.cross(&helper).normalize();
    let v = n.cross(&u);
    Ok(loop_3d.iter().map(|p| [p.dot(&u), p.dot(&v)]).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flatten_empty_errors() {
        let r = flatten_boundary(&[], Vector3::z());
        assert!(matches!(r, Err(MeshPartError::Empty(_))));
    }

    #[test]
    fn flatten_bad_normal_errors() {
        let r = flatten_boundary(&[Vector3::zeros()], Vector3::zeros());
        assert!(matches!(r, Err(MeshPartError::BadParameter { .. })));
    }

    #[test]
    fn flatten_xy_plane_basic() {
        let pts = vec![
            Vector3::new(1.0, 0.0, 5.0),
            Vector3::new(0.0, 2.0, 5.0),
        ];
        let r = flatten_boundary(&pts, Vector3::z()).unwrap();
        assert_eq!(r.len(), 2);
        // Both 2D points should drop the z=5 offset.
        assert!(r[0][0].abs() + r[0][1].abs() > 0.0);
    }
}
