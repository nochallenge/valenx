//! BRep → triangle-mesh tessellation.
//!
//! [`solid_to_mesh`] runs truck-meshalgo's constrained-Delaunay
//! tessellator on the solid's bounding shells, then walks the
//! resulting `truck_polymesh::PolygonMesh` to produce a
//! [`valenx_mesh::Mesh`] of `Tri3` elements. Quads coming out of
//! truck are split into two triangles on the fly; n-gons (n > 4) are
//! fan-triangulated around the first vertex.
//!
//! The resulting mesh is suitable for both the egui viewport
//! renderer (triangle soup) and for STL export — it lives in
//! Valenx's canonical coordinate frame and uses the standard
//! element-block layout (one `Tri3` block, flat connectivity).

use nalgebra::Vector3;
use truck_meshalgo::prelude::MeshableShape;
use truck_meshalgo::tessellation::MeshedShape;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::solid::{CadError, Solid};

/// Tessellate the solid's boundary into a triangle mesh.
///
/// `tolerance` is the chord-error budget passed to truck-meshalgo —
/// roughly the maximum distance any sampled point on the BRep
/// surface deviates from the resulting triangle facet. Smaller
/// values produce denser meshes.
pub fn solid_to_mesh(solid: &Solid, tolerance: f64) -> Result<Mesh, CadError> {
    if !tolerance.is_finite() {
        return Err(CadError::Tessellation(format!(
            "tessellation tolerance must be finite, got {tolerance}"
        )));
    }
    if tolerance <= 0.0 {
        return Err(CadError::Tessellation(format!(
            "tessellation tolerance must be > 0, got {tolerance}"
        )));
    }

    // Mesh-backed solids short-circuit: the cached mesh IS the
    // tessellation output. This is the round-trip that lets fillet /
    // chamfer results flow back through `solid_to_mesh` (which the
    // mesh-toolbox calls before pushing into the viewport) without
    // re-tessellating something that was never a BRep to begin with.
    // The tolerance argument is intentionally ignored for the mesh
    // variant — we have one cached resolution and no way to refine it.
    if let Some(cached) = solid.cached_mesh() {
        return Ok(cached.clone());
    }

    // truck-meshalgo: Solid.triangulation(tol) → triangulated solid,
    // then .to_polygon() merges all the face polygons into a single
    // PolygonMesh. NB: positions inside the polygon are NOT
    // de-duplicated across faces — vertices on shared edges appear
    // multiple times. That's fine for STL / viewport rendering and
    // matches what every other DCC tool does on tessellation. If a
    // downstream consumer needs welded vertices it can run
    // `valenx_mesh::boolean::merge_coincident_nodes` after.
    let polygon = solid
        .try_inner()
        .map_err(|e| CadError::Tessellation(format!("mesh-backed handled above, got {e}")))?
        .triangulation(tolerance)
        .to_polygon();

    let positions = polygon.positions();
    let mut nodes = Vec::with_capacity(positions.len());
    for p in positions.iter() {
        nodes.push(Vector3::new(p.x, p.y, p.z));
    }

    // Build a single Tri3 block, splitting quads / n-gons as we go.
    let mut connectivity: Vec<u32> = Vec::new();
    for tri in polygon.tri_faces() {
        connectivity.push(tri[0].pos as u32);
        connectivity.push(tri[1].pos as u32);
        connectivity.push(tri[2].pos as u32);
    }
    for quad in polygon.quad_faces() {
        // Split a quad into (0,1,2) + (0,2,3).
        let v: [u32; 4] = [
            quad[0].pos as u32,
            quad[1].pos as u32,
            quad[2].pos as u32,
            quad[3].pos as u32,
        ];
        connectivity.extend_from_slice(&[v[0], v[1], v[2]]);
        connectivity.extend_from_slice(&[v[0], v[2], v[3]]);
    }
    for ngon in polygon.other_faces() {
        if ngon.len() < 3 {
            continue;
        }
        // Fan around the first vertex: (v0, v_i, v_{i+1}).
        for i in 1..ngon.len() - 1 {
            connectivity.push(ngon[0].pos as u32);
            connectivity.push(ngon[i].pos as u32);
            connectivity.push(ngon[i + 1].pos as u32);
        }
    }

    if connectivity.is_empty() {
        return Err(CadError::Tessellation(
            "tessellation produced zero triangles — the input solid \
             may have been degenerate"
                .into(),
        ));
    }

    let block = ElementBlock {
        element_type: ElementType::Tri3,
        connectivity,
    };
    let mut mesh = Mesh::new("cad");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{box_solid, sphere};

    #[test]
    fn tessellate_cube_yields_triangles() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let mesh = solid_to_mesh(&cube, 0.1).unwrap();
        assert!(
            mesh.total_elements() >= 12,
            "a cube tessellates to at least 12 triangles (2 per face), got {}",
            mesh.total_elements()
        );
        assert!(!mesh.nodes.is_empty());
        // All elements should be Tri3.
        for block in &mesh.element_blocks {
            assert_eq!(block.element_type, ElementType::Tri3);
        }
    }

    #[test]
    fn tessellate_sphere_with_fine_tolerance() {
        let s = sphere(1.0).unwrap();
        let mesh = solid_to_mesh(&s, 0.1).unwrap();
        assert!(
            mesh.total_elements() > 100,
            "fine sphere tessellation should produce >100 triangles, got {}",
            mesh.total_elements()
        );
    }

    #[test]
    fn tessellate_rejects_bad_tolerance() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        assert!(matches!(
            solid_to_mesh(&cube, -0.1),
            Err(CadError::Tessellation(_))
        ));
        assert!(matches!(
            solid_to_mesh(&cube, 0.0),
            Err(CadError::Tessellation(_))
        ));
        assert!(matches!(
            solid_to_mesh(&cube, f64::INFINITY),
            Err(CadError::Tessellation(_))
        ));
    }

    /// Task 25 regression — `Solid::from_mesh(m).solid_to_mesh()`
    /// round-trips the same mesh without re-tessellating. The
    /// short-circuit in `solid_to_mesh` is the bridge that lets fillet
    /// / chamfer output flow back through the toolbox path.
    #[test]
    fn mesh_backed_solid_round_trips() {
        use valenx_mesh::{ElementBlock, ElementType};
        let mut mesh = Mesh::new("round-trip");
        // Single triangle is the minimum non-empty mesh.
        mesh.nodes.push(Vector3::new(0.0, 0.0, 0.0));
        mesh.nodes.push(Vector3::new(1.0, 0.0, 0.0));
        mesh.nodes.push(Vector3::new(0.0, 1.0, 0.0));
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = vec![0, 1, 2];
        mesh.element_blocks.push(block);
        mesh.recompute_stats();

        let original_node_count = mesh.nodes.len();
        let original_tri_count = mesh.total_elements();

        let solid = crate::Solid::from_mesh(mesh.clone());
        let restored = solid_to_mesh(&solid, 0.5).expect("round-trip succeeds");
        assert_eq!(restored.nodes.len(), original_node_count);
        assert_eq!(restored.total_elements(), original_tri_count);
        // Vertex positions must match exactly — no re-tessellation.
        for (i, n) in mesh.nodes.iter().enumerate() {
            assert!(
                (restored.nodes[i] - n).norm() < 1e-12,
                "round-trip changed node {i}"
            );
        }
    }
}
