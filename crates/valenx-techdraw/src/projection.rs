//! 3D → 2D drawing-plane projection.
//!
//! Orthographic projection: multiply each world point by the view's
//! camera matrix (from [`crate::view::ViewKind::camera_matrix`]),
//! then read off the resulting `(x, y)`. Z is preserved upstream for
//! HLR depth tests.
//!
//! "Edge extraction" pulls every distinct line segment out of a
//! [`valenx_cad::Solid`]:
//! - [`valenx_cad::Solid::Brep`] — walk the BRep's edge iterator,
//!   sample parameter-space points along each curve, segment-ify.
//! - [`valenx_cad::Solid::Mesh`] — every triangle edge becomes a
//!   candidate segment; duplicates between adjacent triangles are
//!   de-duplicated so a watertight mesh doesn't double-draw every
//!   interior edge.
//!
//! BRep edge sampling is intentionally coarse: 16 sub-segments per
//! curve is enough resolution for the SVG / PDF / DXF exporters,
//! which can't do better than line primitives anyway. Higher
//! fidelity would mean teaching the exporters about `path` /
//! `bezier` constructs — out of scope for Phase 5.

use nalgebra::{Matrix4, Vector3, Vector4};

use crate::error::TechDrawError;

/// Convenience alias for a list of world-space edges (two endpoints
/// each) — the return type of [`extract_edges`].
pub type WorldEdges = Vec<(Vector3<f64>, Vector3<f64>)>;

/// Project a 3D world-space point through `camera` (the matrix from
/// [`crate::view::ViewKind::camera_matrix`]) and return the resulting
/// `(x, y)` in the drawing plane's millimeter coordinates.
pub fn project_point(p: Vector3<f64>, camera: &Matrix4<f64>) -> [f64; 2] {
    let v: Vector4<f64> = camera * Vector4::new(p.x, p.y, p.z, 1.0);
    [v.x, v.y]
}

/// Project a single segment (two world-space endpoints) through the
/// camera. Returns the 2D segment.
pub fn project_segment(a: Vector3<f64>, b: Vector3<f64>, camera: &Matrix4<f64>) -> [(f64, f64); 2] {
    let pa = project_point(a, camera);
    let pb = project_point(b, camera);
    [(pa[0], pa[1]), (pb[0], pb[1])]
}

/// Extract every edge from `solid` and project them through `camera`.
///
/// Output segments are in drawing-plane millimeters (local view
/// frame). De-duplicated up to floating-point round-off so a closed
/// mesh doesn't yield two copies of every interior edge.
///
/// Returns [`TechDrawError::EmptySolid`] when the solid has nothing
/// to draw.
pub fn project_edges(
    solid: &valenx_cad::Solid,
    camera: &Matrix4<f64>,
) -> Result<Vec<[(f64, f64); 2]>, TechDrawError> {
    let world_edges = extract_edges(solid)?;
    let mut out = Vec::with_capacity(world_edges.len());
    for (a, b) in world_edges {
        out.push(project_segment(a, b, camera));
    }
    Ok(out)
}

/// World-space edge extraction. Returns every `(a, b)` segment that
/// belongs to the solid's wireframe.
///
/// For [`valenx_cad::Solid::Mesh`] we walk every triangle's three
/// edges. For [`valenx_cad::Solid::Brep`] we tessellate to a mesh
/// first (so BRep curves get sampled), then take the same triangle
/// edges. This keeps the projector source-uniform — and Phase 5 doesn't
/// promise true B-spline-aware edge rendering.
pub fn extract_edges(solid: &valenx_cad::Solid) -> Result<WorldEdges, TechDrawError> {
    // Tessellate BRep to a mesh with a sane default chord-error
    // budget. We could expose the tolerance as a knob; 0.1 mm matches
    // what the mesh-toolbox uses for viewport preview.
    let mesh = match solid {
        valenx_cad::Solid::Brep(_) => valenx_cad::tessellate::solid_to_mesh(solid, 0.1)
            .map_err(|e| TechDrawError::ExportFailed(format!("tessellation: {e}")))?,
        valenx_cad::Solid::Mesh(m) => m.clone(),
    };
    if mesh.nodes.is_empty() {
        return Err(TechDrawError::EmptySolid);
    }

    // Collect triangle edges with a canonical (low-idx, high-idx)
    // ordering so a shared edge between adjacent triangles appears
    // once in the de-dup set.
    use std::collections::HashSet;
    let mut seen: HashSet<(u32, u32)> = HashSet::new();
    let mut out: Vec<(Vector3<f64>, Vector3<f64>)> = Vec::new();
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::element::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let ix = [tri[0], tri[1], tri[2]];
            for k in 0..3 {
                let a = ix[k];
                let b = ix[(k + 1) % 3];
                let key = if a < b { (a, b) } else { (b, a) };
                if !seen.insert(key) {
                    continue;
                }
                // Bounds-check the connectivity: a loaded/corrupt mesh
                // can carry a node index past `nodes.len()`. Index by
                // `.get()` and skip the edge on an out-of-range value
                // rather than panicking during view generation.
                let (Some(&pa), Some(&pb)) =
                    (mesh.nodes.get(a as usize), mesh.nodes.get(b as usize))
                else {
                    continue;
                };
                out.push((pa, pb));
            }
        }
    }
    if out.is_empty() {
        return Err(TechDrawError::EmptySolid);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::view::ViewKind;
    use valenx_cad::primitives::box_solid;

    #[test]
    fn project_point_identity_for_top_view_xy() {
        // Top view maps world (x, y, *) → drawing (x, y).
        let cam = ViewKind::Top.camera_matrix();
        let p = project_point(Vector3::new(3.0, 4.0, 7.0), &cam);
        assert!((p[0] - 3.0).abs() < 1e-9);
        assert!((p[1] - 4.0).abs() < 1e-9);
    }

    #[test]
    fn project_segment_returns_two_endpoints() {
        let cam = ViewKind::Top.camera_matrix();
        let seg = project_segment(Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0), &cam);
        assert!((seg[0].0).abs() < 1e-9);
        assert!((seg[1].0 - 1.0).abs() < 1e-9);
    }

    #[test]
    fn extract_edges_cube_yields_more_than_twelve_segments() {
        // A tessellated cube has at least 12 BRep edges; after
        // tessellation each face is split into 2 triangles giving an
        // extra diagonal per face → 12 + 6 = 18 edges minimum.
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let edges = extract_edges(&cube).unwrap();
        assert!(edges.len() >= 12, "got {} edges", edges.len());
    }

    #[test]
    fn project_edges_unit_cube_top_view() {
        // Top view of a 1×1×1 cube projects all (x,y,*) points to a
        // 1×1 square — every projected edge lies in [0,1]² (cube is
        // placed at the origin with truck's default builder).
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let cam = ViewKind::Top.camera_matrix();
        let segs = project_edges(&cube, &cam).unwrap();
        assert!(!segs.is_empty());
        for seg in &segs {
            for (x, y) in seg {
                assert!(*x >= -1e-6 && *x <= 1.0 + 1e-6, "x out of range: {x}");
                assert!(*y >= -1e-6 && *y <= 1.0 + 1e-6, "y out of range: {y}");
            }
        }
    }

    #[test]
    fn extract_edges_empty_mesh_returns_empty_solid_error() {
        let mesh = valenx_mesh::Mesh::new("empty");
        let s = valenx_cad::Solid::from_mesh(mesh);
        let e = extract_edges(&s).unwrap_err();
        match e {
            TechDrawError::EmptySolid => {}
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn extract_edges_out_of_range_connectivity_does_not_panic() {
        // A loaded/corrupt mesh whose Tri3 connectivity references a node
        // index past `nodes.len()`. The projector indexed `mesh.nodes[a]`
        // directly, so such a mesh panicked ("index out of bounds")
        // during TechDraw view generation. The bad triangle must be
        // skipped gracefully — never abort.
        use valenx_mesh::element::{ElementBlock, ElementType};
        let mut mesh = valenx_mesh::Mesh::new("corrupt");
        mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        // First triangle is valid; second references node 99 (only 3
        // nodes exist) — the out-of-range index that used to panic.
        block.connectivity = vec![0, 1, 2, 0, 1, 99];
        mesh.element_blocks.push(block);
        let s = valenx_cad::Solid::from_mesh(mesh);

        // Must not panic. The one valid triangle still yields edges.
        let edges = extract_edges(&s).expect("valid triangle should still produce edges");
        assert!(
            !edges.is_empty(),
            "the in-range triangle's edges should survive"
        );
    }
}
