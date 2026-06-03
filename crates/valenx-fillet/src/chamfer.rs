//! [`apply_chamfer`] — replace every sharp convex edge of a mesh with
//! a flat bevel strip of width `distance`.
//!
//! Chamfer is the "flat" variant of fillet: each sharp edge becomes
//! a single quad (two triangles) connecting the two face-offset
//! points. Same v1 simplifications as fillet (no endpoint
//! repositioning, no triangle clipping, no corner blending — see
//! [`crate::fillet`] module docs).
//!
//! Geometry: on each face share a point `r` units away from the
//! shared edge along the face bitangent, then connect the two with
//! a straight quad strip.

use std::collections::HashMap;

use nalgebra::Vector3;
use valenx_mesh::{ElementBlock, ElementType, Mesh};

use crate::edge_graph::{EdgeGraph, EdgeKey};
use crate::error::FilletError;

/// Spatial dedup tolerance — same as fillet.
const DEDUP_TOL: f64 = 1e-6;

fn bucket(v: f64) -> i64 {
    (v / DEDUP_TOL).round() as i64
}

/// A single flat-bevel quad strip for one edge.
#[derive(Clone, Debug)]
struct ChamferStrip {
    /// Four corner vertices, world space: `[p0+r*b0, p1+r*b0, p0+r*b1, p1+r*b1]`.
    vertices: [Vector3<f64>; 4],
    /// Two triangles wound outward.
    triangles: [(usize, usize, usize); 2],
}

/// Build the chamfer quad for the given edge.
fn build_chamfer_strip(
    p0: Vector3<f64>,
    p1: Vector3<f64>,
    b0: Vector3<f64>,
    b1: Vector3<f64>,
    distance: f64,
) -> ChamferStrip {
    let v0 = p0 + b0 * distance;
    let v1 = p1 + b0 * distance;
    let v2 = p0 + b1 * distance;
    let v3 = p1 + b1 * distance;
    ChamferStrip {
        vertices: [v0, v1, v2, v3],
        // Triangles: (0,1,3) + (0,3,2). Same winding pattern as
        // cyl_strip's quads.
        triangles: [(0, 1, 3), (0, 3, 2)],
    }
}

/// Apply a constant-distance chamfer to every convex sharp edge of
/// `mesh` whose dihedral angle exceeds `threshold_rad`.
///
/// # Parameters
/// - `mesh`: source triangle mesh.
/// - `distance`: chamfer width in mesh units — how far inward each
///   face is offset before the bevel connects them. Must be > 0.
/// - `threshold_rad`: minimum dihedral angle (radians) for an edge
///   to qualify.
///
/// # Errors
/// Same as [`crate::apply_fillet`].
///
/// # v1 limitations
/// See [`crate::fillet`] module docs.
pub fn apply_chamfer(mesh: &Mesh, distance: f64, threshold_rad: f64) -> Result<Mesh, FilletError> {
    if distance <= 0.0 {
        return Err(FilletError::BadParameter {
            name: "distance",
            reason: format!("must be > 0, got {distance}"),
        });
    }
    if threshold_rad < 0.0 {
        return Err(FilletError::BadParameter {
            name: "threshold_rad",
            reason: format!("must be >= 0, got {threshold_rad}"),
        });
    }

    let graph = EdgeGraph::from_mesh(mesh);
    if graph.triangles.is_empty() {
        return Err(FilletError::EmptyMesh);
    }

    let bevelable = graph.filletable_edges(mesh, threshold_rad);

    let mut out = Mesh::new(format!("{}-chamfered", mesh.id));
    out.nodes = mesh.nodes.clone();
    let mut connectivity: Vec<u32> = Vec::new();
    for tri in &graph.triangles {
        connectivity.push(tri.v[0] as u32);
        connectivity.push(tri.v[1] as u32);
        connectivity.push(tri.v[2] as u32);
    }

    let mut dedup: HashMap<(i64, i64, i64), usize> = HashMap::new();
    for (idx, p) in out.nodes.iter().copied().enumerate() {
        dedup
            .entry((bucket(p.x), bucket(p.y), bucket(p.z)))
            .or_insert(idx);
    }

    for key in bevelable {
        let p0 = mesh.nodes[key.0];
        let p1 = mesh.nodes[key.1];
        if (p1 - p0).norm() < 1e-12 {
            return Err(FilletError::DegenerateEdge {
                from: key.0,
                to: key.1,
            });
        }
        let tris = &graph.adjacency[&key];
        if tris.len() != 2 {
            continue;
        }
        let Some(b0) = face_bitangent(mesh, &graph, tris[0], key) else {
            continue;
        };
        let Some(b1) = face_bitangent(mesh, &graph, tris[1], key) else {
            continue;
        };
        let strip = build_chamfer_strip(p0, p1, b0, b1, distance);
        let mut map = [0usize; 4];
        for (i, &p) in strip.vertices.iter().enumerate() {
            map[i] = intern(&mut out.nodes, &mut dedup, p);
        }
        for &(a, b, c) in &strip.triangles {
            connectivity.push(map[a] as u32);
            connectivity.push(map[b] as u32);
            connectivity.push(map[c] as u32);
        }
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = connectivity;
    out.element_blocks.push(block);
    out.recompute_stats();
    Ok(out)
}

fn intern(
    nodes: &mut Vec<Vector3<f64>>,
    dedup: &mut HashMap<(i64, i64, i64), usize>,
    p: Vector3<f64>,
) -> usize {
    let key = (bucket(p.x), bucket(p.y), bucket(p.z));
    if let Some(&idx) = dedup.get(&key) {
        return idx;
    }
    let idx = nodes.len();
    nodes.push(p);
    dedup.insert(key, idx);
    idx
}

/// Compute the bitangent of the given triangle for the given shared
/// edge: direction in the triangle plane perpendicular to the edge,
/// pointing toward the triangle's opposite vertex. Returns `None`
/// for a degenerate triangle.
///
/// Duplicates the helper in `fillet.rs` rather than importing it to
/// keep the chamfer module self-contained.
fn face_bitangent(
    mesh: &Mesh,
    graph: &EdgeGraph,
    tri_idx: usize,
    key: EdgeKey,
) -> Option<Vector3<f64>> {
    let tri = &graph.triangles[tri_idx];
    let opp = tri.v.iter().copied().find(|&v| v != key.0 && v != key.1)?;
    let q = mesh.nodes[opp];
    let p0 = mesh.nodes[key.0];
    let p1 = mesh.nodes[key.1];
    let edge = p1 - p0;
    let edge_len2 = edge.norm_squared();
    if edge_len2 < 1e-30 {
        return None;
    }
    let v = q - p0;
    let along = edge * (v.dot(&edge) / edge_len2);
    let perp = v - along;
    let len = perp.norm();
    if len < 1e-30 {
        return None;
    }
    Some(perp / len)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_mesh_returns_empty_mesh_error() {
        let m = Mesh::new("empty");
        let err = apply_chamfer(&m, 0.1, 0.5).unwrap_err();
        assert!(matches!(err, FilletError::EmptyMesh));
    }

    #[test]
    fn negative_distance_returns_bad_parameter() {
        let m = crate::edge_graph::tests::unit_cube();
        let err = apply_chamfer(&m, -0.1, 0.5).unwrap_err();
        assert!(matches!(
            err,
            FilletError::BadParameter {
                name: "distance",
                ..
            }
        ));
    }

    #[test]
    fn build_chamfer_strip_quad_shape() {
        let p0 = Vector3::zeros();
        let p1 = Vector3::x();
        let b0 = Vector3::y();
        let b1 = Vector3::z();
        let strip = build_chamfer_strip(p0, p1, b0, b1, 0.1);
        assert_eq!(strip.vertices.len(), 4);
        assert!((strip.vertices[0] - Vector3::new(0.0, 0.1, 0.0)).norm() < 1e-12);
        assert!((strip.vertices[1] - Vector3::new(1.0, 0.1, 0.0)).norm() < 1e-12);
        assert!((strip.vertices[2] - Vector3::new(0.0, 0.0, 0.1)).norm() < 1e-12);
        assert!((strip.vertices[3] - Vector3::new(1.0, 0.0, 0.1)).norm() < 1e-12);
    }

    /// Task 22 regression: full unit cube chamfer.
    /// - distance 0.1, threshold 45°
    /// - 12 edges qualify
    /// - 2 triangles per chamfer × 12 edges = 24 strip triangles
    /// - + 12 original = 36 total
    /// - bounding box unchanged
    #[test]
    fn unit_cube_chamfer_regression() {
        let m = crate::edge_graph::tests::unit_cube();
        let chamfered = apply_chamfer(&m, 0.1, std::f64::consts::FRAC_PI_4).unwrap();
        let tri_count = chamfered.element_blocks[0].count();
        assert_eq!(
            tri_count,
            12 + 12 * 2,
            "expected 36 triangles, got {tri_count}"
        );
        // Bounding box check.
        let (mut lo, mut hi) = (
            Vector3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
            Vector3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
        );
        for n in &chamfered.nodes {
            lo.x = lo.x.min(n.x);
            lo.y = lo.y.min(n.y);
            lo.z = lo.z.min(n.z);
            hi.x = hi.x.max(n.x);
            hi.y = hi.y.max(n.y);
            hi.z = hi.z.max(n.z);
        }
        assert!((lo - Vector3::zeros()).norm() < 1e-12);
        assert!((hi - Vector3::new(1.0, 1.0, 1.0)).norm() < 1e-12);
    }

    #[test]
    fn chamfer_strip_count_matches_dedup() {
        // After dedup, opposite-face chamfer vertices at e.g. the
        // bottom-front edge share with bottom-left and front-left
        // strips at the (x=0) corner. So the actual node count is
        // < 8 + 12 * 4. Just assert it's positive and less than the
        // upper bound.
        let m = crate::edge_graph::tests::unit_cube();
        let chamfered = apply_chamfer(&m, 0.1, std::f64::consts::FRAC_PI_4).unwrap();
        let upper = 8 + 12 * 4;
        assert!(
            chamfered.nodes.len() < upper,
            "expected dedup to reduce nodes below {upper}, got {}",
            chamfered.nodes.len()
        );
    }
}
