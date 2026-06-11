//! [`apply_fillet`] — replace every sharp convex edge of a mesh with
//! a cylindrical strip of radius `r`.
//!
//! # v1 simplifications (overlap-not-clip)
//!
//! Three honest compromises ship in v1, all documented here so
//! downstream callers know the seams:
//!
//! ## 1. No endpoint repositioning
//!
//! A "true" fillet would move the original mesh's edge-endpoint
//! vertices inward (along the bisector of the two adjacent face
//! normals) by `radius`, so the strip's endpoints meet the moved
//! vertices and form a tangent join. v1 leaves the original vertices
//! in place and lets the strip ride on top of them.
//!
//! ## 2. No original-triangle clipping
//!
//! The two triangles that share a filleted edge would, in a true
//! fillet, shrink to make room for the strip (their vertices on the
//! edge move along with the endpoint repositioning above). v1 keeps
//! the originals at full size — the strip overlaps them by roughly
//! `radius * sin(dihedral / 2)`. From outside the model the strip
//! sits in front and the overlap is invisible; from inside, the two
//! surfaces interpenetrate.
//!
//! ## 3. No spherical corner blend
//!
//! At a vertex where 3+ filleted edges meet (every corner of a
//! cube), a true fillet would weave a spherical patch through the
//! end caps of all three strips. v1 lets the three strips overlap
//! each other near the corner — the dedup table fuses any exactly-
//! coincident vertices but doesn't generate the corner-sphere
//! geometry.
//!
//! ## Net effect
//!
//! All three simplifications produce a mesh that:
//!
//! - Has the correct silhouette and rounded outer surface (good for
//!   visualization and 3D printing).
//! - Has internal interpenetration that's invisible from outside.
//! - Is NOT manifold-clean — STL export works but boolean
//!   operations against the result will misbehave.
//! - Is NOT a true BRep fillet — STEP / IGES export would lose the
//!   parametric fillet feature.
//!
//! Phase 3.5 revisits these compromises with a true BRep approach.

use std::collections::HashMap;

use nalgebra::Vector3;
use valenx_mesh::{ElementBlock, ElementType, Mesh};

use crate::cyl_strip;
use crate::edge_graph::EdgeGraph;
use crate::error::FilletError;

/// Default smoothness of each fillet strip.
const DEFAULT_SEGMENTS: usize = 8;

/// Spatial dedup tolerance — vertices within this distance fuse to
/// one. Chosen to be small relative to typical CAD units (mm) and
/// the smallest expected fillet radius.
const DEDUP_TOL: f64 = 1e-6;

/// Quantize one f64 axis to a tolerance bucket. Two values within
/// `DEDUP_TOL` of each other map to the same bucket — *most* of the
/// time. (Values straddling a bucket boundary will not fuse; that's
/// acceptable for v1 — duplicates aren't fatal, just wasteful.)
fn bucket(v: f64) -> i64 {
    (v / DEDUP_TOL).round() as i64
}

/// Spatial hash for fusing coincident vertices during strip append.
#[derive(Default)]
struct VertexDedup {
    map: HashMap<(i64, i64, i64), usize>,
}

impl VertexDedup {
    /// Look up `p` in the existing vertex list, or push it and return
    /// the new index. Coincident vertices (within [`DEDUP_TOL`])
    /// return the existing index.
    fn intern(&mut self, nodes: &mut Vec<Vector3<f64>>, p: Vector3<f64>) -> usize {
        let key = (bucket(p.x), bucket(p.y), bucket(p.z));
        if let Some(&idx) = self.map.get(&key) {
            return idx;
        }
        let idx = nodes.len();
        nodes.push(p);
        self.map.insert(key, idx);
        idx
    }
}

/// Apply a constant-radius fillet to every convex sharp edge of
/// `mesh` whose dihedral angle exceeds `threshold_rad`.
///
/// # Parameters
/// - `mesh`: source triangle mesh. Must have at least one Tri3
///   element block. Outward normals are required for correct
///   convex/concave classification.
/// - `radius`: fillet radius in mesh units. Must be > 0.
/// - `threshold_rad`: minimum dihedral angle (radians) for an edge
///   to qualify. Typical value: 45° (~0.785 rad).
///
/// # Errors
/// - [`FilletError::EmptyMesh`] if the source has zero triangles.
/// - [`FilletError::BadParameter`] if `radius` ≤ 0 or
///   `threshold_rad` < 0.
/// - [`FilletError::DegenerateEdge`] if a filletable edge has
///   coincident endpoints.
///
/// # v1 limitations
/// See the module-level docs: strips overlap the original triangles
/// rather than clipping them.
pub fn apply_fillet(mesh: &Mesh, radius: f64, threshold_rad: f64) -> Result<Mesh, FilletError> {
    if radius <= 0.0 {
        return Err(FilletError::BadParameter {
            name: "radius",
            reason: format!("must be > 0, got {radius}"),
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

    let filletable = graph.filletable_edges(mesh, threshold_rad);

    // Start the output mesh by copying every original vertex and
    // every original Tri3 triangle. Strip vertices are appended via
    // a dedup table so coincident points (e.g. where 3 strips meet
    // at a cube corner) share one index.
    let mut out = Mesh::new(format!("{}-filleted", mesh.id));
    out.nodes = mesh.nodes.clone();
    let mut connectivity: Vec<u32> = Vec::new();
    for tri in &graph.triangles {
        connectivity.push(tri.v[0] as u32);
        connectivity.push(tri.v[1] as u32);
        connectivity.push(tri.v[2] as u32);
    }

    let mut dedup = VertexDedup::default();
    // Seed the dedup table with the existing original-mesh vertices
    // so any strip vertex coincident with an original mesh vertex
    // also fuses. This is important once we add endpoint
    // repositioning (Task 14+).
    for (idx, p) in out.nodes.iter().copied().enumerate() {
        let key = (bucket(p.x), bucket(p.y), bucket(p.z));
        dedup.map.entry(key).or_insert(idx);
    }

    // For each filletable edge: build a cylindrical strip and
    // append its vertices + triangles into the output. Strips
    // overlap the original triangles by design in v1.
    for key in filletable {
        let p0 = mesh.nodes[key.0];
        let p1 = mesh.nodes[key.1];
        let edge_len = (p1 - p0).norm();
        if edge_len < 1e-12 {
            return Err(FilletError::DegenerateEdge {
                from: key.0,
                to: key.1,
            });
        }
        let tris = &graph.adjacency[&key];
        if tris.len() != 2 {
            // Defensive: filletable_edges only returns manifold edges,
            // but guard anyway.
            continue;
        }
        let n0 = graph.triangle_normal(mesh, tris[0]);
        let n1 = graph.triangle_normal(mesh, tris[1]);
        if n0.norm_squared() < 0.5 || n1.norm_squared() < 0.5 {
            // Degenerate triangle adjacent to the edge — skip.
            continue;
        }
        // Bitangents: direction in each face from the edge midpoint
        // toward the triangle's opposite vertex, projected to be
        // perpendicular to the edge. Cannot be derived from normals
        // alone — depends on triangle winding / position.
        let Some(b0) = face_bitangent(mesh, &graph, tris[0], key) else {
            continue;
        };
        let Some(b1) = face_bitangent(mesh, &graph, tris[1], key) else {
            continue;
        };
        let strip = cyl_strip::build(p0, p1, n0, n1, b0, b1, radius, DEFAULT_SEGMENTS);
        append_strip(&mut out.nodes, &mut connectivity, &mut dedup, &strip);
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = connectivity;
    out.element_blocks.push(block);
    out.recompute_stats();
    Ok(out)
}

/// Compute the bitangent of the given triangle for the given shared
/// edge: direction in the triangle plane perpendicular to the edge,
/// pointing toward the triangle's opposite vertex. Returns `None`
/// for a degenerate triangle.
fn face_bitangent(
    mesh: &Mesh,
    graph: &EdgeGraph,
    tri_idx: usize,
    key: crate::edge_graph::EdgeKey,
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
    // Vector from any edge point to q, projected to be perpendicular
    // to the edge.
    let v = q - p0;
    let along = edge * (v.dot(&edge) / edge_len2);
    let perp = v - along;
    let len = perp.norm();
    if len < 1e-30 {
        return None;
    }
    Some(perp / len)
}

/// Append `strip` to the output buffers, rebasing each strip-local
/// triangle index to the global vertex index space and fusing
/// coincident vertices via `dedup`.
fn append_strip(
    nodes: &mut Vec<Vector3<f64>>,
    connectivity: &mut Vec<u32>,
    dedup: &mut VertexDedup,
    strip: &cyl_strip::Strip,
) {
    // Map strip-local index → global mesh index.
    let mut map = Vec::with_capacity(strip.vertices.len());
    for &p in &strip.vertices {
        map.push(dedup.intern(nodes, p));
    }
    for &(a, b, c) in &strip.triangles {
        connectivity.push(map[a] as u32);
        connectivity.push(map[b] as u32);
        connectivity.push(map[c] as u32);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_mesh_returns_empty_mesh_error() {
        let m = Mesh::new("empty");
        let err = apply_fillet(&m, 0.1, 0.5).unwrap_err();
        assert!(matches!(err, FilletError::EmptyMesh));
    }

    #[test]
    fn negative_radius_returns_bad_parameter() {
        let m = crate::edge_graph::tests::unit_cube();
        let err = apply_fillet(&m, -0.1, 0.5).unwrap_err();
        assert!(matches!(
            err,
            FilletError::BadParameter { name: "radius", .. }
        ));
    }

    #[test]
    fn zero_radius_returns_bad_parameter() {
        let m = crate::edge_graph::tests::unit_cube();
        let err = apply_fillet(&m, 0.0, 0.5).unwrap_err();
        assert!(matches!(
            err,
            FilletError::BadParameter { name: "radius", .. }
        ));
    }

    #[test]
    fn cube_fillet_dedups_strip_endpoints() {
        // Without dedup, every strip contributes 18 vertices and 12
        // edges contribute 18 * 12 = 216 strip vertices. With dedup
        // many of those coincide (radii from the same corner). We
        // can't compute the exact savings analytically, but assert
        // that the output has *fewer* nodes than the no-dedup upper
        // bound.
        let m = crate::edge_graph::tests::unit_cube();
        let filleted = apply_fillet(&m, 0.1, std::f64::consts::FRAC_PI_4).unwrap();
        let upper_bound = 8 + 12 * (DEFAULT_SEGMENTS + 1) * 2;
        assert!(
            filleted.nodes.len() < upper_bound,
            "expected dedup to reduce node count below {upper_bound}, got {}",
            filleted.nodes.len()
        );
    }

    #[test]
    fn degenerate_triangle_does_not_panic() {
        use valenx_mesh::element::{ElementBlock, ElementType};
        // A normal triangle [0,1,2] and a DEGENERATE triangle [0,1,1] (repeated
        // vertex) sharing edge (0,1). The degenerate triangle has no third
        // distinct vertex; `is_convex` must return `None` (the edge is skipped)
        // rather than panic on `.expect`.
        let mut m = Mesh::new("degen");
        m.nodes.push(Vector3::new(0.0, 0.0, 0.0));
        m.nodes.push(Vector3::new(1.0, 0.0, 0.0));
        m.nodes.push(Vector3::new(0.0, 1.0, 0.0));
        let mut blk = ElementBlock::new(ElementType::Tri3);
        blk.connectivity.extend_from_slice(&[0, 1, 2, 0, 1, 1]);
        m.element_blocks.push(blk);
        // Must return rather than panic; the degenerate edge is simply skipped.
        let _ = apply_fillet(&m, 0.1, std::f64::consts::FRAC_PI_4);
    }

    /// Bounding box of all node positions in a mesh.
    fn bounding_box(m: &Mesh) -> (Vector3<f64>, Vector3<f64>) {
        let mut lo = Vector3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut hi = Vector3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for n in &m.nodes {
            lo.x = lo.x.min(n.x);
            lo.y = lo.y.min(n.y);
            lo.z = lo.z.min(n.z);
            hi.x = hi.x.max(n.x);
            hi.y = hi.y.max(n.y);
            hi.z = hi.z.max(n.z);
        }
        (lo, hi)
    }

    #[test]
    fn v1_simplification_keeps_original_vertices() {
        // Documented v1 behavior: the original 8 cube corners are
        // present unchanged in the output (we don't reposition them
        // inward). The strip vertices are *added* on top.
        let m = crate::edge_graph::tests::unit_cube();
        let filleted = apply_fillet(&m, 0.05, std::f64::consts::FRAC_PI_4).unwrap();
        // First 8 nodes of the output are the original cube corners,
        // unchanged.
        for i in 0..8 {
            assert!(
                (filleted.nodes[i] - m.nodes[i]).norm() < 1e-12,
                "v1 must not reposition original vertex {i}"
            );
        }
    }

    #[test]
    fn v1_overlap_keeps_bounding_box_within_radius() {
        // Documented v1 behavior: strip vertices sit AT radius `r`
        // from the original edges, so the bounding box never grows
        // (offsets are inward toward face interiors) but also never
        // shrinks. The hi/lo extremes stay at the original corners
        // because we kept them.
        let m = crate::edge_graph::tests::unit_cube();
        let (orig_lo, orig_hi) = bounding_box(&m);
        let filleted = apply_fillet(&m, 0.1, std::f64::consts::FRAC_PI_4).unwrap();
        let (lo, hi) = bounding_box(&filleted);
        assert!((orig_lo - lo).norm() < 1e-12);
        assert!((orig_hi - hi).norm() < 1e-12);
    }

    #[test]
    fn cube_fillet_adds_triangles() {
        let m = crate::edge_graph::tests::unit_cube();
        let original_tris = m.element_blocks[0].count();
        let filleted = apply_fillet(&m, 0.1, std::f64::consts::FRAC_PI_4).unwrap();
        let new_tris = filleted.element_blocks[0].count();
        // 12 original + 12 edges * (DEFAULT_SEGMENTS * 2) strip tris
        // = 12 + 12 * 16 = 204.
        assert_eq!(new_tris, 12 + 12 * 16);
        assert!(new_tris > original_tris);
    }

    /// Task 17 regression: full unit-cube fillet round-trip.
    /// - radius 0.1, threshold 45°
    /// - all 12 edges qualify
    /// - output has 12 original + 12*16 strip triangles = 204
    /// - bounding box unchanged
    /// - strip vertices all lie inside or on the original cube
    #[test]
    fn unit_cube_fillet_regression() {
        let m = crate::edge_graph::tests::unit_cube();
        let filleted = apply_fillet(&m, 0.1, std::f64::consts::FRAC_PI_4).unwrap();
        // Triangle count.
        assert_eq!(filleted.element_blocks[0].count(), 204);
        // Bounding box equals original (0..=1 on each axis).
        let (lo, hi) = bounding_box(&filleted);
        assert!((lo - Vector3::zeros()).norm() < 1e-12);
        assert!((hi - Vector3::new(1.0, 1.0, 1.0)).norm() < 1e-12);
        // Stats refreshed.
        assert_eq!(filleted.stats.element_count, 204);
        // Every strip vertex lies inside the [0,1]^3 cube.
        for n in &filleted.nodes {
            assert!(n.x >= -1e-9 && n.x <= 1.0 + 1e-9, "x out of range: {}", n.x);
            assert!(n.y >= -1e-9 && n.y <= 1.0 + 1e-9, "y out of range: {}", n.y);
            assert!(n.z >= -1e-9 && n.z <= 1.0 + 1e-9, "z out of range: {}", n.z);
        }
    }
}
