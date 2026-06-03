//! Discrete mean-curvature per vertex.
//!
//! Computed from the cotangent Laplacian (Meyer / Desbrun / Schroeder
//! / Barr, 2002). The cotangent-weighted Laplace-Beltrami at a vertex
//! `v` is
//!
//! ```text
//! Lap(v) = (1 / (2*A_v)) * Σ_{j in N(v)} (cot α_ij + cot β_ij) * (v_j - v)
//! ```
//!
//! and the mean curvature normal is `H * n = -Lap(v) / 2`. We return
//! the *magnitude* `H` per vertex so callers can use it as a positive
//! weight.
//!
//! ## Boundary vertices
//!
//! The cotangent formula is the discrete Laplace-Beltrami operator
//! only for a vertex with a **closed** one-ring. A boundary vertex has
//! an open fan, so the cotangent sum does *not* vanish even on a
//! perfectly flat mesh — it would report a spurious non-zero
//! curvature. Boundary-vertex mean curvature is a genuinely different
//! quantity, so this routine detects boundary vertices (any incident
//! edge used by a single triangle) and reports `0.0` for them rather
//! than a meaningless value.

use std::collections::HashMap;

use nalgebra::Vector3;

use valenx_mesh::element::ElementType;
use valenx_mesh::Mesh;

/// Compute per-vertex discrete mean curvature `H` for every node in
/// `mesh`. Returns a vector with `mesh.nodes.len()` entries; vertices
/// touched by no Tri3 — and vertices on a mesh boundary — are 0.0.
pub fn per_vertex(mesh: &Mesh) -> Vec<f64> {
    let n = mesh.nodes.len();
    let mut lap = vec![Vector3::<f64>::zeros(); n];
    let mut area = vec![0.0_f64; n];
    // Directed-edge use counts: an undirected edge with total use 1 is
    // a boundary edge; both its endpoints are boundary vertices.
    let mut edge_use: HashMap<(usize, usize), u32> = HashMap::new();
    let bump = |a: usize, b: usize, m: &mut HashMap<(usize, usize), u32>| {
        let key = if a < b { (a, b) } else { (b, a) };
        *m.entry(key).or_insert(0) += 1;
    };

    for block in &mesh.element_blocks {
        if !matches!(block.element_type, ElementType::Tri3) {
            continue;
        }
        for tri in block.connectivity.chunks(3) {
            if tri.len() < 3 {
                continue;
            }
            let (i, j, k) = (tri[0] as usize, tri[1] as usize, tri[2] as usize);
            if i >= n || j >= n || k >= n {
                continue;
            }
            bump(i, j, &mut edge_use);
            bump(j, k, &mut edge_use);
            bump(k, i, &mut edge_use);
            let vi = mesh.nodes[i];
            let vj = mesh.nodes[j];
            let vk = mesh.nodes[k];
            // Triangle area (1/3 to each vertex barycentric).
            let face_area = 0.5 * (vj - vi).cross(&(vk - vi)).norm();
            area[i] += face_area / 3.0;
            area[j] += face_area / 3.0;
            area[k] += face_area / 3.0;
            // Cotangent of each angle, distributed to the OPPOSITE
            // edge's endpoints (standard cot-Laplacian recipe).
            let cot_i = cot_angle(vi, vj, vk);
            let cot_j = cot_angle(vj, vk, vi);
            let cot_k = cot_angle(vk, vi, vj);
            // Edge (j, k) opposite to i — contributes cot_i to both.
            lap[j] += cot_i * (vk - vj);
            lap[k] += cot_i * (vj - vk);
            // Edge (k, i) opposite to j.
            lap[k] += cot_j * (vi - vk);
            lap[i] += cot_j * (vk - vi);
            // Edge (i, j) opposite to k.
            lap[i] += cot_k * (vj - vi);
            lap[j] += cot_k * (vi - vj);
        }
    }

    // Mark boundary vertices — endpoints of any singly-used edge.
    let mut on_boundary = vec![false; n];
    for (&(a, b), &count) in &edge_use {
        if count == 1 {
            on_boundary[a] = true;
            on_boundary[b] = true;
        }
    }

    let mut h = vec![0.0_f64; n];
    for v in 0..n {
        // The cot-Laplacian H-estimate is only valid for a closed
        // one-ring; a boundary vertex's open fan gives a spurious
        // value, so report 0 there.
        if on_boundary[v] {
            continue;
        }
        let a = area[v].max(1e-12);
        let lap_v = lap[v] / (2.0 * a);
        // |mean curvature| = |Lap| / 2.
        h[v] = 0.5 * lap_v.norm();
    }
    h
}

/// `cot(angle at `a` in triangle a-b-c).
fn cot_angle(a: Vector3<f64>, b: Vector3<f64>, c: Vector3<f64>) -> f64 {
    let u = b - a;
    let v = c - a;
    let dot = u.dot(&v);
    let cross = u.cross(&v).norm();
    if cross < 1e-12 {
        return 0.0;
    }
    dot / cross
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_mesh::element::ElementBlock;

    fn unit_quad() -> Mesh {
        // Two tris forming a unit square in the z = 0 plane.
        let mut m = Mesh::new("q");
        m.nodes.push(Vector3::new(0.0, 0.0, 0.0));
        m.nodes.push(Vector3::new(1.0, 0.0, 0.0));
        m.nodes.push(Vector3::new(1.0, 1.0, 0.0));
        m.nodes.push(Vector3::new(0.0, 1.0, 0.0));
        let mut b = ElementBlock::new(ElementType::Tri3);
        b.connectivity.extend_from_slice(&[0, 1, 2, 0, 2, 3]);
        m.element_blocks.push(b);
        m.recompute_stats();
        m
    }

    #[test]
    fn flat_mesh_has_near_zero_curvature() {
        let m = unit_quad();
        let k = per_vertex(&m);
        assert_eq!(k.len(), 4);
        for v in k {
            assert!(v.abs() < 1e-6, "flat curvature should be ~0, got {v}");
        }
    }

    #[test]
    fn curvature_vector_length_matches_node_count() {
        let m = Mesh::new("empty");
        let k = per_vertex(&m);
        assert!(k.is_empty());
    }
}
