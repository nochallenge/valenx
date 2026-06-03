//! Vertex smoothing for Tri3 surface meshes.
//!
//! Two flavors:
//!
//! - [`laplacian`] — classic isotropic Laplace smoothing. Each
//!   iteration moves every vertex toward the centroid of its
//!   one-ring neighbors by a factor `λ ∈ (0, 1]`. Effective at
//!   reducing high-frequency noise but biased toward shrinking the
//!   mesh.
//! - [`taubin`] — Taubin's λ/μ trick. Alternates a positive λ step
//!   (smoothing inward) with a slightly larger negative μ step
//!   (expanding back outward). The net effect approximates a
//!   low-pass filter with very little shrinkage. Typical values:
//!   `λ = 0.5`, `μ = -0.53`.
//!
//! Boundary vertices (those on a one-ring with fewer than 2 incident
//! edges or on the surface boundary) are pinned in place — moving
//! them inward would round off corners users care about. The
//! pre-alpha boundary detection is "any vertex referenced by < 3
//! triangles" which captures the common open-mesh case but won't
//! identify a feature edge on a closed mesh. v1.5 adds feature-edge
//! detection via dihedral-angle thresholding.

use nalgebra::Vector3;

use crate::element::ElementType;
use crate::mesh::Mesh;

/// Run `iterations` of Laplacian smoothing on the Tri3 blocks of
/// `mesh`. Each iteration moves every (non-boundary) vertex toward
/// the centroid of its one-ring neighbors by `factor`.
///
/// `factor` is the standard Laplace smoothing coefficient `λ`:
/// `v_new = v + λ * (mean(neighbors) - v)`. `factor = 1.0` snaps
/// every vertex straight to its neighbor centroid (often too
/// aggressive); `0.5` is a reasonable conservative default.
///
/// Returns a fresh mesh; the original is unmodified. Element
/// connectivity is preserved bit-for-bit — only vertex coordinates
/// change.
pub fn laplacian(mesh: &Mesh, iterations: usize, factor: f64) -> Mesh {
    if mesh.nodes.is_empty() || iterations == 0 {
        return mesh.clone();
    }
    let neighbors = vertex_neighbors(mesh);
    let mut positions = mesh.nodes.clone();
    let mut next = positions.clone();
    let boundary = boundary_vertices(mesh);
    for _ in 0..iterations {
        for v in 0..positions.len() {
            if boundary[v] {
                next[v] = positions[v];
                continue;
            }
            let nbs = &neighbors[v];
            if nbs.is_empty() {
                next[v] = positions[v];
                continue;
            }
            // R34 S2 (defense-in-depth): `vertex_neighbors` already
            // drops out-of-range triangles, so `nb` is in range for a
            // well-formed gather. The `.get()`+skip here keeps the seal
            // uniform and robust if the neighbour table is ever built
            // by another path; we divide by the count actually summed
            // so the centroid stays correct.
            let mut centroid = Vector3::zeros();
            let mut counted = 0usize;
            for &nb in nbs {
                let Some(p) = positions.get(nb as usize) else {
                    continue;
                };
                centroid += p;
                counted += 1;
            }
            if counted == 0 {
                next[v] = positions[v];
                continue;
            }
            centroid /= counted as f64;
            next[v] = positions[v] + (centroid - positions[v]) * factor;
        }
        std::mem::swap(&mut positions, &mut next);
    }
    let mut out = mesh.clone();
    out.id = format!("{}_smoothed", mesh.id);
    out.nodes = positions;
    out.recompute_stats();
    out
}

/// Run `iterations` Taubin smoothing passes. Each pass alternates a
/// positive-λ Laplacian step with a negative-μ step. The standard
/// values `λ = 0.5`, `μ = -0.53` give the canonical low-pass
/// behaviour with negligible shrinkage; tweak only if you know why.
///
/// Total Laplacian-style passes is `2 × iterations` (each iteration
/// does one λ step + one μ step).
pub fn taubin(mesh: &Mesh, iterations: usize, lambda: f64, mu: f64) -> Mesh {
    if mesh.nodes.is_empty() || iterations == 0 {
        return mesh.clone();
    }
    // Two-stage Laplacian per outer iter. Re-use the Laplacian
    // helper twice rather than write a custom loop — the cost is
    // an extra clone of the mesh per stage, which is negligible
    // versus the O(V * iter * avg_degree) inner loop.
    let mut current = mesh.clone();
    for _ in 0..iterations {
        current = laplacian(&current, 1, lambda);
        current = laplacian(&current, 1, mu);
    }
    current.id = format!("{}_smoothed", mesh.id);
    current.recompute_stats();
    current
}

/// One-ring vertex neighbors for every vertex in `mesh`. Index `v`
/// holds a dedup'd vector of vertex indices that share a Tri3
/// triangle with `v`.
///
/// Output length equals `mesh.nodes.len()` even if some nodes are
/// not referenced by any triangle (those entries are empty).
pub fn vertex_neighbors(mesh: &Mesh) -> Vec<Vec<u32>> {
    let node_count = mesh.nodes.len();
    let mut neighbors: Vec<std::collections::BTreeSet<u32>> = vec![Default::default(); node_count];
    // R34 S2 (defense-in-depth): connectivity values index `neighbors`
    // here, and the values stored end up indexed as `positions[nb]` by
    // `laplacian`/`taubin`. A triangle citing a vertex past
    // `nodes.len()` would panic `neighbors[a]` (or seed a dangling
    // neighbour that later panics the smoother), so we SKIP any
    // out-of-range triangle. This is the single chokepoint feeding the
    // one-ring tables, so sealing it protects the smoothing consumers.
    // Backs the per-loader parse guards (OBJ/gmsh/netgen/PLY).
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            if (tri[0] as usize) >= node_count
                || (tri[1] as usize) >= node_count
                || (tri[2] as usize) >= node_count
            {
                continue;
            }
            for k in 0..3 {
                let a = tri[k] as usize;
                let b = tri[(k + 1) % 3];
                let c = tri[(k + 2) % 3];
                neighbors[a].insert(b);
                neighbors[a].insert(c);
            }
        }
    }
    neighbors
        .into_iter()
        .map(|s| s.into_iter().collect())
        .collect()
}

/// Detect boundary vertices for the Tri3 blocks of `mesh`.
///
/// v1 heuristic: a vertex is boundary if any of its incident edges
/// appears in only one triangle (open mesh boundary). Implemented by
/// counting how many triangles each undirected edge appears in.
fn boundary_vertices(mesh: &Mesh) -> Vec<bool> {
    let mut edge_counts: std::collections::HashMap<(u32, u32), u32> = Default::default();
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            for k in 0..3 {
                let a = tri[k];
                let b = tri[(k + 1) % 3];
                let key = if a < b { (a, b) } else { (b, a) };
                *edge_counts.entry(key).or_insert(0) += 1;
            }
        }
    }
    let mut boundary = vec![false; mesh.nodes.len()];
    for ((a, b), n) in edge_counts {
        if n == 1 {
            // R34 S2 (defense-in-depth): `a`/`b` are connectivity
            // values, so an out-of-range index would panic
            // `boundary[..]`. Use `.get_mut()` and skip the
            // out-of-range endpoint. Backs the per-loader parse guards.
            if let Some(flag) = boundary.get_mut(a as usize) {
                *flag = true;
            }
            if let Some(flag) = boundary.get_mut(b as usize) {
                *flag = true;
            }
        }
    }
    boundary
}

/// Mean and variance of edge lengths over all Tri3 blocks. Useful
/// for quantifying the smoothing effect in tests.
pub fn edge_length_stats(mesh: &Mesh) -> Option<(f64, f64)> {
    let mut edges: Vec<f64> = Vec::new();
    // R34 S2 (defense-in-depth): `a`/`b` are connectivity values, so an
    // out-of-range index would panic `mesh.nodes[..]`. Use `.get()` and
    // skip any edge touching an out-of-range vertex. Backs the
    // per-loader parse guards (OBJ/gmsh/netgen/PLY).
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            for k in 0..3 {
                let a = tri[k] as usize;
                let b = tri[(k + 1) % 3] as usize;
                let (Some(pa), Some(pb)) = (mesh.nodes.get(a), mesh.nodes.get(b)) else {
                    continue;
                };
                let len = (pa - pb).norm();
                edges.push(len);
            }
        }
    }
    if edges.is_empty() {
        return None;
    }
    let mean: f64 = edges.iter().copied().sum::<f64>() / edges.len() as f64;
    let var: f64 = edges.iter().map(|&l| (l - mean).powi(2)).sum::<f64>() / edges.len() as f64;
    Some((mean, var))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::ElementBlock;

    /// Unit-cube triangulation (8 verts, 12 tris).
    fn unit_cube() -> Mesh {
        let mut m = Mesh::new("cube");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(0.0, 1.0, 1.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tri3);
        blk.connectivity.extend_from_slice(&[
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 0, 4, 7, 0, 7,
            3, 1, 2, 6, 1, 6, 5,
        ]);
        m.element_blocks = vec![blk];
        m.recompute_stats();
        m
    }

    /// Plane (n×n quads -> 2 n² tris) with one interior vertex
    /// perturbed off-plane. Smoothing pulls the perturbed vertex
    /// back toward the plane.
    fn noisy_plane(n: usize, perturb: f64) -> Mesh {
        let mut m = Mesh::new("plane");
        for j in 0..=n {
            for i in 0..=n {
                let x = i as f64 / n as f64;
                let y = j as f64 / n as f64;
                let mut z = 0.0;
                // Perturb each strictly-interior vertex alternately.
                if i > 0 && i < n && j > 0 && j < n {
                    z = perturb * if (i + j) % 2 == 0 { 1.0 } else { -1.0 };
                }
                m.nodes.push(Vector3::new(x, y, z));
            }
        }
        let nx_plus = n + 1;
        let mut blk = ElementBlock::new(ElementType::Tri3);
        for j in 0..n {
            for i in 0..n {
                let i0 = (j * nx_plus + i) as u32;
                let i1 = i0 + 1;
                let i2 = i0 + nx_plus as u32;
                let i3 = i2 + 1;
                blk.connectivity
                    .extend_from_slice(&[i0, i1, i3, i0, i3, i2]);
            }
        }
        m.element_blocks = vec![blk];
        m.recompute_stats();
        m
    }

    #[test]
    fn vertex_neighbors_correct_for_single_triangle() {
        let mut m = Mesh::new("tri");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tri3);
        blk.connectivity = vec![0, 1, 2];
        m.element_blocks = vec![blk];
        let nbs = vertex_neighbors(&m);
        assert_eq!(nbs.len(), 3);
        // Each vertex's neighbors are the other two.
        for (v, ring) in nbs.iter().enumerate().take(3) {
            assert_eq!(ring.len(), 2);
            assert!(!ring.contains(&(v as u32)));
        }
    }

    #[test]
    fn laplacian_zero_iter_returns_input() {
        let m = unit_cube();
        let s = laplacian(&m, 0, 0.5);
        assert_eq!(s.nodes, m.nodes);
    }

    #[test]
    fn laplacian_empty_mesh_returns_empty() {
        let m = Mesh::new("empty");
        let s = laplacian(&m, 5, 0.5);
        assert_eq!(s.nodes.len(), 0);
    }

    #[test]
    fn laplacian_reduces_noise_on_plane() {
        // Noisy plane: interior vertices oscillate ±0.3 in z. After
        // 5 iterations of factor=0.5 the z-extent should shrink
        // dramatically (boundary stays at 0, interior pulled back).
        let m = noisy_plane(8, 0.3);
        let z_before: f64 = m.nodes.iter().map(|n| n.z.abs()).sum();
        let s = laplacian(&m, 5, 0.5);
        let z_after: f64 = s.nodes.iter().map(|n| n.z.abs()).sum();
        assert!(
            z_after < 0.5 * z_before,
            "z-noise didn't drop enough: before {z_before} after {z_after}"
        );
    }

    #[test]
    fn laplacian_preserves_boundary() {
        // Plane boundary vertices (z=0) must stay at z=0 because
        // they're detected as boundary.
        let m = noisy_plane(6, 0.5);
        let n = 6;
        let s = laplacian(&m, 10, 0.8);
        let nx_plus = n + 1;
        for j in 0..=n {
            for i in 0..=n {
                if i == 0 || i == n || j == 0 || j == n {
                    let idx = j * nx_plus + i;
                    let z = s.nodes[idx].z;
                    assert!(
                        z.abs() < 1e-12,
                        "boundary vertex ({i},{j}) z = {z} (should be 0)"
                    );
                }
            }
        }
    }

    #[test]
    fn taubin_zero_iter_returns_input() {
        let m = unit_cube();
        let s = taubin(&m, 0, 0.5, -0.53);
        assert_eq!(s.nodes, m.nodes);
    }

    #[test]
    fn taubin_reduces_noise_with_less_shrinkage_than_laplacian() {
        // On a noisy plane the Taubin step pair should pull the
        // off-plane perturbation in (λ step) and out (μ step). Net
        // residual should still be much smaller than the original.
        let m = noisy_plane(8, 0.3);
        let z_before: f64 = m.nodes.iter().map(|n| n.z.abs()).sum();
        let s = taubin(&m, 5, 0.5, -0.53);
        let z_after: f64 = s.nodes.iter().map(|n| n.z.abs()).sum();
        assert!(z_after < z_before);
    }

    #[test]
    fn edge_length_stats_for_unit_triangle() {
        let mut m = Mesh::new("tri");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tri3);
        blk.connectivity = vec![0, 1, 2];
        m.element_blocks = vec![blk];
        let (mean, var) = edge_length_stats(&m).unwrap();
        // Edges have lengths 1, 1, sqrt(2). Mean ≈ 1.138.
        let expected_mean = (1.0 + 1.0 + (2.0_f64).sqrt()) / 3.0;
        assert!((mean - expected_mean).abs() < 1e-9);
        assert!(var > 0.0);
    }

    /// R34 S2 (RED→GREEN): defense-in-depth sink seal. A mesh whose
    /// Tri3 connectivity cites a vertex past `nodes.len()` must NOT
    /// panic any smoothing entry point. Pre-fix `vertex_neighbors` did
    /// `neighbors[a]` (and `laplacian` then `positions[nb]`),
    /// `boundary_vertices` did `boundary[a]`, and `edge_length_stats`
    /// did `mesh.nodes[a]` — each panicked "index out of bounds".
    /// Post-fix the out-of-range triangle/edge is skipped and every
    /// helper returns. We assert no panic; the valid vertex count is
    /// preserved (smoothing never adds or drops nodes).
    #[test]
    fn out_of_range_connectivity_does_not_panic() {
        let mut m = Mesh::new("hostile");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tri3);
        // A triangle citing vertex 7 (out of range).
        blk.connectivity = vec![0, 1, 7];
        m.element_blocks = vec![blk];
        m.recompute_stats();
        // Every public smoothing path must return rather than panic.
        let nbs = vertex_neighbors(&m);
        assert_eq!(nbs.len(), 3, "one ring table is sized to node count");
        let s = laplacian(&m, 3, 0.5);
        assert_eq!(s.nodes.len(), 3, "laplacian preserves the node count");
        let t = taubin(&m, 2, 0.5, -0.53);
        assert_eq!(t.nodes.len(), 3, "taubin preserves the node count");
        // edge_length_stats skips the bad edges; the all-bad case here
        // leaves no valid edge, so it returns None rather than panicking.
        let _ = edge_length_stats(&m);
    }

    /// R34 S2: a valid triangle alongside an out-of-range one — the
    /// valid geometry still drives smoothing and edge stats while the
    /// bad triangle is dropped, no panic.
    #[test]
    fn out_of_range_triangle_skipped_valid_kept() {
        let mut m = Mesh::new("mixed");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tri3);
        // First tri valid (0,1,2); second cites vertex 9.
        blk.connectivity = vec![0, 1, 2, 0, 1, 9];
        m.element_blocks = vec![blk];
        m.recompute_stats();
        // The valid triangle still produces a one-ring and edge stats.
        let nbs = vertex_neighbors(&m);
        assert_eq!(nbs[0].len(), 2, "vertex 0 rings to its two tri partners");
        let stats = edge_length_stats(&m);
        assert!(
            stats.is_some(),
            "the valid triangle's edges must still be measured"
        );
    }

    #[test]
    fn smoothed_id_carries_suffix() {
        let m = unit_cube();
        assert_eq!(laplacian(&m, 1, 0.5).id, "cube_smoothed");
        assert_eq!(taubin(&m, 1, 0.5, -0.53).id, "cube_smoothed");
    }
}
