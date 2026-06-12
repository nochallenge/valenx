//! Mesh deformation — Laplacian smoothing + biharmonic deformation.

use nalgebra::Vector3;

use crate::error::LibiglError;
use crate::triangle::TriMesh;

/// Uniform-Laplacian smoothing (in-place per call). Each iteration
/// moves every vertex toward the average of its 1-ring neighbours by
/// step `lambda ∈ (0, 1]`.
pub fn laplacian_smooth(mesh: &TriMesh, iter: usize, lambda: f64) -> Result<TriMesh, LibiglError> {
    if !(0.0..=1.0).contains(&lambda) {
        return Err(LibiglError::BadParameter {
            name: "lambda",
            reason: format!("must be in [0, 1], got {lambda}"),
        });
    }
    let mut out = mesh.clone();
    let one_ring = mesh.vertex_one_ring();
    for _ in 0..iter {
        let mut new_verts = out.vertices.clone();
        for (i, neigh) in one_ring.iter().enumerate() {
            if neigh.is_empty() {
                continue;
            }
            let mut sum = Vector3::zeros();
            for &j in neigh {
                sum += out.vertices[j];
            }
            let avg = sum / neigh.len() as f64;
            new_verts[i] = out.vertices[i] * (1.0 - lambda) + avg * lambda;
        }
        out.vertices = new_verts;
    }
    Ok(out)
}

/// Biharmonic surface deformation. v1: each "handle" vertex is held
/// fixed; every other vertex is pulled toward the weighted mean of its
/// nearest handle's offset (geodesic distance approximated by 1-ring
/// hop count). This is a substantially-simplified surrogate for the
/// full bi-Laplacian solve; the API surface matches.
///
/// `handles` is a list of `(vertex_id, target_position)` pairs.
pub fn biharmonic(
    mesh: &TriMesh,
    handles: &[(usize, Vector3<f64>)],
) -> Result<TriMesh, LibiglError> {
    if handles.is_empty() {
        return Err(LibiglError::BadParameter {
            name: "handles",
            reason: "need at least one handle".into(),
        });
    }
    let one_ring = mesh.vertex_one_ring();
    // BFS hop distance from each handle.
    let mut hop: Vec<Vec<usize>> = handles
        .iter()
        .map(|(seed, _)| bfs_hops(*seed, &one_ring, mesh.vertices.len()))
        .collect();
    // Compute per-vertex offset = inverse-distance-weighted handle deltas.
    let mut out = mesh.clone();
    for (i, v) in mesh.vertices.iter().enumerate() {
        let mut wsum = 0.0;
        let mut delta_sum = Vector3::zeros();
        for (h_idx, (h_vert, target)) in handles.iter().enumerate() {
            let d = hop[h_idx].get(i).copied().unwrap_or(usize::MAX);
            if d == usize::MAX {
                continue;
            }
            let w = 1.0 / ((d as f64).max(1.0)).powi(2);
            let delta = *target - mesh.vertices[*h_vert];
            wsum += w;
            delta_sum += delta * w;
        }
        if wsum > 0.0 {
            out.vertices[i] = v + delta_sum / wsum;
        }
    }
    // Pin handles exactly.
    for (vid, target) in handles {
        out.vertices[*vid] = *target;
    }
    let _ = hop.drain(..);
    Ok(out)
}

fn bfs_hops(seed: usize, one_ring: &[std::collections::BTreeSet<usize>], n: usize) -> Vec<usize> {
    let mut dist = vec![usize::MAX; n];
    if seed >= n {
        return dist;
    }
    dist[seed] = 0;
    let mut q = std::collections::VecDeque::new();
    q.push_back(seed);
    while let Some(v) = q.pop_front() {
        let d = dist[v];
        for &n2 in &one_ring[v] {
            if dist[n2] == usize::MAX {
                dist[n2] = d + 1;
                q.push_back(n2);
            }
        }
    }
    dist
}
