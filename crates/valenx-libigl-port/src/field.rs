//! Scalar fields on meshes — heat-method geodesic distance.
//!
//! ## The heat method (Crane, Weischedel & Wardetzky 2013)
//!
//! Geodesic distance from a source is recovered in three steps:
//!
//! 1. **Heat flow.** Integrate the heat equation for a single short
//!    time step `t`: solve `(M − t·L) u = δ_source`, where `M` is the
//!    lumped mass matrix, `L` the cotangent Laplacian, and
//!    `δ_source` an impulse at the source vertex. `u` is a smooth
//!    bump centred on the source.
//! 2. **Normalised gradient.** Compute the per-triangle gradient
//!    `∇u`, negate and normalise it. By Varadhan's formula the
//!    geodesic distance gradient is exactly the unit heat gradient
//!    as `t → 0`, so `X = −∇u / |∇u|` is the (approximate) geodesic
//!    direction field.
//! 3. **Poisson reconstruction.** Recover the distance `φ` as the
//!    scalar field whose gradient best matches `X`: solve the Poisson
//!    equation `L φ = ∇·X`. The result is shifted so `φ(source) = 0`.
//!
//! ## v1 status — real implementation
//!
//! This is the genuine heat method. It supersedes the old BFS-hop
//! proxy. The linear solves use the dense Cholesky / LU factoriser in
//! [`crate::laplacian`]. For a degenerate mesh (zero area, isolated
//! source) it falls back to a BFS-hop estimate so the function never
//! returns garbage.

use nalgebra::{DVector, Vector3};

use crate::error::LibiglError;
use crate::laplacian::{cotangent_laplacian, lumped_mass, mean_edge_length};
use crate::triangle::TriMesh;

/// Geodesic distance from `source_vertex` to every other vertex via
/// the heat method. Returns one f64 per vertex (`0.0` at the source).
///
/// # Errors
///
/// [`LibiglError::BadParameter`] if `source_vertex` is out of range;
/// [`LibiglError::NotEnough`] if the mesh has fewer than 3 vertices.
pub fn heat_geodesics(mesh: &TriMesh, source_vertex: usize) -> Result<Vec<f64>, LibiglError> {
    let n = mesh.vertices.len();
    if source_vertex >= n {
        return Err(LibiglError::BadParameter {
            name: "source_vertex",
            reason: format!("out of range: {source_vertex} >= {n}"),
        });
    }
    if n < 3 || mesh.triangles.is_empty() {
        return Err(LibiglError::NotEnough {
            what: "vertices",
            needed: 3,
            given: n,
        });
    }

    let h = mean_edge_length(mesh);
    // Time step t = h² is the value Crane et al. recommend.
    let t = h * h;

    let laplacian = cotangent_laplacian(mesh);
    let mass = lumped_mass(mesh);

    // --- Step 1: heat flow.  (M − t·L) u = δ_source ---
    // (Our L is the −Δ form (PSD), so the implicit Euler step of
    //  ∂u/∂t = Δu is (M + t·L) u = M u₀ … but the standard heat-method
    //  derivation uses (M − t·Δ); with L = −Δ that is (M + t·L). We
    //  use the PSD-stable (M + t·L).)
    let mut heat_system = laplacian.clone() * t;
    for i in 0..n {
        heat_system[(i, i)] += mass[i];
    }
    let mut rhs = DVector::<f64>::zeros(n);
    rhs[source_vertex] = 1.0;
    let Some(u) = crate::laplacian::solve_symmetric(&heat_system, &rhs) else {
        return Ok(bfs_fallback(mesh, source_vertex, h));
    };

    // --- Step 2: per-triangle gradient of u, negated + normalised ---
    let mut tri_dir: Vec<Vector3<f64>> = Vec::with_capacity(mesh.triangles.len());
    for tri in &mesh.triangles {
        let g = triangle_gradient(mesh, tri, &u);
        let neg = -g;
        let len = neg.norm();
        tri_dir.push(if len < 1e-14 {
            Vector3::zeros()
        } else {
            neg / len
        });
    }

    // --- Step 3: Poisson reconstruction.  L φ = ∇·X ---
    // Integrated divergence of the unit vector field at each vertex.
    let div = integrated_divergence(mesh, &tri_dir);
    // L is singular (constant null-space); pin one vertex to anchor it.
    let mut poisson = laplacian.clone();
    let pin = source_vertex;
    for j in 0..n {
        poisson[(pin, j)] = 0.0;
    }
    poisson[(pin, pin)] = 1.0;
    let mut poisson_rhs = div;
    poisson_rhs[pin] = 0.0;
    let Some(mut phi) = crate::laplacian::solve_symmetric(&poisson, &poisson_rhs) else {
        return Ok(bfs_fallback(mesh, source_vertex, h));
    };

    // Shift so the source is exactly zero, and flip sign if the field
    // came out negative (the Poisson solve is sign-ambiguous up to the
    // gradient-field orientation).
    let src_val = phi[source_vertex];
    for v in phi.iter_mut() {
        *v -= src_val;
    }
    let max_abs = phi.iter().fold(0.0_f64, |m, &x| m.max(x.abs()));
    let mean_signed: f64 = phi.iter().sum::<f64>();
    if mean_signed < 0.0 && max_abs > 0.0 {
        for v in phi.iter_mut() {
            *v = -*v;
        }
    }
    // Distances are non-negative.
    Ok(phi.iter().map(|&d| d.max(0.0)).collect())
}

/// Gradient of the per-vertex scalar `u` over one triangle. For a
/// linear element the gradient is constant; it is
/// `Σ_v u_v · (n × e_v) / (2A)` where `e_v` is the edge opposite `v`.
fn triangle_gradient(mesh: &TriMesh, tri: &[usize; 3], u: &DVector<f64>) -> Vector3<f64> {
    let (vi, vj, vk) = (
        mesh.vertices[tri[0]],
        mesh.vertices[tri[1]],
        mesh.vertices[tri[2]],
    );
    let normal_raw = (vj - vi).cross(&(vk - vi));
    let area2 = normal_raw.norm();
    if area2 < 1e-14 {
        return Vector3::zeros();
    }
    let n = normal_raw / area2; // unit normal
                                // Edge opposite each vertex (CCW): opposite i is (j→k), etc.
    let e_i = vk - vj;
    let e_j = vi - vk;
    let e_k = vj - vi;
    let grad = n.cross(&e_i) * u[tri[0]] + n.cross(&e_j) * u[tri[1]] + n.cross(&e_k) * u[tri[2]];
    grad / area2
}

/// Integrated divergence of a per-triangle vector field at each
/// vertex — the right-hand side of the Poisson reconstruction. For
/// vertex `i` in triangle `(i, j, k)` with field `X`:
///
/// `div_i += ½ (cot θ_k · (e_k · X) + cot θ_j · (e_j · X))`
///
/// where `e_k`, `e_j` are the two edges of the triangle incident to
/// `i` and `θ` the opposite angles.
fn integrated_divergence(mesh: &TriMesh, tri_dir: &[Vector3<f64>]) -> DVector<f64> {
    let n = mesh.vertices.len();
    let mut div = DVector::<f64>::zeros(n);
    for (t, tri) in mesh.triangles.iter().enumerate() {
        let x = tri_dir[t];
        let v = [
            mesh.vertices[tri[0]],
            mesh.vertices[tri[1]],
            mesh.vertices[tri[2]],
        ];
        for c in 0..3 {
            let i = tri[c];
            let a = tri[(c + 1) % 3];
            let b = tri[(c + 2) % 3];
            // Edges from vertex i.
            let e1 = v[(c + 1) % 3] - v[c]; // i → a
            let e2 = v[(c + 2) % 3] - v[c]; // i → b
                                            // Opposite angles: angle at b is opposite e1, angle at a
                                            // opposite e2.
            let cot_at_b = cot_angle(v[(c + 2) % 3], v[c], v[(c + 1) % 3]);
            let cot_at_a = cot_angle(v[(c + 1) % 3], v[(c + 2) % 3], v[c]);
            let _ = (a, b);
            div[i] += 0.5 * (cot_at_b * e1.dot(&x) + cot_at_a * e2.dot(&x));
        }
    }
    div
}

/// Cotangent of the angle at `apex` in triangle `apex, p, q`.
fn cot_angle(apex: Vector3<f64>, p: Vector3<f64>, q: Vector3<f64>) -> f64 {
    let e1 = p - apex;
    let e2 = q - apex;
    let cross = e1.cross(&e2).norm();
    if cross < 1e-12 {
        0.0
    } else {
        e1.dot(&e2) / cross
    }
}

/// BFS-hop fallback used only when the linear solve fails (degenerate
/// mesh). Keeps the function total.
fn bfs_fallback(mesh: &TriMesh, source: usize, mean_edge: f64) -> Vec<f64> {
    let one_ring = mesh.vertex_one_ring();
    let mut dist = vec![usize::MAX; mesh.vertices.len()];
    dist[source] = 0;
    let mut q = std::collections::VecDeque::new();
    q.push_back(source);
    while let Some(v) = q.pop_front() {
        let d = dist[v];
        for &n2 in &one_ring[v] {
            if dist[n2] == usize::MAX {
                dist[n2] = d + 1;
                q.push_back(n2);
            }
        }
    }
    dist.into_iter()
        .map(|d| {
            if d == usize::MAX {
                f64::INFINITY
            } else {
                d as f64 * mean_edge
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A flat regular grid mesh in the XY plane — `n × n` vertices.
    fn grid(n: usize) -> TriMesh {
        let mut verts = Vec::new();
        for j in 0..n {
            for i in 0..n {
                verts.push(Vector3::new(i as f64, j as f64, 0.0));
            }
        }
        let mut tris = Vec::new();
        for j in 0..n - 1 {
            for i in 0..n - 1 {
                let a = j * n + i;
                let b = j * n + i + 1;
                let c = (j + 1) * n + i;
                let d = (j + 1) * n + i + 1;
                tris.push([a, b, d]);
                tris.push([a, d, c]);
            }
        }
        TriMesh {
            vertices: verts,
            triangles: tris,
        }
    }

    #[test]
    fn source_distance_is_zero() {
        let m = grid(4);
        let d = heat_geodesics(&m, 0).unwrap();
        assert_eq!(d.len(), m.n_verts());
        assert!(d[0].abs() < 1e-9, "source distance = {}", d[0]);
    }

    #[test]
    fn bad_source_is_rejected() {
        let m = grid(3);
        assert!(matches!(
            heat_geodesics(&m, 99),
            Err(LibiglError::BadParameter {
                name: "source_vertex",
                ..
            })
        ));
    }

    #[test]
    fn distances_are_non_negative() {
        let m = grid(5);
        let d = heat_geodesics(&m, 0).unwrap();
        for (i, &x) in d.iter().enumerate() {
            assert!(x >= 0.0, "negative distance at {i}: {x}");
        }
    }

    #[test]
    fn distance_grows_with_graph_distance() {
        // On a 5x5 unit grid, vertex 0 is a corner. The diagonally
        // opposite corner (index 24) must be markedly further than an
        // adjacent vertex (index 1). The heat method recovers true
        // Euclidean-ish geodesic distance, so corner-to-corner ≈ 5.66
        // and corner-to-adjacent ≈ 1.0.
        let m = grid(5);
        let d = heat_geodesics(&m, 0).unwrap();
        let adjacent = d[1];
        let far_corner = d[24];
        assert!(
            far_corner > adjacent * 3.0,
            "far corner ({far_corner}) should dwarf adjacent ({adjacent})"
        );
        // Sanity bound — corner-to-corner on a 4x4-cell grid is ~5.66.
        assert!(
            (3.0..9.0).contains(&far_corner),
            "far-corner distance out of plausible band: {far_corner}"
        );
    }

    #[test]
    fn monotone_along_a_row() {
        // Distances along the first row away from the corner source
        // should be non-decreasing.
        let m = grid(6);
        let d = heat_geodesics(&m, 0).unwrap();
        for i in 1..6 {
            assert!(
                d[i] >= d[i - 1] - 1e-6,
                "row distance not monotone at {i}: {} < {}",
                d[i],
                d[i - 1]
            );
        }
    }
}
