//! Edge-cutting utilities used for unfolding.

use crate::error::LibiglError;
use crate::triangle::TriMesh;

/// Random edge cuts — returns `n_cuts` sets, each containing one
/// random edge encoded as `(min_vertex_id, max_vertex_id)`.
///
/// Deterministic via a tiny linear-congruential generator seeded by
/// `0x9e3779b97f4a7c15 ^ n_cuts` — so re-running with the same
/// `(mesh, n_cuts)` returns the same cut set; tests can assert on it.
pub fn random_cuts(mesh: &TriMesh, n_cuts: usize) -> Result<Vec<Vec<(usize, usize)>>, LibiglError> {
    if mesh.triangles.is_empty() {
        return Err(LibiglError::NotEnough {
            what: "triangles",
            needed: 1,
            given: 0,
        });
    }
    // Enumerate unique edges.
    let mut edges: std::collections::BTreeSet<(usize, usize)> = std::collections::BTreeSet::new();
    for tri in &mesh.triangles {
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let edge = if a < b { (a, b) } else { (b, a) };
            edges.insert(edge);
        }
    }
    let edges: Vec<(usize, usize)> = edges.into_iter().collect();
    if edges.is_empty() {
        return Err(LibiglError::NotEnough {
            what: "edges",
            needed: 1,
            given: 0,
        });
    }
    let mut rng_state: u64 = 0x9e3779b97f4a7c15u64 ^ (n_cuts as u64);
    let mut out = Vec::with_capacity(n_cuts);
    for _ in 0..n_cuts {
        // Tiny xorshift.
        rng_state ^= rng_state << 13;
        rng_state ^= rng_state >> 7;
        rng_state ^= rng_state << 17;
        let idx = (rng_state as usize) % edges.len();
        out.push(vec![edges[idx]]);
    }
    Ok(out)
}
