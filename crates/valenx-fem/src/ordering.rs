//! **Fill-reducing matrix reordering** for the sparse direct solvers
//! (Phase 24.8).
//!
//! ## Why this exists
//!
//! The native solvers factorise the global stiffness matrix with
//! [`nalgebra_sparse::factorization::CscCholesky`], whose documentation
//! states plainly: *"the current implementation performs no fill-in
//! reduction."* The factor `L` of a sparse SPD matrix `K` is only as
//! sparse as the **bandwidth** of `K` allows — and the bandwidth
//! depends entirely on how the nodes are *numbered*.
//!
//! The structured Tet4 box mesh happens to number its nodes
//! `x`-fastest, which gives a small bandwidth, so the original solvers
//! factorise quickly. But a quadratic-tetrahedron mesh appends its
//! mid-edge nodes *after* all the corner nodes, scattering coupled DOFs
//! far apart in the numbering — a huge bandwidth, a near-dense factor,
//! and a factorisation that is orders of magnitude slower than it needs
//! to be.
//!
//! Every production FEA solver fixes this with a **fill-reducing
//! reordering**. This module ships **Reverse Cuthill-McKee (RCM)** — a
//! classical, robust bandwidth-minimising permutation — and the small
//! amount of glue to apply a node permutation to a CSC system, solve,
//! and unpermute the answer. RCM typically shrinks the Tet10
//! factorisation from minutes to a fraction of a second.

use std::collections::VecDeque;

use nalgebra_sparse::{CooMatrix, CscMatrix};

/// Compute a **Reverse Cuthill-McKee** permutation of an `n`-vertex
/// graph given as an adjacency list.
///
/// Returns `perm` where `perm[new_index] = old_index`: row/column
/// `new_index` of the permuted matrix is row/column `perm[new_index]`
/// of the original. Applying RCM to the sparsity graph of a symmetric
/// matrix produces a numbering with a small bandwidth, which keeps the
/// Cholesky factor sparse.
///
/// The algorithm: repeatedly pick a low-degree start vertex in each
/// remaining connected component, breadth-first traverse it visiting
/// each frontier's neighbours in ascending-degree order (Cuthill-
/// McKee), then **reverse** the whole ordering (the "reverse" of RCM,
/// which further reduces the factor's fill profile).
pub fn reverse_cuthill_mckee(adjacency: &[Vec<usize>]) -> Vec<usize> {
    let n = adjacency.len();
    let mut visited = vec![false; n];
    let mut order: Vec<usize> = Vec::with_capacity(n);
    let degree: Vec<usize> = adjacency.iter().map(|a| a.len()).collect();

    // Process every connected component.
    while order.len() < n {
        // Start vertex: the unvisited vertex of minimum degree (a
        // standard, cheap pseudo-peripheral substitute).
        let start = (0..n)
            .filter(|&v| !visited[v])
            .min_by_key(|&v| degree[v])
            .expect("an unvisited vertex must remain");

        let mut queue: VecDeque<usize> = VecDeque::new();
        visited[start] = true;
        queue.push_back(start);
        while let Some(v) = queue.pop_front() {
            order.push(v);
            // Visit neighbours in ascending-degree order.
            let mut nbrs: Vec<usize> = adjacency[v]
                .iter()
                .copied()
                .filter(|&w| !visited[w])
                .collect();
            nbrs.sort_by_key(|&w| degree[w]);
            for w in nbrs {
                if !visited[w] {
                    visited[w] = true;
                    queue.push_back(w);
                }
            }
        }
    }

    // Reverse — the "R" of RCM.
    order.reverse();
    order
}

/// Build the **node adjacency graph** of a structural FE system from
/// its global DOF coupling.
///
/// `n_nodes` is the node count, `dof_per_node` the DOFs each node
/// carries (3 for a continuum mesh, 6 for a beam frame), and
/// `element_nodes` is one entry per element listing the node indices
/// it connects. Two nodes are adjacent iff some element references
/// both — exactly the coupling pattern of the assembled stiffness
/// matrix's node blocks.
pub fn node_adjacency(n_nodes: usize, element_nodes: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let mut sets: Vec<std::collections::BTreeSet<usize>> =
        vec![std::collections::BTreeSet::new(); n_nodes];
    for elem in element_nodes {
        for &a in elem {
            for &b in elem {
                if a != b && a < n_nodes && b < n_nodes {
                    sets[a].insert(b);
                }
            }
        }
    }
    sets.into_iter().map(|s| s.into_iter().collect()).collect()
}

/// The **bandwidth** of a node adjacency graph under a given
/// numbering — the largest `|i − j|` over every adjacent node pair.
///
/// A small bandwidth means a sparse Cholesky factor. Exposed so a
/// caller / test can confirm RCM actually shrank it.
pub fn graph_bandwidth(adjacency: &[Vec<usize>]) -> usize {
    let mut bw = 0;
    for (i, nbrs) in adjacency.iter().enumerate() {
        for &j in nbrs {
            bw = bw.max(i.abs_diff(j));
        }
    }
    bw
}

/// The bandwidth a node adjacency graph would have *after* applying a
/// node permutation `perm` (`perm[new] = old`).
pub fn permuted_bandwidth(adjacency: &[Vec<usize>], perm: &[usize]) -> usize {
    // inverse[old] = new.
    let mut inverse = vec![0usize; perm.len()];
    for (new, &old) in perm.iter().enumerate() {
        inverse[old] = new;
    }
    let mut bw = 0;
    for (i, nbrs) in adjacency.iter().enumerate() {
        for &j in nbrs {
            bw = bw.max(inverse[i].abs_diff(inverse[j]));
        }
    }
    bw
}

/// Expand a **node** permutation into the matching **DOF** permutation.
///
/// If node `new` maps to old node `perm[new]`, then DOF
/// `dof_per_node·new + c` maps to old DOF `dof_per_node·perm[new] + c`.
/// Returns `dof_perm` with `dof_perm[new_dof] = old_dof`.
pub fn node_perm_to_dof_perm(node_perm: &[usize], dof_per_node: usize) -> Vec<usize> {
    let mut dof_perm = Vec::with_capacity(node_perm.len() * dof_per_node);
    for &old_node in node_perm {
        for c in 0..dof_per_node {
            dof_perm.push(dof_per_node * old_node + c);
        }
    }
    dof_perm
}

/// Symmetrically permute a CSC matrix: return `P·A·Pᵀ` where `P` is
/// the permutation with `perm[new] = old`.
///
/// Entry `(i,j)` of the result is entry `(perm[i], perm[j])` of `a`.
pub fn permute_csc(a: &CscMatrix<f64>, perm: &[usize]) -> CscMatrix<f64> {
    let n = a.nrows();
    debug_assert_eq!(n, perm.len(), "permutation length must match");
    // inverse[old] = new.
    let mut inverse = vec![0usize; n];
    for (new, &old) in perm.iter().enumerate() {
        inverse[old] = new;
    }
    let mut coo = CooMatrix::<f64>::new(n, n);
    for (r, c, v) in a.triplet_iter() {
        coo.push(inverse[r], inverse[c], *v);
    }
    CscMatrix::from(&coo)
}

/// Apply a DOF permutation to a dense right-hand-side vector — return
/// `b_perm` with `b_perm[new] = b[perm[new]]`.
pub fn permute_vector(b: &[f64], perm: &[usize]) -> Vec<f64> {
    perm.iter().map(|&old| b[old]).collect()
}

/// Undo a DOF permutation on a solution vector — return `x` with
/// `x[perm[new]] = x_perm[new]`.
pub fn unpermute_vector(x_perm: &[f64], perm: &[usize]) -> Vec<f64> {
    let mut x = vec![0.0; x_perm.len()];
    for (new, &old) in perm.iter().enumerate() {
        x[old] = x_perm[new];
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rcm_shrinks_a_deliberately_bad_numbering() {
        // A path graph 0-1-2-...-n is already bandwidth 1; but number
        // it so the two ends of every edge are far apart and RCM must
        // recover a small bandwidth.
        let n = 20;
        // Adjacency of a simple chain, but we feed the chain in a
        // shuffled vertex labelling: vertex i is chain position
        // (i*7) mod n (a permutation since gcd(7,20)=1).
        let pos = |i: usize| (i * 7) % n;
        let mut adj = vec![Vec::new(); n];
        for chain in 0..n - 1 {
            let (a, b) = (
                (0..n).find(|&v| pos(v) == chain).unwrap(),
                (0..n).find(|&v| pos(v) == chain + 1).unwrap(),
            );
            adj[a].push(b);
            adj[b].push(a);
        }
        let bad = graph_bandwidth(&adj);
        let perm = reverse_cuthill_mckee(&adj);
        let good = permuted_bandwidth(&adj, &perm);
        // A chain reorders to bandwidth 1; RCM must get close.
        assert!(
            good < bad,
            "RCM did not reduce the bandwidth: {bad} → {good}"
        );
        assert!(good <= 2, "RCM bandwidth {good} should be ~1 for a chain");
    }

    #[test]
    fn rcm_is_a_valid_permutation() {
        let adj = vec![
            vec![1, 2],
            vec![0, 2, 3],
            vec![0, 1],
            vec![1],
            vec![], // an isolated vertex (its own component)
        ];
        let perm = reverse_cuthill_mckee(&adj);
        assert_eq!(perm.len(), 5);
        let mut sorted = perm.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, vec![0, 1, 2, 3, 4], "perm must be a bijection");
    }

    #[test]
    fn node_adjacency_from_elements_is_correct() {
        // Two triangles sharing an edge: 0-1-2 and 1-2-3.
        let elems = vec![vec![0, 1, 2], vec![1, 2, 3]];
        let adj = node_adjacency(4, &elems);
        assert_eq!(adj[0], vec![1, 2]);
        assert_eq!(adj[3], vec![1, 2]);
        // Node 1 touches everything except itself.
        assert_eq!(adj[1], vec![0, 2, 3]);
    }

    #[test]
    fn permute_csc_round_trips_through_a_solve() {
        // P·A·Pᵀ permuted, then a permuted RHS solved and unpermuted,
        // must reproduce the original system's solution.
        use nalgebra_sparse::factorization::CscCholesky;
        // A small SPD matrix.
        let mut coo = CooMatrix::<f64>::new(4, 4);
        for i in 0..4 {
            coo.push(i, i, 4.0);
        }
        coo.push(0, 1, 1.0);
        coo.push(1, 0, 1.0);
        coo.push(2, 3, -1.0);
        coo.push(3, 2, -1.0);
        coo.push(1, 2, 0.5);
        coo.push(2, 1, 0.5);
        let a = CscMatrix::from(&coo);
        let b = [1.0, 2.0, 3.0, 4.0];

        // Direct solve.
        let chol = CscCholesky::factor(&a).unwrap();
        let rhs = nalgebra::DVector::from_row_slice(&b);
        let x_direct = chol.solve(&rhs);

        // Permuted solve.
        let perm = vec![2, 0, 3, 1];
        let a_p = permute_csc(&a, &perm);
        let b_p = permute_vector(&b, &perm);
        let chol_p = CscCholesky::factor(&a_p).unwrap();
        let x_p = chol_p.solve(&nalgebra::DVector::from_row_slice(&b_p));
        let x_back = unpermute_vector(x_p.column(0).as_slice(), &perm);

        for i in 0..4 {
            assert!(
                (x_back[i] - x_direct[(i, 0)]).abs() < 1e-9,
                "permuted solve disagreed at {i}"
            );
        }
    }

    #[test]
    fn dof_perm_expands_node_perm() {
        let node_perm = vec![2, 0, 1];
        let dof_perm = node_perm_to_dof_perm(&node_perm, 3);
        // New node 0 = old node 2 → DOFs 6,7,8.
        assert_eq!(&dof_perm[0..3], &[6, 7, 8]);
        // New node 1 = old node 0 → DOFs 0,1,2.
        assert_eq!(&dof_perm[3..6], &[0, 1, 2]);
    }
}
