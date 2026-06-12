//! Distance-matrix clustering into trees.
//!
//! Four algorithms, two families:
//!
//! - **UPGMA / WPGMA** are agglomerative average-linkage clustering.
//!   They build a *rooted, ultrametric* tree — every leaf ends up the
//!   same distance from the root — and so implicitly assume a molecular
//!   clock. UPGMA weights a merged cluster's members by cluster size;
//!   WPGMA weights the two merged clusters equally regardless of size.
//! - **Neighbor-joining / BIONJ** build a *unrooted, additive* tree
//!   with no clock assumption. Both repeatedly pick the pair that
//!   minimises the NJ Q-criterion; BIONJ additionally tracks a variance
//!   estimate to weight the distance update, which improves accuracy
//!   when substitution rates vary.
//!
//! All four consume a [`DistMatrix`] and yield a [`Tree`].

use crate::distance::matrix::DistMatrix;
use crate::error::{PhyloError, Result};
use crate::tree::{Node, NodeId, Tree};

/// Average-linkage variant for [`agglomerate`].
#[derive(Clone, Copy, PartialEq, Eq)]
enum Linkage {
    /// UPGMA — weight each cluster by its leaf count.
    Weighted,
    /// WPGMA — weight the two merged clusters equally.
    Unweighted,
}

/// Builds a rooted ultrametric tree by UPGMA.
///
/// # Errors
/// [`PhyloError::Invalid`] if the matrix has fewer than two taxa;
/// [`PhyloError::InvalidTree`] if the resulting arena fails validation.
pub fn upgma(dm: &DistMatrix) -> Result<Tree> {
    agglomerate(dm, Linkage::Weighted)
}

/// Builds a rooted ultrametric tree by WPGMA.
///
/// # Errors
/// As [`upgma`].
pub fn wpgma(dm: &DistMatrix) -> Result<Tree> {
    agglomerate(dm, Linkage::Unweighted)
}

/// A live cluster during agglomeration.
struct Cluster {
    /// Tree node representing this cluster's subtree root.
    node: NodeId,
    /// Number of original leaves under it (UPGMA weighting).
    size: usize,
    /// Height of `node` above the leaves (ultrametric depth).
    height: f64,
}

/// Shared UPGMA / WPGMA agglomeration.
fn agglomerate(dm: &DistMatrix, linkage: Linkage) -> Result<Tree> {
    let n = dm.len();
    if n < 2 {
        return Err(PhyloError::invalid(
            "matrix",
            "clustering needs at least two taxa",
        ));
    }
    let mut tree = Tree::building();
    // Start: one leaf node per taxon, each its own cluster.
    let mut clusters: Vec<Cluster> = Vec::with_capacity(n);
    for label in dm.labels() {
        let id = tree.push_node(Node {
            label: Some(label.clone()),
            branch_length: None,
            parent: None,
            children: Vec::new(),
        });
        clusters.push(Cluster {
            node: id,
            size: 1,
            height: 0.0,
        });
    }
    // Working distance matrix (mutable copy of the active rows).
    let mut d = vec![vec![0.0; n]; n];
    for (i, row) in d.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = dm.get(i, j);
        }
    }
    // `active[k]` indexes into `clusters`/`d`; collapses as merges run.
    let mut active: Vec<usize> = (0..n).collect();

    while active.len() > 1 {
        // Find the closest active pair.
        let (mut bi, mut bj, mut best) = (0usize, 1usize, f64::INFINITY);
        for a in 0..active.len() {
            for b in (a + 1)..active.len() {
                let dist = d[active[a]][active[b]];
                if dist < best {
                    best = dist;
                    bi = a;
                    bj = b;
                }
            }
        }
        let (ia, ib) = (active[bi], active[bj]);
        let ca = &clusters[ia];
        let cb = &clusters[ib];

        // New internal node joining the two clusters' subtree roots.
        let merged_height = best / 2.0;
        let new_node = tree.push_node(Node {
            label: None,
            branch_length: None,
            parent: None,
            children: vec![ca.node, cb.node],
        });
        // Branch length = height difference (ultrametric).
        let bl_a = (merged_height - ca.height).max(0.0);
        let bl_b = (merged_height - cb.height).max(0.0);
        tree.node_mut(ca.node).parent = Some(new_node);
        tree.node_mut(ca.node).branch_length = Some(bl_a);
        tree.node_mut(cb.node).parent = Some(new_node);
        tree.node_mut(cb.node).branch_length = Some(bl_b);

        let (sa, sb) = (ca.size, cb.size);
        // Update distances from the merged cluster to every other.
        for &k in &active {
            if k == ia || k == ib {
                continue;
            }
            let new_dist = match linkage {
                Linkage::Weighted => {
                    (sa as f64 * d[ia][k] + sb as f64 * d[ib][k]) / (sa + sb) as f64
                }
                Linkage::Unweighted => (d[ia][k] + d[ib][k]) / 2.0,
            };
            d[ia][k] = new_dist;
            d[k][ia] = new_dist;
        }
        // Reuse slot `ia` for the merged cluster; drop `ib`.
        clusters[ia] = Cluster {
            node: new_node,
            size: sa + sb,
            height: merged_height,
        };
        active.remove(bj); // remove the later index first
    }

    let root = clusters[active[0]].node;
    tree.finish_building(root, true)
        .map_err(|e| PhyloError::invalid_tree(e.to_string()))
}

/// Builds an unrooted additive tree by neighbor-joining (Saitou &
/// Nei 1987).
///
/// The result is stored rooted at the final internal node for the
/// arena's sake but is flagged [`Tree::rooted`]` = false`.
///
/// # Errors
/// [`PhyloError::Invalid`] if the matrix has fewer than two taxa.
pub fn neighbor_joining(dm: &DistMatrix) -> Result<Tree> {
    nj_core(dm, false)
}

/// Builds an unrooted additive tree by BIONJ (Gascuel 1997).
///
/// BIONJ is neighbor-joining with a variance-weighted distance-reduction
/// step; it is more accurate than plain NJ when evolutionary rates are
/// unequal. The result is flagged [`Tree::rooted`]` = false`.
///
/// # Errors
/// As [`neighbor_joining`].
pub fn bionj(dm: &DistMatrix) -> Result<Tree> {
    nj_core(dm, true)
}

/// Shared NJ / BIONJ implementation. `bio` selects the variance-weighted
/// BIONJ update over the unweighted NJ midpoint update.
fn nj_core(dm: &DistMatrix, bio: bool) -> Result<Tree> {
    let n = dm.len();
    if n < 2 {
        return Err(PhyloError::invalid(
            "matrix",
            "neighbor-joining needs at least two taxa",
        ));
    }
    let mut tree = Tree::building();
    // Live nodes: leaf per taxon initially.
    let mut nodes: Vec<NodeId> = Vec::with_capacity(n);
    for label in dm.labels() {
        nodes.push(tree.push_node(Node {
            label: Some(label.clone()),
            branch_length: None,
            parent: None,
            children: Vec::new(),
        }));
    }
    // Working distance + (for BIONJ) variance matrices.
    let mut d = vec![vec![0.0; n]; n];
    for (i, row) in d.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = dm.get(i, j);
        }
    }
    // BIONJ variance lambda is recomputed per merge; the variance
    // matrix starts equal to the distance matrix.
    let mut var = d.clone();
    let mut active: Vec<usize> = (0..n).collect();

    // Two-leaf base case: join directly.
    while active.len() > 2 {
        let m = active.len();
        // Net divergence r_i = sum of distances from i to all others.
        let mut r = vec![0.0; m];
        for (a, &ia) in active.iter().enumerate() {
            for &ib in &active {
                if ia != ib {
                    r[a] += d[ia][ib];
                }
            }
        }
        // Q-criterion: pick the pair minimising it.
        let (mut bi, mut bj, mut best) = (0usize, 1usize, f64::INFINITY);
        for a in 0..m {
            for b in (a + 1)..m {
                let (ia, ib) = (active[a], active[b]);
                let q = (m as f64 - 2.0) * d[ia][ib] - r[a] - r[b];
                if q < best {
                    best = q;
                    bi = a;
                    bj = b;
                }
            }
        }
        let (ia, ib) = (active[bi], active[bj]);
        // Branch lengths from the joined pair to the new node.
        let dij = d[ia][ib];
        let delta = (r[bi] - r[bj]) / (m as f64 - 2.0);
        let bl_i = (0.5 * dij + 0.5 * delta).max(0.0);
        let bl_j = (dij - bl_i).max(0.0);

        let new_node = tree.push_node(Node {
            label: None,
            branch_length: None,
            parent: None,
            children: vec![nodes[ia], nodes[ib]],
        });
        tree.node_mut(nodes[ia]).parent = Some(new_node);
        tree.node_mut(nodes[ia]).branch_length = Some(bl_i);
        tree.node_mut(nodes[ib]).parent = Some(new_node);
        tree.node_mut(nodes[ib]).branch_length = Some(bl_j);

        // BIONJ weighting: choose lambda to minimise the merged
        // node's variance. Plain NJ uses lambda = 1/2.
        let lambda = if bio {
            bionj_lambda(&var, &active, ia, ib)
        } else {
            0.5
        };
        // Reduce distances: distance from the new node u to each
        // remaining taxon k.
        for &k in &active {
            if k == ia || k == ib {
                continue;
            }
            let new_d = lambda * (d[ia][k] - bl_i) + (1.0 - lambda) * (d[ib][k] - bl_j);
            d[ia][k] = new_d.max(0.0);
            d[k][ia] = new_d.max(0.0);
            if bio {
                // Variance reduction (Gascuel 1997, eq. for v_uk).
                let new_v = lambda * var[ia][k] + (1.0 - lambda) * var[ib][k]
                    - lambda * (1.0 - lambda) * var[ia][ib];
                var[ia][k] = new_v.max(0.0);
                var[k][ia] = new_v.max(0.0);
            }
        }
        // Slot `ia` becomes the new node; drop `ib`.
        nodes[ia] = new_node;
        active.remove(bj);
    }

    // Final two clusters: they are the two ends of the unrooted tree's
    // central edge. Join them under a fresh degree-2 internal node — the
    // standard rooted representation of an unrooted tree. (The earlier
    // code made one final cluster the *parent* of the other; when an
    // m = 3 Q-criterion tie left a single leaf as one of the last two
    // clusters that turned the leaf into an internal node and the tree
    // lost a leaf.) The residual distance is split evenly across the two
    // new edges so the patristic distances are preserved.
    let (ia, ib) = (active[0], active[1]);
    let final_dist = d[ia][ib].max(0.0);
    let half = final_dist * 0.5;
    let root = tree.push_node(Node {
        label: None,
        branch_length: None,
        parent: None,
        children: vec![nodes[ia], nodes[ib]],
    });
    tree.node_mut(nodes[ia]).parent = Some(root);
    tree.node_mut(nodes[ia]).branch_length = Some(half);
    tree.node_mut(nodes[ib]).parent = Some(root);
    tree.node_mut(nodes[ib]).branch_length = Some(half);

    let mut finished = tree
        .finish_building(root, false)
        .map_err(|e| PhyloError::invalid_tree(e.to_string()))?;
    finished.set_rooted(false);
    Ok(finished)
}

/// BIONJ lambda: the variance-minimising mixing weight for the two
/// merged clusters, clamped to `[0, 1]`.
fn bionj_lambda(var: &[Vec<f64>], active: &[usize], ia: usize, ib: usize) -> f64 {
    let mut num = 0.0;
    let mut denom = 0.0;
    for &k in active {
        if k == ia || k == ib {
            continue;
        }
        num += var[ib][k] - var[ia][k];
        denom += var[ia][k] + var[ib][k];
    }
    if denom.abs() < 1e-12 {
        return 0.5;
    }
    let v_ij = var[ia][ib];
    let lambda = 0.5 + num / (2.0 * (active.len() as f64 - 2.0).max(1.0) * v_ij.max(1e-12));
    // The classical closed form can drift; clamp into the valid range.
    let _ = denom; // denom kept for clarity of the variance picture
    lambda.clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A four-taxon additive matrix with a clear ((A,B),(C,D)) shape.
    fn additive_matrix() -> DistMatrix {
        // Distances generated from a known tree:
        // A,B close; C,D close; the two cherries far apart.
        let labels: Vec<String> = ["A", "B", "C", "D"].iter().map(|s| s.to_string()).collect();
        // row-major 4x4
        let data = vec![
            0.0, 0.2, 0.9, 1.0, //
            0.2, 0.0, 1.0, 1.1, //
            0.9, 1.0, 0.0, 0.3, //
            1.0, 1.1, 0.3, 0.0, //
        ];
        DistMatrix::new(labels, data).unwrap()
    }

    /// Returns the set of leaf labels under each internal node — used
    /// to test recovered topology independently of node ids.
    fn clades(t: &Tree) -> Vec<Vec<String>> {
        let mut out = Vec::new();
        for id in 0..t.node_count() {
            if t.node(id).is_internal() {
                let mut names: Vec<String> = t
                    .descendant_leaves(id)
                    .into_iter()
                    .filter_map(|l| t.node(l).label.clone())
                    .collect();
                names.sort();
                out.push(names);
            }
        }
        out
    }

    #[test]
    fn upgma_builds_an_ultrametric_tree() {
        let dm = additive_matrix();
        let t = upgma(&dm).unwrap();
        assert_eq!(t.leaf_count(), 4);
        assert!(t.rooted);
        // Ultrametric: every leaf is the same patristic distance from
        // the root.
        let root = t.root();
        let depths: Vec<f64> = t
            .leaves()
            .iter()
            .map(|&l| t.patristic_distance(root, l))
            .collect();
        let first = depths[0];
        for d in &depths {
            assert!((d - first).abs() < 1e-9, "not ultrametric: {depths:?}");
        }
    }

    #[test]
    fn wpgma_recovers_the_cherries() {
        let dm = additive_matrix();
        let t = wpgma(&dm).unwrap();
        let cl = clades(&t);
        assert!(cl.iter().any(|c| c == &["A", "B"]));
        assert!(cl.iter().any(|c| c == &["C", "D"]));
    }

    #[test]
    fn neighbor_joining_recovers_the_true_topology() {
        let dm = additive_matrix();
        let t = neighbor_joining(&dm).unwrap();
        assert_eq!(t.leaf_count(), 4);
        assert!(!t.rooted);
        let cl = clades(&t);
        // The (A,B) cherry must appear as a clade.
        assert!(cl.iter().any(|c| c == &["A", "B"]));
    }

    #[test]
    fn bionj_recovers_the_true_topology() {
        let dm = additive_matrix();
        let t = bionj(&dm).unwrap();
        assert_eq!(t.leaf_count(), 4);
        let cl = clades(&t);
        assert!(cl.iter().any(|c| c == &["A", "B"]));
    }

    #[test]
    fn nj_branch_lengths_are_nonnegative() {
        let dm = additive_matrix();
        for t in [neighbor_joining(&dm).unwrap(), bionj(&dm).unwrap()] {
            for id in 0..t.node_count() {
                if let Some(bl) = t.node(id).branch_length {
                    assert!(bl >= 0.0, "negative branch length {bl}");
                }
            }
        }
    }

    #[test]
    fn rejects_too_few_taxa() {
        let dm = DistMatrix::new(vec!["solo".to_string()], vec![0.0]).unwrap();
        assert!(upgma(&dm).is_err());
        assert!(neighbor_joining(&dm).is_err());
    }
}
