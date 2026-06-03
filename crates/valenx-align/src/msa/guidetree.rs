//! Pairwise-distance matrix and guide-tree construction.
//!
//! Progressive multiple-sequence alignment needs an *order* in which
//! to merge sequences — the guide tree. This module builds one in two
//! steps:
//!
//! 1. [`distance_matrix`] — an all-vs-all pairwise distance matrix,
//!    where the distance of two sequences is `1 −` their global
//!    alignment identity.
//! 2. [`upgma`] / [`neighbor_joining`] — cluster the matrix into a
//!    binary [`GuideTree`]: UPGMA produces an ultrametric (rooted)
//!    tree, neighbor-joining an additive (unrooted, here rooted at the
//!    last join) one. Either works as a merge order.
//!
//! The tree is a flat arena of [`TreeNode`]s; leaves carry the index
//! of their sequence, internal nodes carry their two children.

use crate::error::{AlignError, Result};
use crate::matrix::ScoringScheme;
use crate::pairwise::global::needleman_wunsch;

/// One node of a [`GuideTree`]. Either a leaf (a sequence) or an
/// internal node joining two children.
#[derive(Clone, Debug, PartialEq)]
pub enum TreeNode {
    /// A leaf: the index of the sequence it represents.
    Leaf(usize),
    /// An internal node: indices into the tree arena of its two
    /// children, plus the height (branch context) at which they join.
    Internal {
        /// Arena index of the left child.
        left: usize,
        /// Arena index of the right child.
        right: usize,
        /// Join height — half the cluster distance for UPGMA.
        height: f64,
    },
}

/// A binary guide tree stored as a flat arena. The last node pushed is
/// the root.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct GuideTree {
    /// All nodes; the final element is the root.
    pub nodes: Vec<TreeNode>,
}

impl GuideTree {
    /// Arena index of the root, or `None` for an empty tree.
    pub fn root(&self) -> Option<usize> {
        self.nodes.len().checked_sub(1)
    }

    /// Number of leaves (sequences) in the tree.
    pub fn leaf_count(&self) -> usize {
        self.nodes
            .iter()
            .filter(|n| matches!(n, TreeNode::Leaf(_)))
            .count()
    }

    /// The sequence indices in left-to-right (post-order leaf) order —
    /// a sensible progressive-merge order to fall back on.
    pub fn leaf_order(&self) -> Vec<usize> {
        let mut order = Vec::new();
        if let Some(r) = self.root() {
            self.collect_leaves(r, &mut order);
        }
        order
    }

    fn collect_leaves(&self, idx: usize, out: &mut Vec<usize>) {
        match &self.nodes[idx] {
            TreeNode::Leaf(s) => out.push(*s),
            TreeNode::Internal { left, right, .. } => {
                self.collect_leaves(*left, out);
                self.collect_leaves(*right, out);
            }
        }
    }
}

/// A symmetric pairwise distance matrix.
#[derive(Clone, Debug, PartialEq)]
pub struct DistanceMatrix {
    /// Number of sequences `n`.
    n: usize,
    /// Flattened `n × n` distances, row-major. `d[i][i] == 0`.
    data: Vec<f64>,
}

impl DistanceMatrix {
    /// A zero `n × n` matrix.
    pub fn zeros(n: usize) -> Self {
        DistanceMatrix {
            n,
            data: vec![0.0; n * n],
        }
    }

    /// Number of sequences.
    pub fn len(&self) -> usize {
        self.n
    }

    /// `true` if the matrix is empty.
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// The distance between sequences `i` and `j`.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.n + j]
    }

    /// Sets `d(i, j) = d(j, i) = v`.
    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        self.data[i * self.n + j] = v;
        self.data[j * self.n + i] = v;
    }
}

/// Builds an all-vs-all pairwise distance matrix from sequences.
///
/// The distance of two sequences is `1 − identity`, where identity is
/// the fraction of identical columns in their global
/// (Needleman-Wunsch) alignment. Distances lie in `[0, 1]`.
///
/// Returns [`AlignError::Invalid`] for fewer than one sequence.
pub fn distance_matrix(seqs: &[&[u8]], scheme: &ScoringScheme) -> Result<DistanceMatrix> {
    if seqs.is_empty() {
        return Err(AlignError::invalid("seqs", "need >= 1 sequence"));
    }
    let n = seqs.len();
    let mut d = DistanceMatrix::zeros(n);
    for i in 0..n {
        for j in (i + 1)..n {
            let al = needleman_wunsch(seqs[i], seqs[j], scheme)?;
            let dist = 1.0 - al.percent_identity();
            d.set(i, j, dist);
        }
    }
    Ok(d)
}

/// Builds a guide tree from a distance matrix by **UPGMA**
/// (unweighted pair-group method with arithmetic mean).
///
/// At each step the two closest clusters are merged; the new cluster's
/// distance to every other is the size-weighted average of its
/// children's distances. The join height is half the merge distance,
/// giving an ultrametric tree. Returns [`AlignError::Invalid`] for an
/// empty matrix.
pub fn upgma(dm: &DistanceMatrix) -> Result<GuideTree> {
    let n = dm.len();
    if n == 0 {
        return Err(AlignError::invalid("matrix", "empty distance matrix"));
    }

    let mut tree = GuideTree::default();
    // Active clusters: (arena-node-index, member-count).
    let mut clusters: Vec<(usize, usize)> = Vec::with_capacity(n);
    for s in 0..n {
        tree.nodes.push(TreeNode::Leaf(s));
        clusters.push((tree.nodes.len() - 1, 1));
    }

    // Working distance matrix over the *current* cluster list.
    let mut dist: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| dm.get(i, j)).collect())
        .collect();

    while clusters.len() > 1 {
        // Find the closest pair.
        let (mut bi, mut bj, mut best) = (0usize, 1usize, f64::INFINITY);
        for (i, row) in dist.iter().enumerate() {
            for (j, &d) in row.iter().enumerate().skip(i + 1) {
                if d < best {
                    best = d;
                    bi = i;
                    bj = j;
                }
            }
        }

        let (ci, si) = clusters[bi];
        let (cj, sj) = clusters[bj];
        tree.nodes.push(TreeNode::Internal {
            left: ci,
            right: cj,
            height: best / 2.0,
        });
        let new_node = tree.nodes.len() - 1;
        let new_size = si + sj;

        // New cluster's distance to every surviving cluster k. The
        // UPGMA update is the size-weighted average of the two merged
        // clusters' distances.
        let mut new_row = Vec::with_capacity(clusters.len());
        for k in 0..clusters.len() {
            if k == bi || k == bj {
                continue;
            }
            let d_new = (si as f64 * dist[bi][k] + sj as f64 * dist[bj][k]) / new_size as f64;
            new_row.push((k, d_new));
        }

        // Rebuild the cluster list and distance matrix without bi, bj.
        let keep: Vec<usize> = (0..clusters.len()).filter(|&k| k != bi && k != bj).collect();
        let mut next_clusters: Vec<(usize, usize)> = keep.iter().map(|&k| clusters[k]).collect();
        let mut next_dist: Vec<Vec<f64>> = vec![vec![0.0; keep.len() + 1]; keep.len() + 1];
        for (a, &ka) in keep.iter().enumerate() {
            for (b, &kb) in keep.iter().enumerate() {
                next_dist[a][b] = dist[ka][kb];
            }
        }
        // Append the merged cluster as the last row/column.
        for (a, &ka) in keep.iter().enumerate() {
            let dval = new_row
                .iter()
                .find(|&&(k, _)| k == ka)
                .map(|&(_, v)| v)
                .unwrap_or(0.0);
            next_dist[a][keep.len()] = dval;
            next_dist[keep.len()][a] = dval;
        }
        next_clusters.push((new_node, new_size));

        clusters = next_clusters;
        dist = next_dist;
    }

    Ok(tree)
}

/// Builds a guide tree from a distance matrix by **neighbor-joining**.
///
/// NJ does not assume a molecular clock: at each step it picks the pair
/// minimising the Q-criterion (which corrects pairwise distance for
/// each cluster's average distance to all others), so it recovers the
/// correct topology even with unequal evolutionary rates. The result
/// is returned rooted at the final join. Returns
/// [`AlignError::Invalid`] for an empty matrix.
pub fn neighbor_joining(dm: &DistanceMatrix) -> Result<GuideTree> {
    let n = dm.len();
    if n == 0 {
        return Err(AlignError::invalid("matrix", "empty distance matrix"));
    }
    if n <= 2 {
        // Degenerate: UPGMA gives the same 1- or 2-leaf tree.
        return upgma(dm);
    }

    let mut tree = GuideTree::default();
    let mut clusters: Vec<usize> = Vec::with_capacity(n); // arena node ids
    for s in 0..n {
        tree.nodes.push(TreeNode::Leaf(s));
        clusters.push(tree.nodes.len() - 1);
    }
    let mut dist: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| dm.get(i, j)).collect())
        .collect();

    while clusters.len() > 2 {
        let r = clusters.len();
        // Net divergence of each cluster.
        let div: Vec<f64> = (0..r)
            .map(|i| (0..r).map(|j| dist[i][j]).sum())
            .collect();

        // Q-criterion: minimise (r-2)*d(i,j) - div[i] - div[j].
        let (mut bi, mut bj, mut best) = (0usize, 1usize, f64::INFINITY);
        for i in 0..r {
            for j in (i + 1)..r {
                let q = (r as f64 - 2.0) * dist[i][j] - div[i] - div[j];
                if q < best {
                    best = q;
                    bi = i;
                    bj = j;
                }
            }
        }

        tree.nodes.push(TreeNode::Internal {
            left: clusters[bi],
            right: clusters[bj],
            height: dist[bi][bj] / 2.0,
        });
        let new_node = tree.nodes.len() - 1;

        // Distance of the new cluster u to every other k:
        // d(u,k) = (d(i,k) + d(j,k) - d(i,j)) / 2.
        let keep: Vec<usize> = (0..r).filter(|&k| k != bi && k != bj).collect();
        let mut next_clusters: Vec<usize> = keep.iter().map(|&k| clusters[k]).collect();
        let mut next_dist = vec![vec![0.0; keep.len() + 1]; keep.len() + 1];
        for (a, &ka) in keep.iter().enumerate() {
            for (b, &kb) in keep.iter().enumerate() {
                next_dist[a][b] = dist[ka][kb];
            }
        }
        for (a, &ka) in keep.iter().enumerate() {
            let d_uk = (dist[bi][ka] + dist[bj][ka] - dist[bi][bj]) / 2.0;
            next_dist[a][keep.len()] = d_uk.max(0.0);
            next_dist[keep.len()][a] = d_uk.max(0.0);
        }
        next_clusters.push(new_node);
        clusters = next_clusters;
        dist = next_dist;
    }

    // Final two clusters join at the root.
    tree.nodes.push(TreeNode::Internal {
        left: clusters[0],
        right: clusters[1],
        height: dist[0][1] / 2.0,
    });
    Ok(tree)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::{GapCost, ScoringScheme, SubstitutionMatrix};

    fn dna_scheme() -> ScoringScheme {
        ScoringScheme::new(SubstitutionMatrix::dna_simple(1, -1), GapCost::new(0, 1))
    }

    #[test]
    fn distance_matrix_basic_properties() {
        let seqs: &[&[u8]] = &[b"ACGTACGT", b"ACGTACGT", b"TTTTTTTT"];
        let dm = distance_matrix(seqs, &dna_scheme()).unwrap();
        // Identical sequences -> distance 0.
        assert!(dm.get(0, 1).abs() < 1e-9);
        // Diagonal is 0.
        assert!(dm.get(2, 2).abs() < 1e-9);
        // Dissimilar pair -> large distance.
        assert!(dm.get(0, 2) > 0.5);
        // Symmetric.
        assert!((dm.get(0, 2) - dm.get(2, 0)).abs() < 1e-9);
    }

    #[test]
    fn upgma_builds_binary_tree() {
        let seqs: &[&[u8]] = &[b"ACGTACGT", b"ACGTACGA", b"TTTTTTTT", b"TTTTTTTA"];
        let dm = distance_matrix(seqs, &dna_scheme()).unwrap();
        let tree = upgma(&dm).unwrap();
        assert_eq!(tree.leaf_count(), 4);
        // 4 leaves + 3 internal nodes = 7 arena entries.
        assert_eq!(tree.nodes.len(), 7);
        assert!(tree.root().is_some());
    }

    #[test]
    fn upgma_groups_similar_sequences() {
        // Two tight pairs: {0,1} nearly identical, {2,3} nearly
        // identical, the pairs very different from each other.
        let seqs: &[&[u8]] = &[b"AAAAAAAAAA", b"AAAAAAAAAT", b"CCCCCCCCCC", b"CCCCCCCCCG"];
        let dm = distance_matrix(seqs, &dna_scheme()).unwrap();
        let tree = upgma(&dm).unwrap();
        // The first merge should pair within {0,1} or {2,3}.
        if let TreeNode::Internal { left, right, .. } = &tree.nodes[4] {
            let pair: Vec<&TreeNode> = vec![&tree.nodes[*left], &tree.nodes[*right]];
            let leaves: Vec<usize> = pair
                .iter()
                .filter_map(|n| match n {
                    TreeNode::Leaf(s) => Some(*s),
                    _ => None,
                })
                .collect();
            assert_eq!(leaves.len(), 2);
            let same_group = (leaves.contains(&0) && leaves.contains(&1))
                || (leaves.contains(&2) && leaves.contains(&3));
            assert!(same_group, "first merge should pair a tight group");
        } else {
            panic!("node 4 should be the first internal join");
        }
    }

    #[test]
    fn neighbor_joining_builds_tree() {
        let seqs: &[&[u8]] = &[b"ACGTACGT", b"ACGTACGA", b"TTTTTTTT", b"TTTTTTTA", b"GGGGGGGG"];
        let dm = distance_matrix(seqs, &dna_scheme()).unwrap();
        let tree = neighbor_joining(&dm).unwrap();
        assert_eq!(tree.leaf_count(), 5);
        let order = tree.leaf_order();
        assert_eq!(order.len(), 5);
        // Every sequence appears exactly once.
        let mut sorted = order.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn single_sequence_tree() {
        let seqs: &[&[u8]] = &[b"ACGT"];
        let dm = distance_matrix(seqs, &dna_scheme()).unwrap();
        let tree = upgma(&dm).unwrap();
        assert_eq!(tree.leaf_count(), 1);
        assert_eq!(tree.leaf_order(), vec![0]);
    }

    #[test]
    fn empty_matrix_rejected() {
        assert!(upgma(&DistanceMatrix::zeros(0)).is_err());
        assert!(neighbor_joining(&DistanceMatrix::zeros(0)).is_err());
    }
}
