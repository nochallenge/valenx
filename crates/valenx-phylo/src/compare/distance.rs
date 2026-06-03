//! Topological distances between two trees.
//!
//! - [`robinson_foulds`] — the Robinson-Foulds (1981) symmetric
//!   distance: the number of *bipartitions* (splits of the leaf set
//!   induced by internal edges) present in one tree but not the other.
//!   It is the standard "how different are these two topologies"
//!   metric.
//! - [`quartet_distance`] — the quartet distance: over all
//!   `C(n, 4)` four-leaf subsets, the number resolved differently by
//!   the two trees. Quartet distance is finer-grained than RF — a
//!   single misplaced taxon changes few bipartitions but many
//!   quartets.
//!
//! Both treat the trees as **unrooted** (RF on rooted trees would also
//! count the root bipartition); both require the two trees to share an
//! identical leaf-label set.

use crate::error::{PhyloError, Result};
use crate::tree::{NodeId, Tree};
use std::collections::HashSet;

/// Outcome of a Robinson-Foulds comparison.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RfResult {
    /// The symmetric distance: splits unique to one tree or the other.
    pub distance: usize,
    /// Non-trivial bipartitions in the first tree.
    pub splits_a: usize,
    /// Non-trivial bipartitions in the second tree.
    pub splits_b: usize,
    /// Bipartitions present in both trees.
    pub shared: usize,
}

impl RfResult {
    /// The normalised RF distance in `[0, 1]`: `distance / (splits_a +
    /// splits_b)`. Returns 0 when both trees are stars (no internal
    /// splits).
    pub fn normalized(&self) -> f64 {
        let denom = self.splits_a + self.splits_b;
        if denom == 0 {
            0.0
        } else {
            self.distance as f64 / denom as f64
        }
    }
}

/// A bipartition encoded as a sorted set of leaf-label indices — the
/// "smaller side" of the split, canonicalised so the two trees compare
/// equal regardless of rooting.
type Bipartition = Vec<usize>;

/// Builds a stable index for a tree's leaf labels.
///
/// Returns the sorted label list and a map from label to index.
fn leaf_index(tree: &Tree) -> (Vec<String>, std::collections::HashMap<String, usize>) {
    let mut labels = tree.leaf_labels();
    labels.sort();
    labels.dedup();
    let map = labels
        .iter()
        .cloned()
        .enumerate()
        .map(|(i, l)| (l, i))
        .collect();
    (labels, map)
}

/// Collects every non-trivial bipartition of a tree, each canonicalised
/// against the shared leaf-index space `index`.
fn bipartitions(
    tree: &Tree,
    index: &std::collections::HashMap<String, usize>,
    n_leaves: usize,
) -> Result<HashSet<Bipartition>> {
    let mut splits = HashSet::new();
    // For each internal, non-root node, the leaves below it form one
    // side of a split.
    for id in 0..tree.node_count() {
        let node = tree.node(id);
        if node.is_leaf() || node.parent.is_none() {
            continue;
        }
        let mut side: Vec<usize> = Vec::new();
        for leaf in tree.descendant_leaves(id) {
            let label = tree.node(leaf).label.as_deref().ok_or_else(|| {
                PhyloError::invalid("tree", "leaf without a label")
            })?;
            let idx = *index.get(label).ok_or_else(|| {
                PhyloError::invalid("tree", format!("leaf `{label}` not in both trees"))
            })?;
            side.push(idx);
        }
        side.sort_unstable();
        side.dedup();
        // Trivial splits (a single leaf, or all-but-one) carry no
        // information — skip them.
        if side.len() < 2 || side.len() > n_leaves - 2 {
            continue;
        }
        splits.insert(canonical(side, n_leaves));
    }
    Ok(splits)
}

/// Canonicalises a split: keeps whichever side does **not** contain
/// leaf index 0, so the two unrooted representations of the same split
/// hash equal.
fn canonical(side: Vec<usize>, n_leaves: usize) -> Bipartition {
    if side.contains(&0) {
        // Return the complement.
        let set: HashSet<usize> = side.into_iter().collect();
        (0..n_leaves).filter(|i| !set.contains(i)).collect()
    } else {
        side
    }
}

/// Computes the Robinson-Foulds distance between two trees.
///
/// # Errors
/// [`PhyloError::Invalid`] if the trees do not share an identical leaf
/// set.
pub fn robinson_foulds(a: &Tree, b: &Tree) -> Result<RfResult> {
    let (labels_a, index) = leaf_index(a);
    let (labels_b, _) = leaf_index(b);
    if labels_a != labels_b {
        return Err(PhyloError::invalid(
            "trees",
            "Robinson-Foulds needs identical leaf sets",
        ));
    }
    let n = labels_a.len();
    let sa = bipartitions(a, &index, n)?;
    let sb = bipartitions(b, &index, n)?;
    let shared = sa.intersection(&sb).count();
    let distance = (sa.len() - shared) + (sb.len() - shared);
    Ok(RfResult {
        distance,
        splits_a: sa.len(),
        splits_b: sb.len(),
        shared,
    })
}

/// Computes the quartet distance between two trees.
///
/// For each of the `C(n, 4)` leaf quadruples this checks whether the
/// two trees induce the same resolved quartet topology; the distance is
/// the count of quadruples that disagree. An unresolved quartet (a
/// polytomy in one tree) counts as a disagreement only when the other
/// tree resolves it.
///
/// This is the direct `O(n⁴)` algorithm — fine for the tree sizes a v1
/// handles; the sub-quadratic Brodal algorithm is out of scope.
///
/// # Errors
/// [`PhyloError::Invalid`] if the trees do not share an identical leaf
/// set, or have fewer than four leaves.
pub fn quartet_distance(a: &Tree, b: &Tree) -> Result<usize> {
    let (labels_a, index) = leaf_index(a);
    let (labels_b, _) = leaf_index(b);
    if labels_a != labels_b {
        return Err(PhyloError::invalid(
            "trees",
            "quartet distance needs identical leaf sets",
        ));
    }
    let n = labels_a.len();
    if n < 4 {
        return Err(PhyloError::invalid("trees", "need at least four leaves"));
    }
    // Leaf id by shared index, for both trees.
    let leaf_of = |tree: &Tree| -> Result<Vec<NodeId>> {
        let mut v = vec![usize::MAX; n];
        for leaf in tree.leaves() {
            let label = tree
                .node(leaf)
                .label
                .as_deref()
                .ok_or_else(|| PhyloError::invalid("tree", "leaf without a label"))?;
            v[index[label]] = leaf;
        }
        Ok(v)
    };
    let la = leaf_of(a)?;
    let lb = leaf_of(b)?;

    let mut disagree = 0usize;
    for i in 0..n {
        for j in (i + 1)..n {
            for k in (j + 1)..n {
                for l in (k + 1)..n {
                    let qa = quartet_topology(a, la[i], la[j], la[k], la[l]);
                    let qb = quartet_topology(b, lb[i], lb[j], lb[k], lb[l]);
                    if qa != qb {
                        disagree += 1;
                    }
                }
            }
        }
    }
    Ok(disagree)
}

/// The three possible resolved topologies of four leaves (plus
/// "unresolved").
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Quartet {
    /// `ab | cd`
    AbCd,
    /// `ac | bd`
    AcBd,
    /// `ad | bc`
    AdBc,
    /// A polytomy leaves the quartet unresolved.
    Star,
}

/// Determines which way a tree resolves the quartet `(w, x, y, z)`.
///
/// The quartet is `pq | rs` iff the path connecting `p` and `q` shares
/// no internal node with the path connecting `r` and `s`. Implemented
/// via lowest-common-ancestor depths.
fn quartet_topology(
    tree: &Tree,
    w: NodeId,
    x: NodeId,
    y: NodeId,
    z: NodeId,
) -> Quartet {
    // Depth (edge count from root) of each node.
    let depth = |mut n: NodeId| -> usize {
        let mut d = 0;
        while let Some(p) = tree.node(n).parent {
            d += 1;
            n = p;
        }
        d
    };
    // The "split" of a quartet is decided by which pairing has the
    // deepest internal connection. For each of the three pairings,
    // measure depth(lca(pair1)) + depth(lca(pair2)); the pairing that
    // keeps the two cherries separate maximises it.
    let score = |p: NodeId, q: NodeId, r: NodeId, s: NodeId| -> usize {
        depth(tree.lca(p, q)) + depth(tree.lca(r, s))
    };
    let ab_cd = score(w, x, y, z);
    let ac_bd = score(w, y, x, z);
    let ad_bc = score(w, z, x, y);
    let max = ab_cd.max(ac_bd).max(ad_bc);
    // A unique maximum gives the resolved topology; a tie is an
    // unresolved (star) quartet.
    let winners = [ab_cd, ac_bd, ad_bc]
        .iter()
        .filter(|&&v| v == max)
        .count();
    if winners != 1 {
        Quartet::Star
    } else if ab_cd == max {
        Quartet::AbCd
    } else if ac_bd == max {
        Quartet::AcBd
    } else {
        Quartet::AdBc
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::read_newick;

    #[test]
    fn identical_trees_have_zero_rf() {
        let t = read_newick("((A,B),(C,D));").unwrap();
        let rf = robinson_foulds(&t, &t).unwrap();
        assert_eq!(rf.distance, 0);
        assert_eq!(rf.shared, rf.splits_a);
        assert!((rf.normalized()).abs() < 1e-12);
    }

    #[test]
    fn different_topologies_have_nonzero_rf() {
        let t1 = read_newick("((A,B),(C,D));").unwrap();
        let t2 = read_newick("((A,C),(B,D));").unwrap();
        let rf = robinson_foulds(&t1, &t2).unwrap();
        assert!(rf.distance > 0);
        // For two conflicting four-taxon trees the RF distance is 2.
        assert_eq!(rf.distance, 2);
    }

    #[test]
    fn rf_is_symmetric() {
        let t1 = read_newick("(((A,B),C),(D,E));").unwrap();
        let t2 = read_newick("(((A,C),B),(D,E));").unwrap();
        let ab = robinson_foulds(&t1, &t2).unwrap();
        let ba = robinson_foulds(&t2, &t1).unwrap();
        assert_eq!(ab.distance, ba.distance);
    }

    #[test]
    fn rf_rejects_mismatched_leaf_sets() {
        let t1 = read_newick("((A,B),(C,D));").unwrap();
        let t2 = read_newick("((A,B),(C,E));").unwrap();
        assert!(robinson_foulds(&t1, &t2).is_err());
    }

    #[test]
    fn identical_trees_have_zero_quartet_distance() {
        let t = read_newick("(((A,B),C),(D,E));").unwrap();
        assert_eq!(quartet_distance(&t, &t).unwrap(), 0);
    }

    #[test]
    fn conflicting_trees_have_positive_quartet_distance() {
        let t1 = read_newick("((A,B),(C,D));").unwrap();
        let t2 = read_newick("((A,C),(B,D));").unwrap();
        // The single quartet ABCD is resolved oppositely.
        assert_eq!(quartet_distance(&t1, &t2).unwrap(), 1);
    }

    #[test]
    fn quartet_distance_rejects_small_trees() {
        let t = read_newick("((A,B),C);").unwrap();
        assert!(quartet_distance(&t, &t).is_err());
    }

    #[test]
    fn quartet_distance_grows_with_more_conflict() {
        // Caterpillar vs a re-arranged caterpillar — several quartets
        // flip.
        let t1 = read_newick("((((A,B),C),D),E);").unwrap();
        let t2 = read_newick("((((A,E),C),D),B);").unwrap();
        let qd = quartet_distance(&t1, &t2).unwrap();
        assert!(qd > 0 && qd <= 5);
    }
}
