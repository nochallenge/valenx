//! Consensus trees.
//!
//! A consensus tree summarises a *set* of trees — typically bootstrap
//! replicates or a Bayesian posterior sample — into one tree showing
//! the clades the set agrees on.
//!
//! - **Strict consensus** keeps only the bipartitions present in
//!   *every* input tree. It is the most conservative summary.
//! - **Majority-rule consensus** keeps the bipartitions present in
//!   *more than half* the input trees. Each kept clade is labelled with
//!   its support frequency.
//!
//! Both are built by the same routine: tally every input tree's
//! bipartitions, keep those clearing the threshold, and assemble a tree
//! by inserting the kept clades from largest to smallest (a larger
//! clade must be inserted before the smaller clades it contains).

use crate::error::{PhyloError, Result};
use crate::tree::{Node, NodeId, Tree};
use std::collections::HashMap;

/// Which consensus rule to apply.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConsensusKind {
    /// Keep only clades present in every input tree.
    Strict,
    /// Keep clades present in more than half the input trees.
    MajorityRule,
}

/// Builds a consensus tree from a non-empty slice of trees.
///
/// Every input tree must share an identical leaf-label set. Each
/// retained internal node is labelled with its support frequency
/// (e.g. `"0.85"`); the root is unlabelled.
///
/// # Errors
/// [`PhyloError::Invalid`] if `trees` is empty or the leaf sets differ.
pub fn consensus_tree(trees: &[Tree], kind: ConsensusKind) -> Result<Tree> {
    if trees.is_empty() {
        return Err(PhyloError::invalid("trees", "no trees supplied"));
    }
    // Shared leaf-label index from the first tree.
    let mut labels = trees[0].leaf_labels();
    labels.sort();
    labels.dedup();
    let n = labels.len();
    if n < 2 {
        return Err(PhyloError::invalid("trees", "trees need at least two leaves"));
    }
    let index: HashMap<String, usize> = labels
        .iter()
        .cloned()
        .enumerate()
        .map(|(i, l)| (l, i))
        .collect();
    for t in trees {
        let mut ls = t.leaf_labels();
        ls.sort();
        ls.dedup();
        if ls != labels {
            return Err(PhyloError::invalid(
                "trees",
                "consensus needs identical leaf sets",
            ));
        }
    }

    // Tally bipartitions. A bipartition is a sorted bitmask-ish vector
    // of leaf indices; canonicalised so it never contains index 0.
    let mut counts: HashMap<Vec<usize>, usize> = HashMap::new();
    for t in trees {
        for split in tree_bipartitions(t, &index, n)? {
            *counts.entry(split).or_insert(0) += 1;
        }
    }

    let total = trees.len();
    let threshold = match kind {
        ConsensusKind::Strict => total,
        // "more than half" — strict majority.
        ConsensusKind::MajorityRule => total / 2 + 1,
    };

    // Kept clades, each as (leaf-index set, support frequency).
    let mut kept: Vec<(Vec<usize>, f64)> = counts
        .into_iter()
        .filter(|(_, c)| *c >= threshold)
        .map(|(split, c)| (split, c as f64 / total as f64))
        .collect();
    // Insert largest clades first.
    kept.sort_by(|a, b| b.0.len().cmp(&a.0.len()));

    build_consensus(&labels, &kept)
}

/// Collects a tree's non-trivial **rooted clades** as sorted leaf-index
/// vectors — the descendant-leaf set of every non-root internal node.
///
/// The input trees here are rooted (consensus of rooted trees), so the
/// clades are kept rooted and are NOT folded to unrooted bipartitions:
/// for `((A,B),(C,D))` the clades `{A,B}` and `{C,D}` are distinct and
/// must both be counted (folding by "remove index 0" would collapse the
/// two into one split, double-count it, and lose one of the two clades
/// when the tree is rebuilt). The returned vector is deduplicated so a
/// tree contributes each clade at most once.
fn tree_bipartitions(
    tree: &Tree,
    index: &HashMap<String, usize>,
    n: usize,
) -> Result<Vec<Vec<usize>>> {
    let mut seen: std::collections::HashSet<Vec<usize>> =
        std::collections::HashSet::new();
    for id in 0..tree.node_count() {
        let node = tree.node(id);
        if node.is_leaf() || node.parent.is_none() {
            continue;
        }
        let mut side: Vec<usize> = Vec::new();
        for leaf in tree.descendant_leaves(id) {
            let label = tree
                .node(leaf)
                .label
                .as_deref()
                .ok_or_else(|| PhyloError::invalid("tree", "leaf without a label"))?;
            side.push(*index.get(label).ok_or_else(|| {
                PhyloError::invalid("tree", "leaf not in the shared set")
            })?);
        }
        side.sort_unstable();
        side.dedup();
        // Keep only non-trivial clades (≥ 2 leaves, not the whole set).
        if side.len() < 2 || side.len() >= n {
            continue;
        }
        seen.insert(side);
    }
    Ok(seen.into_iter().collect())
}

/// Assembles a consensus tree from leaf labels and the kept clades
/// (largest-first), each clade carrying its support frequency.
fn build_consensus(labels: &[String], clades: &[(Vec<usize>, f64)]) -> Result<Tree> {
    let mut tree = Tree::building();
    // One leaf node per label.
    let leaf_nodes: Vec<NodeId> = labels
        .iter()
        .map(|l| {
            tree.push_node(Node {
                label: Some(l.clone()),
                branch_length: None,
                parent: None,
                children: Vec::new(),
            })
        })
        .collect();
    // The root, initially the parent of every leaf (a star tree).
    let root = tree.push_node(Node {
        label: None,
        branch_length: None,
        parent: None,
        children: leaf_nodes.clone(),
    });
    for &leaf in &leaf_nodes {
        tree.node_mut(leaf).parent = Some(root);
    }
    // `cluster_node[i]` = the node, a direct child of `i`'s current
    // common parent, that leaf `i` descends through. Initially every
    // leaf is a direct child of the star root, so each leaf's cluster
    // node is the LEAF ITSELF — not the root. (Initialising these to the
    // root collapsed every clade's distinct-cluster set to a single node
    // and the clade was skipped, so the consensus came out as a star.)
    let mut cluster_node: Vec<NodeId> = leaf_nodes.clone();

    // Insert each clade by re-parenting its leaves' current cluster
    // members under a new internal node.
    for (members, support) in clades {
        // The set of distinct current cluster nodes the clade's
        // members belong to.
        let mut parent_of: Option<NodeId> = None;
        let mut ok = true;
        for &m in members {
            let p = tree.node(cluster_node[m]).parent.unwrap_or(cluster_node[m]);
            let _ = p;
        }
        // All members must currently share one parent for the clade
        // to be insertable as a nested group (it does, because clades
        // are inserted largest-first and the kept set is compatible
        // for strict/majority consensus). Their parent is the LCA of
        // their cluster nodes.
        let members_clusters: Vec<NodeId> =
            members.iter().map(|&m| cluster_node[m]).collect();
        // Find their common parent.
        for &c in &members_clusters {
            let par = tree.node(c).parent;
            match (parent_of, par) {
                (None, Some(p)) => parent_of = Some(p),
                (Some(p0), Some(p)) if p0 == p => {}
                (Some(_), Some(_)) => {
                    ok = false;
                    break;
                }
                _ => {
                    ok = false;
                    break;
                }
            }
        }
        let Some(common_parent) = parent_of else {
            continue;
        };
        if !ok {
            // Incompatible clade — skip it (keeps the tree valid even
            // if the caller passed an unfiltered set).
            continue;
        }
        // The children of `common_parent` that belong to this clade.
        let clade_children: Vec<NodeId> = members_clusters
            .iter()
            .copied()
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();
        if clade_children.len() < 2 {
            continue;
        }
        // New internal node grouping those children.
        let new_node = tree.push_node(Node {
            label: Some(format!("{support:.2}")),
            branch_length: None,
            parent: Some(common_parent),
            children: clade_children.clone(),
        });
        // Detach the clade children from `common_parent`, attach the
        // new node in their place.
        tree.node_mut(common_parent)
            .children
            .retain(|c| !clade_children.contains(c));
        tree.node_mut(common_parent).children.push(new_node);
        for &c in &clade_children {
            tree.node_mut(c).parent = Some(new_node);
        }
        // Every leaf in this clade now sits one level deeper.
        for &m in members {
            cluster_node[m] = new_node;
        }
    }

    tree.finish_building(root, false)
        .map_err(|e| PhyloError::invalid_tree(e.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::read_newick;

    #[test]
    fn strict_consensus_of_identical_trees_is_that_tree() {
        let t = read_newick("((A,B),(C,D));").unwrap();
        let c = consensus_tree(&[t.clone(), t.clone(), t], ConsensusKind::Strict).unwrap();
        assert_eq!(c.leaf_count(), 4);
        // The (A,B) and (C,D) clades survive.
        let cl: Vec<Vec<String>> = (0..c.node_count())
            .filter(|&id| c.node(id).is_internal())
            .map(|id| {
                let mut v: Vec<String> = c
                    .descendant_leaves(id)
                    .into_iter()
                    .filter_map(|l| c.node(l).label.clone())
                    .collect();
                v.sort();
                v
            })
            .collect();
        assert!(cl.iter().any(|x| x == &["A", "B"]));
        assert!(cl.iter().any(|x| x == &["C", "D"]));
    }

    #[test]
    fn strict_consensus_drops_conflicting_clades() {
        // Two conflicting trees => the strict consensus is a star.
        let t1 = read_newick("((A,B),(C,D));").unwrap();
        let t2 = read_newick("((A,C),(B,D));").unwrap();
        let c = consensus_tree(&[t1, t2], ConsensusKind::Strict).unwrap();
        // No internal node other than the root.
        let internal = (0..c.node_count())
            .filter(|&id| c.node(id).is_internal())
            .count();
        assert_eq!(internal, 1, "expected a star tree");
    }

    #[test]
    fn majority_rule_keeps_a_clade_seen_in_most_trees() {
        // (A,B) appears in two of three trees => kept by majority rule.
        let t1 = read_newick("((A,B),(C,D));").unwrap();
        let t2 = read_newick("((A,B),(C,D));").unwrap();
        let t3 = read_newick("((A,C),(B,D));").unwrap();
        let c = consensus_tree(&[t1, t2, t3], ConsensusKind::MajorityRule).unwrap();
        let cl: Vec<Vec<String>> = (0..c.node_count())
            .filter(|&id| c.node(id).is_internal())
            .map(|id| {
                let mut v: Vec<String> = c
                    .descendant_leaves(id)
                    .into_iter()
                    .filter_map(|l| c.node(l).label.clone())
                    .collect();
                v.sort();
                v
            })
            .collect();
        assert!(cl.iter().any(|x| x == &["A", "B"]), "clades: {cl:?}");
    }

    #[test]
    fn majority_clades_are_labelled_with_support() {
        let t1 = read_newick("((A,B),(C,D));").unwrap();
        let t2 = read_newick("((A,B),(C,D));").unwrap();
        let t3 = read_newick("((A,B),(C,D));").unwrap();
        let c = consensus_tree(&[t1, t2, t3], ConsensusKind::MajorityRule).unwrap();
        // The (A,B) node should be labelled "1.00".
        let ab = c.find("1.00");
        assert!(ab.is_some(), "expected a 100%-support label");
    }

    #[test]
    fn rejects_empty_or_mismatched_input() {
        assert!(consensus_tree(&[], ConsensusKind::Strict).is_err());
        let t1 = read_newick("((A,B),(C,D));").unwrap();
        let t2 = read_newick("((A,B),(C,E));").unwrap();
        assert!(consensus_tree(&[t1, t2], ConsensusKind::Strict).is_err());
    }
}

