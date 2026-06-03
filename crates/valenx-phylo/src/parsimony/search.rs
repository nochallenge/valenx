//! Large-parsimony heuristic tree search.
//!
//! Finding the *most*-parsimonious tree is NP-hard, so this is a
//! hill-climbing heuristic: start from a tree, generate neighbouring
//! topologies by **NNI** (nearest-neighbour interchange) and **SPR**
//! (subtree prune-and-regraft) rearrangements, keep any move that
//! lowers the [Fitch](super::fitch) parsimony score, and repeat until
//! no move improves — a local optimum.
//!
//! - **NNI** swaps the two subtrees across an internal edge; every
//!   internal edge has two NNI neighbours.
//! - **SPR** detaches a subtree and reattaches it on a different edge;
//!   the SPR neighbourhood is much larger and escapes many NNI optima.
//!
//! The search is deliberately simple — one restart from the supplied
//! starting tree, first-improvement acceptance — but it is a real
//! topology search, not a stub. The starting tree is typically a
//! neighbor-joining tree ([`crate::distance::neighbor_joining`]).

use crate::error::{PhyloError, Result};
use crate::parsimony::fitch::fitch_parsimony;
use crate::tree::{Node, NodeId, Tree};

/// Configuration for [`parsimony_search`].
#[derive(Debug, Clone)]
pub struct ParsimonySearch {
    /// Number of states in the alphabet (4 for nucleotides).
    pub n_states: u8,
    /// Whether to include SPR moves (slower, larger neighbourhood).
    pub use_spr: bool,
    /// Hard cap on hill-climb iterations (a safety bound).
    pub max_iterations: usize,
}

impl Default for ParsimonySearch {
    fn default() -> Self {
        ParsimonySearch {
            n_states: 4,
            use_spr: true,
            max_iterations: 200,
        }
    }
}

/// What [`parsimony_search`] found.
#[derive(Debug, Clone)]
pub struct SearchReport {
    /// The best (lowest-score) tree discovered.
    pub tree: Tree,
    /// Parsimony score of [`tree`](Self::tree).
    pub score: usize,
    /// Parsimony score of the starting tree.
    pub start_score: usize,
    /// Number of accepted improving moves.
    pub moves_accepted: usize,
    /// Number of hill-climb iterations performed.
    pub iterations: usize,
}

/// Runs the NNI / SPR hill-climbing large-parsimony search.
///
/// `start` is the initial topology; `alignment` maps leaf labels to
/// `u8` state rows (see [`fitch_parsimony`]).
///
/// # Errors
/// [`PhyloError`] propagated from [`fitch_parsimony`] (bad alignment,
/// missing leaf row, …).
pub fn parsimony_search(
    start: &Tree,
    alignment: &[(String, Vec<u8>)],
    cfg: &ParsimonySearch,
) -> Result<SearchReport> {
    let score_of = |t: &Tree| -> Result<usize> {
        fitch_parsimony(t, alignment, cfg.n_states).map(|r| r.score)
    };

    let start_score = score_of(start)?;
    let mut best = start.clone();
    let mut best_score = start_score;
    let mut moves_accepted = 0usize;
    let mut iterations = 0usize;

    loop {
        if iterations >= cfg.max_iterations {
            break;
        }
        iterations += 1;
        let mut improved = false;

        // Generate the neighbourhood and accept the first improvement.
        let mut neighbours = nni_neighbours(&best);
        if cfg.use_spr {
            neighbours.extend(spr_neighbours(&best));
        }
        for candidate in neighbours {
            if let Ok(s) = score_of(&candidate) {
                if s < best_score {
                    best = candidate;
                    best_score = s;
                    moves_accepted += 1;
                    improved = true;
                    break;
                }
            }
        }
        if !improved {
            break;
        }
    }

    Ok(SearchReport {
        tree: best,
        score: best_score,
        start_score,
        moves_accepted,
        iterations,
    })
}

// --- Tree-rearrangement move generators -------------------------------

/// All NNI neighbours of `tree`.
///
/// For every internal edge `(parent, child)` where `child` is also
/// internal, the two subtrees of `child` can each swap with a sibling
/// subtree of `parent` — two neighbours per qualifying edge.
pub(crate) fn nni_neighbours(tree: &Tree) -> Vec<Tree> {
    let mut out = Vec::new();
    for child in 0..tree.node_count() {
        let node = tree.node(child);
        let Some(parent) = node.parent else { continue };
        if node.is_leaf() {
            continue;
        }
        // `parent` must have another child to swap with.
        let siblings: Vec<NodeId> = tree
            .node(parent)
            .children
            .iter()
            .copied()
            .filter(|&c| c != child)
            .collect();
        if siblings.is_empty() || node.children.len() < 2 {
            continue;
        }
        let sib = siblings[0];
        // Swap `sib` with each of `child`'s children in turn.
        for &gc in &tree.node(child).children {
            if let Some(t) = swap_subtrees(tree, sib, gc) {
                out.push(t);
            }
        }
    }
    out
}

/// Swaps the subtrees rooted at `a` and `b` (exchanges their parents),
/// returning the rearranged tree, or `None` if the swap is degenerate
/// (one is an ancestor of the other, or they share a parent).
fn swap_subtrees(tree: &Tree, a: NodeId, b: NodeId) -> Option<Tree> {
    if a == b {
        return None;
    }
    let pa = tree.node(a).parent?;
    let pb = tree.node(b).parent?;
    if pa == pb {
        return None;
    }
    // Reject if one is an ancestor of the other.
    if is_ancestor(tree, a, b) || is_ancestor(tree, b, a) {
        return None;
    }
    let mut t = tree.clone();
    // Re-point parents.
    t.node_mut(a).parent = Some(pb);
    t.node_mut(b).parent = Some(pa);
    // Re-point child lists.
    let slot_a = t.node(pa).children.iter().position(|&c| c == a)?;
    let slot_b = t.node(pb).children.iter().position(|&c| c == b)?;
    t.node_mut(pa).children[slot_a] = b;
    t.node_mut(pb).children[slot_b] = a;
    t.validate().ok()?;
    Some(t)
}

/// All SPR neighbours of `tree`.
///
/// A prunable subtree (any non-root, non-root-child node) is detached;
/// its former parent is suppressed; the subtree is regrafted onto every
/// other edge. Self / ancestor regrafts are skipped. To keep the
/// neighbourhood tractable each prunable subtree contributes a regraft
/// onto every legal edge.
pub(crate) fn spr_neighbours(tree: &Tree) -> Vec<Tree> {
    let mut out = Vec::new();
    let root = tree.root();
    for prune in 0..tree.node_count() {
        if prune == root {
            continue;
        }
        let parent = tree.node(prune).parent.expect("non-root has parent");
        if parent == root {
            // Pruning a root child would leave the root unary —
            // handled by NNI elsewhere; skip for SPR simplicity.
            continue;
        }
        // Candidate regraft edges = every node except the pruned
        // subtree, its parent, and its grandparent path.
        let subtree: std::collections::HashSet<NodeId> =
            tree.descendant_set(prune);
        for regraft in 0..tree.node_count() {
            if regraft == root
                || subtree.contains(&regraft)
                || regraft == parent
            {
                continue;
            }
            // Skip the trivial no-op (regraft onto a current sibling
            // edge keeps the same tree).
            if tree.node(regraft).parent == Some(parent) {
                continue;
            }
            if let Some(t) = spr_move(tree, prune, regraft) {
                out.push(t);
            }
        }
    }
    out
}

/// Performs one SPR move: prune the subtree at `prune`, suppress its
/// parent, and regraft it on the edge above `regraft`.
fn spr_move(tree: &Tree, prune: NodeId, regraft: NodeId) -> Option<Tree> {
    let mut t = tree.clone();
    let parent = t.node(prune).parent?;
    let grandparent = t.node(parent).parent?;
    // The pruned subtree's sibling under `parent`.
    let sibling = *t
        .node(parent)
        .children
        .iter()
        .find(|&&c| c != prune)?;

    // --- Suppress `parent`: connect `sibling` directly to
    // `grandparent`, summing the two branch lengths.
    let bl_sibling = t.node(sibling).branch_length.unwrap_or(0.0);
    let bl_parent = t.node(parent).branch_length.unwrap_or(0.0);
    t.node_mut(sibling).parent = Some(grandparent);
    t.node_mut(sibling).branch_length = Some(bl_sibling + bl_parent);
    let gp_slot = t
        .node(grandparent)
        .children
        .iter()
        .position(|&c| c == parent)?;
    t.node_mut(grandparent).children[gp_slot] = sibling;

    // --- Regraft: reuse the now-free `parent` node as the new
    // attachment point on the edge above `regraft`.
    let regraft_old_parent = t.node(regraft).parent?;
    // The regraft target may have been the suppressed `parent` — but
    // that is excluded by the caller. It could be `sibling` though,
    // which is fine.
    let bl_regraft = t.node(regraft).branch_length.unwrap_or(0.0);
    let rp_slot = t
        .node(regraft_old_parent)
        .children
        .iter()
        .position(|&c| c == regraft)?;
    // `parent` now sits between `regraft_old_parent` and `regraft`.
    t.node_mut(regraft_old_parent).children[rp_slot] = parent;
    t.node_mut(parent).parent = Some(regraft_old_parent);
    t.node_mut(parent).branch_length = Some(bl_regraft / 2.0);
    t.node_mut(parent).children = vec![regraft, prune];
    t.node_mut(regraft).parent = Some(parent);
    t.node_mut(regraft).branch_length = Some(bl_regraft / 2.0);
    t.node_mut(prune).parent = Some(parent);

    t.validate().ok()?;
    Some(t)
}

/// `true` if `anc` is a (strict or improper) ancestor of `desc`.
fn is_ancestor(tree: &Tree, anc: NodeId, desc: NodeId) -> bool {
    let mut x = desc;
    loop {
        if x == anc {
            return true;
        }
        match tree.node(x).parent {
            Some(p) => x = p,
            None => return false,
        }
    }
}

// --- A small Tree extension used by the move generators ---------------

impl Tree {
    /// Set of `id` plus all of its descendants.
    pub(crate) fn descendant_set(
        &self,
        id: NodeId,
    ) -> std::collections::HashSet<NodeId> {
        let mut set = std::collections::HashSet::new();
        let mut stack = vec![id];
        while let Some(x) = stack.pop() {
            if set.insert(x) {
                stack.extend(self.node(x).children.iter().copied());
            }
        }
        set
    }
}

/// Builds a balanced caterpillar starting tree from leaf labels — a
/// convenience for callers with no NJ tree handy.
///
/// # Errors
/// [`PhyloError::Invalid`] if fewer than three labels are supplied.
pub fn caterpillar_start(labels: &[String]) -> Result<Tree> {
    if labels.len() < 3 {
        return Err(PhyloError::invalid(
            "labels",
            "need at least three taxa to build a tree",
        ));
    }
    let mut t = Tree::building();
    // First cherry: (l0, l1).
    let l0 = t.push_node(leaf_node(&labels[0]));
    let l1 = t.push_node(leaf_node(&labels[1]));
    let mut current_internal = t.push_node(Node {
        label: None,
        branch_length: None,
        parent: None,
        children: vec![l0, l1],
    });
    t.node_mut(l0).parent = Some(current_internal);
    t.node_mut(l1).parent = Some(current_internal);
    // Chain the rest.
    for label in &labels[2..] {
        let leaf = t.push_node(leaf_node(label));
        let new_internal = t.push_node(Node {
            label: None,
            branch_length: None,
            parent: None,
            children: vec![current_internal, leaf],
        });
        t.node_mut(current_internal).parent = Some(new_internal);
        t.node_mut(leaf).parent = Some(new_internal);
        current_internal = new_internal;
    }
    t.finish_building(current_internal, true)
        .map_err(|e| PhyloError::invalid_tree(e.to_string()))
}

/// A leaf node carrying `label` and a default branch length of 0.1.
fn leaf_node(label: &str) -> Node {
    Node {
        label: Some(label.to_string()),
        branch_length: Some(0.1),
        parent: None,
        children: Vec::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::read_newick;

    fn row(label: &str, states: &[u8]) -> (String, Vec<u8>) {
        (label.to_string(), states.to_vec())
    }

    #[test]
    fn nni_generates_distinct_neighbours() {
        let t = read_newick("((A,B),(C,D));").unwrap();
        let n = nni_neighbours(&t);
        assert!(!n.is_empty(), "expected NNI neighbours");
        for cand in &n {
            assert_eq!(cand.leaf_count(), 4);
            assert!(cand.validate().is_ok());
        }
    }

    #[test]
    fn spr_generates_valid_neighbours() {
        let t = read_newick("(((A,B),C),(D,E));").unwrap();
        let n = spr_neighbours(&t);
        assert!(!n.is_empty(), "expected SPR neighbours");
        for cand in &n {
            assert_eq!(cand.leaf_count(), 5);
            assert!(cand.validate().is_ok(), "SPR produced an invalid tree");
            // Leaf set is preserved.
            assert_eq!(cand.leaf_labels(), t.leaf_labels());
        }
    }

    #[test]
    fn caterpillar_start_is_valid() {
        let labels: Vec<String> = ["A", "B", "C", "D", "E"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let t = caterpillar_start(&labels).unwrap();
        assert_eq!(t.leaf_count(), 5);
        assert!(t.validate().is_ok());
        assert!(caterpillar_start(&labels[..2]).is_err());
    }

    #[test]
    fn search_finds_the_parsimonious_topology() {
        // Data favours ((A,B),(C,D)): A,B share state 0; C,D share 1
        // across many columns.
        let labels: Vec<String> = ["A", "B", "C", "D"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        // Start from the WRONG topology ((A,C),(B,D)).
        let start = read_newick("((A,C),(B,D));").unwrap();
        let aln = vec![
            row("A", &[0, 0, 0, 0, 0]),
            row("B", &[0, 0, 0, 0, 0]),
            row("C", &[1, 1, 1, 1, 1]),
            row("D", &[1, 1, 1, 1, 1]),
        ];
        let _ = labels;
        let cfg = ParsimonySearch::default();
        let report = parsimony_search(&start, &aln, &cfg).unwrap();
        // The correct topology costs 1; the wrong one costs more.
        assert!(report.score <= report.start_score);
        assert_eq!(report.score, 5); // one change per column on the
                                     // best tree
    }

    #[test]
    fn search_on_an_already_optimal_tree_does_nothing() {
        let start = read_newick("((A,B),(C,D));").unwrap();
        let aln = vec![
            row("A", &[0]),
            row("B", &[0]),
            row("C", &[1]),
            row("D", &[1]),
        ];
        let report = parsimony_search(&start, &aln, &ParsimonySearch::default()).unwrap();
        assert_eq!(report.score, report.start_score);
        assert_eq!(report.moves_accepted, 0);
    }
}
