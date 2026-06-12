//! Structure distance — base-pair distance and tree-edit distance.
//!
//! Two metrics for "how different are these two secondary
//! structures":
//!
//! - **Base-pair distance** ([`base_pair_distance`]) — the size of
//!   the symmetric difference of the two pair sets: the number of
//!   pairs in exactly one of the structures. Fast, `O(n)`, the metric
//!   the inverse-folding walk optimises.
//! - **Tree-edit distance** ([`tree_edit_distance`]) — RNA secondary
//!   structures are *trees* (loops nested in loops); the tree-edit
//!   distance is the minimum number of node insert / delete / relabel
//!   operations to turn one structure tree into the other. This is
//!   the Zhang-Shasha (1989) ordered-tree-edit DP, the metric behind
//!   RNAdistance.

use crate::error::{Result, RnaStructError};
use crate::structure::Structure;

/// The base-pair distance between two structures of equal length.
///
/// Equal to `|P₁ △ P₂|` — pairs present in exactly one structure.
///
/// # Errors
/// [`RnaStructError::Structure`] if the two structures differ in
/// length.
pub fn base_pair_distance(a: &Structure, b: &Structure) -> Result<usize> {
    if a.len() != b.len() {
        return Err(RnaStructError::structure(format!(
            "cannot compare structures of length {} and {}",
            a.len(),
            b.len()
        )));
    }
    let mut diff = 0;
    for i in 0..a.len() {
        // count a pair once, at its 5' base
        if let Some(j) = a.partner(i) {
            if i < j && b.partner(i) != Some(j) {
                diff += 1;
            }
        }
        if let Some(j) = b.partner(i) {
            if i < j && a.partner(i) != Some(j) {
                diff += 1;
            }
        }
    }
    Ok(diff)
}

/// The fraction of positions whose pairing state matches between two
/// structures — a normalised similarity in `[0, 1]`.
///
/// # Errors
/// [`RnaStructError::Structure`] on a length mismatch.
pub fn pairing_agreement(a: &Structure, b: &Structure) -> Result<f64> {
    if a.len() != b.len() {
        return Err(RnaStructError::structure("structures differ in length"));
    }
    if a.is_empty() {
        return Ok(1.0);
    }
    let agree = (0..a.len())
        .filter(|&i| a.partner(i) == b.partner(i))
        .count();
    Ok(agree as f64 / a.len() as f64)
}

// ---------------------------------------------------------------------
// Tree-edit distance (Zhang-Shasha)
// ---------------------------------------------------------------------

/// A node of the ordered structure tree.
///
/// Each base pair becomes a `P` node; each unpaired base an `U` node.
/// Children appear in 5′→3′ order. A virtual root holds the exterior
/// loop's elements.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum NodeKind {
    /// A base pair.
    Pair,
    /// An unpaired base.
    Unpaired,
    /// The virtual exterior root.
    Root,
}

/// An ordered tree in left-to-right post-order, the form Zhang-Shasha
/// consumes.
struct OrderedTree {
    /// `kind[v]` — the label of post-order node `v`.
    kind: Vec<NodeKind>,
    /// `lmld[v]` — the leftmost-leaf descendant of node `v`
    /// (post-order index).
    lmld: Vec<usize>,
    /// The post-order indices of the key roots.
    keyroots: Vec<usize>,
}

/// Builds the ordered structure tree of a (nested) structure.
fn build_tree(s: &Structure) -> OrderedTree {
    // Recursive build of a nested representation, then flatten to
    // post-order.
    enum Tmp {
        Node(NodeKind, Vec<Tmp>),
    }

    fn build_region(s: &Structure, lo: usize, hi: usize) -> Vec<Tmp> {
        // elements of the loop spanning [lo, hi)
        let mut out = Vec::new();
        let mut k = lo;
        while k < hi {
            match s.partner(k) {
                Some(p) if p > k && p < hi => {
                    let children = build_region(s, k + 1, p);
                    out.push(Tmp::Node(NodeKind::Pair, children));
                    k = p + 1;
                }
                _ => {
                    out.push(Tmp::Node(NodeKind::Unpaired, Vec::new()));
                    k += 1;
                }
            }
        }
        out
    }

    let root = Tmp::Node(NodeKind::Root, build_region(s, 0, s.len()));

    // Flatten to post-order, recording leftmost-leaf descendants.
    let mut kind = Vec::new();
    let mut lmld = Vec::new();

    fn flatten(node: &Tmp, kind: &mut Vec<NodeKind>, lmld: &mut Vec<usize>) -> usize {
        let Tmp::Node(k, children) = node;
        let mut leftmost: Option<usize> = None;
        for c in children {
            let cl = flatten(c, kind, lmld);
            if leftmost.is_none() {
                // leftmost leaf of the first child subtree
                leftmost = Some(lmld[cl]);
            }
        }
        let idx = kind.len();
        kind.push(*k);
        let my_lmld = leftmost.unwrap_or(idx); // a leaf is its own lmld
        lmld.push(my_lmld);
        idx
    }

    flatten(&root, &mut kind, &mut lmld);

    // Key roots: a node is a key root if it has no parent or it is not
    // the leftmost child of its parent. Equivalently: the set of
    // nodes v such that no node w>v has lmld[w]==lmld[v].
    let n = kind.len();
    let mut keyroots = Vec::new();
    let mut seen_lmld = std::collections::HashSet::new();
    for v in (0..n).rev() {
        if seen_lmld.insert(lmld[v]) {
            keyroots.push(v);
        }
    }
    keyroots.sort_unstable();

    OrderedTree {
        kind,
        lmld,
        keyroots,
    }
}

/// The cost of relabelling node label `a` to label `b` (0 if equal,
/// 1 otherwise). Insert / delete each cost 1.
fn relabel_cost(a: NodeKind, b: NodeKind) -> usize {
    if a == b {
        0
    } else {
        1
    }
}

/// The Zhang-Shasha ordered tree-edit distance between two
/// structures.
///
/// Unlike [`base_pair_distance`], this does **not** require the two
/// structures to have the same length — it is a genuine edit distance
/// over the structure trees.
pub fn tree_edit_distance(a: &Structure, b: &Structure) -> usize {
    let ta = build_tree(a);
    let tb = build_tree(b);
    let na = ta.kind.len();
    let nb = tb.kind.len();
    // treedist[i][j] — final answer accumulates here.
    let mut treedist = vec![vec![0usize; nb]; na];

    for &ki in &ta.keyroots {
        for &kj in &tb.keyroots {
            // forest-distance DP over the subtrees rooted at ki, kj.
            let il = ta.lmld[ki];
            let jl = tb.lmld[kj];
            let rows = ki - il + 2;
            let cols = kj - jl + 2;
            let mut fd = vec![vec![0usize; cols]; rows];
            // First row / column of the forest-distance matrix: pure
            // insertions / deletions.
            for dj in 1..cols {
                fd[0][dj] = dj;
            }
            for (di, row) in fd.iter_mut().enumerate().skip(1) {
                row[0] = di;
            }
            for di in 1..rows {
                for dj in 1..cols {
                    let i = il + di - 1;
                    let j = jl + dj - 1;
                    let del = fd[di - 1][dj] + 1;
                    let ins = fd[di][dj - 1] + 1;
                    if ta.lmld[i] == il && tb.lmld[j] == jl {
                        // both are subtrees rooted at i / j
                        let sub = fd[di - 1][dj - 1] + relabel_cost(ta.kind[i], tb.kind[j]);
                        let m = del.min(ins).min(sub);
                        fd[di][dj] = m;
                        treedist[i][j] = m;
                    } else {
                        let pi = ta.lmld[i] - il;
                        let pj = tb.lmld[j] - jl;
                        let sub = fd[pi][pj] + treedist[i][j];
                        fd[di][dj] = del.min(ins).min(sub);
                    }
                }
            }
        }
    }

    treedist[na - 1][nb - 1]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_structures_have_zero_distance() {
        let s = Structure::from_dot_bracket("(((...)))").unwrap();
        assert_eq!(base_pair_distance(&s, &s).unwrap(), 0);
        assert_eq!(tree_edit_distance(&s, &s), 0);
        assert!((pairing_agreement(&s, &s).unwrap() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn base_pair_distance_counts_symmetric_difference() {
        let a = Structure::from_dot_bracket("(((...)))").unwrap();
        let b = Structure::from_dot_bracket("((.....))").unwrap();
        // a has pairs {0-8,1-7,2-6}; b has {0-8,1-7}. Difference = 1.
        assert_eq!(base_pair_distance(&a, &b).unwrap(), 1);
    }

    #[test]
    fn base_pair_distance_fully_disjoint() {
        // Both structures are 8 nt. `a` pairs (0,7) and (1,6); `b` pairs
        // (2,5) — no shared pair. The base-pair distance is the
        // symmetric difference: 2 + 1 = 3.
        let a = Structure::from_dot_bracket("((....))").unwrap();
        let b = Structure::from_dot_bracket("..(..)..").unwrap();
        let d = base_pair_distance(&a, &b).unwrap();
        assert_eq!(d, 3);
    }

    #[test]
    fn base_pair_distance_rejects_length_mismatch() {
        let a = Structure::from_dot_bracket("(((...)))").unwrap();
        let b = Structure::from_dot_bracket("(...)").unwrap();
        assert!(base_pair_distance(&a, &b).is_err());
    }

    #[test]
    fn tree_edit_distance_is_symmetric() {
        let a = Structure::from_dot_bracket("(((...)))").unwrap();
        let b = Structure::from_dot_bracket("((.....))").unwrap();
        assert_eq!(tree_edit_distance(&a, &b), tree_edit_distance(&b, &a));
    }

    #[test]
    fn tree_edit_distance_grows_with_difference() {
        let a = Structure::from_dot_bracket("(((...)))").unwrap();
        let close = Structure::from_dot_bracket("((.....))").unwrap();
        let far = Structure::from_dot_bracket(".........").unwrap();
        let d_close = tree_edit_distance(&a, &close);
        let d_far = tree_edit_distance(&a, &far);
        assert!(d_close > 0);
        assert!(
            d_far >= d_close,
            "unfolded structure should be at least as far ({d_far} vs {d_close})"
        );
    }

    #[test]
    fn tree_edit_distance_handles_different_lengths() {
        let a = Structure::from_dot_bracket("(((...)))").unwrap();
        let b = Structure::from_dot_bracket("(...)").unwrap();
        // does not panic, returns a positive distance
        let d = tree_edit_distance(&a, &b);
        assert!(d > 0);
    }

    #[test]
    fn pairing_agreement_partial() {
        let a = Structure::from_dot_bracket("(((...)))").unwrap();
        let b = Structure::from_dot_bracket("((.....))").unwrap();
        let agr = pairing_agreement(&a, &b).unwrap();
        assert!((0.0..1.0).contains(&agr));
    }
}
