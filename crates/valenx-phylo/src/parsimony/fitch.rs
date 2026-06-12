//! Fitch small-parsimony (Fitch 1971).
//!
//! Given a *fixed* tree and an aligned character matrix, Fitch's
//! algorithm computes — in two arena traversals — the minimum number
//! of state changes the tree implies, plus a set of optimal states for
//! every internal node.
//!
//! 1. **Postorder (bottom-up).** Each leaf's state set is the singleton
//!    of its observed state. An internal node's set is the
//!    *intersection* of its children's sets if that intersection is
//!    non-empty; otherwise it is their *union*, and that union counts
//!    one change. The total change count is the parsimony score.
//! 2. **Preorder (top-down).** A final pass narrows each internal set
//!    to a single ancestral state consistent with its parent — the
//!    standard Fitch ancestral assignment.
//!
//! State sets are stored as `u32` bitmasks, so up to 32 states (ample
//! for the 4 nucleotides or 20 amino acids) and the intersection /
//! union are single bit-ops. A gap / missing character is the
//! all-states wildcard mask.

use crate::error::{PhyloError, Result};
use crate::tree::{NodeId, Tree};

/// Result of a Fitch small-parsimony analysis over an alignment.
#[derive(Debug, Clone, PartialEq)]
pub struct FitchResult {
    /// Total parsimony score summed over every alignment column.
    pub score: usize,
    /// Per-column parsimony score (length = alignment width).
    pub site_scores: Vec<usize>,
    /// Reconstructed ancestral states, indexed `[node][column]`. Leaf
    /// rows echo the observed states. A state of `u8::MAX` marks a
    /// node/column the algorithm left fully ambiguous.
    pub ancestral: Vec<Vec<u8>>,
}

/// The mask meaning "any state" — used for gaps and missing data.
fn wildcard_mask(n_states: u8) -> u32 {
    if n_states >= 32 {
        u32::MAX
    } else {
        (1u32 << n_states) - 1
    }
}

/// Runs Fitch small-parsimony on a tree and an aligned character
/// matrix.
///
/// `alignment` maps a *leaf label* to its row of `u8` state indices
/// (one entry per column). A state byte of `u8::MAX` is treated as a
/// gap / missing wildcard. Every leaf of `tree` must have an entry; all
/// rows must share one width.
///
/// `n_states` is the alphabet size (4 for nucleotides, 20 for amino
/// acids); it must be `1..=32`.
///
/// # Errors
/// - [`PhyloError::Invalid`] if `n_states` is out of range, the
///   alignment is empty, or a leaf has no row.
/// - [`PhyloError::Dimension`] if the rows differ in width.
pub fn fitch_parsimony(
    tree: &Tree,
    alignment: &[(String, Vec<u8>)],
    n_states: u8,
) -> Result<FitchResult> {
    if !(1..=32).contains(&n_states) {
        return Err(PhyloError::invalid("n_states", "must be in 1..=32"));
    }
    if alignment.is_empty() {
        return Err(PhyloError::invalid("alignment", "no sequences supplied"));
    }
    let width = alignment[0].1.len();
    if width == 0 {
        return Err(PhyloError::invalid("alignment", "zero-width alignment"));
    }
    for (name, row) in alignment {
        if row.len() != width {
            return Err(PhyloError::dimension(width, row.len(), "alignment rows"));
        }
        let _ = name;
    }

    let n = tree.node_count();
    let wild = wildcard_mask(n_states);

    // Resolve each leaf to its alignment row.
    let row_for = |id: NodeId| -> Result<&Vec<u8>> {
        let label = tree
            .node(id)
            .label
            .as_deref()
            .ok_or_else(|| PhyloError::invalid("tree", "leaf without a label"))?;
        alignment
            .iter()
            .find(|(name, _)| name == label)
            .map(|(_, row)| row)
            .ok_or_else(|| PhyloError::invalid("alignment", format!("no row for leaf `{label}`")))
    };

    let post = tree.postorder();
    let pre = tree.preorder();

    let mut site_scores = vec![0usize; width];
    // Final single-state ancestral assignment.
    let mut ancestral = vec![vec![u8::MAX; width]; n];

    // `col` indexes several per-column arrays (site_scores, ancestral
    // rows, the alignment rows) — a range loop is the clearest form.
    #[allow(clippy::needless_range_loop)]
    for col in 0..width {
        // Bottom-up: state-set masks.
        let mut mask = vec![0u32; n];
        for &id in &post {
            let node = tree.node(id);
            if node.is_leaf() {
                let s = row_for(id)?[col];
                mask[id] = if s == u8::MAX || s >= n_states {
                    wild
                } else {
                    1u32 << s
                };
            } else {
                let mut inter = wild;
                let mut union = 0u32;
                for &c in &node.children {
                    inter &= mask[c];
                    union |= mask[c];
                }
                if inter != 0 {
                    mask[id] = inter;
                } else {
                    mask[id] = union;
                    site_scores[col] += 1;
                }
            }
        }
        // Top-down: pick one state per internal node, preferring one
        // shared with the parent.
        for &id in &pre {
            let node = tree.node(id);
            let chosen = if let Some(p) = node.parent {
                let parent_bit = 1u32 << ancestral[p][col].min(31);
                if ancestral[p][col] != u8::MAX && (mask[id] & parent_bit) != 0 {
                    ancestral[p][col]
                } else {
                    lowest_set_bit(mask[id])
                }
            } else {
                lowest_set_bit(mask[id])
            };
            ancestral[id][col] = chosen;
        }
    }

    let score = site_scores.iter().sum();
    Ok(FitchResult {
        score,
        site_scores,
        ancestral,
    })
}

/// Index of the lowest set bit of `mask`, or `u8::MAX` if `mask == 0`.
fn lowest_set_bit(mask: u32) -> u8 {
    if mask == 0 {
        u8::MAX
    } else {
        mask.trailing_zeros() as u8
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::read_newick;

    /// Alignment helper: `(label, states)` from a string of digits.
    fn col(label: &str, states: &[u8]) -> (String, Vec<u8>) {
        (label.to_string(), states.to_vec())
    }

    #[test]
    fn invariant_site_costs_nothing() {
        let tree = read_newick("((A,B),(C,D));").unwrap();
        // Every taxon has state 0 — no changes needed.
        let aln = vec![
            col("A", &[0]),
            col("B", &[0]),
            col("C", &[0]),
            col("D", &[0]),
        ];
        let r = fitch_parsimony(&tree, &aln, 4).unwrap();
        assert_eq!(r.score, 0);
    }

    #[test]
    fn one_clade_differs_costs_one() {
        let tree = read_newick("((A,B),(C,D));").unwrap();
        // (A,B) = state 0, (C,D) = state 1 => one change on the
        // internal edge.
        let aln = vec![
            col("A", &[0]),
            col("B", &[0]),
            col("C", &[1]),
            col("D", &[1]),
        ];
        let r = fitch_parsimony(&tree, &aln, 4).unwrap();
        assert_eq!(r.score, 1);
        assert_eq!(r.site_scores, vec![1]);
    }

    #[test]
    fn fully_homoplasious_site() {
        let tree = read_newick("((A,B),(C,D));").unwrap();
        // A=C=0, B=D=1 — incongruent with the tree => two changes.
        let aln = vec![
            col("A", &[0]),
            col("B", &[1]),
            col("C", &[0]),
            col("D", &[1]),
        ];
        let r = fitch_parsimony(&tree, &aln, 4).unwrap();
        assert_eq!(r.score, 2);
    }

    #[test]
    fn multi_column_score_is_the_sum() {
        let tree = read_newick("((A,B),(C,D));").unwrap();
        let aln = vec![
            col("A", &[0, 0]),
            col("B", &[0, 1]),
            col("C", &[1, 0]),
            col("D", &[1, 1]),
        ];
        let r = fitch_parsimony(&tree, &aln, 4).unwrap();
        // Column 0: one change. Column 1: two changes.
        assert_eq!(r.site_scores, vec![1, 2]);
        assert_eq!(r.score, 3);
    }

    #[test]
    fn ancestral_states_are_assigned() {
        let tree = read_newick("((A,B),(C,D));").unwrap();
        let aln = vec![
            col("A", &[0]),
            col("B", &[0]),
            col("C", &[1]),
            col("D", &[1]),
        ];
        let r = fitch_parsimony(&tree, &aln, 4).unwrap();
        // Every internal node gets a concrete (non-MAX) state.
        for id in 0..tree.node_count() {
            if tree.node(id).is_internal() {
                assert_ne!(r.ancestral[id][0], u8::MAX);
            }
        }
    }

    #[test]
    fn gaps_are_wildcards() {
        let tree = read_newick("((A,B),(C,D));").unwrap();
        // D is a gap — it should never force an extra change.
        let aln = vec![
            col("A", &[0]),
            col("B", &[0]),
            col("C", &[0]),
            col("D", &[u8::MAX]),
        ];
        let r = fitch_parsimony(&tree, &aln, 4).unwrap();
        assert_eq!(r.score, 0);
    }

    #[test]
    fn rejects_bad_input() {
        let tree = read_newick("(A,B);").unwrap();
        assert!(fitch_parsimony(&tree, &[], 4).is_err());
        let ragged = vec![col("A", &[0, 0]), col("B", &[0])];
        assert!(fitch_parsimony(&tree, &ragged, 4).is_err());
    }
}
