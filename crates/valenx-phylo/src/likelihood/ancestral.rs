//! Marginal ancestral-state reconstruction.
//!
//! Given a tree, a substitution model and an alignment, this estimates
//! the probability distribution of nucleotide states at every internal
//! node — the *marginal* posterior, computed independently per node
//! (as opposed to the *joint* most-likely assignment over all nodes at
//! once).
//!
//! The computation combines a downward and an upward likelihood pass:
//!
//! - The **partial (downward) likelihood** `L_v[s]` of a node — its
//!   subtree's data given state `s` — is exactly the
//!   [Felsenstein](super::felsenstein) conditional-likelihood vector.
//! - The **partial (upward) likelihood** `U_v[s]` is the probability of
//!   all data *outside* `v`'s subtree given `v` is in state `s`.
//!
//! The marginal posterior at node `v` is then proportional to
//! `π[s] · L_v[s] · U_v[s]` (for the root, `U` is the equilibrium
//! vector). Normalising over `s` gives the per-node state distribution.

use crate::error::{PhyloError, Result};
use crate::likelihood::model::SubstModel;
use crate::tree::{NodeId, Tree};
use nalgebra::Matrix4;

/// Per-node marginal ancestral-state estimate for one alignment.
#[derive(Debug, Clone, PartialEq)]
pub struct AncestralResult {
    /// Posterior state distribution, indexed `[node][column][state]`.
    /// Each `[state]` 4-vector sums to 1. Leaf rows reflect the
    /// observed (or, for gaps, equilibrium-weighted) state.
    pub posteriors: Vec<Vec<[f64; 4]>>,
    /// Most-probable state per node and column (the marginal MAP
    /// estimate), `[node][column]`.
    pub map_states: Vec<Vec<u8>>,
}

impl AncestralResult {
    /// Posterior distribution at one node / column.
    ///
    /// # Panics
    /// If either index is out of range.
    pub fn posterior(&self, node: NodeId, column: usize) -> [f64; 4] {
        self.posteriors[node][column]
    }
}

/// Maps a nucleotide byte to a 0..3 index, or `None` for a gap.
fn base_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

/// Reconstructs marginal ancestral states for every internal node.
///
/// `alignment` maps leaf labels to nucleotide rows (`A`/`C`/`G`/`T`;
/// other bytes are missing data).
///
/// # Errors
/// - [`PhyloError::Invalid`] on an empty alignment / a leaf with no
///   row.
/// - [`PhyloError::Dimension`] if the rows differ in width.
/// - any error from [`SubstModel::transition_engine`].
pub fn ancestral_states(
    tree: &Tree,
    model: &SubstModel,
    alignment: &[(String, Vec<u8>)],
) -> Result<AncestralResult> {
    if alignment.is_empty() {
        return Err(PhyloError::invalid("alignment", "no sequences supplied"));
    }
    let width = alignment[0].1.len();
    if width == 0 {
        return Err(PhyloError::invalid("alignment", "zero-width alignment"));
    }
    for (_, row) in alignment {
        if row.len() != width {
            return Err(PhyloError::dimension(width, row.len(), "alignment rows"));
        }
    }
    let engine = model.transition_engine()?;
    let freqs = engine.frequencies();
    let n = tree.node_count();

    // Resolve leaf rows.
    let mut leaf_rows: Vec<(NodeId, &[u8])> = Vec::new();
    for &leaf in &tree.leaves() {
        let label = tree
            .node(leaf)
            .label
            .as_deref()
            .ok_or_else(|| PhyloError::invalid("tree", "leaf without a label"))?;
        let row = alignment
            .iter()
            .find(|(name, _)| name == label)
            .map(|(_, r)| r.as_slice())
            .ok_or_else(|| {
                PhyloError::invalid("alignment", format!("no row for leaf `{label}`"))
            })?;
        leaf_rows.push((leaf, row));
    }

    let post = tree.postorder();
    let pre = tree.preorder();

    // Precompute each non-root node's P(t) once — branch length is
    // column-independent.
    let p_matrices: Vec<Option<Matrix4<f64>>> = (0..n)
        .map(|id| {
            tree.node(id)
                .parent
                .map(|_| engine.p(tree.node(id).branch_length.unwrap_or(0.1).max(0.0)))
        })
        .collect();

    let mut posteriors = vec![vec![[0.25_f64; 4]; width]; n];
    let mut map_states = vec![vec![0u8; width]; n];

    for col in 0..width {
        // --- Downward pass: Felsenstein partial likelihoods L.
        let mut down = vec![[1.0_f64; 4]; n];
        for &(leaf, row) in &leaf_rows {
            match base_index(row[col]) {
                Some(b) => {
                    down[leaf] = [0.0; 4];
                    down[leaf][b] = 1.0;
                }
                None => down[leaf] = [1.0; 4],
            }
        }
        for &id in &post {
            let node = tree.node(id);
            if node.is_leaf() {
                continue;
            }
            let mut acc = [1.0_f64; 4];
            for &child in &node.children {
                let p = p_matrices[child].expect("child has a branch");
                let cl = down[child];
                for s in 0..4 {
                    let mut sum = 0.0;
                    for x in 0..4 {
                        sum += p[(s, x)] * cl[x];
                    }
                    acc[s] *= sum;
                }
            }
            down[id] = acc;
        }

        // --- Upward pass: partial likelihoods U for the rest of the
        // tree. Root's U is the equilibrium distribution.
        let mut up = vec![[1.0_f64; 4]; n];
        up[tree.root()] = freqs;
        for &id in &pre {
            let node = tree.node(id);
            if node.parent.is_none() {
                continue;
            }
            let parent = node.parent.unwrap();
            // Sibling product: parent's U times every sibling's
            // downward contribution.
            let mut sib_product = up[parent];
            for &sib in &tree.node(parent).children {
                if sib == id {
                    continue;
                }
                let p = p_matrices[sib].expect("sibling has a branch");
                let cl = down[sib];
                for s in 0..4 {
                    let mut sum = 0.0;
                    for x in 0..4 {
                        sum += p[(s, x)] * cl[x];
                    }
                    sib_product[s] *= sum;
                }
            }
            // Propagate across this node's own branch: U_v[s] =
            // Σ_x P(t)[x][s] · sib_product[x].
            let p = p_matrices[id].expect("node has a branch");
            let mut u = [0.0_f64; 4];
            for s in 0..4 {
                let mut sum = 0.0;
                for x in 0..4 {
                    sum += p[(x, s)] * sib_product[x];
                }
                u[s] = sum;
            }
            up[id] = u;
        }

        // --- Marginal posterior: π[s] · L[s] · U[s], normalised.
        for id in 0..n {
            let mut post_vec = [0.0_f64; 4];
            for s in 0..4 {
                post_vec[s] = freqs[s] * down[id][s] * up[id][s];
            }
            let total: f64 = post_vec.iter().sum();
            if total > 0.0 {
                for v in &mut post_vec {
                    *v /= total;
                }
            } else {
                post_vec = [0.25; 4];
            }
            posteriors[id][col] = post_vec;
            // MAP state = arg max.
            let (best, _) = post_vec.iter().enumerate().fold(
                (0usize, f64::NEG_INFINITY),
                |(bi, bv), (i, &v)| if v > bv { (i, v) } else { (bi, bv) },
            );
            map_states[id][col] = best as u8;
        }
    }

    Ok(AncestralResult {
        posteriors,
        map_states,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::read_newick;

    fn row(label: &str, seq: &str) -> (String, Vec<u8>) {
        (label.to_string(), seq.as_bytes().to_vec())
    }

    #[test]
    fn posteriors_are_normalised() {
        let tree = read_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();
        let aln = vec![
            row("A", "ACGT"),
            row("B", "ACGT"),
            row("C", "AGGT"),
            row("D", "AGGT"),
        ];
        let r = ancestral_states(&tree, &SubstModel::Jc69, &aln).unwrap();
        for node in 0..tree.node_count() {
            for col in 0..4 {
                let p = r.posterior(node, col);
                let sum: f64 = p.iter().sum();
                assert!((sum - 1.0).abs() < 1e-9, "node {node} col {col}");
                assert!(p.iter().all(|&x| x >= -1e-12));
            }
        }
    }

    #[test]
    fn invariant_column_reconstructs_the_shared_state() {
        // Every taxon has A at column 0 => the ancestor is almost
        // certainly A.
        let tree = read_newick("((A:0.05,B:0.05):0.05,(C:0.05,D:0.05):0.05);").unwrap();
        let aln = vec![row("A", "A"), row("B", "A"), row("C", "A"), row("D", "A")];
        let r = ancestral_states(&tree, &SubstModel::Jc69, &aln).unwrap();
        let root = tree.root();
        // Index 0 == 'A'.
        assert_eq!(r.map_states[root][0], 0);
        assert!(r.posterior(root, 0)[0] > 0.9, "{:?}", r.posterior(root, 0));
    }

    #[test]
    fn clade_specific_state_propagates_to_its_ancestor() {
        // (A,B) = G, (C,D) = A. The (A,B) ancestor should favour G.
        let tree = read_newick("((A:0.05,B:0.05):0.1,(C:0.05,D:0.05):0.1);").unwrap();
        let aln = vec![row("A", "G"), row("B", "G"), row("C", "A"), row("D", "A")];
        let r = ancestral_states(&tree, &SubstModel::Jc69, &aln).unwrap();
        let ab = tree.lca(tree.find("A").unwrap(), tree.find("B").unwrap());
        // Index 2 == 'G'.
        assert_eq!(r.map_states[ab][0], 2, "{:?}", r.posterior(ab, 0));
    }

    #[test]
    fn leaf_posteriors_match_observations() {
        let tree = read_newick("((A:0.1,B:0.1):0.1,C:0.1);").unwrap();
        let aln = vec![row("A", "C"), row("B", "C"), row("C", "C")];
        let r = ancestral_states(&tree, &SubstModel::Jc69, &aln).unwrap();
        let a = tree.find("A").unwrap();
        // Leaf A observed C (index 1) => posterior concentrates there.
        assert!(r.posterior(a, 0)[1] > 0.99);
    }

    #[test]
    fn rejects_bad_input() {
        let tree = read_newick("((A,B),C);").unwrap();
        assert!(ancestral_states(&tree, &SubstModel::Jc69, &[]).is_err());
    }
}
