//! Felsenstein's pruning algorithm (Felsenstein 1981).
//!
//! The pruning algorithm computes the likelihood of an alignment given
//! a tree, a [substitution model](super::model) and branch lengths in
//! a single postorder traversal — without it, summing over every
//! possible assignment of states to internal nodes would cost `4^m`.
//!
//! For one column, each node carries a **conditional-likelihood
//! vector** `L[s]` = the probability of the data in that node's
//! subtree, given the node is in state `s`:
//!
//! - **Leaf:** `L` is the indicator of the observed base (a gap is the
//!   all-ones vector — every state equally consistent).
//! - **Internal:** for each child `c` with branch length `t_c` and
//!   transition matrix `P(t_c)`, form
//!   `L_c'[s] = Σ_x P(t_c)[s][x] · L_c[x]`; the node's vector is the
//!   element-wise product of the `L_c'` over its children.
//!
//! The column likelihood is `Σ_s π[s] · L_root[s]`; the tree
//! log-likelihood is the sum of the per-column log-likelihoods. To
//! avoid underflow on long alignments the per-column likelihood is
//! taken in log space immediately.

use crate::error::{PhyloError, Result};
use crate::likelihood::gamma::DiscreteGamma;
use crate::likelihood::model::{SubstModel, TransitionMatrix};
use crate::tree::{NodeId, Tree};

/// A tree paired with a substitution model — the unit a likelihood is
/// computed for. Branch lengths live on the [`Tree`].
#[derive(Debug, Clone)]
pub struct LikelihoodModel {
    /// The (possibly optimised) tree topology and branch lengths.
    pub tree: Tree,
    /// The nucleotide substitution model.
    pub model: SubstModel,
}

/// Maps a nucleotide byte to a 0..3 index, or `None` for a gap /
/// ambiguity / non-ACGT byte (treated as missing data).
fn base_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

/// A leaf id paired with the alignment row that belongs to it.
type LeafRow<'a> = (NodeId, &'a [u8]);

/// Resolves an alignment so row `i` belongs to leaf `leaf_order[i]`.
///
/// Returns, for each leaf id, its sequence as `&[u8]`, plus the shared
/// alignment width. Errors if a leaf has no row or rows differ in
/// width.
fn resolve_rows<'a>(
    tree: &Tree,
    alignment: &'a [(String, Vec<u8>)],
) -> Result<(Vec<LeafRow<'a>>, usize)> {
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
    let mut rows = Vec::new();
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
        rows.push((leaf, row));
    }
    Ok((rows, width))
}

/// Log-likelihood of an alignment given a tree and a substitution
/// model (single rate, no rate heterogeneity).
///
/// `alignment` maps leaf labels to nucleotide byte rows (`A`/`C`/`G`/`T`;
/// any other byte is treated as missing data).
///
/// # Errors
/// - [`PhyloError::Invalid`] on an empty alignment or a leaf with no
///   row.
/// - [`PhyloError::Dimension`] if the rows differ in width.
/// - any error from [`SubstModel::transition_engine`].
pub fn log_likelihood(
    tree: &Tree,
    model: &SubstModel,
    alignment: &[(String, Vec<u8>)],
) -> Result<f64> {
    let engine = model.transition_engine()?;
    let (rows, width) = resolve_rows(tree, alignment)?;
    let mut total = 0.0;
    let mut scratch = Scratch::new(tree.node_count());
    for col in 0..width {
        total += column_log_likelihood(tree, &engine, &rows, col, 1.0, &mut scratch);
    }
    Ok(total)
}

/// Log-likelihood of an alignment under a tree, a model **and**
/// discrete-gamma rate heterogeneity.
///
/// Each site's likelihood is averaged over the gamma rate categories:
/// `L_site = (1/k) Σ_cat L_site(rate_cat)`. Branch lengths are scaled
/// by each category's relative rate.
///
/// # Errors
/// As [`log_likelihood`].
pub fn log_likelihood_gamma(
    tree: &Tree,
    model: &SubstModel,
    alignment: &[(String, Vec<u8>)],
    gamma: &DiscreteGamma,
) -> Result<f64> {
    let engine = model.transition_engine()?;
    let (rows, width) = resolve_rows(tree, alignment)?;
    let cat_prob = gamma.category_probability();
    let mut total = 0.0;
    let mut scratch = Scratch::new(tree.node_count());
    for col in 0..width {
        // Average the per-category column likelihoods, in linear space,
        // then take the log.
        let mut site_like = 0.0;
        for &rate in gamma.rates() {
            let ll = column_log_likelihood(tree, &engine, &rows, col, rate, &mut scratch);
            site_like += cat_prob * ll.exp();
        }
        total += site_like.max(f64::MIN_POSITIVE).ln();
    }
    Ok(total)
}

/// Reusable per-node conditional-likelihood storage.
struct Scratch {
    /// `clv[node]` is the node's 4-vector for the current column.
    clv: Vec<[f64; 4]>,
}

impl Scratch {
    fn new(n: usize) -> Self {
        Scratch {
            clv: vec![[0.0; 4]; n],
        }
    }
}

/// Felsenstein pruning for one column at a given relative `rate`;
/// returns the column's log-likelihood.
fn column_log_likelihood(
    tree: &Tree,
    engine: &TransitionMatrix,
    rows: &[(NodeId, &[u8])],
    col: usize,
    rate: f64,
    scratch: &mut Scratch,
) -> f64 {
    // Leaf vectors: indicator of the observed base.
    for &(leaf, row) in rows {
        let v = &mut scratch.clv[leaf];
        match base_index(row[col]) {
            Some(b) => {
                *v = [0.0; 4];
                v[b] = 1.0;
            }
            None => *v = [1.0; 4], // gap / missing: all states allowed
        }
    }
    // Postorder: internal vectors as products over children.
    for &id in &tree.postorder() {
        let node = tree.node(id);
        if node.is_leaf() {
            continue;
        }
        let mut acc = [1.0_f64; 4];
        for &child in &node.children {
            let bl = node_branch_length(tree, child) * rate;
            let p = engine.p(bl);
            let child_clv = scratch.clv[child];
            // Partial vector for this child: P(t) · L_child.
            for s in 0..4 {
                let mut sum = 0.0;
                for x in 0..4 {
                    sum += p[(s, x)] * child_clv[x];
                }
                acc[s] *= sum;
            }
        }
        scratch.clv[id] = acc;
    }
    // Root: weight by equilibrium frequencies.
    let freqs = engine.frequencies();
    let root_clv = scratch.clv[tree.root()];
    let like: f64 = (0..4).map(|s| freqs[s] * root_clv[s]).sum();
    like.max(f64::MIN_POSITIVE).ln()
}

/// Branch length of the edge above `child`, defaulting to a tiny
/// positive value when unset (a zero-length branch makes `P(t)` the
/// identity and is biologically meaningful, but `None` usually means
/// "unknown" — a small length keeps the likelihood finite).
fn node_branch_length(tree: &Tree, child: NodeId) -> f64 {
    tree.node(child).branch_length.unwrap_or(0.1).max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::read_newick;

    fn row(label: &str, seq: &str) -> (String, Vec<u8>) {
        (label.to_string(), seq.as_bytes().to_vec())
    }

    #[test]
    fn likelihood_is_a_negative_log_probability() {
        let tree = read_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();
        let aln = vec![
            row("A", "ACGT"),
            row("B", "ACGT"),
            row("C", "ACGA"),
            row("D", "ACGA"),
        ];
        let ll = log_likelihood(&tree, &SubstModel::Jc69, &aln).unwrap();
        // A log-likelihood of a probability is finite and <= 0.
        assert!(ll.is_finite() && ll <= 0.0, "ll = {ll}");
    }

    #[test]
    fn identical_sequences_score_better_when_clustered_correctly() {
        // A,B identical; C,D identical. The tree grouping them should
        // out-score the tree that splits the pairs.
        let aln = vec![
            row("A", "AAAAAAAA"),
            row("B", "AAAAAAAA"),
            row("C", "GGGGGGGG"),
            row("D", "GGGGGGGG"),
        ];
        let good = read_newick("((A:0.05,B:0.05):0.2,(C:0.05,D:0.05):0.2);").unwrap();
        let bad = read_newick("((A:0.05,C:0.05):0.2,(B:0.05,D:0.05):0.2);").unwrap();
        let model = SubstModel::Jc69;
        let ll_good = log_likelihood(&good, &model, &aln).unwrap();
        let ll_bad = log_likelihood(&bad, &model, &aln).unwrap();
        assert!(ll_good > ll_bad, "good {ll_good} !> bad {ll_bad}");
    }

    #[test]
    fn longer_alignment_lowers_total_log_likelihood() {
        // Each independent column multiplies the probability => the
        // log-likelihood becomes more negative.
        let tree = read_newick("((A:0.1,B:0.1):0.1,C:0.1);").unwrap();
        let short = vec![row("A", "AC"), row("B", "AC"), row("C", "AC")];
        let long = vec![
            row("A", "ACACACAC"),
            row("B", "ACACACAC"),
            row("C", "ACACACAC"),
        ];
        let model = SubstModel::Jc69;
        let ll_short = log_likelihood(&tree, &model, &short).unwrap();
        let ll_long = log_likelihood(&tree, &model, &long).unwrap();
        assert!(ll_long < ll_short);
    }

    #[test]
    fn gamma_likelihood_is_finite_and_differs_from_uniform() {
        let tree = read_newick("((A:0.2,B:0.2):0.2,(C:0.2,D:0.2):0.2);").unwrap();
        let aln = vec![
            row("A", "ACGTACGT"),
            row("B", "ACGTTCGT"),
            row("C", "AGGTACGA"),
            row("D", "ACGTACGA"),
        ];
        let model = SubstModel::Hky85 {
            kappa: 2.0,
            freqs: [0.3, 0.2, 0.25, 0.25],
        };
        let gamma = DiscreteGamma::new(0.5, 4).unwrap();
        let ll_g = log_likelihood_gamma(&tree, &model, &aln, &gamma).unwrap();
        let ll_u = log_likelihood(&tree, &model, &aln).unwrap();
        assert!(ll_g.is_finite() && ll_g <= 0.0);
        // Strong heterogeneity moves the score away from the uniform
        // model.
        assert!((ll_g - ll_u).abs() > 1e-6);
    }

    #[test]
    fn gaps_are_handled_as_missing_data() {
        let tree = read_newick("((A:0.1,B:0.1):0.1,C:0.1);").unwrap();
        let aln = vec![row("A", "AC-T"), row("B", "ACGT"), row("C", "ACGT")];
        let ll = log_likelihood(&tree, &SubstModel::Jc69, &aln).unwrap();
        assert!(ll.is_finite());
    }

    #[test]
    fn rejects_a_missing_leaf_row() {
        let tree = read_newick("((A,B),C);").unwrap();
        let aln = vec![row("A", "ACGT"), row("B", "ACGT")]; // C missing
        assert!(log_likelihood(&tree, &SubstModel::Jc69, &aln).is_err());
    }
}
