//! Posterior summaries — MAP tree, consensus tree, clade posterior
//! probabilities.
//!
//! Given a tree sample from an MCMC chain ([`super::chain::ChainResult`])
//! the typical reporting is:
//!
//! - the **majority-rule consensus** of the trees, with each retained
//!   internal clade labelled by its posterior probability (= the
//!   fraction of trees in the sample that contained that clade),
//! - the **MAP tree** — the (topology, branch lengths) sample with the
//!   single highest log posterior,
//! - the per-clade **posterior probability** table — every clade ever
//!   sampled, with its frequency.
//!
//! The consensus + clade-probability table here are *the* BEAST 2 /
//! MrBayes summary, so this module is what the wider tool wires into.

use crate::compare::consensus::{consensus_tree, ConsensusKind};
use crate::error::{PhyloError, Result};
use crate::tree::Tree;
use std::collections::HashMap;

/// A clade posterior — a leaf-label set + the fraction of input trees
/// that contained it.
#[derive(Debug, Clone)]
pub struct CladePosterior {
    /// Sorted list of leaf labels that make up the clade.
    pub clade: Vec<String>,
    /// Posterior probability of this clade (between 0 and 1).
    pub probability: f64,
}

/// Posterior summary of a tree sample.
#[derive(Debug, Clone)]
pub struct PosteriorSummary {
    /// Majority-rule consensus tree (clades present in `>50 %` of
    /// trees). Each retained internal node is labelled with its
    /// support frequency (e.g. `"0.85"`).
    pub consensus: Tree,
    /// MAP tree (the highest-posterior sample) — branch lengths
    /// preserved.
    pub map_tree: Tree,
    /// All clades ever sampled, with their posterior probabilities,
    /// sorted by probability (descending).
    pub clade_probabilities: Vec<CladePosterior>,
}

/// Computes a posterior summary from a tree sample plus per-sample log
/// posteriors. The MAP tree is the one with the highest log posterior.
///
/// # Errors
/// [`PhyloError::Invalid`] if `trees` is empty, if the trees / log
/// posteriors disagree in length, or if the trees disagree on their
/// leaf set.
pub fn summarise_posterior(
    trees: &[Tree],
    log_posteriors: &[f64],
) -> Result<PosteriorSummary> {
    if trees.is_empty() {
        return Err(PhyloError::invalid("trees", "empty tree sample"));
    }
    if trees.len() != log_posteriors.len() {
        return Err(PhyloError::dimension(
            trees.len(),
            log_posteriors.len(),
            "tree / log-posterior counts",
        ));
    }
    // MAP: argmax of the log-posterior trace.
    let mut best = 0usize;
    let mut best_lp = log_posteriors[0];
    for (i, &lp) in log_posteriors.iter().enumerate().skip(1) {
        if lp > best_lp {
            best_lp = lp;
            best = i;
        }
    }
    let map_tree = trees[best].clone();
    // Consensus tree on the same sample.
    let consensus = consensus_tree(trees, ConsensusKind::MajorityRule)?;
    let clade_probabilities = clade_posterior_table(trees);
    Ok(PosteriorSummary {
        consensus,
        map_tree,
        clade_probabilities,
    })
}

/// Builds the clade-posterior table: every non-trivial clade that ever
/// appeared in the sample, paired with the fraction of trees that had
/// it.
pub fn clade_posterior_table(trees: &[Tree]) -> Vec<CladePosterior> {
    if trees.is_empty() {
        return Vec::new();
    }
    let mut counts: HashMap<Vec<String>, usize> = HashMap::new();
    for tree in trees {
        // For each tree, collect the set of non-trivial clades it
        // contains, deduplicate (a clade may correspond to multiple
        // internal node ids on a malformed tree — keep one count per
        // distinct leaf set per tree), and tally.
        let mut seen: std::collections::HashSet<Vec<String>> =
            std::collections::HashSet::new();
        for id in 0..tree.node_count() {
            let node = tree.node(id);
            if node.is_leaf() || node.parent.is_none() {
                continue;
            }
            let mut clade: Vec<String> = tree
                .descendant_leaves(id)
                .into_iter()
                .filter_map(|l| tree.node(l).label.clone())
                .collect();
            clade.sort();
            clade.dedup();
            if clade.len() < 2 || clade.len() >= tree.leaf_count() {
                continue;
            }
            seen.insert(clade);
        }
        for clade in seen {
            *counts.entry(clade).or_insert(0) += 1;
        }
    }
    let total = trees.len() as f64;
    let mut table: Vec<CladePosterior> = counts
        .into_iter()
        .map(|(clade, c)| CladePosterior {
            clade,
            probability: c as f64 / total,
        })
        .collect();
    table.sort_by(|a, b| {
        b.probability
            .partial_cmp(&a.probability)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.clade.cmp(&b.clade))
    });
    table
}

/// Posterior probability of a specific named clade (looked up by its
/// sorted leaf-label vector). Returns 0.0 if never sampled.
pub fn clade_probability(trees: &[Tree], clade: &[&str]) -> f64 {
    let mut wanted: Vec<String> = clade.iter().map(|s| s.to_string()).collect();
    wanted.sort();
    let table = clade_posterior_table(trees);
    table
        .into_iter()
        .find(|cp| cp.clade == wanted)
        .map(|cp| cp.probability)
        .unwrap_or(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::read_newick;

    #[test]
    fn map_tree_is_argmax_log_posterior() {
        let trees = vec![
            read_newick("((A,B),(C,D));").unwrap(),
            read_newick("((A,C),(B,D));").unwrap(),
            read_newick("((A,D),(B,C));").unwrap(),
        ];
        let lps = vec![-200.0, -150.0, -300.0];
        let summary = summarise_posterior(&trees, &lps).unwrap();
        // Index 1 has the highest log posterior.
        assert_eq!(summary.map_tree.leaf_count(), 4);
        // The kept clade for the MAP should be {A, C}.
        let cl: Vec<Vec<String>> = (0..summary.map_tree.node_count())
            .filter(|&id| summary.map_tree.node(id).is_internal())
            .map(|id| {
                let mut v: Vec<String> = summary
                    .map_tree
                    .descendant_leaves(id)
                    .into_iter()
                    .filter_map(|l| summary.map_tree.node(l).label.clone())
                    .collect();
                v.sort();
                v
            })
            .collect();
        assert!(cl.iter().any(|c| c == &["A", "C"]), "map clades: {cl:?}");
    }

    #[test]
    fn clade_probabilities_count_correctly() {
        let trees = vec![
            read_newick("((A,B),(C,D));").unwrap(),
            read_newick("((A,B),(C,D));").unwrap(),
            read_newick("((A,C),(B,D));").unwrap(),
        ];
        let p = clade_probability(&trees, &["A", "B"]);
        assert!((p - 2.0 / 3.0).abs() < 1e-9, "p = {p}");
        let q = clade_probability(&trees, &["A", "C"]);
        assert!((q - 1.0 / 3.0).abs() < 1e-9, "q = {q}");
    }

    #[test]
    fn consensus_has_a_kept_clade_label() {
        let trees = vec![
            read_newick("((A,B),(C,D));").unwrap(),
            read_newick("((A,B),(C,D));").unwrap(),
            read_newick("((A,B),(C,D));").unwrap(),
        ];
        let lps = vec![-200.0; 3];
        let summary = summarise_posterior(&trees, &lps).unwrap();
        // The (A,B) consensus node carries the 1.00 label.
        let label = summary
            .consensus
            .nodes()
            .iter()
            .filter_map(|n| n.label.clone())
            .find(|s| s == "1.00");
        assert!(label.is_some(), "no 1.00 support label");
    }

    #[test]
    fn empty_sample_is_rejected() {
        let trees: Vec<Tree> = Vec::new();
        let lps: Vec<f64> = Vec::new();
        assert!(summarise_posterior(&trees, &lps).is_err());
    }

    #[test]
    fn mismatched_lengths_rejected() {
        let trees = vec![read_newick("(A,B);").unwrap()];
        let lps = vec![-1.0, -2.0];
        assert!(summarise_posterior(&trees, &lps).is_err());
    }

    #[test]
    fn clade_table_is_sorted_descending() {
        let trees = vec![
            read_newick("((A,B),(C,D));").unwrap(),
            read_newick("((A,B),(C,D));").unwrap(),
            read_newick("((A,B),(C,D));").unwrap(),
            read_newick("((A,C),(B,D));").unwrap(),
        ];
        let table = clade_posterior_table(&trees);
        assert!(!table.is_empty());
        for w in table.windows(2) {
            assert!(w[0].probability >= w[1].probability);
        }
    }
}
