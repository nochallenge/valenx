//! Non-parametric bootstrap support.
//!
//! Bootstrapping estimates how robust each clade of an inferred tree is
//! to sampling noise in the alignment (Felsenstein 1985):
//!
//! 1. Resample the alignment's columns *with replacement* to the
//!    original width — a bootstrap replicate alignment.
//! 2. Infer a tree from that replicate (here: a distance tree, via
//!    [`crate::distance`]).
//! 3. Repeat for `n_replicates` replicates.
//! 4. For every clade of a **reference tree**, report the fraction of
//!    replicate trees that contain the same bipartition — its bootstrap
//!    support.
//!
//! High support (≳ 70-95 %) means the clade is recovered consistently;
//! low support means the data barely distinguish it from alternatives.
//! The reference tree is returned with each internal node labelled by
//! its support percentage.

use crate::distance::cluster::neighbor_joining;
use crate::distance::matrix::{distance_matrix, DistanceModel};
use crate::error::{PhyloError, Result};
use crate::rng::Rng;
use crate::tree::{NodeId, Tree};
use std::collections::{HashMap, HashSet};
use valenx_align::Msa;

/// Outcome of a bootstrap analysis.
#[derive(Debug, Clone)]
pub struct BootstrapResult {
    /// The reference tree, each internal node relabelled with its
    /// integer support percentage (`"0".."100"`).
    pub tree: Tree,
    /// Support fraction (`0.0..=1.0`) for every internal, non-root node
    /// of the reference tree, indexed by node id. Leaf / root entries
    /// are `0.0`.
    pub support: Vec<f64>,
    /// Number of bootstrap replicates performed.
    pub replicates: usize,
}

/// Runs a non-parametric bootstrap and maps support onto a reference
/// tree.
///
/// `msa` is the original alignment; `labels` names its rows (row
/// order); `reference` is the tree whose clades are scored; `model`
/// selects the distance correction used to rebuild each replicate;
/// `n_replicates` is the replicate count; `seed` seeds the
/// deterministic column resampler.
///
/// # Errors
/// - [`PhyloError::Invalid`] on an empty alignment / zero replicates.
/// - [`PhyloError::Dimension`] if `labels.len() != msa.depth()`.
/// - any error propagated from the distance pipeline.
pub fn bootstrap_support(
    msa: &Msa,
    labels: &[String],
    reference: &Tree,
    model: DistanceModel,
    n_replicates: usize,
    seed: u64,
) -> Result<BootstrapResult> {
    if n_replicates == 0 {
        return Err(PhyloError::invalid("n_replicates", "must be positive"));
    }
    let depth = msa.depth();
    if depth < 3 {
        return Err(PhyloError::invalid(
            "msa",
            "bootstrap needs at least three sequences",
        ));
    }
    if labels.len() != depth {
        return Err(PhyloError::dimension(depth, labels.len(), "alignment rows"));
    }
    let width = msa.width();
    if width == 0 {
        return Err(PhyloError::invalid("msa", "zero-width alignment"));
    }

    // Shared leaf index for canonical bipartitions.
    let mut sorted_labels = reference.leaf_labels();
    sorted_labels.sort();
    let n = sorted_labels.len();
    let index: HashMap<String, usize> = sorted_labels
        .iter()
        .cloned()
        .enumerate()
        .map(|(i, l)| (l, i))
        .collect();

    // Reference bipartitions, by node id.
    let ref_splits: Vec<Option<Vec<usize>>> = (0..reference.node_count())
        .map(|id| node_bipartition(reference, id, &index, n))
        .collect();

    // Count, per replicate, how often each reference split recurs.
    let mut hits = vec![0usize; reference.node_count()];
    let mut rng = Rng::new(seed);

    for _ in 0..n_replicates {
        // Resample columns with replacement.
        let replicate = resample_columns(msa, &mut rng);
        // Rebuild a tree (distance + NJ).
        let dm = distance_matrix(&replicate, labels, model)?;
        let rep_tree = neighbor_joining(&dm)?;
        // Collect the replicate's bipartition set.
        let rep_set: HashSet<Vec<usize>> = (0..rep_tree.node_count())
            .filter_map(|id| node_bipartition(&rep_tree, id, &index, n))
            .collect();
        // Tally agreement with each reference clade.
        for (id, split) in ref_splits.iter().enumerate() {
            if let Some(s) = split {
                if rep_set.contains(s) {
                    hits[id] += 1;
                }
            }
        }
    }

    // Support fractions + a relabelled copy of the reference tree.
    let support: Vec<f64> = hits
        .iter()
        .map(|&h| h as f64 / n_replicates as f64)
        .collect();
    let mut tree = reference.clone();
    for id in 0..tree.node_count() {
        if ref_splits[id].is_some() {
            let pct = (support[id] * 100.0).round() as i64;
            tree.node_mut(id).label = Some(pct.to_string());
        }
    }

    Ok(BootstrapResult {
        tree,
        support,
        replicates: n_replicates,
    })
}

/// Resamples an MSA's columns with replacement to the same width.
fn resample_columns(msa: &Msa, rng: &mut Rng) -> Msa {
    let width = msa.width();
    let depth = msa.depth();
    let picks: Vec<usize> = (0..width).map(|_| rng.below(width)).collect();
    let mut rows = Vec::with_capacity(depth);
    for r in 0..depth {
        let src = &msa.rows[r];
        let row: Vec<u8> = picks.iter().map(|&c| src[c]).collect();
        rows.push(row);
    }
    // Equal lengths by construction — `expect` documents that.
    Msa::new(rows).expect("resampled rows are equal length")
}

/// The canonical bipartition (leaf-index vector, side without index 0)
/// induced by an internal, non-root node — or `None` for a
/// leaf / root / trivial split.
fn node_bipartition(
    tree: &Tree,
    id: NodeId,
    index: &HashMap<String, usize>,
    n: usize,
) -> Option<Vec<usize>> {
    let node = tree.node(id);
    if node.is_leaf() || node.parent.is_none() {
        return None;
    }
    let mut side: Vec<usize> = Vec::new();
    for leaf in tree.descendant_leaves(id) {
        let label = tree.node(leaf).label.as_deref()?;
        side.push(*index.get(label)?);
    }
    side.sort_unstable();
    side.dedup();
    if side.len() < 2 || side.len() > n - 2 {
        return None;
    }
    if side.contains(&0) {
        let set: HashSet<usize> = side.iter().copied().collect();
        Some((0..n).filter(|i| !set.contains(i)).collect())
    } else {
        Some(side)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::read_newick;

    /// A clean four-taxon alignment where (A,B) and (C,D) are strongly
    /// supported (many informative columns).
    fn clean_msa() -> (Msa, Vec<String>) {
        let rows = vec![
            b"AAAAAAAAAACCCCC".to_vec(), // A
            b"AAAAAAAAAACCCCC".to_vec(), // B
            b"GGGGGGGGGGTTTTT".to_vec(), // C
            b"GGGGGGGGGGTTTTT".to_vec(), // D
        ];
        let labels: Vec<String> = ["A", "B", "C", "D"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        (Msa::new(rows).unwrap(), labels)
    }

    #[test]
    fn strong_signal_gives_high_support() {
        let (msa, labels) = clean_msa();
        let reference = read_newick("((A,B),(C,D));").unwrap();
        let result = bootstrap_support(
            &msa,
            &labels,
            &reference,
            DistanceModel::JukesCantor,
            100,
            42,
        )
        .unwrap();
        assert_eq!(result.replicates, 100);
        // The (A,B) clade should be recovered in (nearly) every
        // replicate.
        let ab = reference.lca(
            reference.find("A").unwrap(),
            reference.find("B").unwrap(),
        );
        assert!(result.support[ab] > 0.9, "support = {}", result.support[ab]);
    }

    #[test]
    fn support_values_are_fractions() {
        let (msa, labels) = clean_msa();
        let reference = read_newick("((A,B),(C,D));").unwrap();
        let result = bootstrap_support(
            &msa,
            &labels,
            &reference,
            DistanceModel::PDistance,
            50,
            7,
        )
        .unwrap();
        for &s in &result.support {
            assert!((0.0..=1.0).contains(&s));
        }
    }

    #[test]
    fn internal_nodes_get_a_support_label() {
        let (msa, labels) = clean_msa();
        let reference = read_newick("((A,B),(C,D));").unwrap();
        let result = bootstrap_support(
            &msa,
            &labels,
            &reference,
            DistanceModel::JukesCantor,
            30,
            1,
        )
        .unwrap();
        // The relabelled tree has a numeric label on the (A,B) node.
        let ab = result.tree.lca(
            result.tree.find("A").unwrap(),
            result.tree.find("B").unwrap(),
        );
        let label = result.tree.node(ab).label.as_deref().unwrap_or("");
        assert!(label.parse::<i64>().is_ok(), "label `{label}` not numeric");
    }

    #[test]
    fn is_deterministic_for_a_seed() {
        let (msa, labels) = clean_msa();
        let reference = read_newick("((A,B),(C,D));").unwrap();
        let r1 = bootstrap_support(
            &msa,
            &labels,
            &reference,
            DistanceModel::JukesCantor,
            40,
            99,
        )
        .unwrap();
        let r2 = bootstrap_support(
            &msa,
            &labels,
            &reference,
            DistanceModel::JukesCantor,
            40,
            99,
        )
        .unwrap();
        assert_eq!(r1.support, r2.support);
    }

    #[test]
    fn rejects_bad_input() {
        let (msa, labels) = clean_msa();
        let reference = read_newick("((A,B),(C,D));").unwrap();
        // Zero replicates.
        assert!(bootstrap_support(
            &msa,
            &labels,
            &reference,
            DistanceModel::PDistance,
            0,
            1
        )
        .is_err());
    }
}
