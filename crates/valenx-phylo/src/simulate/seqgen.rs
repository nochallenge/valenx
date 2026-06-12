//! Sequence evolution along a tree (Seq-Gen-class).
//!
//! Given a tree, a [substitution model](crate::likelihood::model) and a
//! sequence length, this simulates an alignment: a root sequence is
//! drawn from the model's equilibrium frequencies, then evolved down
//! every branch.
//!
//! Down a branch of length `t`, each site changes state according to
//! the transition-probability matrix `P(t)`: the descendant base is
//! sampled from the row of `P(t)` indexed by the ancestral base. A
//! preorder traversal therefore evolves the whole tree, and the leaf
//! sequences are the simulated alignment — the inverse of what
//! [`crate::likelihood`] infers, and the standard way to generate test
//! data for inference methods.
//!
//! Optional discrete-gamma rate heterogeneity assigns each site a fixed
//! relative rate (drawn once, at the root) that scales every branch
//! length for that site.

use crate::error::{PhyloError, Result};
use crate::likelihood::gamma::DiscreteGamma;
use crate::likelihood::model::SubstModel;
use crate::rng::Rng;
use crate::tree::{NodeId, Tree};

/// A simulated alignment: one row per tree leaf.
#[derive(Debug, Clone, PartialEq)]
pub struct SimulatedAlignment {
    /// `(leaf label, nucleotide sequence)` pairs, in leaf-id order.
    pub rows: Vec<(String, Vec<u8>)>,
    /// Per-site relative rate actually used (all `1.0` when no gamma
    /// heterogeneity was requested).
    pub site_rates: Vec<f64>,
}

impl SimulatedAlignment {
    /// Sequence length (number of columns).
    pub fn length(&self) -> usize {
        self.rows.first().map(|(_, s)| s.len()).unwrap_or(0)
    }

    /// Number of taxa (rows).
    pub fn n_taxa(&self) -> usize {
        self.rows.len()
    }

    /// Looks up a row by leaf label.
    pub fn row(&self, label: &str) -> Option<&[u8]> {
        self.rows
            .iter()
            .find(|(l, _)| l == label)
            .map(|(_, s)| s.as_slice())
    }
}

/// Index 0..3 → nucleotide byte.
const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Simulates an alignment of length `length` down `tree` under `model`.
///
/// `gamma` is optional discrete-gamma rate heterogeneity; pass `None`
/// for a single rate. `seed` drives the deterministic RNG.
///
/// # Errors
/// - [`PhyloError::Invalid`] if `length` is zero or the tree has no
///   leaves.
/// - any error from [`SubstModel::transition_engine`].
pub fn simulate_sequences(
    tree: &Tree,
    model: &SubstModel,
    length: usize,
    gamma: Option<&DiscreteGamma>,
    seed: u64,
) -> Result<SimulatedAlignment> {
    if length == 0 {
        return Err(PhyloError::invalid("length", "must be positive"));
    }
    if tree.leaf_count() == 0 {
        return Err(PhyloError::invalid("tree", "tree has no leaves"));
    }
    let engine = model.transition_engine()?;
    let freqs = engine.frequencies();
    let mut rng = Rng::new(seed);

    // Per-site relative rates.
    let site_rates: Vec<f64> = match gamma {
        None => vec![1.0; length],
        Some(g) => (0..length)
            .map(|_| g.rates()[rng.below(g.n_categories())])
            .collect(),
    };

    let n = tree.node_count();
    // `seq[node]` — the node's evolved sequence as 0..3 state indices.
    let mut seq: Vec<Vec<u8>> = vec![Vec::new(); n];

    // Root sequence: draw each site from the equilibrium frequencies.
    let root = tree.root();
    seq[root] = (0..length)
        .map(|_| rng.weighted_index(&freqs) as u8)
        .collect();

    // Preorder: evolve each child from its parent.
    for &id in &tree.preorder() {
        if id == root {
            continue;
        }
        let parent = tree.node(id).parent.expect("non-root has a parent");
        let base_bl = tree.node(id).branch_length.unwrap_or(0.1).max(0.0);
        // Group sites by their rate so P(t) is built once per distinct
        // (rate · branch-length) value.
        let parent_seq = seq[parent].clone();
        let mut child_seq = vec![0u8; length];
        // Distinct rates present.
        let mut distinct: Vec<f64> = site_rates.clone();
        distinct.sort_by(|a, b| a.partial_cmp(b).unwrap());
        distinct.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
        for &rate in &distinct {
            let p = engine.p(base_bl * rate);
            for site in 0..length {
                if (site_rates[site] - rate).abs() > 1e-12 {
                    continue;
                }
                let from = parent_seq[site] as usize;
                // Sample the descendant base from P(t)'s `from` row.
                let row = [p[(from, 0)], p[(from, 1)], p[(from, 2)], p[(from, 3)]];
                child_seq[site] = rng.weighted_index(&row) as u8;
            }
        }
        seq[id] = child_seq;
    }

    // Collect leaf rows as nucleotide bytes.
    let mut rows = Vec::new();
    for &leaf in &tree.leaves() {
        let label = tree
            .node(leaf)
            .label
            .clone()
            .unwrap_or_else(|| format!("leaf{leaf}"));
        let bytes: Vec<u8> = seq[leaf].iter().map(|&s| BASES[s as usize]).collect();
        rows.push((label, bytes));
    }

    Ok(SimulatedAlignment { rows, site_rates })
}

/// Internal helper kept for symmetry with the inference side — the
/// nucleotide byte for a 0..3 state index.
#[allow(dead_code)]
fn base_byte(state: u8) -> u8 {
    BASES[state as usize % 4]
}

/// Internal helper — guards a node id against an out-of-range arena.
#[allow(dead_code)]
fn is_valid_node(tree: &Tree, id: NodeId) -> bool {
    id < tree.node_count()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::read_newick;

    #[test]
    fn simulates_an_alignment_of_the_right_shape() {
        let tree = read_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();
        let aln = simulate_sequences(&tree, &SubstModel::Jc69, 50, None, 42).unwrap();
        assert_eq!(aln.n_taxa(), 4);
        assert_eq!(aln.length(), 50);
        for (_, s) in &aln.rows {
            assert_eq!(s.len(), 50);
            assert!(s.iter().all(|&b| b"ACGT".contains(&b)));
        }
    }

    #[test]
    fn is_deterministic_for_a_seed() {
        let tree = read_newick("((A:0.2,B:0.2):0.2,C:0.2);").unwrap();
        let a = simulate_sequences(&tree, &SubstModel::Jc69, 40, None, 7).unwrap();
        let b = simulate_sequences(&tree, &SubstModel::Jc69, 40, None, 7).unwrap();
        assert_eq!(a.rows, b.rows);
    }

    #[test]
    fn short_branches_keep_sequences_similar() {
        // Near-zero branches => sister taxa almost identical.
        let tree = read_newick("((A:0.001,B:0.001):0.001,C:0.001);").unwrap();
        let aln = simulate_sequences(&tree, &SubstModel::Jc69, 200, None, 1).unwrap();
        let a = aln.row("A").unwrap();
        let b = aln.row("B").unwrap();
        let diffs = a.iter().zip(b).filter(|(x, y)| x != y).count();
        assert!(diffs < 20, "expected near-identical, got {diffs} diffs");
    }

    #[test]
    fn long_branches_diverge_sequences() {
        // Long branches => substantial divergence from the root.
        let tree = read_newick("(A:2.0,B:2.0);").unwrap();
        let aln = simulate_sequences(&tree, &SubstModel::Jc69, 300, None, 3).unwrap();
        let a = aln.row("A").unwrap();
        let b = aln.row("B").unwrap();
        let diffs = a.iter().zip(b).filter(|(x, y)| x != y).count();
        // Two long branches => well above zero divergence.
        assert!(diffs > 60, "expected divergence, got {diffs}");
    }

    #[test]
    fn gamma_heterogeneity_assigns_per_site_rates() {
        let tree = read_newick("((A:0.2,B:0.2):0.2,C:0.2);").unwrap();
        let gamma = DiscreteGamma::new(0.5, 4).unwrap();
        let aln = simulate_sequences(&tree, &SubstModel::Jc69, 100, Some(&gamma), 11).unwrap();
        assert_eq!(aln.site_rates.len(), 100);
        // With four categories at least two distinct rates should be
        // present.
        let mut rates = aln.site_rates.clone();
        rates.sort_by(|a, b| a.partial_cmp(b).unwrap());
        rates.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
        assert!(rates.len() >= 2, "expected varied site rates");
    }

    #[test]
    fn hky_model_simulation_runs() {
        let tree = read_newick("((A:0.3,B:0.3):0.2,(C:0.3,D:0.3):0.2);").unwrap();
        let model = SubstModel::Hky85 {
            kappa: 3.0,
            freqs: [0.3, 0.2, 0.25, 0.25],
        };
        let aln = simulate_sequences(&tree, &model, 80, None, 5).unwrap();
        assert_eq!(aln.length(), 80);
        assert_eq!(aln.n_taxa(), 4);
    }

    #[test]
    fn rejects_zero_length() {
        let tree = read_newick("(A,B);").unwrap();
        assert!(simulate_sequences(&tree, &SubstModel::Jc69, 0, None, 1).is_err());
    }
}
