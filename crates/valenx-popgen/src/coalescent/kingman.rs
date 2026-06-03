//! The backward (Kingman) coalescent simulator.
//!
//! Kingman's coalescent (1982) models the genealogy of a sample of `n`
//! lineages by running *backward* in time. With `k` lineages extant,
//! every pair coalesces independently, so the waiting time to the next
//! coalescence is exponential with rate `C(k,2) / N` (in generations,
//! `N` the diploid effective size). At each event two uniformly-chosen
//! lineages join. After `n - 1` coalescences a single lineage — the
//! MRCA — remains and the join history is a rooted ultrametric tree.
//!
//! This module produces the genealogy as a [`valenx_phylo::Tree`],
//! reusing all of `valenx-phylo`'s tree machinery. It supports:
//!
//! - a **constant** effective size,
//! - a **piecewise-constant** size history ([`PopHistory`]), so
//!   bottlenecks and expansions stretch or compress branch lengths,
//! - a **structured** sample drawn from several demes connected by
//!   migration ([`structured_coalescent`]).

use crate::error::{PopgenError, Result};
use crate::rng::Rng;
use valenx_phylo::tree::{Node as PhyloNode, NodeId, Tree};

/// An effective-population-size history played backward in time.
#[derive(Clone, Debug, PartialEq)]
pub enum PopHistory {
    /// A single constant diploid effective size for all time.
    Constant(f64),
    /// A piecewise-constant history: `(duration, size)` segments
    /// ordered from the present backward. The final segment extends to
    /// infinity (its duration is ignored).
    Piecewise(Vec<(f64, f64)>),
}

impl PopHistory {
    /// The effective size at time `t` before the present.
    pub fn size_at(&self, t: f64) -> f64 {
        match self {
            PopHistory::Constant(n) => *n,
            PopHistory::Piecewise(segs) => {
                let mut acc = 0.0;
                for (dur, size) in segs {
                    acc += dur;
                    if t < acc {
                        return *size;
                    }
                }
                segs.last().map(|(_, s)| *s).unwrap_or(1.0)
            }
        }
    }

    /// Validates the history.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] on a non-positive size or empty
    /// piecewise list.
    pub fn validate(&self) -> Result<()> {
        match self {
            PopHistory::Constant(n) => {
                if *n <= 0.0 {
                    return Err(PopgenError::invalid(
                        "effective_size",
                        "must be positive",
                    ));
                }
            }
            PopHistory::Piecewise(segs) => {
                if segs.is_empty() {
                    return Err(PopgenError::invalid(
                        "pop_history",
                        "piecewise history has no segments",
                    ));
                }
                if segs.iter().any(|&(_, s)| s <= 0.0) {
                    return Err(PopgenError::invalid(
                        "pop_history",
                        "every segment size must be positive",
                    ));
                }
            }
        }
        Ok(())
    }
}

/// Simulates a Kingman coalescent genealogy for `n` lineages.
///
/// `labels` names the tips, `pop` is the size history and `seed` fixes
/// the deterministic RNG. The returned tree is rooted, ultrametric and
/// has branch lengths in generations.
///
/// # Errors
/// [`PopgenError::Invalid`] if fewer than two labels are given or the
/// history is invalid; [`PopgenError::Model`] if tree assembly fails.
pub fn coalescent(labels: &[String], pop: &PopHistory, seed: u64) -> Result<Tree> {
    if labels.len() < 2 {
        return Err(PopgenError::invalid(
            "labels",
            "need at least two lineages",
        ));
    }
    pop.validate()?;
    let mut rng = Rng::new(seed);

    // Phylo arena under construction. Tips first.
    let mut nodes: Vec<PhyloNode> = labels
        .iter()
        .map(|l| PhyloNode {
            label: Some(l.clone()),
            branch_length: None,
            parent: None,
            children: Vec::new(),
        })
        .collect();
    // Active lineages: (node id, height at which the lineage starts).
    let mut active: Vec<(NodeId, f64)> = (0..labels.len()).map(|i| (i, 0.0)).collect();
    let mut current_time = 0.0;

    while active.len() > 1 {
        let k = active.len() as f64;
        let pair_rate = k * (k - 1.0) / 2.0;
        let wait = sample_wait(&mut rng, pair_rate, current_time, pop);
        current_time += wait;

        // Pick two distinct lineages.
        let i = rng.below(active.len());
        let mut j = rng.below(active.len() - 1);
        if j >= i {
            j += 1;
        }
        let (ni, hi) = active[i];
        let (nj, hj) = active[j];

        let parent = nodes.len();
        nodes.push(PhyloNode {
            label: None,
            branch_length: None,
            parent: None,
            children: vec![ni, nj],
        });
        nodes[ni].parent = Some(parent);
        nodes[ni].branch_length = Some((current_time - hi).max(0.0));
        nodes[nj].parent = Some(parent);
        nodes[nj].branch_length = Some((current_time - hj).max(0.0));

        let (lo, hi_idx) = if i < j { (i, j) } else { (j, i) };
        active.remove(hi_idx);
        active.remove(lo);
        active.push((parent, current_time));
    }

    let root = active[0].0;
    Tree::new(nodes, root, true).map_err(|e| PopgenError::model(e.to_string()))
}

/// Samples the waiting time to the next coalescence under a (possibly
/// piecewise) size history. For a constant `N` this is `Exp(pair_rate /
/// N)`; for a piecewise history it integrates the time-varying hazard.
fn sample_wait(rng: &mut Rng, pair_rate: f64, start: f64, pop: &PopHistory) -> f64 {
    match pop {
        PopHistory::Constant(n) => rng.exponential(pair_rate / n),
        PopHistory::Piecewise(_) => {
            let mut budget = rng.exponential(1.0);
            let mut t = start;
            let step = 0.5;
            for _ in 0..2_000_000 {
                let n = pop.size_at(t).max(1e-9);
                let hazard = pair_rate / n * step;
                if hazard >= budget {
                    let frac = budget / hazard.max(1e-300);
                    return (t + frac * step) - start;
                }
                budget -= hazard;
                t += step;
            }
            t - start
        }
    }
}

/// Simulates a **structured** coalescent: the sample is split across
/// `deme_sizes.len()` demes, each of its own effective size, with
/// symmetric migration at backward rate `migration_rate` per lineage
/// per generation.
///
/// Two lineages can coalesce only when they are in the *same* deme.
/// Migration events move a single lineage between demes. The result is
/// a single genealogy (a [`valenx_phylo::Tree`]) once every lineage has
/// coalesced.
///
/// `sample_per_deme[d]` is how many tips start in deme `d`; the tip
/// labels are taken in order from `labels`.
///
/// # Errors
/// [`PopgenError::Invalid`] / [`PopgenError::Dimension`] on malformed
/// input; [`PopgenError::Model`] on tree-assembly failure.
pub fn structured_coalescent(
    labels: &[String],
    sample_per_deme: &[usize],
    deme_sizes: &[f64],
    migration_rate: f64,
    seed: u64,
) -> Result<Tree> {
    if deme_sizes.is_empty() {
        return Err(PopgenError::invalid("deme_sizes", "no demes"));
    }
    if sample_per_deme.len() != deme_sizes.len() {
        return Err(PopgenError::dimension(
            deme_sizes.len(),
            sample_per_deme.len(),
            "sample-per-deme vector",
        ));
    }
    if deme_sizes.iter().any(|&n| n <= 0.0) {
        return Err(PopgenError::invalid(
            "deme_sizes",
            "every deme size must be positive",
        ));
    }
    if migration_rate < 0.0 {
        return Err(PopgenError::invalid(
            "migration_rate",
            "must be non-negative",
        ));
    }
    let total: usize = sample_per_deme.iter().sum();
    if total < 2 {
        return Err(PopgenError::invalid(
            "sample_per_deme",
            "need at least two sampled lineages",
        ));
    }
    if labels.len() != total {
        return Err(PopgenError::dimension(
            total,
            labels.len(),
            "tip labels",
        ));
    }
    let mut rng = Rng::new(seed);

    let mut nodes: Vec<PhyloNode> = labels
        .iter()
        .map(|l| PhyloNode {
            label: Some(l.clone()),
            branch_length: None,
            parent: None,
            children: Vec::new(),
        })
        .collect();
    // Active lineage: (node id, height, deme).
    let mut active: Vec<(NodeId, f64, usize)> = Vec::with_capacity(total);
    let mut next_label = 0;
    for (deme, &count) in sample_per_deme.iter().enumerate() {
        for _ in 0..count {
            active.push((next_label, 0.0, deme));
            next_label += 1;
        }
    }
    let mut current_time = 0.0;
    let guard = 10_000_000usize;
    let mut iterations = 0;

    while active.len() > 1 {
        iterations += 1;
        if iterations > guard {
            return Err(PopgenError::model(
                "structured coalescent failed to converge",
            ));
        }
        // Total coalescent hazard: per deme C(k_d, 2) / N_d.
        let mut coal_hazard = 0.0;
        for (deme, &n) in deme_sizes.iter().enumerate() {
            let k = active.iter().filter(|&&(_, _, d)| d == deme).count() as f64;
            coal_hazard += k * (k - 1.0) / 2.0 / n;
        }
        // Total migration hazard: each lineage migrates at
        // migration_rate.
        let mig_hazard = migration_rate * active.len() as f64;
        let total_hazard = coal_hazard + mig_hazard;
        if total_hazard <= 0.0 {
            // No coalescence possible and no migration: lineages are
            // stranded in distinct demes forever.
            return Err(PopgenError::model(
                "structured sample cannot coalesce (zero migration, isolated demes)",
            ));
        }
        let wait = rng.exponential(total_hazard);
        current_time += wait;

        if rng.uniform() * total_hazard < coal_hazard {
            // A coalescence: pick the deme proportional to its hazard.
            let deme = pick_coalescing_deme(&active, deme_sizes, &mut rng);
            let in_deme: Vec<usize> = (0..active.len())
                .filter(|&i| active[i].2 == deme)
                .collect();
            let a = in_deme[rng.below(in_deme.len())];
            let mut b_idx = rng.below(in_deme.len() - 1);
            // pick a distinct second lineage in the deme
            let in_deme_minus: Vec<usize> =
                in_deme.iter().copied().filter(|&i| i != a).collect();
            b_idx %= in_deme_minus.len();
            let b = in_deme_minus[b_idx];

            let (na, ha, _) = active[a];
            let (nb, hb, _) = active[b];
            let parent = nodes.len();
            nodes.push(PhyloNode {
                label: None,
                branch_length: None,
                parent: None,
                children: vec![na, nb],
            });
            nodes[na].parent = Some(parent);
            nodes[na].branch_length = Some((current_time - ha).max(0.0));
            nodes[nb].parent = Some(parent);
            nodes[nb].branch_length = Some((current_time - hb).max(0.0));

            let (lo, hi) = if a < b { (a, b) } else { (b, a) };
            active.remove(hi);
            active.remove(lo);
            active.push((parent, current_time, deme));
        } else {
            // A migration: move one lineage to a uniformly-chosen
            // different deme.
            let idx = rng.below(active.len());
            if deme_sizes.len() > 1 {
                let mut to = rng.below(deme_sizes.len() - 1);
                if to >= active[idx].2 {
                    to += 1;
                }
                active[idx].2 = to;
            }
        }
    }

    let root = active[0].0;
    Tree::new(nodes, root, true).map_err(|e| PopgenError::model(e.to_string()))
}

/// Picks the deme in which the next coalescence happens, weighting by
/// each deme's coalescent hazard.
fn pick_coalescing_deme(
    active: &[(NodeId, f64, usize)],
    deme_sizes: &[f64],
    rng: &mut Rng,
) -> usize {
    let weights: Vec<f64> = deme_sizes
        .iter()
        .enumerate()
        .map(|(deme, &n)| {
            let k = active.iter().filter(|&&(_, _, d)| d == deme).count() as f64;
            k * (k - 1.0) / 2.0 / n
        })
        .collect();
    rng.weighted_index(&weights)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn labels(n: usize) -> Vec<String> {
        (0..n).map(|i| format!("L{i}")).collect()
    }

    fn tree_height(t: &Tree) -> f64 {
        t.patristic_distance(t.root(), t.leaves()[0])
    }

    #[test]
    fn coalescent_is_ultrametric_and_valid() {
        let t = coalescent(&labels(8), &PopHistory::Constant(1000.0), 42).unwrap();
        assert_eq!(t.leaf_count(), 8);
        assert!(t.validate().is_ok());
        let h0 = tree_height(&t);
        for &l in &t.leaves() {
            let h = t.patristic_distance(t.root(), l);
            assert!((h - h0).abs() < 1e-6, "not ultrametric");
        }
    }

    #[test]
    fn coalescent_is_deterministic() {
        let a = coalescent(&labels(6), &PopHistory::Constant(500.0), 7).unwrap();
        let b = coalescent(&labels(6), &PopHistory::Constant(500.0), 7).unwrap();
        assert!((a.total_length() - b.total_length()).abs() < 1e-9);
    }

    #[test]
    fn larger_ne_gives_a_deeper_tree() {
        let small =
            coalescent(&labels(20), &PopHistory::Constant(100.0), 1).unwrap();
        let large =
            coalescent(&labels(20), &PopHistory::Constant(10_000.0), 1).unwrap();
        assert!(tree_height(&large) > tree_height(&small));
    }

    #[test]
    fn expected_tmrca_is_about_right() {
        // E[TMRCA] for n lineages = 2N(1 - 1/n) generations.
        let n_lineages = 10usize;
        let ne = 1000.0;
        let mut acc = 0.0;
        let reps = 300;
        for seed in 0..reps {
            let t = coalescent(&labels(n_lineages), &PopHistory::Constant(ne), seed)
                .unwrap();
            acc += tree_height(&t);
        }
        let mean = acc / reps as f64;
        let expected = 2.0 * ne * (1.0 - 1.0 / n_lineages as f64);
        // Coalescent variance is large; allow a wide tolerance.
        assert!(
            (mean - expected).abs() / expected < 0.2,
            "mean TMRCA {mean} vs expected {expected}"
        );
    }

    #[test]
    fn piecewise_history_runs() {
        let pop = PopHistory::Piecewise(vec![(50.0, 10.0), (100.0, 5000.0)]);
        let t = coalescent(&labels(10), &pop, 99).unwrap();
        assert_eq!(t.leaf_count(), 10);
        assert!(t.validate().is_ok());
    }

    #[test]
    fn coalescent_rejects_bad_input() {
        assert!(coalescent(&labels(1), &PopHistory::Constant(1.0), 1).is_err());
        assert!(coalescent(&labels(4), &PopHistory::Constant(-1.0), 1).is_err());
        assert!(
            coalescent(&labels(4), &PopHistory::Piecewise(vec![]), 1).is_err()
        );
    }

    #[test]
    fn structured_coalescent_builds_a_tree() {
        let t = structured_coalescent(
            &labels(8),
            &[4, 4],
            &[1000.0, 1000.0],
            0.001,
            42,
        )
        .unwrap();
        assert_eq!(t.leaf_count(), 8);
        assert!(t.validate().is_ok());
    }

    #[test]
    fn isolated_demes_cannot_coalesce() {
        // Two demes, lineages in both, zero migration -> failure.
        let r = structured_coalescent(
            &labels(4),
            &[2, 2],
            &[1000.0, 1000.0],
            0.0,
            1,
        );
        assert!(r.is_err());
    }

    #[test]
    fn structured_rejects_dimension_mismatch() {
        let r = structured_coalescent(
            &labels(4),
            &[2, 2],
            &[1000.0], // wrong length
            0.01,
            1,
        );
        assert!(r.is_err());
    }
}
