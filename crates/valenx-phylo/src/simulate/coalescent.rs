//! Coalescent tree simulation.
//!
//! Kingman's coalescent (1982) models the genealogy of a sample of `n`
//! lineages from a population, running *backward* in time. With `k`
//! lineages currently extant, each pair coalesces independently, so the
//! waiting time to the next coalescence is exponential with rate
//! `C(k, 2) / N`, where `N` is the (diploid-scaled) effective
//! population size. At each event two uniformly-chosen lineages join.
//!
//! After `n - 1` coalescences a single lineage — the most recent common
//! ancestor — remains, and the sequence of join events is a rooted,
//! ultrametric tree whose branch lengths are in coalescent time units.
//!
//! [`PopulationSize`] supports a **constant** size and a
//! **piecewise-constant** history (a step function of `N` over time
//! intervals), so bottlenecks and expansions can be simulated.

use crate::error::{PhyloError, Result};
use crate::rng::Rng;
use crate::tree::{Node, NodeId, Tree};

/// An effective-population-size history.
#[derive(Debug, Clone, PartialEq)]
pub enum PopulationSize {
    /// A single constant size `N` for all time.
    Constant(f64),
    /// A piecewise-constant history: `(duration, size)` segments
    /// ordered from the present backward. The final segment's size
    /// extends to infinity (its duration is ignored).
    Piecewise(Vec<(f64, f64)>),
}

impl PopulationSize {
    /// The population size at a given time `t` before the present.
    fn size_at(&self, t: f64) -> f64 {
        match self {
            PopulationSize::Constant(n) => *n,
            PopulationSize::Piecewise(segs) => {
                let mut acc = 0.0;
                for (dur, size) in segs {
                    acc += dur;
                    if t < acc {
                        return *size;
                    }
                }
                // Past the last boundary: the last segment's size.
                segs.last().map(|(_, s)| *s).unwrap_or(1.0)
            }
        }
    }

    /// Validates the size history.
    fn validate(&self) -> Result<()> {
        match self {
            PopulationSize::Constant(n) => {
                if *n <= 0.0 {
                    return Err(PhyloError::invalid("population_size", "must be positive"));
                }
            }
            PopulationSize::Piecewise(segs) => {
                if segs.is_empty() {
                    return Err(PhyloError::invalid(
                        "population_size",
                        "piecewise history has no segments",
                    ));
                }
                if segs.iter().any(|&(_, s)| s <= 0.0) {
                    return Err(PhyloError::invalid(
                        "population_size",
                        "every segment size must be positive",
                    ));
                }
            }
        }
        Ok(())
    }
}

/// Simulates a coalescent genealogy for `n` sampled lineages.
///
/// `labels` names the `n` tips; `pop` is the effective-population-size
/// history; `seed` drives the deterministic RNG. The returned tree is
/// rooted and ultrametric, with branch lengths in coalescent time
/// (units of generations, scaled by `N`).
///
/// For a piecewise history the per-event exponential rate uses the
/// population size sampled at the current time, advancing the time
/// integral segment by segment so changing `N` correctly stretches or
/// compresses waiting times.
///
/// # Errors
/// [`PhyloError::Invalid`] if fewer than two labels are given or the
/// population history is invalid.
pub fn simulate_coalescent(labels: &[String], pop: &PopulationSize, seed: u64) -> Result<Tree> {
    if labels.len() < 2 {
        return Err(PhyloError::invalid("labels", "need at least two lineages"));
    }
    pop.validate()?;
    let mut rng = Rng::new(seed);
    let mut tree = Tree::building();

    // Active lineages: (node id, height at which the lineage starts).
    // Tips all start at height 0.
    let mut active: Vec<(NodeId, f64)> = labels
        .iter()
        .map(|l| {
            let id = tree.push_node(Node {
                label: Some(l.clone()),
                branch_length: None,
                parent: None,
                children: Vec::new(),
            });
            (id, 0.0)
        })
        .collect();

    let mut current_time = 0.0;
    while active.len() > 1 {
        let k = active.len() as f64;
        // Coalescent rate with k lineages: C(k,2) / N.
        // For a piecewise N we integrate forward: draw a unit-rate
        // exponential and convert to real time through the size
        // history (a step-function approximation by small advances).
        let pair_rate = k * (k - 1.0) / 2.0;
        let wait = sample_coalescent_wait(&mut rng, pair_rate, current_time, pop);
        current_time += wait;

        // Pick two distinct lineages to coalesce.
        let i = rng.below(active.len());
        let mut j = rng.below(active.len() - 1);
        if j >= i {
            j += 1;
        }
        let (ni, hi) = active[i];
        let (nj, hj) = active[j];

        // New internal node at `current_time`.
        let parent = tree.push_node(Node {
            label: None,
            branch_length: None,
            parent: None,
            children: vec![ni, nj],
        });
        tree.node_mut(ni).parent = Some(parent);
        tree.node_mut(ni).branch_length = Some((current_time - hi).max(0.0));
        tree.node_mut(nj).parent = Some(parent);
        tree.node_mut(nj).branch_length = Some((current_time - hj).max(0.0));

        // Replace the two coalesced lineages with their ancestor.
        let (lo, hi_idx) = if i < j { (i, j) } else { (j, i) };
        active.remove(hi_idx);
        active.remove(lo);
        active.push((parent, current_time));
    }

    let root = active[0].0;
    tree.finish_building(root, true)
        .map_err(|e| PhyloError::invalid_tree(e.to_string()))
}

/// Samples the waiting time to the next coalescence under a (possibly
/// piecewise) population-size history.
///
/// For a constant `N` this is just `Exp(pair_rate / N)`. For a
/// piecewise history it integrates the hazard: a unit-rate exponential
/// "amount of coalescent intensity" is drawn, then converted to real
/// time by walking the size segments and consuming intensity
/// `pair_rate · Δt / N(t)` until the budget is spent.
fn sample_coalescent_wait(
    rng: &mut Rng,
    pair_rate: f64,
    start_time: f64,
    pop: &PopulationSize,
) -> f64 {
    match pop {
        PopulationSize::Constant(n) => rng.exponential(pair_rate / n),
        PopulationSize::Piecewise(_) => {
            // Budget of cumulative hazard to consume.
            let mut budget = rng.exponential(1.0);
            let mut t = start_time;
            // Advance in small steps; each step's hazard is
            // pair_rate / N(t) · dt.
            let step = 0.01;
            for _ in 0..1_000_000 {
                let n = pop.size_at(t).max(1e-9);
                let hazard = pair_rate / n * step;
                if hazard >= budget {
                    // The event happens partway through this step.
                    let frac = budget / hazard.max(1e-300);
                    return (t + frac * step) - start_time;
                }
                budget -= hazard;
                t += step;
            }
            t - start_time
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn labels(n: usize) -> Vec<String> {
        (0..n).map(|i| format!("L{i}")).collect()
    }

    #[test]
    fn produces_a_valid_ultrametric_tree() {
        let t = simulate_coalescent(&labels(8), &PopulationSize::Constant(1.0), 42).unwrap();
        assert_eq!(t.leaf_count(), 8);
        assert!(t.validate().is_ok());
        assert!(t.rooted);
        // Ultrametric: all tips equidistant from the root.
        let root = t.root();
        let depths: Vec<f64> = t
            .leaves()
            .iter()
            .map(|&l| t.patristic_distance(root, l))
            .collect();
        let d0 = depths[0];
        for d in &depths {
            assert!((d - d0).abs() < 1e-9, "not ultrametric: {depths:?}");
        }
    }

    #[test]
    fn is_deterministic_for_a_seed() {
        let a = simulate_coalescent(&labels(6), &PopulationSize::Constant(2.0), 7).unwrap();
        let b = simulate_coalescent(&labels(6), &PopulationSize::Constant(2.0), 7).unwrap();
        assert_eq!(write_topology(&a), write_topology(&b));
    }

    #[test]
    fn larger_population_gives_a_taller_tree() {
        // Coalescence rate ∝ 1/N, so bigger N => longer waiting times.
        let small = simulate_coalescent(&labels(20), &PopulationSize::Constant(0.5), 1).unwrap();
        let large = simulate_coalescent(&labels(20), &PopulationSize::Constant(50.0), 1).unwrap();
        let height = |t: &Tree| -> f64 { t.patristic_distance(t.root(), t.leaves()[0]) };
        assert!(height(&large) > height(&small));
    }

    #[test]
    fn piecewise_history_runs_and_validates() {
        let pop = PopulationSize::Piecewise(vec![
            (0.5, 0.1), // recent bottleneck
            (1.0, 5.0), // larger ancestral population
        ]);
        let t = simulate_coalescent(&labels(10), &pop, 99).unwrap();
        assert_eq!(t.leaf_count(), 10);
        assert!(t.validate().is_ok());
        assert!(t.total_length() > 0.0);
    }

    #[test]
    fn rejects_bad_input() {
        assert!(simulate_coalescent(&labels(1), &PopulationSize::Constant(1.0), 1).is_err());
        assert!(simulate_coalescent(&labels(4), &PopulationSize::Constant(-1.0), 1).is_err());
        assert!(simulate_coalescent(&labels(4), &PopulationSize::Piecewise(vec![]), 1).is_err());
    }

    /// Topology fingerprint: sorted leaf set under each internal node.
    fn write_topology(t: &Tree) -> Vec<Vec<String>> {
        let mut v: Vec<Vec<String>> = (0..t.node_count())
            .filter(|&id| t.node(id).is_internal())
            .map(|id| {
                let mut names: Vec<String> = t
                    .descendant_leaves(id)
                    .into_iter()
                    .filter_map(|l| t.node(l).label.clone())
                    .collect();
                names.sort();
                names
            })
            .collect();
        v.sort();
        v
    }
}
