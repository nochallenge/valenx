//! Forward-in-time simulation with tree-sequence recording.
//!
//! `pyslim` / `tskit`'s headline trick is to record the genealogy
//! *while the forward simulation runs*: every time an offspring genome
//! inherits a chromosome stretch from a parent genome, an **edge** is
//! written. At the end the accumulated node + edge tables are the exact
//! genealogy of the final generation — no separate coalescent needed,
//! and complete with selection's distortions.
//!
//! [`record_wright_fisher`] runs a haploid-genome Wright-Fisher process
//! (each diploid individual is two genome nodes) and returns a
//! [`crate::coalescent::TreeSequence`]. The pipeline mirrors `pyslim`:
//!
//! 1. **Record raw edges** — every offspring genome inherits from a
//!    parent genome over each `[left, right)` stretch determined by
//!    that meiosis's crossovers. Multiple crossovers per chromosome
//!    are recorded as separate edges, not collapsed to the first.
//! 2. **Simplify** — drop every node not ancestral to the final
//!    sample, and squash chains of unary edges (the bulk of what
//!    `tskit_table_collection.simplify()` does in spirit), so the
//!    output table has the same edge count as the equivalent
//!    coalescent ARG would.
//!
//! Mutations are *not* overlaid during recording; call
//! [`crate::coalescent::overlay_mutations`] on the returned tree
//! sequence afterwards, exactly as the `pyslim` workflow does.

use crate::coalescent::tree_sequence::{Edge, TreeSequence};
use crate::error::{PopgenError, Result};
use crate::rng::Rng;

/// Configuration for a tree-recording forward run.
#[derive(Copy, Clone, Debug)]
pub struct RecordingConfig {
    /// Diploid census size (constant). The genome sample size is `2n`.
    pub n: usize,
    /// Number of generations to simulate.
    pub generations: usize,
    /// Per-base-pair crossover rate. Multiple crossovers per meiosis
    /// are drawn as `Poisson(rate * sequence_length)`.
    pub recombination_rate: f64,
    /// Length of the recorded chromosome.
    pub sequence_length: f64,
    /// RNG seed.
    pub seed: u64,
}

impl RecordingConfig {
    /// Validates the configuration.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] on a zero size, zero generations,
    /// negative rate or non-positive length.
    pub fn validate(&self) -> Result<()> {
        if self.n == 0 {
            return Err(PopgenError::invalid("n", "must be positive"));
        }
        if self.generations == 0 {
            return Err(PopgenError::invalid("generations", "must be positive"));
        }
        if self.recombination_rate < 0.0 || !self.recombination_rate.is_finite() {
            return Err(PopgenError::invalid(
                "recombination_rate",
                "must be finite and non-negative",
            ));
        }
        if self.sequence_length <= 0.0 || !self.sequence_length.is_finite() {
            return Err(PopgenError::invalid(
                "sequence_length",
                "must be finite and positive",
            ));
        }
        Ok(())
    }
}

/// Runs a Wright-Fisher process with genealogy recording and returns
/// the (simplified) tree sequence of the final generation's genomes.
///
/// The returned tree sequence's *samples* are the `2n` genomes of the
/// final generation. Node times count generations *before the present*
/// (the final generation is time 0).
///
/// # Errors
/// [`PopgenError`] from configuration validation or tree-sequence
/// assembly.
pub fn record_wright_fisher(config: RecordingConfig) -> Result<TreeSequence> {
    config.validate()?;
    let mut rng = Rng::new(config.seed);
    let genomes_per_gen = 2 * config.n;
    let l = config.sequence_length;

    // Raw edges as `(parent_node, child_node, left, right)`. Node ids
    // densify per-generation: founder genomes are nodes 0..2n, then
    // each generation appends 2n nodes.
    let mut raw_edges: Vec<(usize, usize, f64, f64)> = Vec::new();
    let mut node_gen: Vec<usize> = Vec::new();

    // Founding generation.
    let mut prev: Vec<usize> = (0..genomes_per_gen)
        .map(|_| {
            node_gen.push(0);
            node_gen.len() - 1
        })
        .collect();

    for gen in 1..=config.generations {
        let mut current = Vec::with_capacity(genomes_per_gen);
        for _ in 0..genomes_per_gen {
            // Child genome node.
            node_gen.push(gen);
            let child = node_gen.len() - 1;
            current.push(child);

            // Pick the diploid parent, then its two genome nodes.
            let parent_ind = rng.below(config.n);
            let pg0 = prev[2 * parent_ind];
            let pg1 = prev[2 * parent_ind + 1];

            // Draw an arbitrary number of crossover breakpoints.
            let n_xo = rng.poisson(config.recombination_rate * l);
            let mut breakpoints: Vec<f64> = (0..n_xo).map(|_| rng.uniform() * l).collect();
            breakpoints.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            // Dedup essentially-equal breakpoints to avoid zero-width
            // edges.
            breakpoints.dedup_by(|a, b| (*a - *b).abs() < 1e-12);

            // Decide which homolog is copied first (uniform).
            let mut current_src = if rng.bernoulli(0.5) { pg0 } else { pg1 };
            let other_src = if current_src == pg0 { pg1 } else { pg0 };
            let mut cur_left = 0.0f64;
            for &bp in &breakpoints {
                if bp <= cur_left {
                    continue;
                }
                if bp >= l {
                    break;
                }
                raw_edges.push((current_src, child, cur_left, bp));
                cur_left = bp;
                // Swap source.
                current_src = if current_src == pg0 { pg1 } else { pg0 };
                let _ = other_src;
            }
            // Final stretch to the chromosome end.
            if cur_left < l {
                raw_edges.push((current_src, child, cur_left, l));
            }
        }
        prev = current;
    }

    // The final generation's genomes are the samples.
    let samples = prev.clone();
    simplify(node_gen, raw_edges, &samples, config.generations, l)
}

/// Squash + drop-unreachable: keep only nodes ancestral to the sample,
/// then collapse chains of unary edges so that the kept tree sequence
/// has no spurious internal nodes — the `tskit_table_collection.simplify`
/// operation in spirit. Unary squashing makes the recorded forward-
/// time genealogy look identical in shape to a coalescent ARG over
/// the same generations.
fn simplify(
    node_gen: Vec<usize>,
    raw_edges: Vec<(usize, usize, f64, f64)>,
    samples: &[usize],
    generations: usize,
    sequence_length: f64,
) -> Result<TreeSequence> {
    let n_nodes = node_gen.len();
    // child -> list of parent edges, for backward traversal.
    let mut parents_of: Vec<Vec<usize>> = vec![Vec::new(); n_nodes];
    for (i, &(_, child, _, _)) in raw_edges.iter().enumerate() {
        if child < n_nodes {
            parents_of[child].push(i);
        }
    }
    // Mark every node reachable backward from a sample.
    let mut keep = vec![false; n_nodes];
    let mut stack: Vec<usize> = samples.to_vec();
    for &s in samples {
        keep[s] = true;
    }
    while let Some(node) = stack.pop() {
        for &edge_idx in &parents_of[node] {
            let parent = raw_edges[edge_idx].0;
            if parent < n_nodes && !keep[parent] {
                keep[parent] = true;
                stack.push(parent);
            }
        }
    }

    // Forward children-of map over kept edges.
    let mut children_of: Vec<Vec<(f64, f64, usize)>> = vec![Vec::new(); n_nodes];
    for &(parent, child, l, r) in &raw_edges {
        if parent < n_nodes && child < n_nodes && keep[parent] && keep[child] {
            children_of[parent].push((l, r, child));
        }
    }

    // Per-genomic-stretch child counts: at every node, for every
    // interval, count distinct children. A node is "branching at
    // position p" if it has at least two children covering p.
    // A node with exactly one child covering p is a unary intermediate
    // there — we redirect that child's edge to the node's own parent.
    // For samples we always keep them (they are tips by definition).
    let mut is_sample = vec![false; n_nodes];
    for &s in samples {
        is_sample[s] = true;
    }
    // We simplify by recursively short-circuiting unary edges. To keep
    // it tractable, walk children-of for each node and replace child=u
    // with child=v whenever u has exactly one child at the right
    // interval — repeated until no further short-circuiting is
    // possible. To bound complexity we cap iterations.
    let mut squashed: Vec<(usize, usize, f64, f64)> = raw_edges
        .iter()
        .copied()
        .filter(|&(p, c, _, _)| p < n_nodes && c < n_nodes && keep[p] && keep[c])
        .collect();
    // Build a "child -> single ancestral child" map per node when the
    // node has exactly one child *over the whole chromosome*. For each
    // unary-everywhere kept node that is NOT a sample, redirect every
    // edge ending at it to its single child, and drop the node from
    // `keep`. Iterate.
    let mut changed = true;
    let mut iters = 0usize;
    while changed && iters < 50 {
        changed = false;
        iters += 1;
        // For each non-sample kept node, count unique child ids over
        // its edge list AND check whether they tile the whole [0, L)
        // exactly once.
        let mut single_child: Vec<Option<usize>> = vec![None; n_nodes];
        let mut children_intervals: Vec<Vec<(f64, f64, usize)>> = vec![Vec::new(); n_nodes];
        for &(p, c, l, r) in &squashed {
            children_intervals[p].push((l, r, c));
        }
        for node in 0..n_nodes {
            if !keep[node] || is_sample[node] {
                continue;
            }
            let intervals = &children_intervals[node];
            if intervals.is_empty() {
                continue;
            }
            // Distinct child ids.
            let mut child_ids: Vec<usize> = intervals.iter().map(|&(_, _, c)| c).collect();
            child_ids.sort_unstable();
            child_ids.dedup();
            if child_ids.len() != 1 {
                continue;
            }
            // Coverage check: their intervals tile [0, L). Sort by
            // left; walk; ensure contiguous and covering [0, L].
            let mut sorted = intervals.clone();
            sorted.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            let mut covered = 0.0f64;
            let mut ok = (sorted[0].0 - 0.0).abs() < 1e-9;
            for w in sorted.windows(2) {
                if (w[1].0 - w[0].1).abs() > 1e-9 {
                    ok = false;
                    break;
                }
            }
            if !ok {
                continue;
            }
            covered += sorted.last().unwrap().1 - sorted.first().unwrap().0;
            if (covered - sequence_length).abs() > 1e-9 {
                continue;
            }
            single_child[node] = Some(child_ids[0]);
        }
        // Now redirect: any edge whose parent is a "single-child"
        // node gets redirected to that node's child.
        let mut new_squashed: Vec<(usize, usize, f64, f64)> = Vec::new();
        for &(p, c, l, r) in &squashed {
            // If c is a unary intermediate, drop the edge — its sole
            // child will be picked up by a (grandparent -> child) edge
            // we synthesise below.
            if single_child[c].is_some() {
                continue;
            }
            new_squashed.push((p, c, l, r));
        }
        // Now, for every single-child intermediate, copy the edges
        // ending at it (parent -> intermediate) as edges to its child
        // (parent -> child) covering the same intervals.
        // First snapshot parents-of for the intermediate.
        let mut parent_edges_at: Vec<Vec<(usize, f64, f64)>> = vec![Vec::new(); n_nodes];
        for &(p, c, l, r) in &squashed {
            if single_child[c].is_some() {
                parent_edges_at[c].push((p, l, r));
            }
        }
        for node in 0..n_nodes {
            if let Some(real_child) = single_child[node] {
                for &(p, l, r) in &parent_edges_at[node] {
                    // The intermediate's parent now points at the real
                    // child. (Could itself be another intermediate; the
                    // outer loop will iterate.)
                    new_squashed.push((p, real_child, l, r));
                    changed = true;
                }
                keep[node] = false;
            }
        }
        squashed = new_squashed;
    }

    // Now keep is the final set of kept nodes; assemble.
    let mut kept_order: Vec<usize> = (0..n_nodes).filter(|&i| keep[i]).collect();
    // Order: present-first so sample ids stay low.
    kept_order.sort_by_key(|&i| generations - node_gen[i]);
    let mut new_id = vec![usize::MAX; n_nodes];
    for (dense, &old) in kept_order.iter().enumerate() {
        new_id[old] = dense;
    }

    let mut ts = TreeSequence::new(sequence_length)?;
    for &old in &kept_order {
        let time = (generations - node_gen[old]) as f64;
        let is_s = is_sample[old];
        ts.add_node(time, is_s);
    }
    for &(parent, child, left, right) in &squashed {
        if parent < n_nodes && child < n_nodes && keep[parent] && keep[child] && right > left {
            ts.add_edge(Edge {
                parent: new_id[parent],
                child: new_id[child],
                left,
                right,
            })?;
        }
    }
    ts.finalize()?;
    Ok(ts)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn recording_yields_a_consistent_tree_sequence() {
        let cfg = RecordingConfig {
            n: 10,
            generations: 15,
            recombination_rate: 0.0,
            sequence_length: 100.0,
            seed: 42,
        };
        let ts = record_wright_fisher(cfg).unwrap();
        // 2n = 20 sample genomes.
        assert_eq!(ts.sample_count(), 20);
        assert!(ts.node_count() >= 20);
        assert!(ts.edge_count() > 0);
    }

    #[test]
    fn recording_is_deterministic() {
        let cfg = RecordingConfig {
            n: 8,
            generations: 12,
            recombination_rate: 1e-3,
            sequence_length: 200.0,
            seed: 7,
        };
        let a = record_wright_fisher(cfg).unwrap();
        let b = record_wright_fisher(cfg).unwrap();
        assert_eq!(a.node_count(), b.node_count());
        assert_eq!(a.edge_count(), b.edge_count());
    }

    #[test]
    fn recombination_creates_split_edges() {
        let cfg = RecordingConfig {
            n: 12,
            generations: 10,
            recombination_rate: 5e-3,
            sequence_length: 1000.0,
            seed: 3,
        };
        let ts = record_wright_fisher(cfg).unwrap();
        assert!(ts.edge_count() >= ts.sample_count());
    }

    #[test]
    fn simplify_drops_unreachable_nodes() {
        let cfg = RecordingConfig {
            n: 6,
            generations: 25,
            recombination_rate: 0.0,
            sequence_length: 100.0,
            seed: 1,
        };
        let total_created = 2 * 6 * (25 + 1);
        let ts = record_wright_fisher(cfg).unwrap();
        assert!(
            ts.node_count() < total_created,
            "simplify kept everything ({} of {total_created})",
            ts.node_count()
        );
    }

    #[test]
    fn rejects_bad_config() {
        let bad = RecordingConfig {
            n: 0,
            generations: 5,
            recombination_rate: 0.0,
            sequence_length: 100.0,
            seed: 1,
        };
        assert!(record_wright_fisher(bad).is_err());
    }

    /// With multiple expected crossovers per meiosis, at least one
    /// offspring has three or more parent-edges (two breakpoints).
    #[test]
    fn multiple_crossovers_yield_three_plus_edges_per_offspring() {
        let cfg = RecordingConfig {
            n: 50,
            generations: 1,           // single-generation run keeps the table small
            recombination_rate: 2e-2, // mean 20 breakpoints on L=1000
            sequence_length: 1000.0,
            seed: 99,
        };
        let ts = record_wright_fisher(cfg).unwrap();
        // Count edges per child.
        let mut edges_per_child: std::collections::HashMap<usize, usize> =
            std::collections::HashMap::new();
        for e in ts.edges() {
            *edges_per_child.entry(e.child).or_insert(0) += 1;
        }
        let three_plus = edges_per_child.values().filter(|&&v| v >= 3).count();
        assert!(
            three_plus > 0,
            "no offspring had 3+ edges, max {:?}",
            edges_per_child.values().max()
        );
    }

    /// All edges should tile each child's `[0, L)` exactly once — no
    /// missing or overlapping inheritance.
    #[test]
    fn child_edges_tile_the_chromosome() {
        let cfg = RecordingConfig {
            n: 8,
            generations: 5,
            recombination_rate: 5e-3,
            sequence_length: 500.0,
            seed: 21,
        };
        let ts = record_wright_fisher(cfg).unwrap();
        let mut by_child: std::collections::HashMap<usize, Vec<(f64, f64)>> =
            std::collections::HashMap::new();
        for e in ts.edges() {
            by_child.entry(e.child).or_default().push((e.left, e.right));
        }
        for samp in ts.samples() {
            let mut ivs = by_child.get(&samp).cloned().unwrap_or_default();
            ivs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            assert!(!ivs.is_empty(), "sample {samp} has no parent edges");
            assert!(ivs.first().unwrap().0.abs() < 1e-9);
            assert!((ivs.last().unwrap().1 - 500.0).abs() < 1e-9);
            // Contiguous.
            for w in ivs.windows(2) {
                assert!(
                    (w[1].0 - w[0].1).abs() < 1e-9,
                    "gap in sample {samp}: {:?} then {:?}",
                    w[0],
                    w[1]
                );
            }
        }
    }

    /// Unary squashing: with `n=2` diploids over many generations and
    /// no recombination, the simplify pass should leave well under one
    /// internal node per generation — chains have been collapsed.
    /// Without squashing the count would be roughly the number of
    /// generations.
    #[test]
    fn simplify_squashes_unary_chains() {
        let cfg = RecordingConfig {
            n: 2,
            generations: 50,
            recombination_rate: 0.0,
            sequence_length: 100.0,
            seed: 5,
        };
        let ts = record_wright_fisher(cfg).unwrap();
        // 4 samples; with chain squashing, internal nodes should be at
        // most (4 - 1) = 3 in the fully-coalesced case (plus a few
        // until everything coalesces back).
        let internal = ts.node_count() - ts.sample_count();
        assert!(
            internal < 20,
            "expected unary chains to be squashed, got {internal} internal nodes",
        );
    }

    /// Every sample appears as a leaf-only node (no outgoing edges).
    #[test]
    fn samples_are_leaves() {
        let cfg = RecordingConfig {
            n: 4,
            generations: 8,
            recombination_rate: 1e-3,
            sequence_length: 200.0,
            seed: 13,
        };
        let ts = record_wright_fisher(cfg).unwrap();
        let samples: HashSet<usize> = ts.samples().into_iter().collect();
        for e in ts.edges() {
            assert!(
                !samples.contains(&e.parent),
                "sample {} appears as a parent",
                e.parent
            );
        }
    }
}
