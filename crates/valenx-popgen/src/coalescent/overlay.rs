//! Mutation overlay: turning a genealogy into genotypes.
//!
//! A coalescent run produces a genealogy with no mutations — branch
//! lengths only. Mutations are dropped on afterward, the way `ms` and
//! `msprime` work: the number of mutations on a branch is
//! `Poisson(mu * branch_length)`, and a mutation on a branch is
//! inherited by *exactly the descendants of that branch's child node*.
//!
//! Two entry points are provided:
//!
//! - [`overlay_on_tree`] — drops infinite-sites mutations on a single
//!   [`valenx_phylo::Tree`] genealogy and returns the resulting
//!   [`GenotypeMatrix`] (one row per leaf).
//! - [`overlay_mutations`] — does the same on a
//!   [`TreeSequence`], correctly handling recombination: a mutation
//!   falls at a genomic position, on an edge covering that position,
//!   and is inherited by the descendants *under the local tree at that
//!   position*.

use crate::coalescent::tree_sequence::{TreeSequence, TsMutation};
use crate::error::{PopgenError, Result};
use crate::infer::GenotypeMatrix;
use crate::rng::Rng;
use valenx_phylo::tree::{NodeId, Tree};

/// Drops infinite-sites mutations on a single genealogy.
///
/// `mutation_rate` is the per-branch-length-unit mutation rate; with a
/// coalescent tree whose branch lengths are in generations and a
/// per-site sequence of length `L`, pass `mu * L` to mutate the whole
/// segment. Each mutation lands on a uniform `[0, 1)` position and is
/// inherited by the leaves descending from the mutated branch.
///
/// The returned [`GenotypeMatrix`] has one row per leaf, ordered by the
/// tree's ascending leaf-id order, and one column per mutation.
///
/// # Errors
/// [`PopgenError::Invalid`] on a negative rate;
/// [`PopgenError::Model`] on genotype-matrix assembly failure.
pub fn overlay_on_tree(tree: &Tree, mutation_rate: f64, seed: u64) -> Result<GenotypeMatrix> {
    if mutation_rate < 0.0 {
        return Err(PopgenError::invalid(
            "mutation_rate",
            "must be non-negative",
        ));
    }
    let mut rng = Rng::new(seed);
    let leaves = tree.leaves();
    let leaf_index: std::collections::HashMap<NodeId, usize> =
        leaves.iter().enumerate().map(|(i, &l)| (l, i)).collect();

    // For each non-root node, draw mutations on its incoming branch.
    // A mutation is (position, set of leaf rows that carry it).
    let mut variants: Vec<(f64, Vec<usize>)> = Vec::new();
    for id in 0..tree.node_count() {
        let node = tree.node(id);
        let bl = match node.branch_length {
            Some(b) => b,
            None => continue, // the root
        };
        let n_mut = rng.poisson(mutation_rate * bl.max(0.0));
        if n_mut == 0 {
            continue;
        }
        // Descendant leaves of this node carry every mutation here.
        let carriers: Vec<usize> = tree
            .descendant_leaves(id)
            .into_iter()
            .filter_map(|l| leaf_index.get(&l).copied())
            .collect();
        for _ in 0..n_mut {
            variants.push((rng.uniform(), carriers.clone()));
        }
    }
    // Sort variants by position for a tidy matrix.
    variants.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let n_rows = leaves.len();
    let n_cols = variants.len();
    let mut rows = vec![vec![0u8; n_cols]; n_rows];
    let mut positions = Vec::with_capacity(n_cols);
    for (col, (pos, carriers)) in variants.into_iter().enumerate() {
        positions.push(pos);
        for r in carriers {
            rows[r][col] = 1;
        }
    }
    GenotypeMatrix::from_rows(rows, positions).map_err(|e| PopgenError::model(e.to_string()))
}

/// Drops infinite-sites mutations onto a [`TreeSequence`], writing the
/// site and mutation tables, and returns the sample genotype matrix.
///
/// Mutations respect recombination: a mutation at genomic `position` is
/// placed on an edge covering `position`, and the carriers are the
/// samples below that edge's child *in the local tree at `position`*.
///
/// `mutation_rate` is the per-base-pair per-generation rate; the
/// expected number of mutations on an edge covering `[l, r)` for time
/// `dt` is `mu * (r - l) * dt`.
///
/// # Errors
/// [`PopgenError::Invalid`] on a negative rate;
/// [`PopgenError::Model`] on table or matrix assembly failure.
pub fn overlay_mutations(
    ts: &mut TreeSequence,
    mutation_rate: f64,
    seed: u64,
) -> Result<GenotypeMatrix> {
    if mutation_rate < 0.0 {
        return Err(PopgenError::invalid(
            "mutation_rate",
            "must be non-negative",
        ));
    }
    let mut rng = Rng::new(seed);
    let samples = ts.samples();
    let sample_index: std::collections::HashMap<usize, usize> =
        samples.iter().enumerate().map(|(i, &s)| (s, i)).collect();

    // For each edge, draw mutations along its branch span.
    // edge gives parent/child and [left, right); branch time is
    // parent.time - child.time.
    let mut placed: Vec<(f64, usize)> = Vec::new(); // (position, child node)
    for edge in ts.edges().to_vec() {
        let span = edge.right - edge.left;
        let dt = ts.nodes()[edge.parent].time - ts.nodes()[edge.child].time;
        let expected = mutation_rate * span.max(0.0) * dt.max(0.0);
        let n_mut = rng.poisson(expected);
        for _ in 0..n_mut {
            let pos = edge.left + rng.uniform() * span;
            placed.push((pos, edge.child));
        }
    }
    placed.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    // Write the site and mutation tables and build the matrix.
    let n_rows = samples.len();
    let mut rows = vec![Vec::<u8>::new(); n_rows];
    let mut positions = Vec::new();
    for (pos, child) in placed {
        let site = ts.add_site(pos, 0);
        ts.add_mutation(TsMutation {
            site,
            node: child,
            derived_state: 1,
        })?;
        positions.push(pos);
        // Carriers: samples descended from `child` in the local tree.
        let carriers = descendants_in_local_tree(ts, pos, child)?;
        for (r, &s) in samples.iter().enumerate() {
            let carries = carriers.contains(&s);
            rows[r].push(u8::from(carries));
        }
        let _ = &sample_index; // index kept for clarity / future use
    }
    GenotypeMatrix::from_rows(rows, positions).map_err(|e| PopgenError::model(e.to_string()))
}

/// Sample-node ids descended from `node` in the local tree at
/// `position`.
fn descendants_in_local_tree(ts: &TreeSequence, position: f64, node: usize) -> Result<Vec<usize>> {
    // child -> parent under the local tree.
    let mut parent_of: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    for e in ts.edges() {
        if position >= e.left && position < e.right {
            parent_of.insert(e.child, e.parent);
        }
    }
    // A sample is a descendant of `node` if walking up reaches `node`.
    let mut out = Vec::new();
    for s in ts.samples() {
        let mut cur = s;
        loop {
            if cur == node {
                out.push(s);
                break;
            }
            match parent_of.get(&cur) {
                Some(&p) => cur = p,
                None => break,
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coalescent::kingman::{coalescent, PopHistory};

    fn labels(n: usize) -> Vec<String> {
        (0..n).map(|i| format!("L{i}")).collect()
    }

    #[test]
    fn overlay_on_tree_produces_a_matrix() {
        let tree = coalescent(&labels(10), &PopHistory::Constant(1000.0), 42).unwrap();
        // Rate scaled so a decent number of mutations land.
        let gm = overlay_on_tree(&tree, 1e-3, 7).unwrap();
        assert_eq!(gm.n_samples(), 10);
        assert!(gm.n_sites() > 0, "no mutations placed");
    }

    #[test]
    fn higher_rate_yields_more_sites() {
        let tree = coalescent(&labels(12), &PopHistory::Constant(1000.0), 1).unwrap();
        let low = overlay_on_tree(&tree, 1e-4, 3).unwrap();
        let high = overlay_on_tree(&tree, 1e-2, 3).unwrap();
        assert!(high.n_sites() > low.n_sites());
    }

    #[test]
    fn overlay_is_deterministic() {
        let tree = coalescent(&labels(8), &PopHistory::Constant(1000.0), 5).unwrap();
        let a = overlay_on_tree(&tree, 1e-3, 9).unwrap();
        let b = overlay_on_tree(&tree, 1e-3, 9).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn singleton_branch_makes_a_singleton_variant() {
        // A mutation on a leaf's terminal branch is carried by exactly
        // one sample. Over many mutations at least one column must be a
        // singleton (derived count 1).
        let tree = coalescent(&labels(15), &PopHistory::Constant(1000.0), 2).unwrap();
        let gm = overlay_on_tree(&tree, 5e-3, 4).unwrap();
        let has_singleton = (0..gm.n_sites()).any(|c| gm.derived_count(c).unwrap() == 1);
        assert!(has_singleton, "no singleton variant appeared");
    }

    #[test]
    fn overlay_on_tree_sequence_writes_tables() {
        use crate::coalescent::arg::{simulate_arg, ArgParams};
        let mut ts =
            simulate_arg(ArgParams::uniform(8, 1000.0, 3e-4, 5000.0, 11).unwrap()).unwrap();
        let gm = overlay_mutations(&mut ts, 1e-4, 13).unwrap();
        assert_eq!(gm.n_samples(), 8);
        // Site and mutation tables grew in lockstep with the matrix.
        assert_eq!(ts.site_count(), gm.n_sites());
        assert_eq!(ts.mutation_count(), gm.n_sites());
    }

    #[test]
    fn overlay_rejects_negative_rate() {
        let tree = coalescent(&labels(4), &PopHistory::Constant(1000.0), 1).unwrap();
        assert!(overlay_on_tree(&tree, -1.0, 1).is_err());
    }
}
