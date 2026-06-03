//! Tree-sequence summary statistics — windowed, branch and site modes.
//!
//! This is the `tskit` statistics-framework approach. Every summary
//! statistic over a sample has three equivalent computational modes:
//!
//! - **site mode** — count the segregating sites in the mutation
//!   table that contribute to the statistic (a single mutation at a
//!   site separating sample-sets `A` and `B` contributes `|A|·|B|` to
//!   site divergence, for instance). This is the *finite-mutation*
//!   estimate.
//! - **branch mode** — sum branch lengths in the local trees,
//!   weighted by the analogous combinatorial factor. Branch-mode is
//!   the *expected* statistic under infinite-sites mutation at rate
//!   `mu`: `site(mu) = mu · branch + O(mu^2)` (Ralph, Thornton & Kelleher
//!   2020, "Efficiently summarizing relationships in large samples").
//! - **per-window** — restrict either mode to a contiguous chromosomal
//!   window `[a, b)`, computing site and edge contributions only on
//!   their overlap with `[a, b)`.
//!
//! The functions in this module operate directly on a
//! [`crate::coalescent::TreeSequence`]. Windowed site π and divergence
//! also take a [`GenotypeMatrix`] for the sample alleles; branch
//! variants need only the tree sequence.

use crate::coalescent::tree_sequence::TreeSequence;
use crate::error::{PopgenError, Result};
use crate::infer::GenotypeMatrix;

/// Distinct genomic break-points of a tree sequence: every unique
/// edge endpoint, sorted ascending and de-duplicated, padded with
/// `0` and `sequence_length`.
fn breakpoints(ts: &TreeSequence) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::with_capacity(2 + 2 * ts.edge_count());
    out.push(0.0);
    out.push(ts.sequence_length());
    for e in ts.edges() {
        out.push(e.left);
        out.push(e.right);
    }
    out.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    out.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
    out
}

/// For position `p`, returns a parent-of array such that
/// `parent_of[c] = Some(parent)` iff there is an edge covering `p`
/// from `parent` to `c`. Length is `ts.node_count()`.
fn parent_of_at(ts: &TreeSequence, position: f64) -> Vec<Option<usize>> {
    let mut p = vec![None; ts.node_count()];
    for e in ts.edges() {
        if position >= e.left && position < e.right {
            p[e.child] = Some(e.parent);
        }
    }
    p
}

/// Walks downward from `node` under the local tree given by
/// `parent_of`, counting the sample nodes (those flagged
/// `is_sample`) reachable.
fn count_descendant_samples(
    ts: &TreeSequence,
    parent_of: &[Option<usize>],
    node: usize,
) -> usize {
    // Build a children-of view by inverting parent_of.
    let mut children_of: Vec<Vec<usize>> = vec![Vec::new(); ts.node_count()];
    for (c, &p) in parent_of.iter().enumerate() {
        if let Some(pp) = p {
            children_of[pp].push(c);
        }
    }
    let nodes = ts.nodes();
    let mut count = 0usize;
    let mut stack = vec![node];
    while let Some(n) = stack.pop() {
        if nodes[n].is_sample {
            count += 1;
        }
        for &c in &children_of[n] {
            stack.push(c);
        }
    }
    count
}

/// Windowed statistic result: one value per requested window.
#[derive(Clone, Debug, PartialEq)]
pub struct WindowedStats {
    /// Per-window values.
    pub values: Vec<f64>,
    /// Per-window `[left, right)` boundaries (length = `values.len()`).
    pub windows: Vec<(f64, f64)>,
}

/// Validates a list of `[left, right)` windows against a tree
/// sequence's chromosome length: each must lie in `[0, L]`, be
/// strictly positive, and the list must be sorted and non-overlapping.
fn validate_windows(ts: &TreeSequence, windows: &[(f64, f64)]) -> Result<()> {
    let l = ts.sequence_length();
    if windows.is_empty() {
        return Err(PopgenError::invalid(
            "windows",
            "need at least one window",
        ));
    }
    let mut last_right = 0.0;
    for (i, &(a, b)) in windows.iter().enumerate() {
        if !a.is_finite() || !b.is_finite() {
            return Err(PopgenError::invalid(
                "windows",
                "window endpoints must be finite",
            ));
        }
        if a < 0.0 - 1e-12 || b > l + 1e-9 {
            return Err(PopgenError::invalid(
                "windows",
                format!(
                    "window {i} [{a}, {b}) is outside the chromosome [0, {l})"
                ),
            ));
        }
        if b <= a {
            return Err(PopgenError::invalid(
                "windows",
                format!("window {i} is empty"),
            ));
        }
        if i > 0 && a < last_right - 1e-9 {
            return Err(PopgenError::invalid(
                "windows",
                "windows must be sorted and non-overlapping",
            ));
        }
        last_right = b;
    }
    Ok(())
}

/// Per-site contribution to nucleotide diversity π:
/// `2 * d * (n - d) / (n * (n - 1))`.
fn site_pi_contrib(d: usize, n: usize) -> f64 {
    if n < 2 {
        return 0.0;
    }
    let d = d as f64;
    let nn = n as f64;
    2.0 * d * (nn - d) / (nn * (nn - 1.0))
}

/// **Site-mode** windowed nucleotide diversity π: for each window,
/// sums per-site π contributions over the sites whose `positions`
/// fall in `[left, right)`. The result is *per-site* (the mean
/// pairwise-difference per site over the window).
///
/// `positions` are taken from `matrix.positions()` and must match the
/// tree sequence's mutation positions in order (the genotype matrix
/// produced by [`crate::coalescent::overlay_mutations`] satisfies this
/// by construction).
///
/// # Errors
/// [`PopgenError::Invalid`] on empty / out-of-range windows.
pub fn windowed_site_diversity(
    ts: &TreeSequence,
    matrix: &GenotypeMatrix,
    windows: &[(f64, f64)],
) -> Result<WindowedStats> {
    validate_windows(ts, windows)?;
    let n = matrix.n_samples();
    let mut values = vec![0.0; windows.len()];
    for (col, &pos) in matrix.positions().iter().enumerate() {
        for (w_idx, &(a, b)) in windows.iter().enumerate() {
            if pos >= a && pos < b {
                let d = matrix.derived_count(col)?;
                values[w_idx] += site_pi_contrib(d, n);
                break;
            }
        }
    }
    // Per-site π is the running sum divided by window width if the
    // caller wants a per-bp rate. We expose the raw "summed
    // contribution over sites in window" — equivalent to S-weighted
    // π — which is the natural quantity for site/branch comparison
    // (msprime returns the same).
    Ok(WindowedStats {
        values,
        windows: windows.to_vec(),
    })
}

/// **Site-mode** windowed segregating-site count: the number of
/// sites with derived count strictly between `0` and `n`.
///
/// # Errors
/// [`PopgenError::Invalid`] on empty / out-of-range windows.
pub fn windowed_segregating_sites(
    ts: &TreeSequence,
    matrix: &GenotypeMatrix,
    windows: &[(f64, f64)],
) -> Result<WindowedStats> {
    validate_windows(ts, windows)?;
    let n = matrix.n_samples();
    let mut values = vec![0.0; windows.len()];
    for (col, &pos) in matrix.positions().iter().enumerate() {
        let d = matrix.derived_count(col)?;
        if d == 0 || d == n {
            continue;
        }
        for (w_idx, &(a, b)) in windows.iter().enumerate() {
            if pos >= a && pos < b {
                values[w_idx] += 1.0;
                break;
            }
        }
    }
    Ok(WindowedStats {
        values,
        windows: windows.to_vec(),
    })
}

/// **Site-mode** divergence (`d_xy`) between two sample sets over
/// windows. At each site contributes `count_a_derived *
/// (size_b - count_b_derived) + count_b_derived *
/// (size_a - count_a_derived)` divided by `size_a * size_b` — the
/// probability a randomly-chosen pair `(a, b)` with `a ∈ A, b ∈ B`
/// differs at the site.
///
/// # Errors
/// [`PopgenError::Invalid`] on empty windows or empty sample sets.
pub fn windowed_site_divergence(
    ts: &TreeSequence,
    matrix: &GenotypeMatrix,
    set_a: &[usize],
    set_b: &[usize],
    windows: &[(f64, f64)],
) -> Result<WindowedStats> {
    validate_windows(ts, windows)?;
    if set_a.is_empty() || set_b.is_empty() {
        return Err(PopgenError::invalid(
            "sample_set",
            "both sample sets must be non-empty",
        ));
    }
    let na = set_a.len() as f64;
    let nb = set_b.len() as f64;
    let mut values = vec![0.0; windows.len()];
    for (col, &pos) in matrix.positions().iter().enumerate() {
        let mut da = 0usize;
        for &s in set_a {
            if s < matrix.n_samples() {
                da += matrix.get(s, col) as usize;
            }
        }
        let mut db = 0usize;
        for &s in set_b {
            if s < matrix.n_samples() {
                db += matrix.get(s, col) as usize;
            }
        }
        let dxy = (da as f64 * (nb - db as f64)
            + db as f64 * (na - da as f64))
            / (na * nb);
        for (w_idx, &(a, b)) in windows.iter().enumerate() {
            if pos >= a && pos < b {
                values[w_idx] += dxy;
                break;
            }
        }
    }
    Ok(WindowedStats {
        values,
        windows: windows.to_vec(),
    })
}

/// **Branch-mode** total nucleotide diversity π over the whole
/// chromosome. The contribution of every edge `e = (parent, child,
/// left, right)` over each tree interval it spans is `branch_length ·
/// width · 2·k·(n - k) / (n·(n-1))` where `k` is the number of sample
/// descendants of `child` in the local tree at that interval and `n`
/// is the total sample count.
///
/// Under infinite-sites mutation at per-base-pair per-generation rate
/// `mu`, the expected site-mode π equals `mu · branch_pi_total`
/// (Ralph, Thornton & Kelleher 2020).
///
/// # Errors
/// [`PopgenError::Invalid`] if the tree sequence has fewer than two
/// samples.
pub fn branch_diversity(ts: &TreeSequence) -> Result<f64> {
    let n = ts.sample_count();
    if n < 2 {
        return Err(PopgenError::invalid(
            "ts",
            "branch diversity needs at least two samples",
        ));
    }
    let bps = breakpoints(ts);
    let nodes = ts.nodes();
    let mut total = 0.0f64;
    for w in bps.windows(2) {
        let (lo, hi) = (w[0], w[1]);
        if hi - lo <= 1e-12 {
            continue;
        }
        let midpoint = 0.5 * (lo + hi);
        let parent_of = parent_of_at(ts, midpoint);
        for e in ts.edges() {
            if midpoint >= e.left && midpoint < e.right {
                let dt = (nodes[e.parent].time - nodes[e.child].time).max(0.0);
                let k = count_descendant_samples(ts, &parent_of, e.child);
                let factor = site_pi_contrib(k, n);
                total += (hi - lo) * dt * factor;
            }
        }
    }
    Ok(total)
}

/// **Branch-mode** windowed π: same as [`branch_diversity`] but
/// restricted to each window in turn.
///
/// # Errors
/// [`PopgenError::Invalid`] on empty / out-of-range windows or too
/// few samples.
pub fn windowed_branch_diversity(
    ts: &TreeSequence,
    windows: &[(f64, f64)],
) -> Result<WindowedStats> {
    validate_windows(ts, windows)?;
    let n = ts.sample_count();
    if n < 2 {
        return Err(PopgenError::invalid(
            "ts",
            "branch diversity needs at least two samples",
        ));
    }
    let bps = breakpoints(ts);
    let nodes = ts.nodes();
    let mut values = vec![0.0; windows.len()];
    for w in bps.windows(2) {
        let (lo, hi) = (w[0], w[1]);
        if hi - lo <= 1e-12 {
            continue;
        }
        let midpoint = 0.5 * (lo + hi);
        let parent_of = parent_of_at(ts, midpoint);
        // Pre-compute descendant counts for every (edge, this
        // sub-interval) pair we'll need below.
        let mut edge_contrib_per_unit = Vec::with_capacity(ts.edge_count());
        for e in ts.edges() {
            if midpoint >= e.left && midpoint < e.right {
                let dt = (nodes[e.parent].time - nodes[e.child].time).max(0.0);
                let k = count_descendant_samples(ts, &parent_of, e.child);
                let factor = site_pi_contrib(k, n);
                edge_contrib_per_unit.push(dt * factor);
            } else {
                edge_contrib_per_unit.push(0.0);
            }
        }
        let total_per_unit: f64 = edge_contrib_per_unit.iter().sum();
        // Allocate by overlap with each window.
        for (w_idx, &(a, b)) in windows.iter().enumerate() {
            let overlap_lo = lo.max(a);
            let overlap_hi = hi.min(b);
            if overlap_hi > overlap_lo {
                values[w_idx] += (overlap_hi - overlap_lo) * total_per_unit;
            }
        }
    }
    Ok(WindowedStats {
        values,
        windows: windows.to_vec(),
    })
}

/// **Branch-mode** total divergence (`d_xy`) between two sample sets.
/// At each tree interval and edge, the contribution is `dt · width ·
/// (k_a · (n_b - k_b) + k_b · (n_a - k_a)) / (n_a · n_b)` where
/// `k_a` (respectively `k_b`) counts descendants of the edge's child
/// in sample set `A` (respectively `B`).
///
/// # Errors
/// [`PopgenError::Invalid`] on empty sample sets or out-of-range
/// indices.
pub fn branch_divergence(
    ts: &TreeSequence,
    set_a: &[usize],
    set_b: &[usize],
) -> Result<f64> {
    if set_a.is_empty() || set_b.is_empty() {
        return Err(PopgenError::invalid(
            "sample_set",
            "both sample sets must be non-empty",
        ));
    }
    let nodes = ts.nodes();
    for &s in set_a.iter().chain(set_b.iter()) {
        if s >= nodes.len() || !nodes[s].is_sample {
            return Err(PopgenError::invalid(
                "sample_set",
                format!("node {s} is not a sample"),
            ));
        }
    }
    let na = set_a.len() as f64;
    let nb = set_b.len() as f64;
    let bps = breakpoints(ts);
    let mut total = 0.0f64;
    for w in bps.windows(2) {
        let (lo, hi) = (w[0], w[1]);
        if hi - lo <= 1e-12 {
            continue;
        }
        let midpoint = 0.5 * (lo + hi);
        let parent_of = parent_of_at(ts, midpoint);
        for e in ts.edges() {
            if midpoint >= e.left && midpoint < e.right {
                let dt = (nodes[e.parent].time - nodes[e.child].time).max(0.0);
                let (ka, kb) =
                    count_set_descendants(ts, &parent_of, e.child, set_a, set_b);
                let dxy = (ka * (nb - kb) + kb * (na - ka)) / (na * nb);
                total += (hi - lo) * dt * dxy;
            }
        }
    }
    Ok(total)
}

/// Number of descendants of `node` in the local tree given by
/// `parent_of`, restricted to two given sample sets — returned as
/// `(count_in_a, count_in_b)`.
fn count_set_descendants(
    ts: &TreeSequence,
    parent_of: &[Option<usize>],
    node: usize,
    set_a: &[usize],
    set_b: &[usize],
) -> (f64, f64) {
    let mut children_of: Vec<Vec<usize>> = vec![Vec::new(); ts.node_count()];
    for (c, &p) in parent_of.iter().enumerate() {
        if let Some(pp) = p {
            children_of[pp].push(c);
        }
    }
    let set_a_set: std::collections::HashSet<usize> = set_a.iter().copied().collect();
    let set_b_set: std::collections::HashSet<usize> = set_b.iter().copied().collect();
    let mut ka = 0.0;
    let mut kb = 0.0;
    let mut stack = vec![node];
    while let Some(n) = stack.pop() {
        if set_a_set.contains(&n) {
            ka += 1.0;
        }
        if set_b_set.contains(&n) {
            kb += 1.0;
        }
        for &c in &children_of[n] {
            stack.push(c);
        }
    }
    (ka, kb)
}

/// Returns evenly-spaced windows partitioning the chromosome into
/// `n_windows` consecutive segments of equal width. Convenience.
pub fn equal_windows(ts: &TreeSequence, n_windows: usize) -> Vec<(f64, f64)> {
    let l = ts.sequence_length();
    let w = l / n_windows.max(1) as f64;
    (0..n_windows.max(1))
        .map(|i| {
            let a = i as f64 * w;
            let b = if i + 1 == n_windows.max(1) {
                l
            } else {
                (i + 1) as f64 * w
            };
            (a, b)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coalescent::arg::{simulate_arg, ArgParams};
    use crate::coalescent::overlay::overlay_mutations;
    use crate::coalescent::tree_sequence::Edge;

    fn build_tiny_ts() -> TreeSequence {
        // 3 samples, single tree, no recombination. a, b coalesce at
        // ab (time 1); ab coalesces with c at root (time 2).
        let mut ts = TreeSequence::new(100.0).unwrap();
        let _a = ts.add_node(0.0, true); // 0
        let _b = ts.add_node(0.0, true); // 1
        let _c = ts.add_node(0.0, true); // 2
        let ab = ts.add_node(1.0, false); // 3
        let root = ts.add_node(2.0, false); // 4
        for &(p, c) in &[(ab, 0), (ab, 1), (root, ab), (root, 2)] {
            ts.add_edge(Edge { parent: p, child: c, left: 0.0, right: 100.0 })
                .unwrap();
        }
        ts.finalize().unwrap();
        ts
    }

    #[test]
    fn branch_diversity_on_a_balanced_tiny_tree() {
        // Expected branch π:
        // tree intervals = single [0, 100).
        // Edges:
        //   ab->a, ab->b: each width 100, dt 1, k=1 → 2*1*2/(3*2)=2/3.
        //     contribution = 100 * 1 * 2/3 = 200/3 each; sum 400/3.
        //   root->ab: width 100, dt 1, k=2 (a,b under ab) →
        //     2*2*1/(3*2)=2/3. contribution = 100 * 1 * 2/3 = 200/3.
        //   root->c: width 100, dt 2, k=1 (c) → 2/3.
        //     contribution = 100 * 2 * 2/3 = 400/3.
        // Total = 400/3 + 200/3 + 400/3 = 1000/3.
        let ts = build_tiny_ts();
        let pi = branch_diversity(&ts).unwrap();
        assert!(
            (pi - 1000.0 / 3.0).abs() < 1e-9,
            "got {pi} vs 333.33...",
        );
    }

    #[test]
    fn windowed_branch_diversity_partitions_total() {
        let ts = build_tiny_ts();
        let total = branch_diversity(&ts).unwrap();
        let windows = equal_windows(&ts, 4);
        let ws = windowed_branch_diversity(&ts, &windows).unwrap();
        let sum: f64 = ws.values.iter().sum();
        assert!(
            (sum - total).abs() < 1e-9,
            "windowed sum {sum} vs total {total}"
        );
        // Each window has equal width and a flat tree → equal value.
        assert!(
            (ws.values[0] - ws.values[3]).abs() < 1e-9,
            "non-uniform values for uniform tree: {:?}",
            ws.values
        );
    }

    #[test]
    fn windowed_branch_diversity_rejects_bad_windows() {
        let ts = build_tiny_ts();
        // Empty list.
        assert!(windowed_branch_diversity(&ts, &[]).is_err());
        // Out of range.
        assert!(windowed_branch_diversity(&ts, &[(0.0, 200.0)]).is_err());
        // Empty interval.
        assert!(windowed_branch_diversity(&ts, &[(50.0, 50.0)]).is_err());
        // Overlapping.
        assert!(windowed_branch_diversity(&ts, &[(0.0, 60.0), (40.0, 80.0)]).is_err());
    }

    #[test]
    fn branch_divergence_on_tiny_tree() {
        // A={0}, B={2}. Walking from each edge's child:
        // ab->a: ka=1, kb=0 → dxy = (1*1 + 0*0) / 1 = 1; width 100, dt 1.
        //   contribution = 100.
        // ab->b: ka=0, kb=0 → 0.
        // root->ab: ka=1, kb=0 → 1; width 100, dt 1 → 100.
        // root->c: ka=0, kb=1 → 1; width 100, dt 2 → 200.
        // Total = 400.
        let ts = build_tiny_ts();
        let dxy = branch_divergence(&ts, &[0], &[2]).unwrap();
        assert!((dxy - 400.0).abs() < 1e-9, "got {dxy}");
    }

    #[test]
    fn branch_divergence_rejects_non_samples() {
        let ts = build_tiny_ts();
        assert!(branch_divergence(&ts, &[3], &[2]).is_err());
        assert!(branch_divergence(&ts, &[], &[2]).is_err());
    }

    /// On an ARG with mutation rate `mu`, `E[site_pi] ≈ mu *
    /// branch_pi`. With 0 recombination and a single big sample, the
    /// ratio should be roughly mu over many seeds.
    #[test]
    fn site_pi_tracks_branch_pi_via_mu() {
        let mu = 5e-3;
        let mut acc_site = 0.0;
        let mut acc_branch = 0.0;
        let reps = 20;
        for seed in 0..reps {
            let mut ts = simulate_arg(
                ArgParams::uniform(12, 1000.0, 0.0, 1000.0, seed).unwrap(),
            )
            .unwrap();
            let bpi = branch_diversity(&ts).unwrap();
            let gm = overlay_mutations(&mut ts, mu, seed + 1000).unwrap();
            // site-mode total π = sum over sites of 2 d (n-d) / (n(n-1)).
            let n = gm.n_samples();
            let site_pi: f64 = (0..gm.n_sites())
                .map(|c| site_pi_contrib(gm.derived_count(c).unwrap(), n))
                .sum();
            acc_site += site_pi;
            acc_branch += bpi;
        }
        let mean_site = acc_site / reps as f64;
        let mean_branch = acc_branch / reps as f64;
        let expected = mu * mean_branch;
        let ratio = mean_site / expected;
        // Tolerance is generous (Poisson noise over 20 reps).
        assert!(
            ratio > 0.5 && ratio < 1.7,
            "mean site π {mean_site} vs expected {expected} (ratio {ratio})",
        );
    }

    #[test]
    fn site_windowed_diversity_partitions_total() {
        let mut ts = simulate_arg(
            ArgParams::uniform(10, 1000.0, 1e-4, 1000.0, 7).unwrap(),
        )
        .unwrap();
        let gm = overlay_mutations(&mut ts, 1e-3, 13).unwrap();
        // Total site π via the matrix.
        let n = gm.n_samples();
        let total: f64 = (0..gm.n_sites())
            .map(|c| site_pi_contrib(gm.derived_count(c).unwrap(), n))
            .sum();
        let windows = equal_windows(&ts, 5);
        let ws = windowed_site_diversity(&ts, &gm, &windows).unwrap();
        let sum: f64 = ws.values.iter().sum();
        assert!((sum - total).abs() < 1e-9, "windowed sum {sum} vs total {total}");
    }

    #[test]
    fn site_segregating_count_matches_matrix() {
        let mut ts = simulate_arg(
            ArgParams::uniform(8, 1000.0, 1e-4, 500.0, 3).unwrap(),
        )
        .unwrap();
        let gm = overlay_mutations(&mut ts, 1e-3, 4).unwrap();
        let windows = equal_windows(&ts, 3);
        let ws = windowed_segregating_sites(&ts, &gm, &windows).unwrap();
        let sum: f64 = ws.values.iter().sum();
        assert!(
            (sum - gm.segregating_sites() as f64).abs() < 1e-9,
            "windowed S sum {sum} vs matrix S {}",
            gm.segregating_sites(),
        );
    }

    #[test]
    fn site_divergence_runs_and_is_non_negative() {
        let mut ts = simulate_arg(
            ArgParams::uniform(10, 1000.0, 1e-4, 500.0, 5).unwrap(),
        )
        .unwrap();
        let gm = overlay_mutations(&mut ts, 1e-3, 6).unwrap();
        let set_a: Vec<usize> = (0..5).collect();
        let set_b: Vec<usize> = (5..10).collect();
        let windows = equal_windows(&ts, 4);
        let ws =
            windowed_site_divergence(&ts, &gm, &set_a, &set_b, &windows).unwrap();
        for &v in &ws.values {
            assert!(v >= 0.0, "divergence negative? {v}");
        }
    }

    /// Branch divergence equals site divergence in expectation:
    /// `E[site_dxy] ≈ mu * branch_dxy`.
    #[test]
    fn site_and_branch_divergence_are_proportional_to_mu() {
        let mu = 5e-3;
        let mut acc_site = 0.0;
        let mut acc_branch = 0.0;
        let reps = 20;
        let set_a: Vec<usize> = (0..6).collect();
        let set_b: Vec<usize> = (6..12).collect();
        for seed in 0..reps {
            let mut ts = simulate_arg(
                ArgParams::uniform(12, 1000.0, 0.0, 1000.0, seed).unwrap(),
            )
            .unwrap();
            let bdxy = branch_divergence(&ts, &set_a, &set_b).unwrap();
            let gm = overlay_mutations(&mut ts, mu, seed + 1000).unwrap();
            // Site-mode total divergence via the matrix.
            let na = set_a.len() as f64;
            let nb = set_b.len() as f64;
            let mut total = 0.0;
            for col in 0..gm.n_sites() {
                let mut da = 0usize;
                for &s in &set_a {
                    da += gm.get(s, col) as usize;
                }
                let mut db = 0usize;
                for &s in &set_b {
                    db += gm.get(s, col) as usize;
                }
                total += (da as f64 * (nb - db as f64)
                    + db as f64 * (na - da as f64))
                    / (na * nb);
            }
            acc_site += total;
            acc_branch += bdxy;
        }
        let mean_site = acc_site / reps as f64;
        let mean_branch = acc_branch / reps as f64;
        let expected = mu * mean_branch;
        let ratio = mean_site / expected;
        assert!(
            ratio > 0.5 && ratio < 1.7,
            "mean site dxy {mean_site} vs expected {expected} (ratio {ratio})",
        );
    }

    #[test]
    fn equal_windows_tile_the_chromosome() {
        let ts = build_tiny_ts();
        let ws = equal_windows(&ts, 5);
        assert_eq!(ws.len(), 5);
        assert!(ws[0].0.abs() < 1e-9);
        assert!((ws[4].1 - ts.sequence_length()).abs() < 1e-9);
        for w in ws.windows(2) {
            assert!((w[0].1 - w[1].0).abs() < 1e-9);
        }
    }
}
