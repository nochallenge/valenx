//! The coalescent with recombination — Hudson's (1983) ARG.
//!
//! ## Hudson-canonical algorithm
//!
//! Each *lineage* is represented as a sorted list of segments
//! `(left, right, node)` — the lineage carries ancestral material over
//! `[left, right)`, and the "current" node id in the tree-sequence
//! tables at that genomic stretch is `node`. Initially every sample
//! lineage is a single segment over the whole chromosome, labelled with
//! the sample node id.
//!
//! Two kinds of event occur, the Hudson 1983 algorithm:
//!
//! - **Coalescence** — two lineages `x` and `y` merge into a single
//!   ancestor lineage `z`. Over every genomic position where *both*
//!   carry material, a new internal node `p` is created and edges
//!   `p → x.node`, `p → y.node` are written over the overlapping
//!   interval; the resulting segment in `z` carries node label `p`.
//!   Over stretches where only `x` (or only `y`) had material, the
//!   segment is carried into `z` unchanged — no edge is written,
//!   because no coalescence event actually happened there. (This is
//!   the canonical msprime / `tskit` simplification — it keeps the
//!   edge table sparse and avoids the spurious unary edges a naive
//!   implementation produces.)
//! - **Recombination** — a single lineage splits at a breakpoint into
//!   two: one carrying the segments left of the breakpoint, one the
//!   segments right of it. No new node is created.
//!
//! Each coalescence segment in `z` is *retired* if no other active
//! lineage still carries that genomic stretch — its grand MRCA has
//! been reached. The process stops when every genomic position has
//! retired (only the trivial all-empty lineage(s) remain).
//!
//! Recombination rate is provided as a piecewise-constant
//! [`RecombinationMap`] — the standard `msprime` / SLiM data structure
//! that lets hot-spot and cold-spot regions be modelled along the
//! chromosome. A flat-rate ARG is the [`RecombinationMap::uniform`]
//! special case.

use crate::coalescent::tree_sequence::{Edge, TreeSequence};
use crate::error::{PopgenError, Result};
use crate::rng::Rng;

/// A piecewise-constant recombination-rate map.
///
/// The chromosome `[0, L)` is partitioned into contiguous windows
/// `[boundaries[i], boundaries[i+1])` carrying per-base-pair
/// per-generation rate `rates[i]`. The final window extends to
/// `L = sequence_length`. A uniform map is the `rates.len() == 1`
/// case, available via [`RecombinationMap::uniform`].
///
/// The map's *cumulative* rate `M(x) = integral_0^x r(u) du` is what
/// the simulator uses: per-window rates are stored together with
/// prefix sums so that integrating over any sub-interval is
/// `O(log windows)` and drawing a breakpoint from the rate is an
/// inverse-CDF lookup.
#[derive(Clone, Debug, PartialEq)]
pub struct RecombinationMap {
    boundaries: Vec<f64>,
    rates: Vec<f64>,
    cum: Vec<f64>,
}

impl RecombinationMap {
    /// Builds a uniform-rate map over `[0, length)`.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] on a non-positive length or a negative
    /// rate.
    pub fn uniform(length: f64, rate: f64) -> Result<Self> {
        if length <= 0.0 || !length.is_finite() {
            return Err(PopgenError::invalid(
                "length",
                "must be finite and positive",
            ));
        }
        if rate < 0.0 || !rate.is_finite() {
            return Err(PopgenError::invalid(
                "rate",
                "must be finite and non-negative",
            ));
        }
        Ok(RecombinationMap {
            boundaries: vec![0.0, length],
            rates: vec![rate],
            cum: vec![0.0, rate * length],
        })
    }

    /// Builds a piecewise-constant map from a vector of boundaries (in
    /// ascending order, starting at `0`) and one rate per window.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] on fewer than two boundaries,
    /// non-ascending boundaries, a leading boundary not at `0`, or a
    /// negative / non-finite rate;
    /// [`PopgenError::Dimension`] if the number of rates is not
    /// `boundaries.len() - 1`.
    pub fn piecewise(boundaries: Vec<f64>, rates: Vec<f64>) -> Result<Self> {
        if boundaries.len() < 2 {
            return Err(PopgenError::invalid(
                "boundaries",
                "need at least one window (two boundaries)",
            ));
        }
        if boundaries[0] != 0.0 {
            return Err(PopgenError::invalid(
                "boundaries",
                "first boundary must be 0",
            ));
        }
        if rates.len() + 1 != boundaries.len() {
            return Err(PopgenError::dimension(
                boundaries.len() - 1,
                rates.len(),
                "rate-map windows",
            ));
        }
        for w in boundaries.windows(2) {
            if w[1] <= w[0] {
                return Err(PopgenError::invalid(
                    "boundaries",
                    "must be strictly ascending",
                ));
            }
        }
        for &r in &rates {
            if r < 0.0 || !r.is_finite() {
                return Err(PopgenError::invalid(
                    "rates",
                    "every rate must be finite and non-negative",
                ));
            }
        }
        let mut cum = Vec::with_capacity(boundaries.len());
        cum.push(0.0);
        for i in 0..rates.len() {
            let span = boundaries[i + 1] - boundaries[i];
            cum.push(cum[i] + rates[i] * span);
        }
        Ok(RecombinationMap {
            boundaries,
            rates,
            cum,
        })
    }

    /// Total chromosome length.
    pub fn sequence_length(&self) -> f64 {
        *self.boundaries.last().expect("validated non-empty")
    }

    /// Cumulative rate from `0` to `x`. Clamps `x` to `[0, L]`.
    pub fn cumulative(&self, x: f64) -> f64 {
        let l = self.sequence_length();
        let x = x.clamp(0.0, l);
        let i = self.window_of(x);
        self.cum[i] + self.rates[i] * (x - self.boundaries[i])
    }

    /// Integrated rate over `[left, right)`. Returns `0` if the
    /// interval is empty.
    pub fn integrate(&self, left: f64, right: f64) -> f64 {
        if right <= left {
            return 0.0;
        }
        self.cumulative(right) - self.cumulative(left)
    }

    /// Inverts the CDF: returns the genomic position `x` for which
    /// `cumulative(x) == target`. Used to sample a breakpoint from the
    /// rate-map's intensity. `target` is clamped to `[0, total_rate]`.
    pub fn position_at_cumulative(&self, target: f64) -> f64 {
        let total = self.total_rate();
        let target = target.clamp(0.0, total);
        // Find i with cum[i] <= target <= cum[i+1].
        let mut lo = 0usize;
        let mut hi = self.rates.len();
        while lo + 1 < hi {
            let mid = (lo + hi) / 2;
            if self.cum[mid] <= target {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let i = lo;
        let remaining = target - self.cum[i];
        let r = self.rates[i];
        if r > 0.0 {
            self.boundaries[i] + remaining / r
        } else {
            // Zero-rate window: target is at its left boundary.
            self.boundaries[i]
        }
    }

    /// Total integrated rate over the whole chromosome.
    pub fn total_rate(&self) -> f64 {
        *self.cum.last().expect("validated non-empty")
    }

    /// Mean per-base-pair rate (`total_rate / sequence_length`).
    pub fn mean_rate(&self) -> f64 {
        self.total_rate() / self.sequence_length()
    }

    /// Window boundaries.
    pub fn boundaries(&self) -> &[f64] {
        &self.boundaries
    }

    /// Per-window rates.
    pub fn rates(&self) -> &[f64] {
        &self.rates
    }

    /// Index of the window containing `x` (clamped).
    fn window_of(&self, x: f64) -> usize {
        let mut lo = 0usize;
        let mut hi = self.boundaries.len() - 1;
        while lo + 1 < hi {
            let mid = (lo + hi) / 2;
            if self.boundaries[mid] <= x {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        lo
    }
}

/// Parameters for an ARG simulation.
#[derive(Clone, Debug)]
pub struct ArgParams {
    /// Number of sampled chromosomes.
    pub sample_size: usize,
    /// Diploid effective population size.
    pub effective_size: f64,
    /// Recombination-rate map along the chromosome.
    pub recombination_map: RecombinationMap,
    /// RNG seed.
    pub seed: u64,
}

impl ArgParams {
    /// Convenience constructor for a uniform-rate run.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] on bad parameters.
    pub fn uniform(
        sample_size: usize,
        effective_size: f64,
        recombination_rate: f64,
        sequence_length: f64,
        seed: u64,
    ) -> Result<Self> {
        let map = RecombinationMap::uniform(sequence_length, recombination_rate)?;
        Ok(ArgParams {
            sample_size,
            effective_size,
            recombination_map: map,
            seed,
        })
    }

    /// Validates the parameters.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] on a sample below 2 or a non-positive
    /// effective size; the map's own validation is enforced at
    /// construction.
    pub fn validate(&self) -> Result<()> {
        if self.sample_size < 2 {
            return Err(PopgenError::invalid(
                "sample_size",
                "need at least two chromosomes",
            ));
        }
        if self.effective_size <= 0.0 || !self.effective_size.is_finite() {
            return Err(PopgenError::invalid(
                "effective_size",
                "must be finite and positive",
            ));
        }
        Ok(())
    }

    /// Chromosome length, from the map.
    pub fn sequence_length(&self) -> f64 {
        self.recombination_map.sequence_length()
    }
}

/// A single segment of a lineage: `[left, right)` carries the ancestral
/// material of the sample, and at that genomic stretch the current
/// tree-sequence node id for this lineage is `node`.
#[derive(Copy, Clone, Debug, PartialEq)]
struct Segment {
    left: f64,
    right: f64,
    node: usize,
}

/// A lineage is a list of segments, kept sorted by `left` and
/// non-overlapping. A lineage's total ancestral span is the sum of
/// segment widths. An empty lineage carries no material.
type Lineage = Vec<Segment>;

/// Total ancestral-material span of a lineage. Available to tests and
/// debug-introspection callers.
#[cfg(test)]
fn span(lin: &Lineage) -> f64 {
    lin.iter().map(|s| s.right - s.left).sum()
}

fn rate_integral(lin: &Lineage, map: &RecombinationMap) -> f64 {
    lin.iter().map(|s| map.integrate(s.left, s.right)).sum()
}

/// Splits the lineage at `breakpoint` into `(left_lineage,
/// right_lineage)`. A segment straddling the breakpoint is split into
/// two segments with the same node id.
fn split_lineage(lin: &Lineage, breakpoint: f64) -> (Lineage, Lineage) {
    let mut left = Vec::new();
    let mut right = Vec::new();
    for &s in lin {
        if s.right <= breakpoint {
            left.push(s);
        } else if s.left >= breakpoint {
            right.push(s);
        } else {
            left.push(Segment {
                left: s.left,
                right: breakpoint,
                node: s.node,
            });
            right.push(Segment {
                left: breakpoint,
                right: s.right,
                node: s.node,
            });
        }
    }
    (left, right)
}

/// Simulates the coalescent with recombination and returns the ARG as
/// a [`TreeSequence`].
///
/// The sample chromosomes are the tree sequence's sample nodes (ids
/// `0..sample_size`). Recombination produces multiple local trees;
/// query them with [`TreeSequence::local_tree`].
///
/// # Errors
/// [`PopgenError::Invalid`] on bad parameters; [`PopgenError::Model`]
/// if the ARG fails to reach a grand MRCA within the internal event
/// cap, or on tree-sequence assembly failure.
pub fn simulate_arg(params: ArgParams) -> Result<TreeSequence> {
    params.validate()?;
    let mut rng = Rng::new(params.seed);
    let l = params.sequence_length();
    let map = &params.recombination_map;
    let mut ts = TreeSequence::new(l)?;

    // Active lineages.
    let mut active: Vec<Lineage> = Vec::with_capacity(params.sample_size);
    for _ in 0..params.sample_size {
        let id = ts.add_node(0.0, true);
        active.push(vec![Segment {
            left: 0.0,
            right: l,
            node: id,
        }]);
    }

    let mut time = 0.0;
    let cap = 5_000_000usize;
    let mut events = 0usize;

    // Loop until at most one lineage remains AND nothing has any
    // un-MRCA'd material (single-lineage state means the sample is
    // fully coalesced).
    while active.len() > 1 {
        events += 1;
        if events > cap {
            return Err(PopgenError::model(
                "ARG simulation exceeded its event cap",
            ));
        }
        let k = active.len() as f64;
        let coal_rate = k * (k - 1.0) / 2.0 / params.effective_size;
        // Per-lineage recombination intensity.
        let lineage_rate: Vec<f64> = active
            .iter()
            .map(|lin| rate_integral(lin, map))
            .collect();
        let recomb_rate: f64 = lineage_rate.iter().sum();
        let total_rate = coal_rate + recomb_rate;
        if total_rate <= 0.0 {
            // No further events can fire. With k > 1 still active and
            // no rate, the configuration is stuck — but since
            // coal_rate is always positive for k >= 2 and a finite N
            // (validated), this branch is unreachable in practice.
            break;
        }
        time += rng.exponential(total_rate);

        if rng.uniform() * total_rate < coal_rate {
            coalesce(&mut active, &mut ts, &mut rng, time)?;
        } else {
            recombine(&mut active, &mut rng, &lineage_rate, map);
        }
    }

    ts.finalize()?;
    Ok(ts)
}

/// Merges two active lineages into one, creating coalescent nodes and
/// edges only over genomic stretches where both carry material.
fn coalesce(
    active: &mut Vec<Lineage>,
    ts: &mut TreeSequence,
    rng: &mut Rng,
    time: f64,
) -> Result<()> {
    let n = active.len();
    let i = rng.below(n);
    let mut j = rng.below(n - 1);
    if j >= i {
        j += 1;
    }
    let (lo, hi) = if i < j { (i, j) } else { (j, i) };
    let y = active.remove(hi);
    let x = active.remove(lo);

    // Merge x and y's segment lists by genomic position. Over an
    // overlap, create one new internal node and emit two edges.
    let mut merged: Lineage = Vec::new();
    let mut ix = 0usize;
    let mut iy = 0usize;
    // We need to advance through both sorted segment lists, splitting
    // segments at the union of their boundary points.
    let mut breakpoints = collect_breakpoints(&x, &y);
    // Iterate adjacent breakpoints as intervals [b_i, b_{i+1}).
    let mut sx = x.iter().copied().peekable();
    let mut sy = y.iter().copied().peekable();
    let _ = (&mut sx, &mut sy, &mut ix, &mut iy);
    breakpoints.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
    // Pre-create per-interval node references (segments) from x and y.
    let x_at = |pos: f64| -> Option<usize> {
        x.iter().find(|s| s.left <= pos && pos < s.right).map(|s| s.node)
    };
    let y_at = |pos: f64| -> Option<usize> {
        y.iter().find(|s| s.left <= pos && pos < s.right).map(|s| s.node)
    };

    // Process each contiguous sub-interval.
    for w in breakpoints.windows(2) {
        let (lo_bp, hi_bp) = (w[0], w[1]);
        if hi_bp - lo_bp < 1e-12 {
            continue;
        }
        let mid = 0.5 * (lo_bp + hi_bp);
        let xn = x_at(mid);
        let yn = y_at(mid);
        match (xn, yn) {
            (Some(xnode), Some(ynode)) => {
                // Overlap: create a new internal node and write edges.
                let parent = ts.add_node(time, false);
                // Edges only if the child differs from the parent
                // (children are samples or earlier ancestors; here we
                // always have distinct labels because parent is new).
                ts.add_edge(Edge {
                    parent,
                    child: xnode,
                    left: lo_bp,
                    right: hi_bp,
                })?;
                ts.add_edge(Edge {
                    parent,
                    child: ynode,
                    left: lo_bp,
                    right: hi_bp,
                })?;
                // The merged lineage's segment carries the new node.
                merged.push(Segment {
                    left: lo_bp,
                    right: hi_bp,
                    node: parent,
                });
            }
            (Some(xnode), None) => {
                merged.push(Segment {
                    left: lo_bp,
                    right: hi_bp,
                    node: xnode,
                });
            }
            (None, Some(ynode)) => {
                merged.push(Segment {
                    left: lo_bp,
                    right: hi_bp,
                    node: ynode,
                });
            }
            (None, None) => {}
        }
    }
    // Coalesce adjacent segments with the same node id.
    let merged = squash(merged);
    // Retire any overlap segment (segment whose node was just created
    // at `time`) that no other active lineage still carries — its
    // grand MRCA was just reached.
    let kept: Lineage = merged
        .into_iter()
        .filter(|s| {
            let new_internal = ts.nodes()[s.node].time >= time - 1e-12
                && ts.nodes()[s.node].time <= time + 1e-12
                && !ts.nodes()[s.node].is_sample;
            if !new_internal {
                // Not a freshly-coalesced segment — always keep.
                true
            } else {
                // Keep only if some other lineage still carries some of
                // this segment's interval — otherwise it has MRCA'd.
                others_overlap(active, s.left, s.right)
            }
        })
        .collect();
    if !kept.is_empty() {
        active.push(kept);
    }
    Ok(())
}

/// Returns true if any active lineage covers any of `[l, r)`.
fn others_overlap(active: &[Lineage], l: f64, r: f64) -> bool {
    for lin in active {
        for s in lin {
            let lo = s.left.max(l);
            let hi = s.right.min(r);
            if hi > lo {
                return true;
            }
        }
    }
    false
}

/// Splits a lineage at a recombination breakpoint.
fn recombine(
    active: &mut Vec<Lineage>,
    rng: &mut Rng,
    lineage_rate: &[f64],
    map: &RecombinationMap,
) {
    let idx = rng.weighted_index(lineage_rate);
    let lin = active[idx].clone();
    let breakpoint = sample_breakpoint_in_lineage(rng, &lin, map);
    let (left_lin, right_lin) = split_lineage(&lin, breakpoint);
    if left_lin.is_empty() || right_lin.is_empty() {
        // The breakpoint landed exactly at a segment boundary or
        // outside the material — treat as a no-op event.
        return;
    }
    active.remove(idx);
    active.push(left_lin);
    active.push(right_lin);
}

/// Samples a recombination breakpoint within a lineage's material
/// according to the rate-map intensity restricted to that material.
fn sample_breakpoint_in_lineage(
    rng: &mut Rng,
    lin: &Lineage,
    map: &RecombinationMap,
) -> f64 {
    let weights: Vec<f64> = lin
        .iter()
        .map(|s| map.integrate(s.left, s.right))
        .collect();
    let total: f64 = weights.iter().sum();
    if total <= 0.0 {
        let s = lin[0];
        return 0.5 * (s.left + s.right);
    }
    let target = rng.uniform() * total;
    let mut acc = 0.0;
    for (i, &w) in weights.iter().enumerate() {
        if acc + w >= target {
            let s = lin[i];
            let into = target - acc;
            let global_target = map.cumulative(s.left) + into;
            return map
                .position_at_cumulative(global_target)
                .clamp(s.left, s.right - 1e-12);
        }
        acc += w;
    }
    let s = *lin.last().expect("non-empty");
    0.5 * (s.left + s.right)
}

/// All boundary points from x and y's segments, sorted.
fn collect_breakpoints(x: &Lineage, y: &Lineage) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::with_capacity(2 * (x.len() + y.len()));
    for s in x.iter().chain(y.iter()) {
        out.push(s.left);
        out.push(s.right);
    }
    out.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    out
}

/// Coalesces adjacent segments sharing the same node id.
fn squash(mut seg: Lineage) -> Lineage {
    seg.sort_by(|a, b| {
        a.left.partial_cmp(&b.left).unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut out: Lineage = Vec::with_capacity(seg.len());
    for s in seg {
        if let Some(last) = out.last_mut() {
            if last.node == s.node && (s.left - last.right).abs() < 1e-12 {
                last.right = s.right;
                continue;
            }
        }
        out.push(s);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_phylo::tree::Tree;

    fn tree_height(t: &Tree) -> f64 {
        t.patristic_distance(t.root(), t.leaves()[0])
    }

    #[test]
    fn rate_map_uniform_round_trips() {
        let m = RecombinationMap::uniform(1000.0, 2e-4).unwrap();
        assert!((m.total_rate() - 0.2).abs() < 1e-12);
        assert!((m.integrate(0.0, 500.0) - 0.1).abs() < 1e-12);
        // CDF inversion is exact for a uniform map: cum(750) = 0.15.
        let x = m.position_at_cumulative(0.15);
        assert!((x - 750.0).abs() < 1e-9, "got {x}");
    }

    #[test]
    fn rate_map_piecewise_integrates_correctly() {
        // [0, 100) rate 1, [100, 200) rate 5, [200, 300) rate 0.
        let m = RecombinationMap::piecewise(
            vec![0.0, 100.0, 200.0, 300.0],
            vec![1.0, 5.0, 0.0],
        )
        .unwrap();
        assert!((m.integrate(0.0, 100.0) - 100.0).abs() < 1e-9);
        // integrate 50..250: 50 (rate 1) * 1 + 100 (rate 5) * 1 + 50
        // (rate 0) * 1 = 50 + 500 + 0 = 550.
        assert!((m.integrate(50.0, 250.0) - 550.0).abs() < 1e-9);
        assert!((m.total_rate() - 600.0).abs() < 1e-9);
        // CDF inversion: cum(100) == 100 -> pos 100.
        assert!((m.position_at_cumulative(100.0) - 100.0).abs() < 1e-9);
        // Half-way through window 2: cum = 100 + 250 = 350 -> pos 150.
        assert!((m.position_at_cumulative(350.0) - 150.0).abs() < 1e-9);
        // At total: pos 200 (the start of the zero-rate window — the
        // zero-rate window contributes nothing and CDF saturates there).
        let p = m.position_at_cumulative(600.0);
        assert!(p >= 200.0 - 1e-9, "expected >= 200, got {p}");
    }

    #[test]
    fn rate_map_rejects_bad_input() {
        assert!(RecombinationMap::uniform(0.0, 1.0).is_err());
        assert!(RecombinationMap::uniform(100.0, -1.0).is_err());
        assert!(RecombinationMap::piecewise(vec![0.0], vec![]).is_err());
        assert!(RecombinationMap::piecewise(vec![1.0, 2.0], vec![1.0]).is_err());
        assert!(
            RecombinationMap::piecewise(vec![0.0, 5.0, 3.0], vec![1.0, 1.0]).is_err()
        );
        assert!(
            RecombinationMap::piecewise(vec![0.0, 5.0], vec![1.0, 2.0]).is_err()
        );
    }

    #[test]
    fn arg_with_no_recombination_is_a_single_tree() {
        let ts = simulate_arg(
            ArgParams::uniform(6, 1000.0, 0.0, 100.0, 42).unwrap(),
        )
        .unwrap();
        assert_eq!(ts.sample_count(), 6);
        assert_eq!(ts.tree_count(), 1);
        let tree = ts.local_tree(50.0).unwrap();
        assert_eq!(tree.leaf_count(), 6);
    }

    #[test]
    fn arg_with_recombination_has_multiple_trees() {
        let ts = simulate_arg(
            ArgParams::uniform(8, 1000.0, 5e-4, 5000.0, 7).unwrap(),
        )
        .unwrap();
        assert_eq!(ts.sample_count(), 8);
        assert!(
            ts.tree_count() > 1,
            "recombination produced only {} tree(s)",
            ts.tree_count()
        );
    }

    #[test]
    fn arg_is_deterministic() {
        let p = ArgParams::uniform(6, 800.0, 1e-4, 2000.0, 3).unwrap();
        let a = simulate_arg(p.clone()).unwrap();
        let b = simulate_arg(p).unwrap();
        assert_eq!(a.node_count(), b.node_count());
        assert_eq!(a.edge_count(), b.edge_count());
        assert_eq!(a.tree_count(), b.tree_count());
    }

    #[test]
    fn every_local_tree_is_extractable() {
        let ts = simulate_arg(
            ArgParams::uniform(5, 1000.0, 3e-4, 3000.0, 11).unwrap(),
        )
        .unwrap();
        for pos in [10.0, 750.0, 1500.0, 2900.0] {
            let tree = ts.local_tree(pos).unwrap();
            assert_eq!(tree.leaf_count(), 5);
            assert!(tree.validate().is_ok());
        }
    }

    #[test]
    fn arg_rejects_bad_params() {
        assert!(
            ArgParams::uniform(1, 1000.0, 0.0, 100.0, 1)
                .and_then(simulate_arg)
                .is_err()
        );
        // Negative N is rejected by validate.
        let bad = ArgParams::uniform(4, -1.0, 0.0, 100.0, 1);
        assert!(bad.is_err() || bad.unwrap().validate().is_err());
    }

    /// Hudson's marginal tree at any position is still Kingman, so
    /// E[TMRCA] for n lineages is `2N(1 - 1/n)`.
    #[test]
    fn arg_marginal_tree_height_matches_kingman_expectation() {
        let n = 8usize;
        let ne = 1000.0;
        let reps = 100;
        let mut acc = 0.0;
        for seed in 0..reps {
            let ts = simulate_arg(
                ArgParams::uniform(n, ne, 0.0, 100.0, seed).unwrap(),
            )
            .unwrap();
            let tree = ts.local_tree(50.0).unwrap();
            acc += tree_height(&tree);
        }
        let mean = acc / reps as f64;
        let expected = 2.0 * ne * (1.0 - 1.0 / n as f64);
        assert!(
            (mean - expected).abs() / expected < 0.25,
            "mean TMRCA {mean} vs expected {expected}"
        );
    }

    /// With a single hotspot, breakpoints concentrate in the hot
    /// region. We use a very large rate ratio for a deterministic
    /// signal at modest sample size.
    #[test]
    fn hotspot_region_creates_more_recombinations() {
        // [0, 1000) cold rate 1e-7, [1000, 1500) hot rate 1e-3,
        // [1500, 2500) cold rate 1e-7.
        let map = RecombinationMap::piecewise(
            vec![0.0, 1000.0, 1500.0, 2500.0],
            vec![1e-7, 1e-3, 1e-7],
        )
        .unwrap();
        let params = ArgParams {
            sample_size: 6,
            effective_size: 5000.0,
            recombination_map: map,
            seed: 42,
        };
        let ts = simulate_arg(params).unwrap();
        // Count interior edge endpoints (breakpoints).
        let mut hot = 0usize;
        let mut cold = 0usize;
        for e in ts.edges() {
            for &x in &[e.left, e.right] {
                // Skip the chromosome ends (always edge endpoints).
                if x > 1e-9 && x < 2500.0 - 1e-9 {
                    if (1000.0..1500.0).contains(&x) {
                        hot += 1;
                    } else {
                        cold += 1;
                    }
                }
            }
        }
        assert!(
            hot > cold,
            "hotspot breakpoints {hot} not greater than cold {cold}",
        );
    }

    /// A zero-rate map produces a single local tree.
    #[test]
    fn zero_rate_map_yields_single_tree() {
        let map = RecombinationMap::piecewise(
            vec![0.0, 500.0, 1000.0],
            vec![0.0, 0.0],
        )
        .unwrap();
        let ts = simulate_arg(ArgParams {
            sample_size: 4,
            effective_size: 1000.0,
            recombination_map: map,
            seed: 1,
        })
        .unwrap();
        assert_eq!(ts.tree_count(), 1);
    }

    #[test]
    fn split_lineage_round_trips() {
        let lin = vec![
            Segment { left: 0.0, right: 50.0, node: 0 },
            Segment { left: 50.0, right: 100.0, node: 1 },
        ];
        let (l, r) = split_lineage(&lin, 75.0);
        assert!((span(&l) - 75.0).abs() < 1e-9);
        assert!((span(&r) - 25.0).abs() < 1e-9);
        let (l, r) = split_lineage(&lin, 25.0);
        assert!((span(&l) - 25.0).abs() < 1e-9);
        assert!((span(&r) - 75.0).abs() < 1e-9);
    }

    #[test]
    fn squash_merges_adjacent_same_node() {
        let v = vec![
            Segment { left: 0.0, right: 10.0, node: 5 },
            Segment { left: 10.0, right: 30.0, node: 5 },
            Segment { left: 30.0, right: 40.0, node: 7 },
        ];
        let s = squash(v);
        assert_eq!(s.len(), 2);
        assert!((s[0].right - 30.0).abs() < 1e-9);
    }

    /// More samples -> a deeper tree on average.
    #[test]
    fn larger_sample_gives_a_deeper_tree() {
        let ne = 1000.0;
        let mut acc_small = 0.0;
        let mut acc_large = 0.0;
        for seed in 0..40 {
            let s = simulate_arg(
                ArgParams::uniform(4, ne, 0.0, 100.0, seed).unwrap(),
            )
            .unwrap();
            let t = simulate_arg(
                ArgParams::uniform(16, ne, 0.0, 100.0, seed).unwrap(),
            )
            .unwrap();
            acc_small += tree_height(&s.local_tree(50.0).unwrap());
            acc_large += tree_height(&t.local_tree(50.0).unwrap());
        }
        assert!(acc_large > acc_small);
    }
}
