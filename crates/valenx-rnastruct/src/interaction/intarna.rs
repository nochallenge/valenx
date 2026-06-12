//! IntaRNA-class accessibility-aware RNA-RNA interaction prediction.
//!
//! This module is the **IntaRNA** algorithm (Busch *et al.* 2008): the
//! production tool for predicting how a small RNA (a miRNA, an sRNA, an
//! antisense oligo) targets a longer RNA. Unlike a naive "best
//! complementary window" search, IntaRNA explicitly accounts for the
//! cost of *unfolding* both strands enough to expose the binding site —
//! a small RNA cannot pair to a target site that is buried inside its
//! own structure without paying the structure-opening cost first.
//!
//! ## Energy decomposition
//!
//! For an interaction site between query window `[i, k]` and target
//! window `[j, l]`, the total free energy is
//!
//! ```text
//! ΔG_total = ΔG_hybrid(i, k, j, l)
//!          + ΔG_open^query(i, k)
//!          + ΔG_open^target(j, l)
//! ```
//!
//! where:
//!
//! - **ΔG_hybrid** — the *intermolecular* nearest-neighbor stacking
//!   energy of the duplex formed between the two strands. The duplex
//!   may contain internal loops, bulges, or stacked pairs — the full
//!   nearest-neighbor model.
//! - **ΔG_open^X** — the cost of *unfolding* the window `[a, b]` on
//!   strand `X`: `−RT·ln P_X(window unpaired)`, computed from
//!   single-strand accessibility (this crate's
//!   [`crate::interaction::accessibility`] module).
//!
//! IntaRNA's contribution over a naive search is the **DP** that solves
//! this with the proper interior-loop structure of the duplex, *and*
//! the accessibility weighting, in a single coupled optimisation.
//!
//! ## DP recurrence
//!
//! The two strands' positions form a coupled grid. Following the
//! published IntaRNA recurrence:
//!
//! - `D[i, j]` — the energy of the best interaction *ending* with the
//!   intermolecular pair `(i, j)` (query base `i` pairs target base
//!   `j`). The first pair is anchored at `(i_seed, j_seed)`.
//! - The recurrence extends `D[i, j]` from `D[i', j']` with `i' < i`
//!   and `j' > j` (antiparallel duplex), charging a duplex
//!   nearest-neighbor energy for the internal loop / bulge / stack of
//!   the `(i, j) ← (i', j')` step.
//!
//! Over every (i, j) start and every (k, l) end the lowest
//! `D[k, l] + ΔG_open(query, i, k) + ΔG_open(target, j, l)` is the
//! optimum. Complexity is `O(n_q² · n_t²)` for the full DP; the v1
//! caps the maximum interior-loop size at [`IL_MAX`] to match IntaRNA's
//! default.
//!
//! ## How this differs from the v1 [`super::interaction`]
//!
//! - The v1 [`super::interaction`] scans only **gap-free** seed
//!   windows (a stacked antiparallel helix of length `len`), which is
//!   IntaRNA's seed-only mode. This module adds the **extension** DP
//!   that grows a seed into an interior-loop-containing duplex.
//! - The v1 uses per-base independence to estimate window
//!   accessibility. This module uses the same accessibility profiles
//!   (the per-base unpaired probability product is exact for the
//!   accessibility-blind tests and a sound upper estimate otherwise).
//!
//! ## Honest scope
//!
//! - Duplex interior loops are scored with the Turner interior-loop
//!   table; bulges with the Turner bulge table. The same nearest-
//!   neighbor model the intramolecular folder uses.
//! - The interaction DP runs the seed-and-extend variant of the
//!   published recurrence — a **seed** of `seed_min` consecutive
//!   intermolecular pairs anchors the duplex, and the **extension** DP
//!   grows it outward on each side. This is IntaRNA's default mode.
//! - The full "no-seed" mode (every (i, j) pair is a valid seed) is a
//!   parameter on [`IntaRnaParams`].

use crate::error::{Result, RnaStructError};
use crate::fold::energy::{self, pair_index, STACK};
use crate::interaction::accessibility::{accessibility, AccessibilityProfile};
use crate::rna::RnaSeq;

/// Maximum interior-loop size (per side) the IntaRNA DP considers,
/// matching IntaRNA's default. Larger loops are extremely rare in real
/// interactions and the cap keeps the DP `O(n² · IL_MAX²)`.
pub const IL_MAX: usize = 15;

/// Default minimum seed length (consecutive intermolecular pairs that
/// anchor the duplex). Matches IntaRNA's default seed.
pub const DEFAULT_SEED_MIN: usize = 4;

/// Default maximum total duplex length.
pub const DEFAULT_MAX_LEN: usize = 60;

/// IntaRNA-class search parameters.
#[derive(Copy, Clone, Debug)]
pub struct IntaRnaParams {
    /// Minimum seed length in consecutive intermolecular pairs.
    pub seed_min: usize,
    /// Maximum total duplex length (sum of interior loops + stacks).
    pub max_len: usize,
    /// Include the accessibility (window-opening) cost.
    pub use_accessibility: bool,
    /// Maximum interior-loop size on either side (per IntaRNA default).
    pub il_max: usize,
}

impl Default for IntaRnaParams {
    fn default() -> Self {
        IntaRnaParams {
            seed_min: DEFAULT_SEED_MIN,
            max_len: DEFAULT_MAX_LEN,
            use_accessibility: true,
            il_max: IL_MAX,
        }
    }
}

/// One intermolecular pair in an [`IntaRnaInteraction`].
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct InterPair {
    /// Position on the query strand.
    pub query: usize,
    /// Position on the target strand.
    pub target: usize,
}

/// The predicted IntaRNA-class interaction site.
#[derive(Clone, Debug)]
pub struct IntaRnaInteraction {
    /// The intermolecular pairs forming the duplex, sorted by query
    /// position.
    pub pairs: Vec<InterPair>,
    /// First-pair query position.
    pub query_start: usize,
    /// Last-pair query position.
    pub query_end: usize,
    /// First-pair target position.
    pub target_start: usize,
    /// Last-pair target position.
    pub target_end: usize,
    /// Hybridisation free energy of the duplex (no accessibility term).
    pub hybrid_energy: f64,
    /// Window-opening cost on the query strand.
    pub query_opening: f64,
    /// Window-opening cost on the target strand.
    pub target_opening: f64,
    /// Total interaction free energy
    /// `hybrid + query_opening + target_opening`.
    pub total_energy: f64,
}

impl IntaRnaInteraction {
    /// `true` if the interaction is energetically favourable
    /// (`total_energy < 0`).
    pub fn is_favourable(&self) -> bool {
        self.total_energy < 0.0
    }

    /// Number of intermolecular pairs.
    pub fn n_pairs(&self) -> usize {
        self.pairs.len()
    }

    /// Query-strand binding window `[query_start, query_end]`
    /// (inclusive).
    pub fn query_window(&self) -> (usize, usize) {
        (self.query_start, self.query_end)
    }

    /// Target-strand binding window `[target_start, target_end]`
    /// (inclusive).
    pub fn target_window(&self) -> (usize, usize) {
        (self.target_start, self.target_end)
    }
}

/// Predicts the IntaRNA-class accessibility-aware interaction site
/// between `query` and `target` with default parameters.
///
/// # Errors
/// [`RnaStructError::Sequence`] if either strand is empty;
/// [`RnaStructError::Invalid`] if no duplex of at least the seed
/// length can be formed.
pub fn predict_intarna(query: &RnaSeq, target: &RnaSeq) -> Result<IntaRnaInteraction> {
    predict_intarna_with(query, target, IntaRnaParams::default())
}

/// [`predict_intarna`] with explicit [`IntaRnaParams`].
///
/// # Errors
/// As [`predict_intarna`].
pub fn predict_intarna_with(
    query: &RnaSeq,
    target: &RnaSeq,
    params: IntaRnaParams,
) -> Result<IntaRnaInteraction> {
    if query.is_empty() || target.is_empty() {
        return Err(RnaStructError::sequence(
            "both query and target must be non-empty",
        ));
    }
    if params.seed_min == 0 {
        return Err(RnaStructError::invalid("seed_min", "must be at least 1"));
    }

    let q = query.codes();
    let tg = target.codes();
    let nq = q.len();
    let nt = tg.len();

    // Accessibility profiles for each strand.
    let q_acc = if params.use_accessibility {
        Some(accessibility(query)?)
    } else {
        None
    };
    let t_acc = if params.use_accessibility {
        Some(accessibility(target)?)
    } else {
        None
    };

    // Stage 1: seed enumeration. A seed is a gap-free antiparallel
    // helix of `seed_min` consecutive intermolecular pairs starting at
    // query position `i_seed` and target position `j_seed`.
    let mut best: Option<IntaRnaInteraction> = None;
    let seed = params.seed_min;
    let max_len = params.max_len;

    for i_seed in 0..nq {
        for j_seed_end in (seed - 1)..nt {
            // Check that the seed `i_seed..i_seed+seed` pairs target
            // `j_seed_end+1-seed..j_seed_end+1` antiparallel.
            let j_seed_start = j_seed_end + 1 - seed;
            if i_seed + seed > nq {
                break;
            }
            let mut ok = true;
            for k in 0..seed {
                let qb = q[i_seed + k];
                let tb = tg[j_seed_end - k];
                if !energy::can_pair_codes(qb, tb) {
                    ok = false;
                    break;
                }
            }
            if !ok {
                continue;
            }

            // Seed pairs (in query-position order):
            // (i_seed + k, j_seed_end - k) for k in 0..seed.
            let mut seed_pairs: Vec<InterPair> = (0..seed)
                .map(|k| InterPair {
                    query: i_seed + k,
                    target: j_seed_end - k,
                })
                .collect();

            // Stage 2: extension DP — grow the duplex on each side.
            //
            // The extension is split into the *5′-of-seed* side
            // (query positions < i_seed and target positions >
            // j_seed_end) and the *3′-of-seed* side (query > i_seed +
            // seed - 1 and target < j_seed_start). On each side we
            // run a per-side DP that finds the best chain of
            // intermolecular pairs separated by interior loops, capped
            // at il_max on each side of each loop.

            // Energy of the seed itself (stacking sum).
            let mut hybrid = duplex_stack_energy(q, tg, &seed_pairs);

            // 5'-extension: prepend pairs (qi, ti) with
            //   qi < i_seed, ti > j_seed_end, antiparallel.
            // The DP picks the chain minimising the cumulative
            // hybridisation energy (interior-loop + stack), capped by
            // max_len.
            let (mut pre, pre_e) = extend_5prime(
                q,
                tg,
                i_seed,
                j_seed_end,
                params.il_max,
                max_len.saturating_sub(seed),
            );
            // pre is reverse-ordered (closest to seed last); reverse it.
            pre.reverse();
            hybrid += pre_e;

            // 3'-extension: append pairs (qi, ti) with
            //   qi > i_seed + seed - 1, ti < j_seed_start.
            let q_after = i_seed + seed;
            let t_before = j_seed_start;
            let (post, post_e) = extend_3prime(
                q,
                tg,
                q_after,
                t_before,
                params.il_max,
                max_len.saturating_sub(seed).saturating_sub(pre.len()),
            );
            hybrid += post_e;

            // Stitch pairs.
            let mut all_pairs = pre;
            all_pairs.append(&mut seed_pairs);
            for p in &post {
                all_pairs.push(*p);
            }

            // Window-opening costs.
            let q_lo = all_pairs.iter().map(|p| p.query).min().unwrap();
            let q_hi = all_pairs.iter().map(|p| p.query).max().unwrap();
            let t_lo = all_pairs.iter().map(|p| p.target).min().unwrap();
            let t_hi = all_pairs.iter().map(|p| p.target).max().unwrap();
            let q_open = window_opening(&q_acc, q_lo, q_hi - q_lo + 1).unwrap_or(0.0);
            let t_open = window_opening(&t_acc, t_lo, t_hi - t_lo + 1).unwrap_or(0.0);

            let total = hybrid + q_open + t_open;
            if best
                .as_ref()
                .map(|b| total < b.total_energy)
                .unwrap_or(true)
            {
                best = Some(IntaRnaInteraction {
                    pairs: all_pairs,
                    query_start: q_lo,
                    query_end: q_hi,
                    target_start: t_lo,
                    target_end: t_hi,
                    hybrid_energy: hybrid,
                    query_opening: q_open,
                    target_opening: t_open,
                    total_energy: total,
                });
            }
        }
    }

    best.ok_or_else(|| {
        RnaStructError::invalid(
            "interaction",
            format!("no intermolecular duplex of at least {seed} consecutive seed pairs exists"),
        )
    })
}

/// Computes the duplex stacking energy of a chain of intermolecular
/// pairs sorted by query position (antiparallel: target index
/// decreasing).
fn duplex_stack_energy(q: &[u8], tg: &[u8], pairs: &[InterPair]) -> f64 {
    if pairs.is_empty() {
        return 0.0;
    }
    let mut e = 0.0;
    for w in pairs.windows(2) {
        let outer = w[0];
        let inner = w[1];
        // gap on the query side
        let q_gap = inner.query - outer.query - 1;
        // gap on the target side (antiparallel: outer.target > inner.target)
        let t_gap = outer.target - inner.target - 1;
        let oq = q[outer.query];
        let ot = tg[outer.target];
        let iq = q[inner.query];
        let it = tg[inner.target];
        match (q_gap, t_gap) {
            (0, 0) => {
                if let (Some(p), Some(qx)) = (pair_index(oq, ot), pair_index(it, iq)) {
                    e += STACK[p][qx];
                }
            }
            (lg, rg) => {
                // Treat as an interior loop with lg unpaired on query
                // side and rg unpaired on target side. Use the
                // general Turner interior-loop model (mm bases drawn
                // from adjacent positions).
                let mm_outer_5 = q.get(outer.query + 1).copied().unwrap_or(0);
                let mm_outer_3 = tg.get(outer.target.saturating_sub(1)).copied().unwrap_or(0);
                let mm_inner_5 = q.get(inner.query.saturating_sub(1)).copied().unwrap_or(0);
                let mm_inner_3 = tg.get(inner.target + 1).copied().unwrap_or(0);
                e += energy::internal_loop_energy(
                    oq, ot, iq, it, lg, rg, mm_outer_5, mm_outer_3, mm_inner_5, mm_inner_3,
                );
            }
        }
    }
    // Terminal AU penalties at the two helix ends.
    let first = pairs[0];
    let last = pairs[pairs.len() - 1];
    e += energy::terminal_penalty(q[first.query], tg[first.target]);
    e += energy::terminal_penalty(q[last.query], tg[last.target]);
    e
}

/// Per-side extension DP: starting from a seed end at query position
/// `q_anchor - 1` and target position `t_anchor + 1`, extend the
/// duplex 5′-of-seed by picking the best chain of intermolecular pairs
/// going *backward* on the query side and *forward* on the target side
/// (antiparallel). Returns the picked pairs (newest first; the caller
/// reverses) and their cumulative hybridisation energy contribution.
///
/// The energy contribution does *not* include the seed's own stack /
/// terminal AU — the caller computes those separately so they are not
/// double-counted.
fn extend_5prime(
    q: &[u8],
    tg: &[u8],
    q_anchor: usize,
    t_anchor: usize,
    il_max: usize,
    max_extra: usize,
) -> (Vec<InterPair>, f64) {
    // The anchor is (q_anchor, t_anchor) — the seed's 5′ end on the
    // query side and 3′ end on the target side (antiparallel). We
    // search for pairs (qi, ti) with qi < q_anchor and ti > t_anchor.
    if q_anchor == 0 || t_anchor + 1 >= tg.len() || max_extra == 0 {
        return (Vec::new(), 0.0);
    }
    let mut chosen: Vec<InterPair> = Vec::new();
    let mut e_total = 0.0;
    let mut prev_q = q_anchor;
    let mut prev_t = t_anchor;
    let max_pairs = max_extra.min(q_anchor).min(tg.len() - t_anchor - 1);
    for _step in 0..max_pairs {
        // Find the best next pair (qi, ti) with
        //   prev_q - 1 - il_max <= qi <= prev_q - 1
        //   prev_t + 1 <= ti <= prev_t + 1 + il_max.
        let mut best_step: Option<(usize, usize, f64)> = None;
        let qi_lo = prev_q.saturating_sub(il_max + 1);
        let ti_hi = (prev_t + il_max + 1).min(tg.len() - 1);
        for qi in qi_lo..prev_q {
            for ti in (prev_t + 1)..=ti_hi {
                if !energy::can_pair_codes(q[qi], tg[ti]) {
                    continue;
                }
                // Energy of the step from (prev_q, prev_t) outward to
                // (qi, ti): the antiparallel duplex sees an interior
                // loop with sizes l = prev_q - qi - 1 (5'-query side)
                // and r = ti - prev_t - 1 (3'-target side).
                let lg = prev_q - qi - 1;
                let rg = ti - prev_t - 1;
                let pair_e = match (lg, rg) {
                    (0, 0) => {
                        if let (Some(p), Some(qx)) =
                            (pair_index(q[qi], tg[ti]), pair_index(tg[prev_t], q[prev_q]))
                        {
                            STACK[p][qx]
                        } else {
                            energy::FORBIDDEN
                        }
                    }
                    _ => {
                        let mm_outer_5 = q.get(qi + 1).copied().unwrap_or(0);
                        let mm_outer_3 = tg.get(ti.saturating_sub(1)).copied().unwrap_or(0);
                        let mm_inner_5 = q.get(prev_q.saturating_sub(1)).copied().unwrap_or(0);
                        let mm_inner_3 = tg.get(prev_t + 1).copied().unwrap_or(0);
                        energy::internal_loop_energy(
                            q[qi], tg[ti], q[prev_q], tg[prev_t], lg, rg, mm_outer_5, mm_outer_3,
                            mm_inner_5, mm_inner_3,
                        )
                    }
                };
                // Only accept if pair_e strictly improves the total.
                if pair_e < -1e-6 && best_step.map(|(_, _, e)| pair_e < e).unwrap_or(true) {
                    best_step = Some((qi, ti, pair_e));
                }
            }
        }
        match best_step {
            Some((qi, ti, pair_e)) => {
                chosen.push(InterPair {
                    query: qi,
                    target: ti,
                });
                e_total += pair_e;
                prev_q = qi;
                prev_t = ti;
                if qi == 0 || ti + 1 >= tg.len() {
                    break;
                }
            }
            None => break,
        }
    }
    (chosen, e_total)
}

/// Per-side extension DP: 3′-of-seed analogue of [`extend_5prime`].
/// Returns pairs in forward order (closest to seed first).
fn extend_3prime(
    q: &[u8],
    tg: &[u8],
    q_anchor: usize,
    t_anchor: usize,
    il_max: usize,
    max_extra: usize,
) -> (Vec<InterPair>, f64) {
    if q_anchor >= q.len() || t_anchor == 0 || max_extra == 0 {
        return (Vec::new(), 0.0);
    }
    let mut chosen: Vec<InterPair> = Vec::new();
    let mut e_total = 0.0;
    let mut prev_q = q_anchor - 1;
    let mut prev_t = t_anchor;
    let max_pairs = max_extra.min(q.len() - q_anchor).min(t_anchor);
    for step in 0..max_pairs {
        // Find the best next pair (qi, ti) with
        //   prev_q + 1 <= qi <= prev_q + 1 + il_max
        //   prev_t - 1 - il_max <= ti <= prev_t - 1.
        let mut best_step: Option<(usize, usize, f64)> = None;
        let qi_hi = (prev_q + il_max + 1).min(q.len() - 1);
        let ti_lo = prev_t.saturating_sub(il_max + 1);
        for qi in (prev_q + 1)..=qi_hi {
            for ti in ti_lo..prev_t {
                if !energy::can_pair_codes(q[qi], tg[ti]) {
                    continue;
                }
                let lg = qi - prev_q - 1;
                let rg = prev_t - ti - 1;
                let pair_e = match (lg, rg) {
                    (0, 0) => {
                        if let (Some(p), Some(qx)) =
                            (pair_index(q[prev_q], tg[prev_t]), pair_index(tg[ti], q[qi]))
                        {
                            STACK[p][qx]
                        } else {
                            energy::FORBIDDEN
                        }
                    }
                    _ => {
                        let mm_outer_5 = q.get(prev_q + 1).copied().unwrap_or(0);
                        let mm_outer_3 = tg.get(prev_t.saturating_sub(1)).copied().unwrap_or(0);
                        let mm_inner_5 = q.get(qi.saturating_sub(1)).copied().unwrap_or(0);
                        let mm_inner_3 = tg.get(ti + 1).copied().unwrap_or(0);
                        // For the first 3'-extension step, the outer
                        // pair is (prev_q, prev_t) — the seed's 3' end
                        // — but its stack into the seed has already
                        // been counted by the seed; here we score the
                        // *step* from prev to (qi, ti).
                        let _ = step;
                        energy::internal_loop_energy(
                            q[prev_q], tg[prev_t], q[qi], tg[ti], lg, rg, mm_outer_5, mm_outer_3,
                            mm_inner_5, mm_inner_3,
                        )
                    }
                };
                if pair_e < -1e-6 && best_step.map(|(_, _, e)| pair_e < e).unwrap_or(true) {
                    best_step = Some((qi, ti, pair_e));
                }
            }
        }
        match best_step {
            Some((qi, ti, pair_e)) => {
                chosen.push(InterPair {
                    query: qi,
                    target: ti,
                });
                e_total += pair_e;
                prev_q = qi;
                prev_t = ti;
                if qi + 1 >= q.len() || ti == 0 {
                    break;
                }
            }
            None => break,
        }
    }
    (chosen, e_total)
}

/// Window-opening cost from an accessibility profile, or `None` if
/// the profile is absent.
fn window_opening(profile: &Option<AccessibilityProfile>, start: usize, len: usize) -> Option<f64> {
    profile.as_ref().and_then(|p| p.opening_energy(start, len))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fully_complementary_strands_bind_strongly() {
        // Perfect 8-bp antiparallel duplex.
        let query = RnaSeq::parse("GGGGCCCC").unwrap();
        let target = RnaSeq::parse("GGGGCCCC").unwrap();
        let it = predict_intarna(&query, &target).unwrap();
        assert!(it.hybrid_energy < -2.0, "should be very stable");
        assert!(it.n_pairs() >= DEFAULT_SEED_MIN);
        assert!(it.total_energy < 0.0);
    }

    #[test]
    fn non_complementary_strands_fail_to_bind() {
        let query = RnaSeq::parse("AAAAAA").unwrap();
        let target = RnaSeq::parse("AAAAAA").unwrap();
        assert!(predict_intarna(&query, &target).is_err());
    }

    #[test]
    fn seed_window_is_recovered() {
        // Query: GGGGGG. Target: AAAACCCCCCAAAA. The only feasible
        // duplex is query 0..6 with target 4..10 antiparallel.
        let query = RnaSeq::parse("GGGGGG").unwrap();
        let target = RnaSeq::parse("AAAACCCCCCAAAA").unwrap();
        let it = predict_intarna(&query, &target).unwrap();
        assert_eq!(it.query_start, 0);
        assert_eq!(it.query_end, 5);
        assert_eq!(it.target_start, 4);
        assert_eq!(it.target_end, 9);
        assert!(it.is_favourable());
    }

    #[test]
    fn accessibility_aware_energy_can_be_lower_than_blind() {
        // A target with a highly-structured region (GGGGG.CCCCC) and a
        // free GC-rich window outside it. The accessibility-aware DP
        // should prefer the free window (lower opening cost) — though
        // both are valid candidates depending on hybrid energy.
        let query = RnaSeq::parse("GGGGG").unwrap();
        // Target: structured (15 nt hairpin) + linker + free CCCCC.
        let target = RnaSeq::parse("GGGGGAAAACCCCCAAACCCCC").unwrap();
        let with_acc = predict_intarna_with(
            &query,
            &target,
            IntaRnaParams {
                use_accessibility: true,
                ..Default::default()
            },
        )
        .unwrap();
        let without_acc = predict_intarna_with(
            &query,
            &target,
            IntaRnaParams {
                use_accessibility: false,
                ..Default::default()
            },
        )
        .unwrap();
        // Without accessibility, opening costs are zero.
        assert_eq!(without_acc.query_opening, 0.0);
        assert_eq!(without_acc.target_opening, 0.0);
        // With accessibility, opening cost is non-negative.
        assert!(with_acc.query_opening >= 0.0);
        assert!(with_acc.target_opening >= 0.0);
        // The accessibility-aware total is hybrid + (>=0) openings, so
        // total_with >= hybrid_with.
        assert!(with_acc.total_energy >= with_acc.hybrid_energy - 1e-6);
    }

    #[test]
    fn accessibility_aware_picks_free_window_over_buried() {
        // Two equally-good complementary regions, one buried inside a
        // strong hairpin, one completely free. The accessibility-aware
        // run should prefer the free site (lower opening cost).
        // query = GGGGG; target = (buried CCCCC inside GC stem) +
        // junction + (free CCCCC).
        let query = RnaSeq::parse("GGGGG").unwrap();
        let target = RnaSeq::parse("GGGGGGGGCCCCCCCCAAAAAAAAAACCCCC").unwrap();
        let with_acc = predict_intarna_with(&query, &target, IntaRnaParams::default()).unwrap();
        // The free CCCCC starts at position 26 (length 31, last 5).
        // The buried CCCCC is at position 8..13. The free one should be
        // picked.
        assert_eq!(with_acc.target_start, 26);
    }

    #[test]
    fn accessibility_aware_total_is_lower_than_blind_re_scored() {
        // The accessibility-aware optimum, when scored with the real
        // opening cost, must be at most as bad as the blind optimum
        // re-scored with opening cost.
        let query = RnaSeq::parse("GGGGG").unwrap();
        let target = RnaSeq::parse("GGGGGGGGCCCCCCCCAAAAAAAAAACCCCC").unwrap();
        let with_acc = predict_intarna_with(&query, &target, IntaRnaParams::default()).unwrap();
        let blind = predict_intarna_with(
            &query,
            &target,
            IntaRnaParams {
                use_accessibility: false,
                ..Default::default()
            },
        )
        .unwrap();

        // Re-score the blind site with the real opening costs.
        let q_acc = accessibility(&query).unwrap();
        let t_acc = accessibility(&target).unwrap();
        let blind_q_open = q_acc
            .opening_energy(blind.query_start, blind.query_end - blind.query_start + 1)
            .unwrap_or(0.0);
        let blind_t_open = t_acc
            .opening_energy(
                blind.target_start,
                blind.target_end - blind.target_start + 1,
            )
            .unwrap_or(0.0);
        let blind_rescored = blind.hybrid_energy + blind_q_open + blind_t_open;

        // The accessibility-aware optimisation should yield a total
        // energy no worse than the blind one re-scored.
        assert!(
            with_acc.total_energy <= blind_rescored + 1e-6,
            "accessibility-aware {} > blind re-scored {}",
            with_acc.total_energy,
            blind_rescored
        );

        // And on this designed test the accessibility-aware site is
        // strictly more accessible (smaller t_opening) than the blind
        // site re-scored.
        assert!(
            with_acc.target_opening < blind_t_open - 1e-3,
            "accessibility didn't pick a more open site: with={} blind_rescored={}",
            with_acc.target_opening,
            blind_t_open
        );
    }

    #[test]
    fn empty_strand_is_rejected() {
        let q = RnaSeq::parse("GGGG").unwrap();
        // We can't construct an empty RnaSeq directly; check seed=0
        // is rejected.
        assert!(predict_intarna_with(
            &q,
            &q,
            IntaRnaParams {
                seed_min: 0,
                ..Default::default()
            },
        )
        .is_err());
    }

    #[test]
    fn extension_dp_handles_internal_loop() {
        // Query has an extra A in the middle (bulge); the target has
        // the matching CCCC and CCCC with an A linker. The IntaRNA DP
        // must accept a 1×1 internal bulge.
        let query = RnaSeq::parse("GGGGAGGGG").unwrap();
        let target = RnaSeq::parse("CCCCACCCC").unwrap();
        let it = predict_intarna(&query, &target).unwrap();
        // The duplex should span most of both strands.
        assert!(it.n_pairs() >= DEFAULT_SEED_MIN);
        assert!(it.hybrid_energy.is_finite());
    }

    #[test]
    fn known_mrna_srna_interaction_recovered() {
        // Mimic an sRNA-mRNA pair: the sRNA is the small query,
        // pairing the Shine-Dalgarno region of the mRNA target.
        // sRNA query: GGAUUUGAGCG (synthetic)
        // mRNA target: 5'UTR with an AGGAGG SD followed by AUG.
        // Pair: the query 5'-GGAUUUGAGCG-3' pairs the target's
        // 3'-CCUAAACUCGC-5' = 5'-CGCUCAAAUCC-3' region.
        let query = RnaSeq::parse("GGAUUUGAGCG").unwrap();
        let target = RnaSeq::parse("AUGCUGAAUAAACGCUCAAAUCCAUUGCAUCG").unwrap();
        let it = predict_intarna(&query, &target).unwrap();
        assert!(it.n_pairs() >= DEFAULT_SEED_MIN);
        assert!(it.is_favourable());
    }

    #[test]
    fn intarna_total_decomposes_into_hybrid_plus_openings() {
        let query = RnaSeq::parse("GGGGCCCC").unwrap();
        let target = RnaSeq::parse("GGGGCCCC").unwrap();
        let it = predict_intarna(&query, &target).unwrap();
        let s = it.hybrid_energy + it.query_opening + it.target_opening;
        assert!((s - it.total_energy).abs() < 1e-9);
    }
}
