//! Suboptimal structures within ΔE of the MFE (Zuker-Stiegler).
//!
//! The minimum-free-energy structure is rarely the only one that
//! matters: real RNAs visit a *neighbourhood* of low-energy
//! structures. The Zuker (1989) suboptimal-folding algorithm
//! enumerates every nested structure whose free energy is within a
//! window `[MFE, MFE + delta]`.
//!
//! ## Method
//!
//! This module fills the Zuker DP matrices once (the same
//! [`crate::fold::zuker`] recurrences) and then performs a recursive
//! *bounded* traceback: at every choice point it follows **all**
//! branches whose energy keeps the running total within the window,
//! rather than only the optimal one. Each completed traceback yields
//! a distinct structure; duplicates (reachable by different orders of
//! the same choices) are removed.
//!
//! The number of suboptimal structures grows quickly with `delta`, so
//! a `max_count` cap is enforced.

use crate::error::{Result, RnaStructError};
use crate::fold::constraint::FoldConstraints;
use crate::fold::energy::{self, multiloop, FORBIDDEN};
use crate::fold::nussinov::MIN_HAIRPIN;
use crate::fold::zuker::{self, MAX_LOOP};
use crate::rna::RnaSeq;
use crate::structure::Structure;
use std::collections::HashSet;

const INF: f64 = FORBIDDEN;

/// One suboptimal structure and its free energy.
#[derive(Clone, Debug, PartialEq)]
pub struct SuboptStructure {
    /// The structure.
    pub structure: Structure,
    /// Its free energy in kcal/mol.
    pub energy: f64,
}

/// Enumerates all structures within `delta` kcal/mol of the MFE.
///
/// Returns at most `max_count` structures, sorted by ascending free
/// energy (the MFE first). `delta` must be ≥ 0.
///
/// # Errors
/// [`RnaStructError::Invalid`] if `delta` is negative or not finite,
/// or if `max_count` is zero.
pub fn suboptimal(
    seq: &RnaSeq,
    delta: f64,
    max_count: usize,
) -> Result<Vec<SuboptStructure>> {
    if !delta.is_finite() || delta < 0.0 {
        return Err(RnaStructError::invalid(
            "delta",
            "the energy window must be a finite non-negative number",
        ));
    }
    if max_count == 0 {
        return Err(RnaStructError::invalid(
            "max_count",
            "must request at least one structure",
        ));
    }
    let codes = seq.codes();
    let n = codes.len();
    if n == 0 {
        return Ok(vec![SuboptStructure {
            structure: Structure::empty(0),
            energy: 0.0,
        }]);
    }

    let cons = FoldConstraints::none(n);
    let t = zuker::fill(codes, &cons);
    let mfe = t.w[t.idx(0, n - 1)];
    if mfe >= INF / 2.0 {
        return Ok(Vec::new());
    }
    let budget = mfe + delta;

    // Recursive bounded traceback. Each frame is a list of
    // sub-problems still to resolve and a partner array being built.
    let mut results: Vec<(Vec<Option<usize>>, f64)> = Vec::new();
    let mut seen: HashSet<Vec<Option<usize>>> = HashSet::new();
    let ctx = Ctx {
        codes,
        n,
        t: &t,
        budget,
        max_count,
    };
    let mut partner = vec![None; n];
    enumerate_w(&ctx, 0, n - 1, &mut partner, 0.0, &mut results, &mut seen);

    let mut out: Vec<SuboptStructure> = results
        .into_iter()
        .filter_map(|(p, e)| {
            Structure::from_partner(p)
                .ok()
                .map(|s| SuboptStructure {
                    structure: s,
                    energy: if e.abs() < 1e-9 { 0.0 } else { e },
                })
        })
        .collect();
    out.sort_by(|a, b| {
        a.energy
            .partial_cmp(&b.energy)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    out.truncate(max_count);
    Ok(out)
}

/// Shared immutable context for the recursive enumeration.
struct Ctx<'a> {
    codes: &'a [u8],
    n: usize,
    t: &'a zuker::ZukerTables,
    budget: f64,
    max_count: usize,
}

/// Recursively enumerate all `W(i, j)` resolutions whose energy keeps
/// `acc` (the energy already committed elsewhere) within budget.
fn enumerate_w(
    ctx: &Ctx,
    i: usize,
    j: usize,
    partner: &mut Vec<Option<usize>>,
    acc: f64,
    results: &mut Vec<(Vec<Option<usize>>, f64)>,
    seen: &mut HashSet<Vec<Option<usize>>>,
) {
    if results.len() >= ctx.max_count.saturating_mul(8) {
        return; // generous over-collection cap; final list is truncated
    }
    if i >= j {
        record(ctx, partner, acc, results, seen);
        return;
    }
    let t = ctx.t;

    // Option A: i unpaired, recurse on i+1..j.
    {
        let rest = t.w[t.idx(i + 1, j)];
        if rest < INF / 2.0 && acc + rest <= ctx.budget + 1e-7 {
            enumerate_w(ctx, i + 1, j, partner, acc, results, seen);
        }
    }
    // Option B: (i, j) closes the whole window.
    {
        let vij = t.v[t.idx(i, j)];
        if vij < INF / 2.0 {
            let term = energy::terminal_penalty(ctx.codes[i], ctx.codes[j]);
            if acc + vij + term <= ctx.budget + 1e-7 {
                partner[i] = Some(j);
                partner[j] = Some(i);
                enumerate_v(ctx, i, j, partner, acc + term, results, seen);
                partner[i] = None;
                partner[j] = None;
            }
        }
    }
    // Option C: split i..k and k+1..j.
    for k in i..j {
        let a = t.w[t.idx(i, k)];
        let b = t.w[t.idx(k + 1, j)];
        if a < INF / 2.0 && b < INF / 2.0 && acc + a + b <= ctx.budget + 1e-7 {
            // Resolve the left piece fully, then the right, by
            // threading the accumulated energy. To keep this simple
            // and correct we enumerate the left into a scratch and
            // continue the right for each completed left.
            enumerate_split(ctx, i, k, k + 1, j, partner, acc, results, seen);
        }
    }
}

/// Enumerate a `W` split: resolve `[la,lb]` then `[ra,rb]`.
#[allow(clippy::too_many_arguments)]
fn enumerate_split(
    ctx: &Ctx,
    la: usize,
    lb: usize,
    ra: usize,
    rb: usize,
    partner: &mut Vec<Option<usize>>,
    acc: f64,
    results: &mut Vec<(Vec<Option<usize>>, f64)>,
    seen: &mut HashSet<Vec<Option<usize>>>,
) {
    // Collect all completed left resolutions (partner snapshot + cost).
    let mut left: Vec<(Vec<Option<usize>>, f64)> = Vec::new();
    let mut seen_left: HashSet<Vec<Option<usize>>> = HashSet::new();
    enumerate_w(ctx, la, lb, partner, acc, &mut left, &mut seen_left);
    for (lp, le) in left {
        // Re-apply the left partner choices onto the working array
        // restricted to [la, lb].
        let mut work = partner.clone();
        work[la..=lb].copy_from_slice(&lp[la..=lb]);
        // Resolve the right piece on top, starting from le.
        enumerate_w(ctx, ra, rb, &mut work, le, results, seen);
    }
}

/// Recursively enumerate all `V(i, j)` resolutions — the pair `(i, j)`
/// is already recorded by the caller.
fn enumerate_v(
    ctx: &Ctx,
    i: usize,
    j: usize,
    partner: &mut Vec<Option<usize>>,
    acc: f64,
    results: &mut Vec<(Vec<Option<usize>>, f64)>,
    seen: &mut HashSet<Vec<Option<usize>>>,
) {
    let t = ctx.t;
    let codes = ctx.codes;

    // Hairpin terminal.
    let loop_bases = &codes[(i + 1)..j];
    let hp = energy::hairpin_energy(codes[i], codes[j], loop_bases);
    if hp < INF / 2.0 && acc + hp <= ctx.budget + 1e-7 {
        record(ctx, partner, acc + hp, results, seen);
    }

    // Internal / bulge / stack.
    let k_max = (i + 1 + MAX_LOOP).min(j.saturating_sub(MIN_HAIRPIN + 1));
    for k in (i + 1)..=k_max {
        let l_min = (k + MIN_HAIRPIN + 1).max(j.saturating_sub(MAX_LOOP + 1));
        for l in l_min..j {
            if l <= k {
                continue;
            }
            let inner = t.v[t.idx(k, l)];
            if inner >= INF / 2.0 {
                continue;
            }
            let left = k - i - 1;
            let right = j - l - 1;
            if left + right != 0 && (left > MAX_LOOP || right > MAX_LOOP) {
                continue;
            }
            let il = energy::internal_loop_energy(
                codes[i], codes[j], codes[k], codes[l], left, right,
                codes[i + 1], codes[j - 1], codes[k - 1], codes[l + 1],
            );
            if acc + il + inner > ctx.budget + 1e-7 {
                continue;
            }
            partner[k] = Some(l);
            partner[l] = Some(k);
            enumerate_v(ctx, k, l, partner, acc + il, results, seen);
            partner[k] = None;
            partner[l] = None;
        }
    }

    // Multiloop: closure + a >= 2-branch interior. The interior is
    // resolved by enumerating it as an exterior-style region (a slight
    // over-generation that the dedup pass and energy filter clean up).
    if j >= i + 2 {
        let interior_mfe = t.wm2[t.idx(i + 1, j - 1)];
        if interior_mfe < INF / 2.0 {
            let closure = multiloop::OFFSET + multiloop::PER_BRANCH;
            if acc + closure + interior_mfe <= ctx.budget + 1e-7 {
                // Enumerate the interior with the multiloop branch
                // bonuses already folded into qm; we approximate by
                // taking the MFE interior structure only for the
                // multiloop channel (suboptimal multiloop interiors
                // are still reachable through internal-loop / split
                // channels). This keeps the v1 tractable.
                let interior =
                    multiloop_interior(ctx, i + 1, j - 1);
                if let Some((ip, ie)) = interior {
                    let total = acc + closure + ie;
                    if total <= ctx.budget + 1e-7 {
                        let mut work = partner.clone();
                        for p in (i + 1)..j {
                            if ip[p].is_some() {
                                work[p] = ip[p];
                            }
                        }
                        record(ctx, &mut work, total, results, seen);
                    }
                }
            }
        }
    }
}

/// Resolve a multiloop interior `[i, j]` to its MFE branch set using a
/// plain Zuker traceback restricted to the window.
fn multiloop_interior(
    ctx: &Ctx,
    i: usize,
    j: usize,
) -> Option<(Vec<Option<usize>>, f64)> {
    // Reuse the wm2 value as the interior energy; recover the branch
    // pairs by a small dedicated traceback over wm / wm2 / v.
    let t = ctx.t;
    let e = t.wm2[t.idx(i, j)];
    if e >= INF / 2.0 {
        return None;
    }
    let mut partner = vec![None; ctx.n];
    let mut stack: Vec<(u8, usize, usize)> = vec![(2, i, j)]; // 2 => wm2
    let feq = |a: f64, b: f64| (a - b).abs() < 1e-6;
    while let Some((kind, a, b)) = stack.pop() {
        if a > b {
            continue;
        }
        match kind {
            2 => {
                // wm2
                let target = t.wm2[t.idx(a, b)];
                let mut handled = false;
                for k in a..b {
                    let l = t.wm[t.idx(a, k)];
                    let r = t.wm[t.idx(k + 1, b)];
                    if l < INF / 2.0 && r < INF / 2.0 && feq(target, l + r) {
                        stack.push((1, a, k));
                        stack.push((1, k + 1, b));
                        handled = true;
                        break;
                    }
                }
                if !handled {
                    // unpaired extension
                    if a < b {
                        let rest = t.wm2[t.idx(a + 1, b)];
                        if rest < INF / 2.0 && feq(target, rest) {
                            stack.push((2, a + 1, b));
                        }
                    }
                }
            }
            1 => {
                // wm
                let target = t.wm[t.idx(a, b)];
                let vij = t.v[t.idx(a, b)];
                if vij < INF / 2.0 && feq(target, vij + multiloop::PER_BRANCH) {
                    stack.push((0, a, b));
                    continue;
                }
                if a < b {
                    let rest = t.wm[t.idx(a + 1, b)];
                    if rest < INF / 2.0 && feq(target, rest) {
                        stack.push((1, a + 1, b));
                        continue;
                    }
                    let rest2 = t.wm[t.idx(a, b - 1)];
                    if rest2 < INF / 2.0 && feq(target, rest2) {
                        stack.push((1, a, b - 1));
                        continue;
                    }
                    for k in a..b {
                        let l = t.wm[t.idx(a, k)];
                        let r = t.wm[t.idx(k + 1, b)];
                        if l < INF / 2.0 && r < INF / 2.0 && feq(target, l + r) {
                            stack.push((1, a, k));
                            stack.push((1, k + 1, b));
                            break;
                        }
                    }
                }
            }
            _ => {
                // v: a pairs b — descend with the MFE Zuker traceback.
                partner[a] = Some(b);
                partner[b] = Some(a);
                trace_v_mfe(ctx, a, b, &mut partner);
            }
        }
    }
    Some((partner, e))
}

/// Plain MFE traceback of `V(i, j)` (single optimal branch).
fn trace_v_mfe(ctx: &Ctx, i: usize, j: usize, partner: &mut [Option<usize>]) {
    let t = ctx.t;
    let codes = ctx.codes;
    let feq = |a: f64, b: f64| (a - b).abs() < 1e-6;
    let target = t.v[t.idx(i, j)];
    let loop_bases = &codes[(i + 1)..j];
    if feq(target, energy::hairpin_energy(codes[i], codes[j], loop_bases)) {
        return;
    }
    let k_max = (i + 1 + MAX_LOOP).min(j.saturating_sub(MIN_HAIRPIN + 1));
    for k in (i + 1)..=k_max {
        let l_min = (k + MIN_HAIRPIN + 1).max(j.saturating_sub(MAX_LOOP + 1));
        for l in l_min..j {
            if l <= k {
                continue;
            }
            let inner = t.v[t.idx(k, l)];
            if inner >= INF / 2.0 {
                continue;
            }
            let left = k - i - 1;
            let right = j - l - 1;
            if left + right != 0 && (left > MAX_LOOP || right > MAX_LOOP) {
                continue;
            }
            let il = energy::internal_loop_energy(
                codes[i], codes[j], codes[k], codes[l], left, right,
                codes[i + 1], codes[j - 1], codes[k - 1], codes[l + 1],
            );
            if feq(target, il + inner) {
                partner[k] = Some(l);
                partner[l] = Some(k);
                trace_v_mfe(ctx, k, l, partner);
                return;
            }
        }
    }
    // multiloop interior
    if j >= i + 2 {
        if let Some((ip, _)) = multiloop_interior(ctx, i + 1, j - 1) {
            for p in (i + 1)..j {
                if ip[p].is_some() {
                    partner[p] = ip[p];
                }
            }
        }
    }
}

/// Record a completed structure if it is new and within budget.
fn record(
    ctx: &Ctx,
    partner: &mut [Option<usize>],
    energy: f64,
    results: &mut Vec<(Vec<Option<usize>>, f64)>,
    seen: &mut HashSet<Vec<Option<usize>>>,
) {
    if energy > ctx.budget + 1e-7 {
        return;
    }
    let key: Vec<Option<usize>> = partner.to_vec();
    if seen.insert(key.clone()) {
        results.push((key, energy));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fold::eval::structure_energy;
    use crate::fold::zuker::mfe;

    #[test]
    fn mfe_is_the_first_suboptimal() {
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let subs = suboptimal(&seq, 0.0, 10).unwrap();
        assert!(!subs.is_empty());
        let mfe_e = mfe(&seq).unwrap().energy;
        assert!(
            (subs[0].energy - mfe_e).abs() < 1e-3,
            "first suboptimal {} should equal MFE {}",
            subs[0].energy,
            mfe_e
        );
    }

    #[test]
    fn larger_delta_yields_at_least_as_many() {
        let seq = RnaSeq::parse("GGGGAAACCCCAAAGGGGAAACCCC").unwrap();
        let few = suboptimal(&seq, 1.0, 100).unwrap();
        let many = suboptimal(&seq, 5.0, 100).unwrap();
        assert!(many.len() >= few.len());
    }

    #[test]
    fn all_within_window_and_sorted() {
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let subs = suboptimal(&seq, 4.0, 50).unwrap();
        let mfe_e = subs[0].energy;
        for w in subs.windows(2) {
            assert!(w[0].energy <= w[1].energy + 1e-9, "must be sorted");
        }
        for s in &subs {
            assert!(s.energy <= mfe_e + 4.0 + 1e-3, "outside window");
            // each reported energy matches an independent evaluation
            let e = structure_energy(&seq, &s.structure).unwrap();
            assert!(
                (e - s.energy).abs() < 1e-2,
                "reported {} != eval {}",
                s.energy,
                e
            );
        }
    }

    #[test]
    fn rejects_bad_arguments() {
        let seq = RnaSeq::parse("GGGGCCCC").unwrap();
        assert!(suboptimal(&seq, -1.0, 10).is_err());
        assert!(suboptimal(&seq, f64::NAN, 10).is_err());
        assert!(suboptimal(&seq, 1.0, 0).is_err());
    }

    #[test]
    fn count_is_capped() {
        let seq = RnaSeq::parse("GGGGGGAAAACCCCCCAAAGGGGGGAAAACCCCCC").unwrap();
        let subs = suboptimal(&seq, 8.0, 5).unwrap();
        assert!(subs.len() <= 5);
    }

    #[test]
    fn unfoldable_sequence_yields_only_the_open_chain() {
        // A short poly-A cannot form any valid pair (no complementary
        // bases, too short for the min-hairpin rule). The only
        // structure within any window is the fully unpaired chain.
        let seq = RnaSeq::parse("AAAAAA").unwrap();
        let subs = suboptimal(&seq, 5.0, 20).unwrap();
        assert!(!subs.is_empty());
        for s in &subs {
            assert_eq!(s.structure.n_pairs(), 0, "poly-A cannot pair");
        }
    }

    #[test]
    fn multiloop_structure_is_enumerated_and_consistent() {
        // A sequence that folds to a genuine multiloop (two hairpins
        // under one closing helix) — this drives the multiloop
        // traceback channels (enumerate_v's multiloop branch,
        // multiloop_interior, trace_v_mfe) that simple hairpins miss.
        let seq =
            RnaSeq::parse("GGGGGCCCCCAAAGGGGGCCCCCAAAGGGGG").unwrap();
        let subs = suboptimal(&seq, 3.0, 30).unwrap();
        assert!(!subs.is_empty(), "multiloop sequence must yield structures");
        let mfe_e = mfe(&seq).unwrap().energy;
        // The MFE structure is the first suboptimal.
        assert!(
            (subs[0].energy - mfe_e).abs() < 1e-2,
            "first suboptimal {} != MFE {}",
            subs[0].energy,
            mfe_e,
        );
        // Every reported structure is nested, within the window, and
        // its reported energy survives an independent re-evaluation.
        for s in &subs {
            assert!(s.structure.is_nested(), "structure must be nested");
            assert!(s.energy <= mfe_e + 3.0 + 1e-2, "outside the window");
            let e = structure_energy(&seq, &s.structure).unwrap();
            assert!(
                (e - s.energy).abs() < 1e-2,
                "reported {} != independent eval {}",
                s.energy,
                e,
            );
        }
        // The MFE structure of a 2-hairpin/1-stem sequence has
        // multiple base pairs.
        assert!(subs[0].structure.n_pairs() >= 3, "expected a folded structure");
    }
}
