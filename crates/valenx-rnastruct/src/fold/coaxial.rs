//! Coaxial-stacking energy correction.
//!
//! In the bare nearest-neighbor model a multiloop (and the exterior
//! loop) is scored as a sum of *independent* helices: a linear
//! `a + b·branches + c·unpaired` term, plus a terminal penalty on each
//! helix end. That ignores a real, large stabilising effect — when two
//! helices in the loop lie **end to end** their terminal pairs stack on
//! each other almost exactly like an interior stacked pair.
//!
//! ViennaRNA's default `-d2` model adds this **coaxial-stacking** term
//! explicitly; it is the single largest contribution missing from a
//! dangle-only multiloop treatment, and the residual that previously
//! kept multi-helix folds from matching `RNAfold -d2` exactly.
//!
//! ## Model
//!
//! This module computes the coaxial-stacking correction for a *given*
//! structure (so [`crate::fold::eval`] is exact) and exposes a
//! per-loop helper the Zuker MFE recurrence uses to make the same term
//! part of the optimisation, not just the post-hoc score.
//!
//! For each loop with ≥ 2 helices we walk the loop boundary 5′→3′,
//! collect the ordered list of helix ends, and look at every pair of
//! **adjacent** helix ends:
//!
//! - **0 unpaired bases** between them → a *flush* coaxial stack,
//!   scored by [`t04::coaxial_flush`] (the `stack` table);
//! - **1 unpaired base** between them → a *mismatch-mediated* coaxial
//!   stack, scored by [`t04::coaxial_mismatch`];
//! - **≥ 2 unpaired bases** → no coaxial interaction.
//!
//! A helix end can take part in **at most one** coaxial stack, so the
//! best assignment is a maximum-weight matching of adjacencies around
//! the loop. The loop boundary is a cycle, so the matching is found by
//! a small linear DP over the cyclic sequence of gaps (the standard
//! "no two chosen edges share a vertex on a cycle" DP).
//!
//! Because every number comes from tables already in
//! [`crate::fold::turner2004`], the coaxial term adds **no new fitted
//! parameters** — it is the published `stack` / mismatch numbers
//! re-used in the geometry ViennaRNA scores them in.

use crate::fold::turner2004 as t04;

/// One helix end on a loop boundary, as seen walking the loop 5′→3′.
///
/// `left` / `right` are the encoded bases of the helix-end pair, with
/// `left` met first walking the loop boundary 5′→3′ and `right` met
/// second. The base met *second* (`right`) is the one immediately 5′ of
/// whatever follows this helix end around the loop — so two flush helix
/// ends `a`, `b` stack as the pair `(a.left, a.right)` against the pair
/// `(b.left, b.right)`.
#[derive(Copy, Clone, Debug)]
pub struct HelixEnd {
    /// Encoded base of the helix-end pair met first walking 5′→3′.
    pub left: u8,
    /// Encoded base of the helix-end pair met second — the base
    /// immediately 5′ of the next loop element.
    pub right: u8,
}

/// The coaxial-stacking energy gained when helix end `a` is immediately
/// followed by helix end `b` around a loop, separated by `gap` unpaired
/// bases of which `bridge` is the single bridging base when `gap == 1`.
///
/// Returns `0.0` (or a small value) when the two ends are too far apart
/// to interact. Never positive for a canonical adjacency.
fn adjacency_bonus(a: HelixEnd, b: HelixEnd, gap: usize, bridge: u8) -> f64 {
    match gap {
        0 => {
            // Flush: `a`'s end pair stacks directly on `b`'s end pair.
            // `a.right` is immediately 5' of `b.left`, so the two pairs
            // (a.left, a.right) and (b.left, b.right) stack as a
            // nearest-neighbor stacked pair.
            let e = t04::coaxial_flush(a.left, a.right, b.left, b.right);
            if e >= t04::INF / 2.0 {
                0.0
            } else {
                e.min(0.0)
            }
        }
        1 => {
            // Mismatch-mediated: one bridging base.
            t04::coaxial_mismatch(a.left, a.right, b.left, b.right, bridge)
                .min(0.0)
        }
        _ => 0.0,
    }
}

/// Best total coaxial-stacking energy for a loop whose helix ends, in
/// 5′→3′ cyclic order, are `ends`, with `gaps[k]` unpaired bases (and
/// bridging base `bridges[k]` when that gap is 1) between `ends[k]` and
/// `ends[(k+1) % n]`.
///
/// Returns a value ≤ 0 (the stabilising correction to add to the loop
/// energy). Each helix end is used in at most one coaxial stack — the
/// result is a maximum-weight matching over the cyclic adjacency list.
///
/// `is_cycle` is `true` for a multiloop (the boundary closes on
/// itself) and `false` for the exterior loop (a path: the last and
/// first helix ends are *not* adjacent through the free ends).
pub fn best_coaxial(
    ends: &[HelixEnd],
    gaps: &[usize],
    bridges: &[u8],
    is_cycle: bool,
) -> f64 {
    let n = ends.len();
    if n < 2 {
        return 0.0;
    }
    // Weight of choosing the coaxial stack across gap k (between end k
    // and end k+1). A chosen gap "uses" both its incident helix ends.
    let weight = |k: usize| -> f64 {
        let a = ends[k];
        let b = ends[(k + 1) % n];
        adjacency_bonus(a, b, gaps[k], bridges[k])
    };

    // Maximum-weight set of gaps such that no two chosen gaps are
    // adjacent (share a helix end). On a path this is the classic
    // "house robber" DP; on a cycle, run it twice — once forbidding the
    // wrap edge, once forcing it — and take the better.
    //
    // We minimise energy, so "weight" here is ≤ 0 and we want the most
    // negative achievable sum.
    let path_best = |edges: &[f64]| -> f64 {
        // edges[k] is the bonus of edge k; choosing edge k forbids
        // k-1 and k+1. dp_take/dp_skip track best (most negative) sum.
        let mut take = f64::INFINITY; // best with edge k chosen
        let mut skip = 0.0_f64; // best with edge k not chosen
        for &e in edges {
            let new_take = skip + e; // must have skipped k-1
            let new_skip = take.min(skip); // free to take or skip k-1
            take = new_take;
            skip = new_skip;
        }
        take.min(skip).min(0.0)
    };

    if !is_cycle || n == 2 {
        // Exterior loop: the gaps form a path of n-1 internal gaps
        // (gap n-1 wraps through the free 5'/3' ends and never stacks).
        // n == 2 on a cycle still has only one usable adjacency pairing.
        let edges: Vec<f64> = (0..n.saturating_sub(1)).map(weight).collect();
        if is_cycle && n == 2 {
            // Two helices in a multiloop: they meet at two gaps; pick
            // the better single coaxial stack (they cannot both stack —
            // each end is shared).
            let g0 = weight(0);
            let g1 = weight(1);
            return g0.min(g1).min(0.0);
        }
        return path_best(&edges);
    }

    // Multiloop: cyclic. Case A — wrap edge (n-1) not chosen: a path
    // over edges 0..n-1. Case B — wrap edge chosen: edges 1..n-2 form
    // the remaining path (0 and n-2 forbidden), plus weight(n-1).
    let all_edges: Vec<f64> = (0..n).map(weight).collect();
    let case_a = path_best(&all_edges[..n - 1]);
    let case_b = if n >= 3 {
        let inner = if n >= 4 {
            path_best(&all_edges[1..n - 2])
        } else {
            0.0
        };
        (all_edges[n - 1] + inner).min(all_edges[n - 1])
    } else {
        all_edges[n - 1]
    };
    case_a.min(case_b).min(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Encoded bases: A=0 C=1 G=2 U=3.
    fn gc() -> HelixEnd {
        // a G-C helix end: G met first, C second.
        HelixEnd { left: 2, right: 1 }
    }

    #[test]
    fn no_helices_no_bonus() {
        assert_eq!(best_coaxial(&[], &[], &[], true), 0.0);
        assert_eq!(best_coaxial(&[gc()], &[3], &[0], true), 0.0);
    }

    #[test]
    fn two_flush_helices_get_a_coaxial_stack() {
        // Two G-C helices, 0 gap on one side, big gap on the other.
        let ends = [gc(), gc()];
        let gaps = [0usize, 5];
        let bridges = [0u8, 0];
        let e = best_coaxial(&ends, &gaps, &bridges, true);
        assert!(e < 0.0, "flush adjacency should stabilise, got {e}");
        // It equals exactly one flush G-C/G-C stack: the helix end
        // pair (left=G, right=C) stacked on (left=G, right=C).
        let flush = t04::coaxial_flush(2, 1, 2, 1);
        assert!((e - flush).abs() < 1e-9, "{e} vs {flush}");
    }

    #[test]
    fn a_helix_end_is_used_at_most_once() {
        // Three flush-adjacent G-C helices in a multiloop: the middle
        // helix can stack with only one neighbour, so at most 1 stack
        // (3 ends, 3 gaps all flush — max matching on a 3-cycle = 1).
        let ends = [gc(), gc(), gc()];
        let gaps = [0usize, 0, 0];
        let bridges = [0u8, 0, 0];
        let e = best_coaxial(&ends, &gaps, &bridges, true);
        let flush = t04::coaxial_flush(2, 1, 2, 1);
        assert!(
            (e - flush).abs() < 1e-9,
            "3-cycle of flush helices must give exactly one stack: {e}"
        );
    }

    #[test]
    fn four_helices_can_take_two_disjoint_stacks() {
        // Four flush-adjacent helices: gaps 0 and 2 are disjoint, so
        // two coaxial stacks are possible.
        let ends = [gc(), gc(), gc(), gc()];
        let gaps = [0usize, 0, 0, 0];
        let bridges = [0u8; 4];
        let e = best_coaxial(&ends, &gaps, &bridges, true);
        let flush = t04::coaxial_flush(2, 1, 2, 1);
        assert!(
            (e - 2.0 * flush).abs() < 1e-9,
            "4-cycle of flush helices should give two stacks: {e}"
        );
    }

    #[test]
    fn exterior_loop_is_a_path_not_a_cycle() {
        // Two helices on the exterior loop with a 0 gap between them and
        // the rest wrapping through the free ends: one coaxial stack.
        let ends = [gc(), gc()];
        let gaps = [0usize, 9];
        let bridges = [0u8, 0];
        let cyc = best_coaxial(&ends, &gaps, &bridges, true);
        let ext = best_coaxial(&ends, &gaps, &bridges, false);
        // With one flush gap both see exactly one stack here.
        assert!(ext < 0.0 && cyc < 0.0);
        assert!((ext - cyc).abs() < 1e-9);
    }

    #[test]
    fn far_apart_helices_do_not_stack() {
        let ends = [gc(), gc()];
        let gaps = [4usize, 4];
        let bridges = [0u8, 0];
        assert_eq!(best_coaxial(&ends, &gaps, &bridges, true), 0.0);
    }
}
