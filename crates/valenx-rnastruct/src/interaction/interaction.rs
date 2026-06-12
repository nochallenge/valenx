//! RNA-RNA interaction prediction (IntaRNA-class).
//!
//! When a small RNA (a miRNA, an sRNA, an antisense oligo) targets a
//! longer RNA, the favourable interaction is the one that combines a
//! strong *intermolecular duplex* with low *accessibility cost* — the
//! target site must be cleared of its own structure first. IntaRNA
//! (Busch *et al.* 2008) formalised this as
//!
//! ```text
//! ΔG_total(site) = ΔG_hybrid(site) + ΔG_open(query) + ΔG_open(target)
//! ```
//!
//! and finds the site minimising it.
//!
//! ## Method (v1)
//!
//! - **Hybridisation** — for every pair of windows, one on the query
//!   and one on the (reverse-complementary-aligned) target, the
//!   duplex energy is computed as a stacked antiparallel helix using
//!   the Turner stacking table ([`crate::fold::energy::STACK`]). This
//!   is a *gap-free* duplex — the standard IntaRNA "seed + extension"
//!   simplification for a v1.
//! - **Accessibility** — the opening free energy of each window from
//!   the single-strand [`crate::interaction::accessibility`] profile.
//!
//! The reported interaction is the window pair minimising `ΔG_total`.

use crate::error::{Result, RnaStructError};
use crate::fold::energy::{self, pair_index, STACK};
use crate::interaction::accessibility::accessibility;
use crate::rna::RnaSeq;

/// A predicted RNA-RNA interaction site.
#[derive(Clone, Debug, PartialEq)]
pub struct Interaction {
    /// Start index of the bound window on the query strand.
    pub query_start: usize,
    /// Start index of the bound window on the target strand.
    pub target_start: usize,
    /// Number of consecutive intermolecular base pairs.
    pub length: usize,
    /// Hybridisation free energy of the duplex, kcal/mol.
    pub hybrid_energy: f64,
    /// Opening cost on the query strand, kcal/mol.
    pub query_opening: f64,
    /// Opening cost on the target strand, kcal/mol.
    pub target_opening: f64,
    /// Total interaction free energy (hybrid + both openings).
    pub total_energy: f64,
}

impl Interaction {
    /// `true` if the interaction is energetically favourable
    /// (`total_energy < 0`).
    pub fn is_favourable(&self) -> bool {
        self.total_energy < 0.0
    }
}

/// Parameters controlling the interaction search.
#[derive(Copy, Clone, Debug)]
pub struct InteractionParams {
    /// Minimum number of consecutive intermolecular pairs a reported
    /// interaction must have (the "seed length").
    pub min_seed: usize,
    /// Maximum duplex length to consider.
    pub max_length: usize,
    /// If `true`, the accessibility (opening-energy) term is included;
    /// if `false`, only the raw hybridisation energy is scored.
    pub use_accessibility: bool,
}

impl Default for InteractionParams {
    fn default() -> Self {
        InteractionParams {
            min_seed: 4,
            max_length: 20,
            use_accessibility: true,
        }
    }
}

/// Predicts the best RNA-RNA interaction site between `query` and
/// `target` with default parameters.
///
/// # Errors
/// [`RnaStructError::Sequence`] if either strand is empty;
/// [`RnaStructError::Invalid`] if no duplex of at least the seed
/// length can be formed.
pub fn predict_interaction(query: &RnaSeq, target: &RnaSeq) -> Result<Interaction> {
    predict_interaction_with(query, target, InteractionParams::default())
}

/// [`predict_interaction`] with explicit [`InteractionParams`].
///
/// # Errors
/// As [`predict_interaction`].
pub fn predict_interaction_with(
    query: &RnaSeq,
    target: &RnaSeq,
    params: InteractionParams,
) -> Result<Interaction> {
    if query.is_empty() || target.is_empty() {
        return Err(RnaStructError::sequence(
            "both query and target must be non-empty",
        ));
    }
    if params.min_seed == 0 {
        return Err(RnaStructError::invalid("min_seed", "must be at least 1"));
    }

    let q = query.codes();
    let tg = target.codes();

    // Accessibility profiles (opening energy uses these).
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

    let mut best: Option<Interaction> = None;

    // A gap-free antiparallel duplex: query position qs+k pairs target
    // position ts+len-1-k. We slide both windows and grow the length.
    for qs in 0..q.len() {
        for ts in 0..tg.len() {
            let max_len = params.max_length.min(q.len() - qs).min(tg.len() - ts);
            for len in params.min_seed..=max_len {
                // Check that every position in the duplex pairs.
                let mut ok = true;
                for k in 0..len {
                    let qb = q[qs + k];
                    let tb = tg[ts + len - 1 - k];
                    if !energy::can_pair_codes(qb, tb) {
                        ok = false;
                        break;
                    }
                }
                if !ok {
                    continue;
                }
                let hybrid = duplex_energy(q, tg, qs, ts, len);
                let q_open = q_acc
                    .as_ref()
                    .and_then(|a| a.opening_energy(qs, len))
                    .unwrap_or(0.0);
                let t_open = t_acc
                    .as_ref()
                    .and_then(|a| a.opening_energy(ts, len))
                    .unwrap_or(0.0);
                let total = hybrid + q_open + t_open;
                if best
                    .as_ref()
                    .map(|b| total < b.total_energy)
                    .unwrap_or(true)
                {
                    best = Some(Interaction {
                        query_start: qs,
                        target_start: ts,
                        length: len,
                        hybrid_energy: hybrid,
                        query_opening: q_open,
                        target_opening: t_open,
                        total_energy: total,
                    });
                }
            }
        }
    }

    best.ok_or_else(|| {
        RnaStructError::invalid(
            "interaction",
            format!(
                "no intermolecular duplex of at least {} consecutive pairs exists",
                params.min_seed
            ),
        )
    })
}

/// Free energy of a gap-free antiparallel intermolecular duplex.
///
/// The duplex pairs query `qs..qs+len` with target `ts..ts+len`
/// antiparallel. The energy is the sum of the `len-1` nearest-neighbor
/// stacks plus a terminal-AU penalty at each end — a standard duplex
/// nearest-neighbor evaluation.
fn duplex_energy(q: &[u8], tg: &[u8], qs: usize, ts: usize, len: usize) -> f64 {
    if len < 2 {
        // a single pair: just the terminal penalty
        return energy::terminal_penalty(q[qs], tg[ts + len - 1]);
    }
    let mut e = 0.0;
    for k in 0..(len - 1) {
        // outer pair: query qs+k  —  target ts+len-1-k
        // inner pair: query qs+k+1 — target ts+len-2-k
        let o_q = q[qs + k];
        let o_t = tg[ts + len - 1 - k];
        let i_q = q[qs + k + 1];
        let i_t = tg[ts + len - 2 - k];
        if let (Some(p), Some(qq)) = (pair_index(o_q, o_t), pair_index(i_t, i_q)) {
            e += STACK[p][qq];
        }
    }
    // terminal penalties at the two helix ends
    e += energy::terminal_penalty(q[qs], tg[ts + len - 1]);
    e += energy::terminal_penalty(q[qs + len - 1], tg[ts]);
    e
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complementary_rnas_interact_favourably() {
        // query GGGGGG pairs target CCCCCC antiparallel
        let query = RnaSeq::parse("GGGGGG").unwrap();
        let target = RnaSeq::parse("CCCCCC").unwrap();
        let it = predict_interaction(&query, &target).unwrap();
        assert!(it.is_favourable(), "complementary RNAs should bind");
        assert!(it.length >= 4);
        assert!(it.hybrid_energy < 0.0);
    }

    #[test]
    fn non_complementary_rnas_have_no_duplex() {
        let query = RnaSeq::parse("AAAAAA").unwrap();
        let target = RnaSeq::parse("AAAAAA").unwrap();
        // A-A cannot pair, so no seed exists.
        assert!(predict_interaction(&query, &target).is_err());
    }

    #[test]
    fn finds_the_right_window() {
        // target has a CCCCCC stretch only at positions 4..10
        let query = RnaSeq::parse("GGGGGG").unwrap();
        let target = RnaSeq::parse("AAAACCCCCCAAAA").unwrap();
        let it = predict_interaction(&query, &target).unwrap();
        assert_eq!(it.target_start, 4);
        assert_eq!(it.length, 6);
    }

    #[test]
    fn empty_sequence_rejected() {
        let q = RnaSeq::parse("GGGG").unwrap();
        // we cannot construct an empty RnaSeq directly; check seed=0
        assert!(predict_interaction_with(
            &q,
            &q,
            InteractionParams {
                min_seed: 0,
                ..Default::default()
            }
        )
        .is_err());
    }

    #[test]
    fn accessibility_toggle_changes_energy() {
        let query = RnaSeq::parse("GGGGGG").unwrap();
        let target = RnaSeq::parse("CCCCCCAAAAGGGGGG").unwrap();
        let with = predict_interaction_with(
            &query,
            &target,
            InteractionParams {
                use_accessibility: true,
                ..Default::default()
            },
        )
        .unwrap();
        let without = predict_interaction_with(
            &query,
            &target,
            InteractionParams {
                use_accessibility: false,
                ..Default::default()
            },
        )
        .unwrap();
        // without accessibility, the openings are zero
        assert_eq!(without.query_opening, 0.0);
        assert_eq!(without.target_opening, 0.0);
        // with accessibility, total includes a (non-negative) cost
        assert!(with.total_energy >= with.hybrid_energy - 1e-6);
    }

    #[test]
    fn min_seed_is_respected() {
        // query GGGG.AAA  vs  target AAA.CCCC — the GGGG/CCCC ends
        // form a contiguous 4-bp antiparallel duplex, but the A·A
        // mismatch in the middle means no contiguous *7*-bp duplex
        // exists anywhere (both strands are length 7, so qs=ts=0 is
        // the only length-7 window).
        let query = RnaSeq::parse("GGGGAAA").unwrap();
        let target = RnaSeq::parse("AAACCCC").unwrap();
        // A 4-pair seed is satisfiable.
        assert!(predict_interaction_with(
            &query,
            &target,
            InteractionParams {
                min_seed: 4,
                ..Default::default()
            },
        )
        .is_ok());
        // A 7-pair seed cannot be met.
        let r = predict_interaction_with(
            &query,
            &target,
            InteractionParams {
                min_seed: 7,
                ..Default::default()
            },
        );
        assert!(r.is_err());
    }
}
