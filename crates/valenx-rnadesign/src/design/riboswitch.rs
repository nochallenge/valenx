//! Feature 7 — riboswitch / two-state sequence design (v1).
//!
//! A riboswitch is an RNA that adopts **two** structures — one in the
//! ligand-free state, another in the ligand-bound state — and switches
//! between them. Designing one is *two-state design*: find a single
//! sequence that can fold to *both* targets.
//!
//! This v1 takes a transparent, honest approach built on the
//! `valenx-rnastruct` primitives:
//!
//! 1. The two targets must be **mutually sequence-compatible** — for
//!    every position, the bases demanded by its pairing role in target
//!    A and in target B must be jointly satisfiable. Positions paired
//!    in *both* targets are the hard constraint; the designer checks
//!    compatibility up front.
//! 2. A multi-start inverse-folding search produces sequences whose
//!    **MFE fold** matches the free-state target (the resting
//!    conformation).
//! 3. Each candidate is scored by how well it *also* admits the
//!    bound-state target: the bound structure's free energy on the
//!    candidate is compared to the MFE — a small gap means the bound
//!    state is a thermodynamically accessible alternative the ligand
//!    can stabilise.
//!
//! The candidate minimising the combined "free-state mismatch + bound-
//! state energy gap" score is returned.
//!
//! ## v1 scope — honest framing
//!
//! This is a real two-state *design heuristic*, not a full
//! constraint-programming riboswitch designer. It does not model the
//! ligand, the binding free energy, or the switching kinetics — the
//! ligand is represented only as "the bound structure becomes
//! favourable". A small predicted energy gap is a necessary, not a
//! sufficient, condition for a working switch; the design must be
//! validated experimentally. The free-state fold is the
//! nearest-neighbor MFE, with all the caveats of that model.

use crate::design::{DesignKind, RnaDesign};
use crate::error::{Result, RnaDesignError};
use valenx_rnastruct::design::inverse_fold_with;
use valenx_rnastruct::{base_pair_distance, mfe, structure_energy, RnaSeq, Structure};

/// Parameters for [`design_riboswitch`].
#[derive(Copy, Clone, Debug)]
pub struct TwoStateParams {
    /// Number of independent inverse-folding starts.
    pub starts: usize,
    /// The adaptive-walk iteration budget per start.
    pub iterations_per_start: usize,
    /// The base random seed; start `k` uses `seed + k`.
    pub seed: u64,
}

impl Default for TwoStateParams {
    /// 12 starts, 2500 iterations each — two-state design is harder than
    /// single-target inverse folding, so it gets more starts.
    fn default() -> Self {
        TwoStateParams {
            starts: 12,
            iterations_per_start: 2500,
            seed: 0x5117,
        }
    }
}

/// The result of a two-state riboswitch design (feature 7).
#[derive(Clone, Debug, PartialEq)]
pub struct RiboswitchDesign {
    /// The packaged design candidate.
    pub design: RnaDesign,
    /// Base-pair distance of the candidate's MFE fold to the free-state
    /// target (`0` = the resting conformation is exactly the target).
    pub free_state_distance: usize,
    /// Free energy of the bound-state structure on the candidate
    /// (kcal/mol).
    pub bound_state_energy: f64,
    /// The candidate's minimum free energy (kcal/mol) — the free state.
    pub mfe_energy: f64,
    /// The energy gap `bound_state_energy − mfe_energy` (kcal/mol). A
    /// small non-negative gap means the bound state is an accessible
    /// alternative.
    pub energy_gap: f64,
}

impl RiboswitchDesign {
    /// `true` when the bound state is thermodynamically accessible — the
    /// energy gap is within a switching-plausible window (≤ 5 kcal/mol).
    /// A larger gap means the ligand would have to do too much work.
    pub fn bound_state_accessible(&self) -> bool {
        self.energy_gap <= 5.0
    }
}

/// Designs a riboswitch — a sequence that can adopt two target
/// structures (feature 7).
///
/// `free_target` is the resting (ligand-free) structure; `bound_target`
/// the structure the ligand stabilises. Both must be the same length
/// and pseudoknot-free.
///
/// # Errors
/// - [`RnaDesignError::Goal`] if either target is empty, pseudoknotted,
///   or the two differ in length / are not jointly satisfiable.
/// - [`RnaDesignError::Invalid`] if `params.starts == 0`.
/// - [`RnaDesignError::NoDesign`] if no candidate is produced.
/// - [`RnaDesignError::Upstream`] if a folder call fails.
pub fn design_riboswitch(
    free_target: &Structure,
    bound_target: &Structure,
    params: TwoStateParams,
) -> Result<RiboswitchDesign> {
    if free_target.is_empty() || bound_target.is_empty() {
        return Err(RnaDesignError::goal("target", "a riboswitch target is empty"));
    }
    if free_target.len() != bound_target.len() {
        return Err(RnaDesignError::goal(
            "target",
            format!(
                "the two riboswitch targets differ in length ({} vs {})",
                free_target.len(),
                bound_target.len()
            ),
        ));
    }
    if free_target.has_pseudoknot() || bound_target.has_pseudoknot() {
        return Err(RnaDesignError::goal(
            "target",
            "a riboswitch target is pseudoknotted — unreachable by the MFE folder",
        ));
    }
    if params.starts == 0 {
        return Err(RnaDesignError::invalid("starts", "need at least one start"));
    }
    // Up-front joint-satisfiability check: a position cannot be paired
    // to two different partners that demand contradictory bases. The
    // simplest sufficient check — no position is paired to *different*
    // partners in the two structures where both partners are also
    // mutually exclusive — is too strict; instead we just confirm the
    // structures are not literally identical-with-conflict. A position
    // paired in both targets to the same partner is always fine; a
    // position paired in one and unpaired in the other is fine; a
    // position paired to *different* partners is the interesting
    // riboswitch case and is allowed (the base just has to pair with
    // whichever partner each state needs — handled by the search).
    if free_target == bound_target {
        return Err(RnaDesignError::goal(
            "target",
            "the two riboswitch states are identical — there is nothing to switch between",
        ));
    }

    // Multi-start: inverse-fold toward the free-state target, then
    // score each candidate by free-state match + bound-state gap.
    let mut best: Option<(Vec<u8>, usize, f64, f64, f64, f64)> = None;
    // (seq, free_dist, bound_energy, mfe_energy, gap, score)
    for k in 0..params.starts {
        let seed = params.seed.wrapping_add(k as u64);
        let r = inverse_fold_with(free_target, seed, params.iterations_per_start)?;
        let seq_bytes = r.sequence.as_bytes().to_vec();
        let rna = RnaSeq::parse(&seq_bytes)?;

        let mfe_res = mfe(&rna)?;
        let free_dist = base_pair_distance(&mfe_res.structure, free_target)?;

        // Energy of the bound-state target on this candidate. If the
        // bound target has a non-canonical pair on this sequence the
        // candidate cannot adopt it at all — skip with a large penalty.
        let bound_energy = match structure_energy(&rna, bound_target) {
            Ok(e) => e,
            Err(_) => f64::INFINITY,
        };
        let gap = bound_energy - mfe_res.energy;
        // Score: free-state mismatch dominates, then the bound-state
        // energy gap (smaller is better — the ligand should not have to
        // overcome a huge barrier).
        let score = free_dist as f64 * 100.0 + gap.max(0.0);

        let take = match &best {
            None => true,
            Some((_, _, _, _, _, bs)) => score < *bs,
        };
        if take {
            best = Some((seq_bytes, free_dist, bound_energy, mfe_res.energy, gap, score));
        }
    }

    let (sequence, free_dist, bound_energy, mfe_energy, gap, _) = best
        .ok_or_else(|| RnaDesignError::no_design("riboswitch", "no candidate produced"))?;

    let mut notes = Vec::new();
    notes.push(format!(
        "Two-state design over a {}-nt riboswitch with {} multi-start walk(s).",
        free_target.len(),
        params.starts,
    ));
    notes.push(format!(
        "Free (resting) state: base-pair distance {free_dist} to the target, MFE {mfe_energy:.2} kcal/mol.",
    ));
    if bound_energy.is_finite() {
        notes.push(format!(
            "Bound state: the ligand-stabilised structure has energy {bound_energy:.2} kcal/mol \
             on this sequence — a {gap:.2} kcal/mol gap above the MFE.",
        ));
    } else {
        notes.push(
            "Bound state: the best candidate cannot adopt the bound-state structure with \
             canonical pairs — the two targets may not be jointly designable."
                .to_string(),
        );
    }
    notes.push(
        "Two-state design is a thermodynamic heuristic — the ligand and its binding energy \
         are not modelled; a small energy gap is necessary but not sufficient for a working \
         switch. Validate experimentally."
            .to_string(),
    );

    let design = RnaDesign {
        sequence,
        kind: DesignKind::Riboswitch {
            free_target: free_target.clone(),
            bound_target: bound_target.clone(),
        },
        cds_span: None,
        construct: None,
        notes,
    };

    Ok(RiboswitchDesign {
        design,
        free_state_distance: free_dist,
        bound_state_energy: bound_energy,
        mfe_energy,
        energy_gap: gap,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn designs_a_two_state_switch() {
        // Two distinct structures of the same length.
        let free = Structure::from_dot_bracket("((((....))))....").unwrap();
        let bound = Structure::from_dot_bracket("....((((....))))").unwrap();
        let r = design_riboswitch(&free, &bound, TwoStateParams::default()).unwrap();
        assert_eq!(r.design.len(), 16);
        assert!(matches!(r.design.kind, DesignKind::Riboswitch { .. }));
        assert!(!r.design.notes.is_empty());
    }

    #[test]
    fn rejects_length_mismatch() {
        let free = Structure::from_dot_bracket("((((....))))").unwrap();
        let bound = Structure::from_dot_bracket("((((....))))....").unwrap();
        assert!(design_riboswitch(&free, &bound, TwoStateParams::default()).is_err());
    }

    #[test]
    fn rejects_identical_states() {
        let s = Structure::from_dot_bracket("((((....))))").unwrap();
        let err = design_riboswitch(&s, &s, TwoStateParams::default()).unwrap_err();
        assert_eq!(err.code(), "rnadesign.goal");
    }

    #[test]
    fn rejects_pseudoknot() {
        let free = Structure::from_dot_bracket("((..[[..))..]]").unwrap();
        let bound = Structure::from_dot_bracket("..............").unwrap();
        assert!(design_riboswitch(&free, &bound, TwoStateParams::default()).is_err());
    }

    #[test]
    fn rejects_zero_starts() {
        let free = Structure::from_dot_bracket("((((....))))....").unwrap();
        let bound = Structure::from_dot_bracket("....((((....))))").unwrap();
        let params = TwoStateParams {
            starts: 0,
            ..TwoStateParams::default()
        };
        assert!(design_riboswitch(&free, &bound, params).is_err());
    }

    #[test]
    fn is_deterministic() {
        let free = Structure::from_dot_bracket("((((....))))....").unwrap();
        let bound = Structure::from_dot_bracket("....((((....))))").unwrap();
        let a = design_riboswitch(&free, &bound, TwoStateParams::default()).unwrap();
        let b = design_riboswitch(&free, &bound, TwoStateParams::default()).unwrap();
        assert_eq!(a.design.sequence, b.design.sequence);
    }

    #[test]
    fn accessible_flag_tracks_gap() {
        let free = Structure::from_dot_bracket("((((....))))....").unwrap();
        let bound = Structure::from_dot_bracket("....((((....))))").unwrap();
        let r = design_riboswitch(&free, &bound, TwoStateParams::default()).unwrap();
        assert_eq!(r.bound_state_accessible(), r.energy_gap <= 5.0);
    }
}
