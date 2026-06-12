//! Feature 5 — structural-RNA sequence design (inverse folding).
//!
//! Given a target secondary structure, find an RNA *sequence* that
//! folds to it — the Eterna problem. This module is a thin
//! orchestration layer over [`valenx_rnastruct::inverse_fold`]: the
//! actual inverse-folding adaptive walk and the Zuker MFE folder are
//! *that* crate's code, never re-implemented here.
//!
//! What this module adds is a **multi-start** wrapper: a single
//! adaptive walk can get stuck; running several walks from different
//! random seeds and keeping the best candidate is the standard way to
//! lift the success rate. [`design_structural`] runs `n` walks (each a
//! call to `valenx_rnastruct::design::inverse_fold_with`) and returns
//! the candidate with the smallest base-pair distance to the target.
//!
//! ## v1 scope
//!
//! Multi-start is independent restarts of the underlying walk — there
//! is no crossover or shared population between starts. The achieved
//! structure is the nearest-neighbor MFE; an inverse-fold that "solves"
//! the target solves it *for that energy model*, which is a strong
//! prediction, not a measured fold.

use crate::design::{DesignKind, RnaDesign};
use crate::error::{Result, RnaDesignError};
use valenx_rnastruct::design::inverse_fold_with;
use valenx_rnastruct::{base_pair_distance, Structure};

/// Parameters for [`design_structural`].
#[derive(Copy, Clone, Debug)]
pub struct StructuralDesignParams {
    /// Number of independent inverse-folding starts. More starts raise
    /// the chance of hitting a sequence that folds exactly to the
    /// target. Must be at least 1.
    pub starts: usize,
    /// The adaptive-walk iteration budget per start (passed straight to
    /// the `valenx-rnastruct` inverse folder).
    pub iterations_per_start: usize,
    /// The base random seed; start `k` uses `seed + k`.
    pub seed: u64,
}

impl Default for StructuralDesignParams {
    /// 8 starts, 2000 iterations each — a solid default for a v1.
    fn default() -> Self {
        StructuralDesignParams {
            starts: 8,
            iterations_per_start: 2000,
            seed: 0x5EED, // a fixed, reproducible seed
        }
    }
}

/// Designs a structural RNA that folds to `target` (feature 5).
///
/// Runs [`StructuralDesignParams::starts`] independent inverse-folding
/// walks and returns the one whose MFE fold is closest (by base-pair
/// distance) to `target`, packaged as an [`RnaDesign`].
///
/// # Errors
/// - [`RnaDesignError::Goal`] if `target` is empty or pseudoknotted.
/// - [`RnaDesignError::Invalid`] if `params.starts == 0`.
/// - [`RnaDesignError::Upstream`] if the underlying folder fails.
pub fn design_structural(target: &Structure, params: StructuralDesignParams) -> Result<RnaDesign> {
    if target.is_empty() {
        return Err(RnaDesignError::goal("target", "target structure is empty"));
    }
    if target.has_pseudoknot() {
        return Err(RnaDesignError::goal(
            "target",
            "target is pseudoknotted — the MFE-based inverse folder cannot reach it",
        ));
    }
    if params.starts == 0 {
        return Err(RnaDesignError::invalid(
            "starts",
            "need at least one inverse-folding start",
        ));
    }

    // Run the multi-start: keep the candidate with the smallest
    // base-pair distance, breaking ties toward the earlier (lower-seed)
    // start for determinism.
    let mut best: Option<(Vec<u8>, Structure, usize, usize)> = None; // (seq, achieved, dist, start)
    for k in 0..params.starts {
        let seed = params.seed.wrapping_add(k as u64);
        let r = inverse_fold_with(target, seed, params.iterations_per_start)?;
        let dist = r.distance;
        let take = match &best {
            None => true,
            Some((_, _, bd, _)) => dist < *bd,
        };
        if take {
            best = Some((r.sequence.as_bytes().to_vec(), r.achieved, dist, k));
        }
        // An exact solve cannot be beaten — stop early.
        if dist == 0 {
            break;
        }
    }

    let (sequence, achieved, distance, winning_start) =
        best.ok_or_else(|| RnaDesignError::no_design("structural", "no candidate produced"))?;

    let match_pct = structure_match_percent(&achieved, target)?;
    let mut notes = Vec::new();
    notes.push(format!(
        "Inverse-folded to a {}-nt target structure over {} multi-start walk(s).",
        target.len(),
        params.starts,
    ));
    notes.push(format!(
        "Best candidate came from start #{winning_start}: base-pair distance {distance} \
         ({match_pct:.0}% of target pairs recovered).",
    ));
    if distance == 0 {
        notes.push(
            "The design folds exactly to the target under the nearest-neighbor MFE model \
             — a strong in-silico prediction; confirm the fold experimentally."
                .to_string(),
        );
    } else {
        notes.push(format!(
            "The design does not fold exactly to the target ({distance} pair(s) differ); \
             the optimisation step can refine it further.",
        ));
    }

    Ok(RnaDesign {
        sequence,
        kind: DesignKind::Structural {
            target: target.clone(),
        },
        cds_span: None,
        construct: None,
        notes,
    })
}

/// The percentage of `target`'s base pairs that `achieved` also
/// contains — a `[0, 100]` recovery score.
fn structure_match_percent(achieved: &Structure, target: &Structure) -> Result<f64> {
    let target_pairs = target.n_pairs();
    if target_pairs == 0 {
        return Ok(100.0);
    }
    // base_pair_distance is the symmetric difference; recovered pairs =
    // target pairs not in the difference.
    let dist = base_pair_distance(achieved, target)?;
    let achieved_pairs = achieved.n_pairs();
    // pairs in exactly one structure = dist; shared = (|A|+|B|-dist)/2.
    let shared = (achieved_pairs + target_pairs).saturating_sub(dist) / 2;
    Ok(100.0 * shared as f64 / target_pairs as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn designs_a_hairpin() {
        let target = Structure::from_dot_bracket("((((....))))").unwrap();
        let d = design_structural(&target, StructuralDesignParams::default()).unwrap();
        assert_eq!(d.len(), 12);
        assert!(!d.is_coding());
        assert!(matches!(d.kind, DesignKind::Structural { .. }));
        assert!(!d.notes.is_empty());
    }

    #[test]
    fn unstructured_target_solves_trivially() {
        let target = Structure::empty(15);
        let d = design_structural(&target, StructuralDesignParams::default()).unwrap();
        // An all-unpaired target is always solvable; the design folds to it.
        let achieved = valenx_rnastruct::mfe(&d.to_rna_seq().unwrap())
            .unwrap()
            .structure;
        assert_eq!(base_pair_distance(&achieved, &target).unwrap(), 0);
    }

    #[test]
    fn rejects_pseudoknot() {
        let pk = Structure::from_dot_bracket("((..[[..))..]]").unwrap();
        assert!(design_structural(&pk, StructuralDesignParams::default()).is_err());
    }

    #[test]
    fn rejects_empty_target() {
        assert!(
            design_structural(&Structure::empty(0), StructuralDesignParams::default()).is_err()
        );
    }

    #[test]
    fn rejects_zero_starts() {
        let target = Structure::from_dot_bracket("(((...)))").unwrap();
        let params = StructuralDesignParams {
            starts: 0,
            ..StructuralDesignParams::default()
        };
        assert!(design_structural(&target, params).is_err());
    }

    #[test]
    fn is_deterministic() {
        let target = Structure::from_dot_bracket("(((...)))").unwrap();
        let a = design_structural(&target, StructuralDesignParams::default()).unwrap();
        let b = design_structural(&target, StructuralDesignParams::default()).unwrap();
        assert_eq!(a.sequence, b.sequence);
    }

    #[test]
    fn structure_match_percent_full_for_unpaired_target() {
        let target = Structure::empty(8);
        let achieved = Structure::empty(8);
        assert!((structure_match_percent(&achieved, &target).unwrap() - 100.0).abs() < 1e-9);
    }
}
