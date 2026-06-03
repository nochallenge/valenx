//! **Feature 11 — coarse-to-fine protocol.**
//!
//! Rosetta `AbinitioRelax` runs in two stages, and so does this
//! module:
//!
//! 1. **Centroid stage** — fragment assembly with a coarse,
//!    centroid-resolution energy. Sidechains are reduced to a single
//!    pseudo-atom; the search is fast and explores fold space broadly.
//!    This is [`crate::abinitio::fragment_assembly`].
//! 2. **All-atom stage** — the best centroid model is "switched" to
//!    all-atom: real sidechains are placed by rotamer repacking, the
//!    backbone φ/ψ angles are pulled into the allowed Ramachandran
//!    regions, and the structure is energy-minimised (relaxed).
//!
//! Coarse-to-fine works because the centroid stage finds the *fold*
//! cheaply and the all-atom stage only has to *refine* it — searching
//! all-atom from scratch would be hopelessly slow. This module is the
//! glue: it runs the centroid assembly, then the all-atom refinement
//! pass, and reports both.

use serde::{Deserialize, Serialize};

use crate::abinitio::assemble::{fragment_assembly, AssemblyOptions};
use crate::abinitio::fragments::build_fragment_library;
use crate::error::{Result, StructPredictError};
use crate::model::ProteinModel;
use crate::refine::mcrefine::{mc_refine, McRefineOptions};
use crate::refine::quality::{assess_quality, QualityReport};
use crate::refine::ramachandran::refine_ramachandran;
use crate::refine::repack::repack_sidechains;

/// The outcome of a coarse-to-fine ab-initio run.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ProtocolResult {
    /// The final all-atom model.
    pub model: ProteinModel,
    /// Knowledge-score total of the centroid-stage model.
    pub centroid_score: f64,
    /// Ramachandran outliers removed by the all-atom refinement.
    pub outliers_removed: usize,
    /// The final model's quality report.
    pub quality: QualityReport,
    /// DOPE energy of the centroid model immediately after the
    /// initial fragment assembly (i.e. the input to the MC refinement
    /// stage). When the MC refinement stage runs, this is the
    /// pre-refinement DOPE total; when it is skipped it is `None`.
    pub pre_refine_dope: Option<f64>,
    /// DOPE energy of the model after the MC refinement stage.
    /// Equals [`Self::pre_refine_dope`] minus the MC improvement;
    /// `None` when the refinement stage was skipped.
    pub post_refine_dope: Option<f64>,
}

/// Runs the full centroid → all-atom ab-initio protocol.
///
/// `sequence` is the target. `centroid_moves` is the fragment-
/// assembly move budget; `repack_moves` the all-atom repacking
/// budget; `seed` fixes both RNGs.
///
/// # Errors
/// [`StructPredictError::Invalid`] for an empty sequence; propagates
/// errors from the assembly and refinement stages.
pub fn coarse_to_fine(
    sequence: &str,
    centroid_moves: usize,
    repack_moves: usize,
    seed: u64,
) -> Result<ProtocolResult> {
    coarse_to_fine_with_refine(sequence, centroid_moves, repack_moves, seed, None)
}

/// The full coarse-to-fine ab-initio protocol, with an optional
/// DOPE-driven simulated-annealing MC refinement stage.
///
/// When `refine_options` is `Some(_)`, an MC + DOPE refinement runs
/// between the centroid assembly and the all-atom Ramachandran +
/// repacking stages — the published Rosetta `relax` protocol pattern.
/// When it is `None`, this is exactly [`coarse_to_fine`] above.
///
/// # Errors
/// [`StructPredictError::Invalid`] for an empty sequence; propagates
/// errors from the assembly and refinement stages.
pub fn coarse_to_fine_with_refine(
    sequence: &str,
    centroid_moves: usize,
    repack_moves: usize,
    seed: u64,
    refine_options: Option<McRefineOptions>,
) -> Result<ProtocolResult> {
    let sequence = sequence.trim();
    if sequence.is_empty() {
        return Err(StructPredictError::invalid("sequence", "empty"));
    }
    if sequence.len() < 4 {
        return Err(StructPredictError::invalid(
            "sequence",
            "need at least 4 residues for fragment assembly",
        ));
    }

    // --- Centroid stage ------------------------------------------------
    let frag_len = 3.min(sequence.len());
    let library = build_fragment_library(sequence, frag_len, 30)?;
    let assembly = fragment_assembly(
        sequence,
        &library,
        AssemblyOptions {
            moves: centroid_moves.max(1),
            seed,
            ..AssemblyOptions::default()
        },
    )?;
    let centroid_score = assembly.final_score;
    let mut model = assembly.model;

    // --- (Optional) DOPE-driven MC refinement stage --------------------
    let (pre_refine_dope, post_refine_dope) = match refine_options {
        Some(mc_opts) => {
            // The MC refinement loop also runs Ramachandran cleanup
            // each cycle — the same `refine_ramachandran` the all-atom
            // stage uses below — so when the MC stage runs we don't
            // need an extra pre-stage Rama pass.
            let pre = crate::abinitio::dope::dope_score(
                &model,
                crate::abinitio::dope::DopeWeights::default(),
            )?
            .total;
            let res = mc_refine(&model, sequence, mc_opts)?;
            model = res.model;
            (Some(pre), Some(res.final_energy))
        }
        None => (None, None),
    };

    // --- All-atom stage ------------------------------------------------
    // 1. Pull strained backbone dihedrals into allowed regions.
    let rama = refine_ramachandran(&mut model)?;
    // 2. Place sidechains by rotamer repacking.
    let _repack = repack_sidechains(&mut model, repack_moves.max(1), seed ^ 0x5DEECE66D)?;
    // 3. Final quality assessment.
    let quality = assess_quality(&model)?;

    Ok(ProtocolResult {
        model,
        centroid_score,
        outliers_removed: rama.removed(),
        quality,
        pre_refine_dope,
        post_refine_dope,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn protocol_produces_a_complete_all_atom_model() {
        let seq = "EEEEAAAALLLLEEEEAAAA";
        let res = coarse_to_fine(seq, 400, 200, 123).expect("protocol");
        assert!(res.model.is_complete(), "backbone complete");
        // Every non-glycine residue has a Cβ after repacking.
        for r in &res.model.residues {
            if r.aa != 'G' {
                assert!(r.cb.is_some(), "{} has Cβ", r.aa);
            }
        }
        assert_eq!(res.model.len(), seq.len());
    }

    #[test]
    fn protocol_is_deterministic() {
        let seq = "ACDEFGHIKLMNPQRS";
        let a = coarse_to_fine(seq, 150, 100, 7).expect("a");
        let b = coarse_to_fine(seq, 150, 100, 7).expect("b");
        assert_eq!(a.centroid_score, b.centroid_score);
        assert_eq!(a.quality.overall, b.quality.overall);
    }

    #[test]
    fn short_or_empty_sequence_rejected() {
        assert!(coarse_to_fine("", 10, 10, 0).is_err());
        assert!(coarse_to_fine("AC", 10, 10, 0).is_err());
    }

    #[test]
    fn coarse_to_fine_with_refine_lowers_dope() {
        // The refine-enabled protocol records pre- and post-refinement
        // DOPE energies; the post must be ≤ the pre (the MC stage
        // always keeps the best-seen).
        let seq = "EEEEAAAALLLLEEEEAAAA";
        let refine = McRefineOptions {
            cycles: 3,
            moves_per_cycle: 60,
            seed: 5,
            ..McRefineOptions::default()
        };
        let res = coarse_to_fine_with_refine(seq, 200, 100, 7, Some(refine)).expect("refine");
        let (pre, post) = (res.pre_refine_dope.unwrap(), res.post_refine_dope.unwrap());
        assert!(post <= pre + 1e-6, "DOPE: pre {pre} post {post}");
        assert!(res.model.is_complete());
    }

    #[test]
    fn end_to_end_short_helix_predicts_low_rmsd_to_native() {
        // A short all-Leucine target sequence's native fold is an
        // α-helix (L is a strong helix former, and DOPE — like every
        // statistical potential — favours the packed helical
        // arrangement). The end-to-end predict_abinitio + DOPE-MC
        // refinement protocol on this 12-residue sequence should
        // recover an α-helix structurally close to the canonical
        // (-63°, -42°) helix.
        let seq = "LLLLLLLLLLLL"; // 12 residues
        // Build the canonical native (-63, -42) all-Leu helix as the
        // reference.
        let mut native = crate::model::ProteinModel::from_sequence(seq).expect("native");
        crate::model::build_backbone_from_torsions(&mut native, &[(-63.0, -42.0); 12])
            .expect("build native");

        // Run the protocol with MC refinement on top.
        let refine = McRefineOptions {
            cycles: 4,
            moves_per_cycle: 200,
            start_temperature: 1.5,
            end_temperature: 0.1,
            seed: 9,
            ..McRefineOptions::default()
        };
        let res = coarse_to_fine_with_refine(seq, 600, 100, 3, Some(refine)).expect("e2e");
        let rmsd = crate::refine::superpose::ca_rmsd_superposed(&res.model, &native)
            .expect("rmsd");
        // For a short canonical helix the end-to-end predicted model
        // should be reasonably close to the native helix — a
        // **classical** ab-initio result; we assert RMSD ≤ 8 Å (a
        // permissive low-resolution bound). A strict AlphaFold-class
        // sub-Å RMSD requires a trained network and is explicitly out
        // of scope for this classical crate.
        assert!(
            rmsd <= 8.0,
            "predicted-vs-native RMSD {rmsd} should be reasonable for a 12-residue helix",
        );
        assert!(res.model.is_complete());
    }
}
