//! # valenx-rnastruct — RNA secondary-structure prediction
//!
//! Round 6 Block 5 of the Valenx roadmap. A native-Rust replacement
//! for the secondary-structure tools ViennaRNA, RNAstructure,
//! mfold / UNAFold, NUPACK, ContraFold and IPknot — pure
//! dynamic-programming algorithms, no neural-network weights and no
//! external processes.
//!
//! It builds on [`valenx_bioseq`] (Block 6.1): every folding entry
//! point accepts a [`valenx_bioseq::Seq`] (transcribing DNA to RNA as
//! needed) through the [`RnaSeq`] wrapper.
//!
//! ## What it does
//!
//! - **Model & I/O** ([`structure`], [`io`]) — the [`Structure`]
//!   base-pair model with dot-bracket I/O (including the pseudoknot
//!   bracket families `[] {} <>`), and connectivity-table (ct) and
//!   bpseq readers / writers.
//! - **Folding** ([`fold`]) — Nussinov maximum-base-pairing, a
//!   Turner-2004 nearest-neighbor energy model with explicit
//!   coaxial-stacking, Zuker minimum-free-energy folding, LinearFold
//!   linear-time beam-search folding for long sequences, free-energy
//!   evaluation, and hard- / SHAPE-soft-constrained folding.
//! - **Ensemble** ([`ensemble`]) — the McCaskill partition function
//!   and base-pair probability matrix, LinearPartition linear-time
//!   partition function for long sequences, centroid and
//!   maximum-expected-accuracy structures, Zuker-Stiegler suboptimal
//!   structures, Boltzmann stochastic sampling, melting curves, and
//!   Kinfold-class stochastic kinetic folding (Monte-Carlo
//!   Metropolis / Kawasaki walks).
//! - **Interaction** ([`interaction`]) — RNA-RNA cofolding,
//!   accessibility profiles, v1 seed-window interaction prediction,
//!   and the full IntaRNA-class accessibility-aware interaction DP
//!   (seed + extension with internal loops / bulges).
//! - **Design** ([`design`]) — Eterna-class inverse folding and
//!   G-quadruplex prediction.
//! - **Comparison** ([`compare`]) — base-pair and tree-edit
//!   structure distance, RNAforester-class structure alignment,
//!   RNAalifold-class consensus folding, restricted (H-type)
//!   pseudoknot folding, and pknotsRG-class pseudoknot folding
//!   (H-type plus kissing-hairpin).
//! - **Specialized** ([`specialized`], [`mod@layout`]) — tRNA
//!   cloverleaf detection, structure statistics, mountain plots,
//!   energy dot plots, 2-D drawing layout, and batch folding reports.
//!
//! ## Errors
//!
//! Every fallible function returns
//! [`Result<_, RnaStructError>`](error::RnaStructError). The error
//! type carries stable [`code`](error::RnaStructError::code) and
//! [`category`](error::RnaStructError::category) accessors for
//! telemetry.
//!
//! ## v1 scope
//!
//! This is a real working v1, not production parity with the 30-year
//! reference tools. Each module documents its own simplifications
//! inline; the notable ones are:
//!
//! - **The Turner-2004 parameter set is the complete published set.**
//!   [`fold::turner2004`] is a verbatim transcription of the published
//!   Turner-2004 nearest-neighbor parameters — the same numbers
//!   ViennaRNA ships in `rna_turner2004.par`: the full 4×4 stacking
//!   table; the complete hairpin / bulge / internal-loop length tables
//!   with the Jacobson-Stockmayer logarithmic extrapolation; the
//!   published triloop / tetraloop small-loop special cases; the full
//!   per-closing-pair terminal-mismatch tables for hairpins and
//!   interior loops; the explicit 1×1 internal-loop energies; the
//!   `dangle5` / `dangle3` dangling-end tables; the linear multiloop
//!   model; the terminal-AU/GU helix-end penalty; and the explicit
//!   **coaxial-stacking** term ([`fold::coaxial`]). Folding energies
//!   reproduce the analytic Turner sum exactly;
//!   [`fold::eval::structure_energy_d2`] adds the coaxial-stacking
//!   correction and reproduces ViennaRNA's default `-d2` `RNAeval`
//!   exactly-to-rounding for any given structure.
//! - **LinearFold / LinearPartition** ([`fold::linear`],
//!   [`ensemble::linear_partition`]) — linear-time beam-search folding
//!   and partition function for long sequences. Beam search is an
//!   *approximate* algorithm: with a wide enough beam it reproduces
//!   the exact Zuker MFE / McCaskill result, but a narrow beam may
//!   miss the global optimum. The exact `O(n³)` Zuker / McCaskill DPs
//!   stay the default for short RNA.
//! - The default Zuker MFE structure search optimises the
//!   dangle-folded multiloop model; coaxial stacking is applied as an
//!   exact energy *re-scoring* ([`fold::zuker::mfe_d2`]). Folding the coaxial
//!   term into the recurrence itself is a documented follow-up.
//! - The Zuker / McCaskill DPs are pseudoknot-free; the
//!   [`compare::pseudoknot`] folder adds only the restricted H-type
//!   class.
//! - The McCaskill base-pair probabilities fold the multiloop channel
//!   into the exterior / internal-loop outside recursion as a
//!   first-order treatment.
//! - Inverse folding is an adaptive walk (not exhaustive); the
//!   2-D layout omits the force-directed overlap-removal pass.

#![forbid(unsafe_code)]

pub mod compare;
pub mod design;
pub mod ensemble;
pub mod error;
pub mod fold;
pub mod interaction;
pub mod io;
pub mod layout;
pub mod rna;
pub mod specialized;
pub mod structure;

// --- Convenience re-exports of the most-used types --------------------

pub use error::{ErrorCategory, Result, RnaStructError};
pub use rna::RnaSeq;
pub use structure::{BasePair, Structure};

pub use fold::constraint::FoldConstraints;
pub use fold::coaxial::{best_coaxial, HelixEnd};
pub use fold::eval::{coaxial_correction, structure_energy, structure_energy_d2};
pub use fold::linear::{
    fold_linear, fold_linear_exact, fold_linear_with_beam, LinearFoldResult,
    DEFAULT_BEAM_SIZE,
};
pub use fold::nussinov::{fold as nussinov_fold, NussinovResult};
pub use fold::shape::fold_with_shape;
pub use fold::zuker::{mfe, mfe_constrained, mfe_d2, MfeResult};

pub use ensemble::centroid::{centroid, mea, MeaResult};
pub use ensemble::kinetics::{
    fold_kinetics, simulate_trajectory, KineticEnsemble, KineticParams, RateModel,
    Trajectory, TrajectoryStep,
};
pub use ensemble::linear_partition::{
    linear_partition, linear_partition_exact, linear_partition_with_beam,
    LinearPartitionResult,
};
pub use ensemble::melting::{melting_curve, MeltingCurve};
pub use ensemble::partition::{partition_function, PartitionFunction};
pub use ensemble::sampling::{sample, sample_with_counts};
pub use ensemble::suboptimal::{suboptimal, SuboptStructure};

pub use interaction::accessibility::{accessibility, AccessibilityProfile};
pub use interaction::cofold::{cofold, CofoldResult};
pub use interaction::intarna::{
    predict_intarna, predict_intarna_with, InterPair, IntaRnaInteraction,
    IntaRnaParams, DEFAULT_MAX_LEN, DEFAULT_SEED_MIN, IL_MAX,
};
pub use interaction::interaction::{predict_interaction, Interaction};

pub use compare::align::{align_structures, StructureAlignment};
pub use compare::consensus::{consensus_structure, ConsensusResult};
pub use compare::distance::{base_pair_distance, tree_edit_distance};
pub use compare::pknots_rg::{
    fold_pknots_rg, fold_pknots_rg_with, PknotsRgParams, PknotsRgResult,
    PseudoknotClass, KISSING_HAIRPIN_PENALTY,
};
pub use compare::pseudoknot::{fold_pseudoknot, PseudoknotResult};

pub use design::{inverse_fold, predict_gquadruplex, GQuadruplex, InverseFoldResult};
pub use io::{read_bpseq, read_ct, write_bpseq, write_ct, StructureRecord};
pub use layout::{layout, Layout};
pub use specialized::dotplot::{mfe_dot_plot, probability_dot_plot, DotPlot};
pub use specialized::mountain::{mountain_plot, MountainPlot};
pub use specialized::report::{batch_fold, folding_report, FoldingReport};
pub use specialized::stats::{structure_stats, StructureStats};
pub use specialized::trna::{scan_trna, TrnaScan};

#[cfg(test)]
mod tests {
    use super::*;

    /// End-to-end: fold a sequence, characterise its ensemble, and
    /// confirm the derived structures are mutually consistent.
    #[test]
    fn fold_and_ensemble_end_to_end() {
        let seq = RnaSeq::parse("GGGGGGGAAAACCCCCCC").unwrap();

        // MFE fold.
        let mfe_result = mfe(&seq).unwrap();
        assert!(mfe_result.structure.is_nested());
        assert!(mfe_result.energy < 0.0);

        // The reported MFE equals an independent evaluation.
        let re = structure_energy(&seq, &mfe_result.structure).unwrap();
        assert!((re - mfe_result.energy).abs() < 1e-4);

        // Ensemble: G <= E_mfe.
        let pf = partition_function(&seq).unwrap();
        assert!(pf.ensemble_free_energy() <= mfe_result.energy + 1e-6);

        // Centroid / MEA are valid structures of the right length.
        let cen = centroid(&pf).unwrap();
        let m = mea(&pf).unwrap();
        assert_eq!(cen.len(), seq.len());
        assert_eq!(m.structure.len(), seq.len());

        // The full report ties it together.
        let report = folding_report(&seq).unwrap();
        assert_eq!(report.length, seq.len());
        assert!((0.0..=1.0).contains(&report.mfe_frequency));
    }

    /// End-to-end: round-trip a structure through dot-bracket and ct.
    #[test]
    fn structure_io_round_trip() {
        let seq = RnaSeq::parse("GGGAAACCC").unwrap();
        let s = mfe(&seq).unwrap().structure;

        // dot-bracket round-trip
        let db = s.to_dot_bracket();
        let back = Structure::from_dot_bracket(&db).unwrap();
        assert_eq!(s.pairs(), back.pairs());

        // ct round-trip
        let rec = StructureRecord::new("demo", seq.clone(), s.clone()).unwrap();
        let ct = write_ct(&rec);
        let parsed = read_ct(&ct).unwrap();
        assert_eq!(parsed.structure, s);
    }

    /// End-to-end: design a sequence for a target then fold it back.
    #[test]
    fn inverse_fold_then_fold_back() {
        let target = Structure::from_dot_bracket("((((....))))").unwrap();
        let designed = inverse_fold(&target, 1).unwrap();
        // the designed sequence folds to something close to the target
        let refold = mfe(&designed.sequence).unwrap().structure;
        let d = base_pair_distance(&refold, &target).unwrap();
        assert!(d <= 4, "inverse fold + refold drifted by {d}");
    }

    #[test]
    fn re_exports_are_wired() {
        // Touch a sampling of the convenience re-exports.
        let seq = RnaSeq::parse("GGGGCCCC").unwrap();
        let _ = nussinov_fold(&seq).unwrap();
        let _ = FoldConstraints::none(seq.len());
        let e = RnaStructError::not_yet("x");
        assert_eq!(e.category_enum(), ErrorCategory::Capability);
        let s = Structure::from_dot_bracket("((....))").unwrap();
        let _ = mfe_dot_plot(&s);
        let _ = mountain_plot(&s);
        let _ = structure_stats(&s);
        let _ = layout(&s);
        let _ = BasePair::new(0, 7).unwrap();
    }
}
