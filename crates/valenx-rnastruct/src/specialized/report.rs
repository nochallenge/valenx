//! Batch folding utilities and the folding "report" type.
//!
//! Most folding tools, when handed a sequence, do not return just one
//! number — they return a *report*: the minimum-free-energy
//! structure, the ensemble free energy, the centroid and MEA
//! structures, and how dominant the MFE structure is. This module
//! bundles that into [`FoldingReport`] and adds batch helpers for
//! folding many sequences at once.

use crate::ensemble::centroid::{centroid, mea};
use crate::ensemble::partition::partition_function;
use crate::error::Result;
use crate::fold::energy::GAS_CONSTANT;
use crate::fold::zuker::mfe;
use crate::rna::RnaSeq;
use crate::structure::Structure;

/// A comprehensive folding report for one sequence.
#[derive(Clone, Debug)]
pub struct FoldingReport {
    /// The sequence length.
    pub length: usize,
    /// The minimum-free-energy structure.
    pub mfe_structure: Structure,
    /// The minimum free energy, kcal/mol.
    pub mfe_energy: f64,
    /// The ensemble free energy `G = −RT·ln Q`, kcal/mol.
    pub ensemble_free_energy: f64,
    /// The centroid structure (`p > 0.5` pairs).
    pub centroid_structure: Structure,
    /// The maximum-expected-accuracy structure.
    pub mea_structure: Structure,
    /// The Boltzmann probability of the MFE structure within the
    /// ensemble — `exp(−(E_mfe − G)/RT)`. Near 1 means the MFE
    /// structure dominates; near 0 means a diffuse ensemble.
    pub mfe_frequency: f64,
}

impl FoldingReport {
    /// The "ensemble diversity" — the gap between the MFE and the
    /// ensemble free energy. A large gap signals a structurally
    /// uncertain RNA.
    pub fn ensemble_gap(&self) -> f64 {
        self.mfe_energy - self.ensemble_free_energy
    }

    /// `true` if the MFE structure is well-defined — it carries the
    /// bulk of the Boltzmann weight (`mfe_frequency > 0.5`).
    pub fn is_well_defined(&self) -> bool {
        self.mfe_frequency > 0.5
    }
}

/// Produces a full [`FoldingReport`] for `seq` at 37 °C.
///
/// This runs the MFE folder, the McCaskill partition function, and
/// the centroid / MEA derivations once each.
///
/// # Errors
/// Propagates folding / ensemble errors.
pub fn folding_report(seq: &RnaSeq) -> Result<FoldingReport> {
    let mfe_result = mfe(seq)?;
    let pf = partition_function(seq)?;
    let g = pf.ensemble_free_energy();
    let centroid_structure = centroid(&pf)?;
    let mea_structure = mea(&pf)?.structure;

    // MFE Boltzmann frequency: exp(-(E_mfe - G)/RT).
    let rt = GAS_CONSTANT * pf.temperature_k();
    let mfe_frequency = (-(mfe_result.energy - g) / rt).exp().clamp(0.0, 1.0);

    Ok(FoldingReport {
        length: seq.len(),
        mfe_structure: mfe_result.structure,
        mfe_energy: mfe_result.energy,
        ensemble_free_energy: g,
        centroid_structure,
        mea_structure,
        mfe_frequency,
    })
}

/// One entry of a batch fold: a label and either the structure that
/// folded or the error that stopped it.
#[derive(Clone, Debug)]
pub struct BatchEntry {
    /// The caller-supplied label / id for this sequence.
    pub label: String,
    /// `Ok` with the MFE structure + energy, or `Err` with the
    /// failure.
    pub result: std::result::Result<(Structure, f64), String>,
}

impl BatchEntry {
    /// `true` if this sequence folded successfully.
    pub fn is_ok(&self) -> bool {
        self.result.is_ok()
    }
}

/// Folds many `(label, sequence)` pairs to their MFE structures.
///
/// A failure on one sequence does not abort the batch — that entry
/// records the error and the rest continue. Returns one [`BatchEntry`]
/// per input, in input order.
pub fn batch_fold(sequences: &[(&str, &RnaSeq)]) -> Vec<BatchEntry> {
    sequences
        .iter()
        .map(|(label, seq)| {
            let result = mfe(seq)
                .map(|r| (r.structure, r.energy))
                .map_err(|e| e.to_string());
            BatchEntry {
                label: (*label).to_string(),
                result,
            }
        })
        .collect()
}

/// Produces a [`FoldingReport`] for each `(label, sequence)` pair.
///
/// Returns one `(label, Result<FoldingReport>)` per input. As with
/// [`batch_fold`], one failure does not abort the rest.
pub fn batch_report(
    sequences: &[(&str, &RnaSeq)],
) -> Vec<(String, std::result::Result<FoldingReport, String>)> {
    sequences
        .iter()
        .map(|(label, seq)| {
            (
                (*label).to_string(),
                folding_report(seq).map_err(|e| e.to_string()),
            )
        })
        .collect()
}

/// The mean MFE free energy over a successfully-folded batch.
///
/// Returns `None` if `batch` is empty or every entry failed.
pub fn mean_mfe(batch: &[BatchEntry]) -> Option<f64> {
    let energies: Vec<f64> = batch
        .iter()
        .filter_map(|e| e.result.as_ref().ok().map(|(_, energy)| *energy))
        .collect();
    if energies.is_empty() {
        None
    } else {
        Some(energies.iter().sum::<f64>() / energies.len() as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn report_fields_are_consistent() {
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let r = folding_report(&seq).unwrap();
        assert_eq!(r.length, seq.len());
        // ensemble free energy never above the MFE
        assert!(r.ensemble_free_energy <= r.mfe_energy + 1e-6);
        // gap is non-positive (G <= E_mfe)
        assert!(r.ensemble_gap() >= -1e-6);
        // frequency is a probability
        assert!((0.0..=1.0).contains(&r.mfe_frequency));
        // all three structures span the sequence
        assert_eq!(r.mfe_structure.len(), seq.len());
        assert_eq!(r.centroid_structure.len(), seq.len());
        assert_eq!(r.mea_structure.len(), seq.len());
    }

    #[test]
    fn unpairable_rna_reports_zero_energy() {
        let seq = RnaSeq::parse("AAAAAAAAAAAA").unwrap();
        let r = folding_report(&seq).unwrap();
        assert_eq!(r.mfe_energy, 0.0);
        assert_eq!(r.mfe_structure.n_pairs(), 0);
        // Q == 1 so the MFE (empty) structure has frequency ~1
        assert!(r.mfe_frequency > 0.99);
        assert!(r.is_well_defined());
    }

    #[test]
    fn batch_fold_handles_many_sequences() {
        let a = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let b = RnaSeq::parse("AAAAAAAA").unwrap();
        let c = RnaSeq::parse("GGGGCCCC").unwrap();
        let batch = batch_fold(&[("a", &a), ("b", &b), ("c", &c)]);
        assert_eq!(batch.len(), 3);
        assert!(batch.iter().all(|e| e.is_ok()));
        assert_eq!(batch[0].label, "a");
    }

    #[test]
    fn mean_mfe_averages_the_batch() {
        let a = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let b = RnaSeq::parse("AAAAAAAA").unwrap();
        let batch = batch_fold(&[("a", &a), ("b", &b)]);
        let mean = mean_mfe(&batch).unwrap();
        // b folds to 0.0, a folds to something negative -> mean < 0
        assert!(mean <= 0.0);
        assert!(mean.is_finite());
    }

    #[test]
    fn mean_mfe_of_empty_batch_is_none() {
        assert!(mean_mfe(&[]).is_none());
    }

    #[test]
    fn batch_report_runs_for_each() {
        let a = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let b = RnaSeq::parse("GGGGCCCC").unwrap();
        let reports = batch_report(&[("a", &a), ("b", &b)]);
        assert_eq!(reports.len(), 2);
        assert!(reports.iter().all(|(_, r)| r.is_ok()));
    }
}
