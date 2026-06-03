//! Feature 20 — mRNA secondary-structure check and minimisation.
//!
//! Stable secondary structure across the **start-codon region** blocks
//! the scanning ribosome and depresses translation; a structured
//! 5′UTR is a classic cause of poor expression. Conversely, *some*
//! structure in the CDS body raises mRNA stability. A good mRNA design
//! therefore keeps the start-codon neighbourhood **open** while not
//! over-destabilising the rest.
//!
//! This module reuses [`valenx_rnastruct`]'s Zuker minimum-free-energy
//! folder ([`valenx_rnastruct::mfe`]) — the folding algorithm and the
//! Turner-2004 energy model are *its* code, not re-implemented here.
//! This module adds the mRNA-design framing:
//!
//! - [`check_structure`] folds a window around the start codon and
//!   reports whether the start codon is *occluded* (base-paired);
//! - [`minimize_start_structure`] scans **synonymous** CDS variants
//!   near the start codon and picks the one whose start region is the
//!   least structured.
//!
//! ## v1 scope
//!
//! The MFE folder's Turner-2004 parameters are a faithful subset (see
//! the `valenx-rnastruct` crate note). The minimisation is a bounded
//! synonymous-codon scan over the first few CDS codons — it is not an
//! exhaustive whole-mRNA sequence optimiser, and it does not co-fold
//! the poly-A tail. The "openness" score is a transparent function of
//! the predicted pairing, not a trained ribosome-loading model.

use crate::error::{GeneditingError, Result};
use crate::sequtil::{is_acgu, reverse_transcribe, transcribe};
use serde::{Deserialize, Serialize};
use valenx_bioseq::ops::translate::GeneticCode;
use valenx_rnastruct::{mfe, RnaSeq};

/// The result of an mRNA start-region structure check.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StructureCheck {
    /// Predicted minimum free energy of the folded window (kcal/mol;
    /// more negative = more stable structure).
    pub mfe: f64,
    /// Dot-bracket structure of the folded window.
    pub dot_bracket: String,
    /// Number of the start codon's three bases that are base-paired
    /// (`0` = fully open, `3` = fully occluded).
    pub start_codon_paired: usize,
    /// `true` when the start codon is occluded enough to impair
    /// scanning (≥ 2 of its 3 bases paired).
    pub start_occluded: bool,
    /// An "openness" score in `[0, 1]` — `1.0` = the start region is
    /// completely unpaired, falling as pairing accumulates.
    pub openness: f64,
}

/// The window, in nucleotides, folded around the start codon: a few
/// bases of 5′UTR plus the first stretch of CDS. A short window keeps
/// folding fast and focuses on the region that gates initiation.
const START_WINDOW: usize = 45;

/// Folds the region around an mRNA start codon and reports occlusion
/// (feature 20).
///
/// `utr5` is the 5′UTR and `cds` the coding sequence (either may be DNA
/// or RNA). A fixed-size window centred on the start codon (the last
/// bases of the 5′UTR + the first CDS bases) is folded with the Zuker
/// MFE algorithm; the three start-codon bases are checked for
/// pairing.
///
/// # Errors
/// - [`GeneditingError::InvalidTarget`] for a non-ACGU 5′UTR / CDS or a
///   CDS shorter than 3 nt.
/// - [`GeneditingError::Invalid`] if the underlying folder rejects the
///   window.
pub fn check_structure(utr5: &[u8], cds: &[u8]) -> Result<StructureCheck> {
    let utr = transcribe(utr5);
    let cds_rna = transcribe(cds);
    if !utr.is_empty() && !is_acgu(&utr) {
        return Err(GeneditingError::invalid_target(
            "region",
            "5'UTR must be A/C/G/U",
        ));
    }
    if cds_rna.len() < 3 || !is_acgu(&cds_rna) {
        return Err(GeneditingError::invalid_target(
            "cds",
            "CDS must be at least one A/C/G/U codon",
        ));
    }
    // Window: up to half the budget of 5'UTR tail + the rest of CDS.
    let utr_take = (START_WINDOW / 3).min(utr.len());
    let cds_take = (START_WINDOW - utr_take).min(cds_rna.len());
    let mut window: Vec<u8> = Vec::with_capacity(utr_take + cds_take);
    window.extend_from_slice(&utr[utr.len() - utr_take..]);
    window.extend_from_slice(&cds_rna[..cds_take]);
    // The start codon's first base sits at index `utr_take` in `window`.
    let start_idx = utr_take;

    let folded = fold_window(&window)?;
    let pairs = paired_in(&folded.1, start_idx, 3.min(cds_take));
    let openness = 1.0 - pairs as f64 / 3.0;

    Ok(StructureCheck {
        mfe: folded.0,
        dot_bracket: folded.1,
        start_codon_paired: pairs,
        start_occluded: pairs >= 2,
        openness: openness.clamp(0.0, 1.0),
    })
}

/// Folds an RNA window, returning `(mfe, dot_bracket)`.
fn fold_window(window: &[u8]) -> Result<(f64, String)> {
    let seq = RnaSeq::parse(window)
        .map_err(|e| GeneditingError::invalid("structure_window", e.to_string()))?;
    let result = mfe(&seq).map_err(|e| GeneditingError::invalid("fold", e.to_string()))?;
    let db = result.structure.to_dot_bracket();
    Ok((result.energy, db))
}

/// Counts the paired positions in a dot-bracket string over the
/// `len` bases starting at `start`.
fn paired_in(dot_bracket: &str, start: usize, len: usize) -> usize {
    dot_bracket
        .as_bytes()
        .iter()
        .skip(start)
        .take(len)
        .filter(|&&c| c == b'(' || c == b')')
        .count()
}

/// The result of a start-region structure minimisation.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StructureMinimization {
    /// The CDS with the least-structured start region (RNA), start and
    /// stop preserved, protein unchanged.
    pub optimized_cds: Vec<u8>,
    /// Start-region openness of the input CDS.
    pub openness_before: f64,
    /// Start-region openness of the optimised CDS.
    pub openness_after: f64,
    /// Number of synonymous CDS variants evaluated.
    pub variants_tried: usize,
}

impl StructureMinimization {
    /// The openness improvement (`after - before`).
    pub fn openness_gain(&self) -> f64 {
        self.openness_after - self.openness_before
    }
}

/// Minimises start-codon-region structure by a synonymous-codon scan
/// (feature 20).
///
/// Holds the 5′UTR and the start codon fixed, then tries each
/// synonymous codon for the first `codons_to_scan` *sense* codons of
/// the CDS, folds the start window for each single-codon swap, and
/// keeps the CDS whose start region is the most open. The protein is
/// never changed.
///
/// # Errors
/// - [`GeneditingError::InvalidTarget`] for a CDS that is empty, not a
///   multiple of three, or non-ACGU.
/// - [`GeneditingError::Invalid`] for `codons_to_scan == 0`.
pub fn minimize_start_structure(
    utr5: &[u8],
    cds: &[u8],
    codons_to_scan: usize,
) -> Result<StructureMinimization> {
    if codons_to_scan == 0 {
        return Err(GeneditingError::invalid(
            "codons_to_scan",
            "must scan at least one codon",
        ));
    }
    let mut cds_rna = transcribe(cds);
    if cds_rna.is_empty() || cds_rna.len() % 3 != 0 || !is_acgu(&cds_rna) {
        return Err(GeneditingError::invalid_target(
            "cds",
            "CDS must be a non-empty A/C/G/U sequence of length divisible by 3",
        ));
    }
    let code = GeneticCode::standard();
    let before = check_structure(utr5, &cds_rna)?.openness;

    let mut best_openness = before;
    let mut tried = 0usize;
    let n_codons = cds_rna.len() / 3;
    // Codon 0 is the start codon — never touch it. Scan codons 1..=k.
    let scan_end = (1 + codons_to_scan).min(n_codons);
    for ci in 1..scan_end {
        let orig: [u8; 3] = [
            cds_rna[ci * 3],
            cds_rna[ci * 3 + 1],
            cds_rna[ci * 3 + 2],
        ];
        let aa = code.translate_codon(&reverse_transcribe(&orig));
        for syn in synonymous_codons(aa, &code) {
            let syn_rna = transcribe(&syn);
            if syn_rna[..] == orig[..] {
                continue;
            }
            // Apply the single-codon swap on a trial copy.
            let mut trial = cds_rna.clone();
            trial[ci * 3] = syn_rna[0];
            trial[ci * 3 + 1] = syn_rna[1];
            trial[ci * 3 + 2] = syn_rna[2];
            tried += 1;
            let open = check_structure(utr5, &trial)?.openness;
            if open > best_openness {
                best_openness = open;
                cds_rna = trial;
            }
        }
    }

    Ok(StructureMinimization {
        optimized_cds: cds_rna,
        openness_before: before,
        openness_after: best_openness,
        variants_tried: tried,
    })
}

/// Every DNA codon synonymous with amino acid `aa` under `code`.
fn synonymous_codons(aa: u8, code: &GeneticCode) -> Vec<[u8; 3]> {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut out = Vec::new();
    for &b0 in &BASES {
        for &b1 in &BASES {
            for &b2 in &BASES {
                let codon = [b0, b1, b2];
                if code.translate_codon(&codon) == aa {
                    out.push(codon);
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn check_structure_reports_openness() {
        // An AU-rich 5'UTR + a simple CDS — little structure expected.
        let c = check_structure(b"AAAAAAAAAAAAAAAAAAAA", b"ATGGCCGCCGCCTAA").unwrap();
        assert!((0.0..=1.0).contains(&c.openness));
        assert!(c.start_codon_paired <= 3);
        assert_eq!(c.dot_bracket.len(), c.dot_bracket.len()); // present
    }

    #[test]
    fn open_start_is_not_occluded() {
        // An unstructured context — the start codon should be open.
        let c = check_structure(b"AAAAAAAAAAAAAAAAAAAA", b"ATGAAAAAAAAATAA").unwrap();
        assert!(!c.start_occluded);
        assert!(c.openness > 0.5);
    }

    #[test]
    fn rejects_non_acgu_inputs() {
        assert!(check_structure(b"NNNN", b"ATGTAA").is_err());
        assert!(check_structure(b"AAAA", b"NN").is_err());
    }

    #[test]
    fn minimization_keeps_protein() {
        let cds = b"ATGCTGCTGCTGCTGTAA";
        let m = minimize_start_structure(b"AAAAAAAAAAAAAAAAAAAA", cds, 3).unwrap();
        let code = GeneticCode::standard();
        let p_before = {
            let dna = reverse_transcribe(&transcribe(cds));
            dna.chunks(3)
                .map(|c| code.translate_codon(c))
                .collect::<Vec<_>>()
        };
        let p_after = {
            let dna = reverse_transcribe(&m.optimized_cds);
            dna.chunks(3)
                .map(|c| code.translate_codon(c))
                .collect::<Vec<_>>()
        };
        assert_eq!(p_before, p_after);
    }

    #[test]
    fn minimization_does_not_reduce_openness() {
        let cds = b"ATGCTGCTGCTGCTGTAA";
        let m = minimize_start_structure(b"GGGCCCGGGCCCGGGCCC", cds, 4).unwrap();
        assert!(m.openness_after >= m.openness_before - 1e-9);
        assert!(m.openness_gain() >= -1e-9);
    }

    #[test]
    fn minimization_rejects_zero_scan() {
        assert!(minimize_start_structure(b"AAAA", b"ATGTAA", 0).is_err());
    }

    #[test]
    fn minimization_rejects_bad_cds() {
        assert!(minimize_start_structure(b"AAAA", b"ATGCT", 2).is_err());
    }

    #[test]
    fn synonymous_codons_translate_identically() {
        let code = GeneticCode::standard();
        // Leucine has 6 codons.
        let leu = synonymous_codons(b'L', &code);
        assert_eq!(leu.len(), 6);
        for c in leu {
            assert_eq!(code.translate_codon(&c), b'L');
        }
    }

    #[test]
    fn paired_in_counts_brackets() {
        assert_eq!(paired_in("..((..))..", 2, 4), 2);
        assert_eq!(paired_in("..........", 0, 5), 0);
    }
}
