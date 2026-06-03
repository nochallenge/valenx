//! Feature 17 — CDS codon optimisation for expression.
//!
//! Synonymous codons are translated identically but expressed very
//! differently — a coding sequence whose codons match the host's
//! preferred set is translated faster and to higher protein levels.
//! The standard metric is the **codon adaptation index (CAI)**: the
//! geometric mean of each codon's relative adaptiveness, `(0, 1]`.
//!
//! This module reuses [`valenx_bioseq`]'s codon machinery
//! ([`valenx_bioseq::cloning::codon_opt`]) — the codon-usage tables,
//! the protein→DNA optimiser and the CAI calculator are *its* code,
//! not re-implemented here. This module adds the mRNA-design framing:
//! it works from a *CDS* (not a bare protein), preserves the start /
//! stop codons, and reports a before / after CAI.
//!
//! ## v1 scope
//!
//! Optimisation here is the simple "pick the host-optimal synonymous
//! codon for every residue" strategy `valenx-bioseq` implements — it
//! maximises CAI. It does not also balance GC content, avoid cryptic
//! splice sites or smooth ramp / rare-codon usage; those are layered
//! on by [`crate::mrna::structure`] (structure) and
//! [`crate::mrna::uridine`] (uridine) as separate passes.

use crate::error::{GeneditingError, Result};
use crate::sequtil::{is_acgu, reverse_transcribe, transcribe};
use serde::{Deserialize, Serialize};
use valenx_bioseq::cloning::codon_opt::{codon_adaptation_index, codon_optimize, Host};
use valenx_bioseq::ops::translate::{translate_default, GeneticCode};
use valenx_bioseq::{Seq, SeqKind};

/// The expression host a CDS is optimised for.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ExpressionHost {
    /// *Homo sapiens* — the host for human therapeutic mRNA.
    Human,
    /// *Escherichia coli* — a common recombinant-protein host.
    EColi,
}

impl ExpressionHost {
    /// Maps to the [`valenx_bioseq`] `Host` enum.
    fn to_bioseq(self) -> Host {
        match self {
            ExpressionHost::Human => Host::Human,
            ExpressionHost::EColi => Host::EColi,
        }
    }

    /// Human-readable name.
    pub fn name(self) -> &'static str {
        match self {
            ExpressionHost::Human => "human",
            ExpressionHost::EColi => "E. coli",
        }
    }
}

/// The result of a CDS codon-optimisation pass.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CodonOptimization {
    /// The optimised CDS, RNA (`A C G U`), start and stop preserved.
    pub optimized_cds: Vec<u8>,
    /// The CAI of the *input* CDS for the host.
    pub cai_before: f64,
    /// The CAI of the optimised CDS for the host.
    pub cai_after: f64,
    /// The host the CDS was optimised for.
    pub host: ExpressionHost,
}

impl CodonOptimization {
    /// The CAI improvement (`cai_after - cai_before`); may be `0` (or
    /// slightly negative due to the geometric mean) if the input was
    /// already optimal.
    pub fn cai_gain(&self) -> f64 {
        self.cai_after - self.cai_before
    }
}

/// Optimises a CDS for expression in a host (feature 17).
///
/// Translates the CDS to protein, calls [`valenx_bioseq`]'s
/// `codon_optimize` to pick host-optimal codons, transcribes the
/// result to RNA and reports the before / after CAI. The original stop
/// codon's amino acid (`*`) is carried through the protein so the
/// optimiser emits a host-preferred stop.
///
/// # Errors
/// - [`GeneditingError::InvalidTarget`] for a CDS that is empty, not a
///   multiple of three, or contains a non-ACGU base.
/// - [`GeneditingError::Invalid`] if the underlying `valenx-bioseq`
///   optimiser rejects the translated protein.
pub fn optimize_cds(cds: &[u8], host: ExpressionHost) -> Result<CodonOptimization> {
    let rna = transcribe(cds);
    if rna.is_empty() || !is_acgu(&rna) {
        return Err(GeneditingError::invalid_target(
            "cds",
            "CDS must be a non-empty A/C/G/U sequence",
        ));
    }
    if rna.len() % 3 != 0 {
        return Err(GeneditingError::invalid_target(
            "cds",
            "CDS length must be a multiple of 3",
        ));
    }
    // valenx-bioseq works on DNA `Seq`s — reverse-transcribe.
    let dna_bytes = reverse_transcribe(&rna);
    let dna = Seq::new(SeqKind::Dna, &dna_bytes)
        .map_err(|e| GeneditingError::invalid_target("cds", e.to_string()))?;

    let code = GeneticCode::standard();
    // Translate to protein, *keeping* the terminal stop as `*` so the
    // optimiser re-emits a stop codon.
    let protein = translate_default(&dna, &code)
        .map_err(|e| GeneditingError::invalid_target("cds", e.to_string()))?;

    let cai_before = codon_adaptation_index(&dna, host.to_bioseq())
        .map_err(|e| GeneditingError::invalid("cds", e.to_string()))?;

    let optimized_dna = codon_optimize(&protein, host.to_bioseq())
        .map_err(|e| GeneditingError::invalid("protein", e.to_string()))?;
    let cai_after = codon_adaptation_index(&optimized_dna, host.to_bioseq())
        .map_err(|e| GeneditingError::invalid("cds", e.to_string()))?;

    let optimized_cds = transcribe(optimized_dna.as_bytes());
    Ok(CodonOptimization {
        optimized_cds,
        cai_before,
        cai_after,
        host,
    })
}

/// The codon adaptation index of a CDS for a host (a thin,
/// mRNA-framed wrapper over [`valenx_bioseq`]'s calculator).
///
/// # Errors
/// [`GeneditingError::InvalidTarget`] for an invalid CDS.
pub fn cds_cai(cds: &[u8], host: ExpressionHost) -> Result<f64> {
    let rna = transcribe(cds);
    if rna.is_empty() || !is_acgu(&rna) || rna.len() % 3 != 0 {
        return Err(GeneditingError::invalid_target(
            "cds",
            "CDS must be a non-empty A/C/G/U sequence of length divisible by 3",
        ));
    }
    let dna_bytes = reverse_transcribe(&rna);
    let dna = Seq::new(SeqKind::Dna, &dna_bytes)
        .map_err(|e| GeneditingError::invalid_target("cds", e.to_string()))?;
    codon_adaptation_index(&dna, host.to_bioseq())
        .map_err(|e| GeneditingError::invalid("cds", e.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;

    // A CDS using deliberately rare human codons so optimisation has
    // room to improve: ATG (Met) + rare Leu codons + stop.
    fn rare_cds() -> Vec<u8> {
        b"ATGCTACTACTACTATAA".to_vec()
    }

    #[test]
    fn optimization_does_not_change_protein() {
        let opt = optimize_cds(&rare_cds(), ExpressionHost::Human).unwrap();
        let code = GeneticCode::standard();
        let before_dna =
            Seq::new(SeqKind::Dna, reverse_transcribe(&transcribe(&rare_cds()))).unwrap();
        let after_dna =
            Seq::new(SeqKind::Dna, reverse_transcribe(&opt.optimized_cds)).unwrap();
        let p1 = translate_default(&before_dna, &code).unwrap();
        let p2 = translate_default(&after_dna, &code).unwrap();
        assert_eq!(p1.as_bytes(), p2.as_bytes());
    }

    #[test]
    fn optimization_does_not_reduce_cai() {
        let opt = optimize_cds(&rare_cds(), ExpressionHost::Human).unwrap();
        // Optimising to host-optimal codons cannot lower CAI.
        assert!(opt.cai_after >= opt.cai_before - 1e-9);
        assert!(opt.cai_gain() >= -1e-9);
    }

    #[test]
    fn optimized_cds_is_rna() {
        let opt = optimize_cds(&rare_cds(), ExpressionHost::Human).unwrap();
        assert!(!opt.optimized_cds.contains(&b'T'));
        assert!(is_acgu(&opt.optimized_cds));
    }

    #[test]
    fn optimized_cds_keeps_length() {
        let opt = optimize_cds(&rare_cds(), ExpressionHost::Human).unwrap();
        assert_eq!(opt.optimized_cds.len(), rare_cds().len());
    }

    #[test]
    fn rejects_non_multiple_of_three() {
        let err = optimize_cds(b"ATGCT", ExpressionHost::Human).unwrap_err();
        assert_eq!(err.code(), "genediting.invalid_target");
    }

    #[test]
    fn rejects_empty_cds() {
        assert!(optimize_cds(b"", ExpressionHost::Human).is_err());
    }

    #[test]
    fn cds_cai_in_unit_interval() {
        let cai = cds_cai(&rare_cds(), ExpressionHost::Human).unwrap();
        assert!(cai > 0.0 && cai <= 1.0);
    }

    #[test]
    fn ecoli_and_human_can_differ() {
        let cds = rare_cds();
        let human = cds_cai(&cds, ExpressionHost::Human).unwrap();
        let ecoli = cds_cai(&cds, ExpressionHost::EColi).unwrap();
        // Both valid CAIs; the host choice is honoured (values may or
        // may not differ for this short CDS, but both are in range).
        assert!(human > 0.0 && ecoli > 0.0);
    }

    #[test]
    fn accepts_rna_cds_input() {
        let opt = optimize_cds(b"AUGCUACUACUACUAUAA", ExpressionHost::Human).unwrap();
        assert!(opt.cai_after > 0.0);
    }
}
