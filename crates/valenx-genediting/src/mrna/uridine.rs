//! Feature 21 — uridine-content and modified-nucleoside planning.
//!
//! Therapeutic mRNA is almost always made with a **modified uridine**
//! (N1-methylpseudouridine, `m1Ψ`) in place of every `U` — it slashes
//! innate-immune activation and raises protein output. Two design
//! consequences follow:
//!
//! - the lower the mRNA's **uridine content**, the less modified
//!   nucleoside is consumed and the lower the residual immunogenicity
//!   even before modification;
//! - among synonymous codons, some carry fewer `U`s than others, so a
//!   **uridine-depleting** synonymous-codon pass can lower the `U`
//!   count without changing the protein.
//!
//! This module measures uridine content ([`uridine_content`]),
//! plans a modified-nucleoside substitution ([`plan_modification`]),
//! and runs a uridine-minimising synonymous-codon optimisation of a
//! CDS ([`minimize_uridine`]).
//!
//! ## v1 scope
//!
//! Uridine minimisation is a per-codon greedy synonymous swap — for
//! each residue it picks the synonymous codon with the fewest `U`s
//! (ties broken by host codon frequency via [`crate::mrna::codon`]
//! when a host is given, else lexicographically). It does not jointly
//! optimise structure or CAI; those are separate passes. It changes
//! only the *sequence* — actual `m1Ψ` incorporation is a wet-lab IVT
//! choice the [`ModificationPlan`] documents.

use crate::error::{GeneditingError, Result};
use crate::mrna::codon::ExpressionHost;
use crate::sequtil::{is_acgu, reverse_transcribe, transcribe};
use serde::{Deserialize, Serialize};
use valenx_bioseq::cloning::codon_opt::{codon_usage_table, CodonUsageTable, Host};
use valenx_bioseq::ops::translate::GeneticCode;

/// A modified nucleoside used in place of standard uridine.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ModifiedNucleoside {
    /// Standard, unmodified uridine.
    Uridine,
    /// Pseudouridine (`Ψ`).
    Pseudouridine,
    /// N1-methylpseudouridine (`m1Ψ`) — the modern therapeutic
    /// standard (used in the approved COVID-19 mRNA vaccines).
    N1MethylPseudouridine,
    /// 5-methoxyuridine (`5moU`).
    Methoxyuridine,
}

impl ModifiedNucleoside {
    /// Human-readable name.
    pub fn name(self) -> &'static str {
        match self {
            ModifiedNucleoside::Uridine => "uridine (unmodified)",
            ModifiedNucleoside::Pseudouridine => "pseudouridine (Ψ)",
            ModifiedNucleoside::N1MethylPseudouridine => "N1-methylpseudouridine (m1Ψ)",
            ModifiedNucleoside::Methoxyuridine => "5-methoxyuridine (5moU)",
        }
    }

    /// `true` when the modification meaningfully reduces innate-immune
    /// activation (everything except unmodified uridine).
    pub fn reduces_immunogenicity(self) -> bool {
        !matches!(self, ModifiedNucleoside::Uridine)
    }
}

/// Uridine fraction of a sequence in `[0, 1]` — the fraction of bases
/// that are `U` (or `T`, treated as `U`).
///
/// An empty sequence returns `0`. This is the metric a modified-
/// nucleoside plan minimises.
pub fn uridine_content(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let u = seq
        .iter()
        .filter(|&&b| matches!(b.to_ascii_uppercase(), b'U' | b'T'))
        .count();
    u as f64 / seq.len() as f64
}

/// A modified-nucleoside substitution plan for an mRNA construct.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ModificationPlan {
    /// The modified nucleoside that will replace every `U`.
    pub nucleoside: ModifiedNucleoside,
    /// Total uridine positions in the transcript that get modified.
    pub uridine_positions: usize,
    /// Uridine fraction of the transcript.
    pub uridine_fraction: f64,
    /// `true` when the chosen modification reduces immunogenicity.
    pub reduces_immunogenicity: bool,
    /// A one-line rationale.
    pub rationale: String,
}

/// Plans a modified-nucleoside substitution for a transcript
/// (feature 21).
///
/// Counts the uridine positions a full `U → modified` substitution
/// would touch and reports the plan. The transcript may be DNA or RNA.
///
/// # Errors
/// [`GeneditingError::InvalidTarget`] for a non-ACGU transcript.
pub fn plan_modification(
    transcript: &[u8],
    nucleoside: ModifiedNucleoside,
) -> Result<ModificationPlan> {
    let rna = transcribe(transcript);
    if rna.is_empty() || !is_acgu(&rna) {
        return Err(GeneditingError::invalid_target(
            "region",
            "transcript must be a non-empty A/C/G/U sequence",
        ));
    }
    let positions = rna.iter().filter(|&&b| b == b'U').count();
    let fraction = uridine_content(&rna);
    let rationale = if nucleoside.reduces_immunogenicity() {
        format!(
            "Replace all {positions} uridines with {} — reduces RIG-I / TLR \
             innate sensing and typically raises protein output.",
            nucleoside.name()
        )
    } else {
        format!(
            "{positions} uridines, unmodified — expect stronger innate-immune \
             activation; consider m1Ψ for a therapeutic construct."
        )
    };
    Ok(ModificationPlan {
        nucleoside,
        uridine_positions: positions,
        uridine_fraction: fraction,
        reduces_immunogenicity: nucleoside.reduces_immunogenicity(),
        rationale,
    })
}

/// The result of a uridine-minimising CDS optimisation.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct UridineMinimization {
    /// The uridine-depleted CDS (RNA), protein and start / stop
    /// preserved.
    pub optimized_cds: Vec<u8>,
    /// Uridine fraction of the input CDS.
    pub uridine_before: f64,
    /// Uridine fraction of the optimised CDS.
    pub uridine_after: f64,
}

impl UridineMinimization {
    /// The uridine reduction (`before - after`); `>= 0`.
    pub fn uridine_reduction(&self) -> f64 {
        self.uridine_before - self.uridine_after
    }
}

/// Minimises the uridine content of a CDS by synonymous-codon choice
/// (feature 21).
///
/// For every sense codon, replaces it with the synonymous codon
/// carrying the fewest `U`s. When several synonymous codons tie on
/// `U`-count, the host's most-frequent codon is chosen if `host` is
/// `Some` (so uridine depletion does not wreck CAI), otherwise the
/// lexicographically smallest. The protein and the start / stop codons
/// are preserved.
///
/// # Errors
/// [`GeneditingError::InvalidTarget`] for a CDS that is empty, not a
/// multiple of three, or non-ACGU.
pub fn minimize_uridine(cds: &[u8], host: Option<ExpressionHost>) -> Result<UridineMinimization> {
    let cds_rna = transcribe(cds);
    if cds_rna.is_empty() || cds_rna.len() % 3 != 0 || !is_acgu(&cds_rna) {
        return Err(GeneditingError::invalid_target(
            "cds",
            "CDS must be a non-empty A/C/G/U sequence of length divisible by 3",
        ));
    }
    let code = GeneticCode::standard();
    let table: Option<CodonUsageTable> = host.map(|h| {
        codon_usage_table(match h {
            ExpressionHost::Human => Host::Human,
            ExpressionHost::EColi => Host::EColi,
        })
    });

    let before = uridine_content(&cds_rna);
    let n_codons = cds_rna.len() / 3;
    let mut out = cds_rna.clone();

    for ci in 0..n_codons {
        let orig: [u8; 3] = [
            cds_rna[ci * 3],
            cds_rna[ci * 3 + 1],
            cds_rna[ci * 3 + 2],
        ];
        let orig_dna = reverse_transcribe(&orig);
        let aa = code.translate_codon(&orig_dna);
        // Keep the start codon and stop codon fixed.
        if ci == 0 || aa == b'*' {
            continue;
        }
        let best = pick_low_u_codon(aa, &code, table.as_ref());
        if let Some(syn) = best {
            let syn_rna = transcribe(&syn);
            out[ci * 3] = syn_rna[0];
            out[ci * 3 + 1] = syn_rna[1];
            out[ci * 3 + 2] = syn_rna[2];
        }
    }
    let after = uridine_content(&out);
    Ok(UridineMinimization {
        optimized_cds: out,
        uridine_before: before,
        uridine_after: after,
    })
}

/// Picks the synonymous DNA codon for `aa` with the fewest `T`/`U`s,
/// breaking ties by host frequency (if a table is given) then
/// lexicographically.
fn pick_low_u_codon(
    aa: u8,
    code: &GeneticCode,
    table: Option<&CodonUsageTable>,
) -> Option<[u8; 3]> {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut best: Option<([u8; 3], usize, f64)> = None; // (codon, u_count, freq)
    for &b0 in &BASES {
        for &b1 in &BASES {
            for &b2 in &BASES {
                let codon = [b0, b1, b2];
                if code.translate_codon(&codon) != aa {
                    continue;
                }
                let u = codon.iter().filter(|&&b| b == b'T').count();
                let freq = table.map(|t| t.frequency(&codon)).unwrap_or(0.0);
                let better = match best {
                    None => true,
                    Some((bc, bu, bf)) => {
                        u < bu
                            || (u == bu && freq > bf)
                            || (u == bu && (freq - bf).abs() < 1e-12 && codon < bc)
                    }
                };
                if better {
                    best = Some((codon, u, freq));
                }
            }
        }
    }
    best.map(|(c, _, _)| c)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uridine_content_counts_u_and_t() {
        assert!((uridine_content(b"UUUU") - 1.0).abs() < 1e-9);
        assert!((uridine_content(b"TTTT") - 1.0).abs() < 1e-9);
        assert!((uridine_content(b"ACGU") - 0.25).abs() < 1e-9);
        assert_eq!(uridine_content(b""), 0.0);
    }

    #[test]
    fn plan_counts_uridines() {
        let p = plan_modification(b"AUGUUUTAA", ModifiedNucleoside::N1MethylPseudouridine)
            .unwrap();
        // AUGUUUUAA after transcription: count the U's.
        assert!(p.uridine_positions > 0);
        assert!(p.reduces_immunogenicity);
        assert!(p.rationale.contains("m1"));
    }

    #[test]
    fn unmodified_plan_flags_immunogenicity() {
        let p = plan_modification(b"AUGUUUUAA", ModifiedNucleoside::Uridine).unwrap();
        assert!(!p.reduces_immunogenicity);
    }

    #[test]
    fn rejects_non_acgu_transcript() {
        assert!(plan_modification(b"NNNN", ModifiedNucleoside::Pseudouridine).is_err());
    }

    #[test]
    fn modified_nucleoside_immunogenicity() {
        assert!(ModifiedNucleoside::N1MethylPseudouridine.reduces_immunogenicity());
        assert!(ModifiedNucleoside::Pseudouridine.reduces_immunogenicity());
        assert!(!ModifiedNucleoside::Uridine.reduces_immunogenicity());
    }

    #[test]
    fn minimize_uridine_does_not_increase_u() {
        // A CDS using U-heavy synonymous codons: Phe UUU, Leu UUA.
        let cds = b"ATGTTTTTATTTTAA";
        let m = minimize_uridine(cds, Some(ExpressionHost::Human)).unwrap();
        assert!(m.uridine_after <= m.uridine_before + 1e-9);
        assert!(m.uridine_reduction() >= -1e-9);
    }

    #[test]
    fn minimize_uridine_keeps_protein() {
        let cds = b"ATGTTTTTATTTTAA";
        let m = minimize_uridine(cds, None).unwrap();
        let code = GeneticCode::standard();
        let p_before: Vec<u8> = reverse_transcribe(&transcribe(cds))
            .chunks(3)
            .map(|c| code.translate_codon(c))
            .collect();
        let p_after: Vec<u8> = reverse_transcribe(&m.optimized_cds)
            .chunks(3)
            .map(|c| code.translate_codon(c))
            .collect();
        assert_eq!(p_before, p_after);
    }

    #[test]
    fn minimize_uridine_rejects_bad_cds() {
        assert!(minimize_uridine(b"ATGCT", None).is_err());
        assert!(minimize_uridine(b"", None).is_err());
    }

    #[test]
    fn pick_low_u_codon_minimises_t() {
        let code = GeneticCode::standard();
        // Glycine: GGN — none contain T, so any is fine; check it
        // returns a Gly codon with zero T.
        let g = pick_low_u_codon(b'G', &code, None).unwrap();
        assert_eq!(code.translate_codon(&g), b'G');
        assert_eq!(g.iter().filter(|&&b| b == b'T').count(), 0);
    }

    #[test]
    fn optimized_cds_is_rna() {
        let m = minimize_uridine(b"ATGTTTTTATTTTAA", None).unwrap();
        assert!(!m.optimized_cds.contains(&b'T'));
    }
}
