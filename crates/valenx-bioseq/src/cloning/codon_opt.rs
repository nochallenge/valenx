//! Codon-usage tables and codon optimization.
//!
//! Provides relative-synonymous-codon-usage tables for *E. coli* and
//! *Homo sapiens*, a protein→DNA codon optimizer that picks the most
//! frequent codon for each amino acid in the chosen host, and the
//! codon adaptation index (CAI) — a measure of how well a coding
//! sequence matches a host's codon preferences.

use crate::error::{BioseqError, Result};
use crate::ops::translate::GeneticCode;
use crate::seq::{Seq, SeqKind, Topology};
use std::collections::HashMap;

/// A target host for codon optimization.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum Host {
    /// *Escherichia coli* K-12.
    EColi,
    /// *Homo sapiens*.
    Human,
}

impl Host {
    /// Human-readable name.
    pub fn name(self) -> &'static str {
        match self {
            Host::EColi => "E. coli",
            Host::Human => "human",
        }
    }
}

/// A codon-usage table: per-codon relative frequency (fraction of that
/// amino acid's total usage, so the synonymous codons for one amino
/// acid sum to 1.0).
#[derive(Clone, Debug)]
pub struct CodonUsageTable {
    /// The host this table describes.
    pub host: Host,
    /// `codon → relative frequency within its synonymous family`.
    freq: HashMap<[u8; 3], f64>,
}

impl CodonUsageTable {
    /// Relative frequency of a codon within its amino-acid family
    /// (`0.0` for an unknown codon).
    pub fn frequency(&self, codon: &[u8]) -> f64 {
        if codon.len() != 3 {
            return 0.0;
        }
        let key = [
            codon[0].to_ascii_uppercase(),
            codon[1].to_ascii_uppercase(),
            codon[2].to_ascii_uppercase(),
        ];
        self.freq.get(&key).copied().unwrap_or(0.0)
    }

    /// The most-frequent (optimal) codon for an amino acid in this
    /// host, or `None` if the amino acid is unknown.
    pub fn optimal_codon(&self, amino_acid: u8, code: &GeneticCode) -> Option<[u8; 3]> {
        let aa = amino_acid.to_ascii_uppercase();
        let mut best: Option<([u8; 3], f64)> = None;
        for (&codon, &f) in &self.freq {
            if code.translate_codon(&codon) == aa {
                match best {
                    None => best = Some((codon, f)),
                    Some((_, bf)) if f > bf => best = Some((codon, f)),
                    // Deterministic tie-break: lexicographically
                    // smallest codon.
                    Some((bc, bf)) if (f - bf).abs() < 1e-12 && codon < bc => {
                        best = Some((codon, f))
                    }
                    Some(_) => {}
                }
            }
        }
        best.map(|(c, _)| c)
    }
}

/// Builds the codon-usage table for a host.
pub fn codon_usage_table(host: Host) -> CodonUsageTable {
    // Relative synonymous codon usage (within-family fractions).
    // Source: Kazusa codon-usage database, rounded; values per family
    // sum to ~1.0. These are representative v1 figures.
    let raw: &[(&str, f64)] = match host {
        Host::EColi => &[
            // Phe
            ("TTT", 0.58), ("TTC", 0.42),
            // Leu
            ("TTA", 0.14), ("TTG", 0.13), ("CTT", 0.12),
            ("CTC", 0.10), ("CTA", 0.04), ("CTG", 0.47),
            // Ile
            ("ATT", 0.49), ("ATC", 0.39), ("ATA", 0.11),
            // Met / Trp
            ("ATG", 1.00), ("TGG", 1.00),
            // Val
            ("GTT", 0.28), ("GTC", 0.20), ("GTA", 0.17), ("GTG", 0.35),
            // Ser
            ("TCT", 0.17), ("TCC", 0.15), ("TCA", 0.14),
            ("TCG", 0.14), ("AGT", 0.16), ("AGC", 0.25),
            // Pro
            ("CCT", 0.18), ("CCC", 0.13), ("CCA", 0.20), ("CCG", 0.49),
            // Thr
            ("ACT", 0.19), ("ACC", 0.40), ("ACA", 0.17), ("ACG", 0.25),
            // Ala
            ("GCT", 0.18), ("GCC", 0.26), ("GCA", 0.23), ("GCG", 0.33),
            // Tyr
            ("TAT", 0.59), ("TAC", 0.41),
            // His
            ("CAT", 0.57), ("CAC", 0.43),
            // Gln
            ("CAA", 0.34), ("CAG", 0.66),
            // Asn
            ("AAT", 0.49), ("AAC", 0.51),
            // Lys
            ("AAA", 0.74), ("AAG", 0.26),
            // Asp
            ("GAT", 0.63), ("GAC", 0.37),
            // Glu
            ("GAA", 0.68), ("GAG", 0.32),
            // Cys
            ("TGT", 0.46), ("TGC", 0.54),
            // Arg
            ("CGT", 0.36), ("CGC", 0.36), ("CGA", 0.07),
            ("CGG", 0.11), ("AGA", 0.07), ("AGG", 0.04),
            // Gly
            ("GGT", 0.35), ("GGC", 0.37), ("GGA", 0.13), ("GGG", 0.15),
            // Stops
            ("TAA", 0.61), ("TAG", 0.09), ("TGA", 0.30),
        ],
        Host::Human => &[
            // Phe
            ("TTT", 0.45), ("TTC", 0.55),
            // Leu
            ("TTA", 0.07), ("TTG", 0.13), ("CTT", 0.13),
            ("CTC", 0.20), ("CTA", 0.07), ("CTG", 0.40),
            // Ile
            ("ATT", 0.36), ("ATC", 0.48), ("ATA", 0.16),
            // Met / Trp
            ("ATG", 1.00), ("TGG", 1.00),
            // Val
            ("GTT", 0.18), ("GTC", 0.24), ("GTA", 0.11), ("GTG", 0.47),
            // Ser
            ("TCT", 0.18), ("TCC", 0.22), ("TCA", 0.15),
            ("TCG", 0.06), ("AGT", 0.15), ("AGC", 0.24),
            // Pro
            ("CCT", 0.28), ("CCC", 0.33), ("CCA", 0.27), ("CCG", 0.11),
            // Thr
            ("ACT", 0.24), ("ACC", 0.36), ("ACA", 0.28), ("ACG", 0.12),
            // Ala
            ("GCT", 0.26), ("GCC", 0.40), ("GCA", 0.23), ("GCG", 0.11),
            // Tyr
            ("TAT", 0.43), ("TAC", 0.57),
            // His
            ("CAT", 0.41), ("CAC", 0.59),
            // Gln
            ("CAA", 0.25), ("CAG", 0.75),
            // Asn
            ("AAT", 0.46), ("AAC", 0.54),
            // Lys
            ("AAA", 0.42), ("AAG", 0.58),
            // Asp
            ("GAT", 0.46), ("GAC", 0.54),
            // Glu
            ("GAA", 0.42), ("GAG", 0.58),
            // Cys
            ("TGT", 0.45), ("TGC", 0.55),
            // Arg
            ("CGT", 0.08), ("CGC", 0.19), ("CGA", 0.11),
            ("CGG", 0.21), ("AGA", 0.20), ("AGG", 0.21),
            // Gly
            ("GGT", 0.16), ("GGC", 0.34), ("GGA", 0.25), ("GGG", 0.25),
            // Stops
            ("TAA", 0.28), ("TAG", 0.20), ("TGA", 0.52),
        ],
    };
    let mut freq = HashMap::new();
    for &(codon, f) in raw {
        let b = codon.as_bytes();
        freq.insert([b[0], b[1], b[2]], f);
    }
    CodonUsageTable { host, freq }
}

/// Codon-optimizes a protein sequence for a host.
///
/// Replaces each amino acid with the host's most-frequent synonymous
/// codon. The standard genetic code (NCBI table 1) is used. A trailing
/// stop is emitted as the host's most-frequent stop codon if the
/// protein ends in `*`. Returns [`BioseqError::Invalid`] for a
/// non-protein input or a residue with no codon.
pub fn codon_optimize(protein: &Seq, host: Host) -> Result<Seq> {
    if protein.kind() != SeqKind::Protein {
        return Err(BioseqError::invalid(
            "kind",
            "codon optimization needs a protein sequence",
        ));
    }
    let table = codon_usage_table(host);
    let code = GeneticCode::standard();
    let mut dna: Vec<u8> = Vec::with_capacity(protein.len() * 3);
    for aa in protein.iter() {
        let codon = table.optimal_codon(aa, &code).ok_or_else(|| {
            BioseqError::invalid(
                "sequence",
                format!("amino acid `{}` has no codon in the {} table", aa as char, host.name()),
            )
        })?;
        dna.extend_from_slice(&codon);
    }
    Ok(Seq::new_unchecked(SeqKind::Dna, dna, Topology::Linear))
}

/// Codon adaptation index (CAI) of a coding sequence for a host.
///
/// CAI is the geometric mean of the *relative adaptiveness* `wᵢ` of
/// each codon, where `wᵢ = freq(codon) / freq(optimal synonym)`. CAI
/// ranges `(0, 1]`; `1.0` means every codon is the host-optimal one.
/// Codons translating to `Met` / `Trp` (single-codon families) are
/// skipped, as is any trailing stop. Returns
/// [`BioseqError::Invalid`] for a non-DNA input or a length not
/// divisible by 3.
pub fn codon_adaptation_index(cds: &Seq, host: Host) -> Result<f64> {
    if cds.kind() != SeqKind::Dna {
        return Err(BioseqError::invalid("kind", "CAI needs a DNA coding sequence"));
    }
    let bytes = cds.as_bytes();
    if bytes.is_empty() || bytes.len() % 3 != 0 {
        return Err(BioseqError::invalid(
            "sequence",
            "CDS length must be a positive multiple of 3",
        ));
    }
    let table = codon_usage_table(host);
    let code = GeneticCode::standard();

    let mut log_sum = 0.0f64;
    let mut counted = 0usize;
    for c in 0..bytes.len() / 3 {
        let codon = &bytes[c * 3..c * 3 + 3];
        let aa = code.translate_codon(codon);
        // Skip single-codon families and stop codons — they carry no
        // adaptive information.
        if matches!(aa, b'M' | b'W' | b'*' | b'X') {
            continue;
        }
        let f = table.frequency(codon);
        let optimal = match table.optimal_codon(aa, &code) {
            Some(o) => table.frequency(&o),
            None => continue,
        };
        if f <= 0.0 || optimal <= 0.0 {
            continue;
        }
        let w = f / optimal;
        log_sum += w.ln();
        counted += 1;
    }
    if counted == 0 {
        return Ok(1.0);
    }
    Ok((log_sum / counted as f64).exp())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tables_have_all_amino_acids() {
        let code = GeneticCode::standard();
        for host in [Host::EColi, Host::Human] {
            let t = codon_usage_table(host);
            for &aa in b"ACDEFGHIKLMNPQRSTVWY" {
                assert!(
                    t.optimal_codon(aa, &code).is_some(),
                    "{} missing codon for {}",
                    host.name(),
                    aa as char
                );
            }
        }
    }

    #[test]
    fn synonymous_families_sum_near_one() {
        let code = GeneticCode::standard();
        let t = codon_usage_table(Host::EColi);
        // Leucine has 6 codons; their frequencies should sum to ~1.
        let mut total = 0.0;
        for c0 in [b'T', b'C', b'A', b'G'] {
            for c1 in [b'T', b'C', b'A', b'G'] {
                for c2 in [b'T', b'C', b'A', b'G'] {
                    let codon = [c0, c1, c2];
                    if code.translate_codon(&codon) == b'L' {
                        total += t.frequency(&codon);
                    }
                }
            }
        }
        assert!((total - 1.0).abs() < 0.05, "Leu family sums to {total}");
    }

    #[test]
    fn ecoli_optimal_leucine_is_ctg() {
        let code = GeneticCode::standard();
        let t = codon_usage_table(Host::EColi);
        assert_eq!(t.optimal_codon(b'L', &code), Some(*b"CTG"));
    }

    #[test]
    fn optimize_produces_translatable_dna() {
        let protein = Seq::new(SeqKind::Protein, "MKVLAAG").unwrap();
        let dna = codon_optimize(&protein, Host::EColi).unwrap();
        assert_eq!(dna.kind(), SeqKind::Dna);
        assert_eq!(dna.len(), protein.len() * 3);
        // Translating it back recovers the protein.
        let code = GeneticCode::standard();
        let back = crate::ops::translate::translate_default(&dna, &code).unwrap();
        assert_eq!(back.as_str(), protein.as_str());
    }

    #[test]
    fn optimize_uses_host_preferences() {
        // For E. coli the optimal Leu codon is CTG.
        let protein = Seq::new(SeqKind::Protein, "L").unwrap();
        let dna = codon_optimize(&protein, Host::EColi).unwrap();
        assert_eq!(dna.as_str(), "CTG");
    }

    #[test]
    fn cai_of_optimized_sequence_is_high() {
        let protein = Seq::new(SeqKind::Protein, "MKVLAAGGRRSSLL").unwrap();
        let optimized = codon_optimize(&protein, Host::EColi).unwrap();
        let cai = codon_adaptation_index(&optimized, Host::EColi).unwrap();
        // An all-optimal sequence should score CAI ~ 1.0.
        assert!(cai > 0.99, "optimized CAI should be ~1, got {cai}");
    }

    #[test]
    fn cai_of_suboptimal_sequence_is_lower() {
        // Build a CDS using deliberately rare E. coli codons.
        // Rare Leu codon CTA, rare Arg codon AGG.
        let cds = Seq::new(SeqKind::Dna, "CTACTAAGGAGG").unwrap();
        let cai = codon_adaptation_index(&cds, Host::EColi).unwrap();
        assert!(cai < 0.5, "rare-codon CAI should be low, got {cai}");
    }

    #[test]
    fn cai_in_unit_interval() {
        let cds = Seq::new(SeqKind::Dna, "ATGAAAGTGCTGGCAGCAGGT").unwrap();
        let cai = codon_adaptation_index(&cds, Host::Human).unwrap();
        assert!(cai > 0.0 && cai <= 1.0, "CAI out of range: {cai}");
    }

    #[test]
    fn wrong_kinds_rejected() {
        let dna = Seq::new(SeqKind::Dna, "ATGAAA").unwrap();
        assert!(codon_optimize(&dna, Host::EColi).is_err());
        let protein = Seq::new(SeqKind::Protein, "MK").unwrap();
        assert!(codon_adaptation_index(&protein, Host::EColi).is_err());
        // Length not a multiple of 3.
        let bad = Seq::new(SeqKind::Dna, "ATGAA").unwrap();
        assert!(codon_adaptation_index(&bad, Host::EColi).is_err());
    }
}
