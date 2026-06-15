//! Reverse translation, Codon Adaptation Index, and GC content.

use crate::code::{is_single_codon, synonymous_codons, translate_codon};
use crate::error::CodonError;
use crate::usage::CodonUsage;

/// Reverse-translate a protein to a DNA coding sequence, choosing the
/// highest-adaptiveness synonymous codon per residue from `usage`.
///
/// Rejects an empty protein and any residue that is not a sense amino acid or
/// has no weighted codon in `usage`.
pub fn reverse_translate_optimal(protein: &str, usage: &CodonUsage) -> Result<String, CodonError> {
    if protein.is_empty() {
        return Err(CodonError::Empty { what: "protein" });
    }
    let mut dna = String::with_capacity(protein.len() * 3);
    for (pos, aa) in protein.chars().enumerate() {
        if aa == '*' || synonymous_codons(aa).is_empty() {
            return Err(CodonError::InvalidResidue { residue: aa, pos });
        }
        let codon = usage
            .optimal_codon(aa)
            .ok_or(CodonError::InvalidResidue { residue: aa, pos })?;
        dna.push_str(&codon);
    }
    Ok(dna)
}

/// Codon Adaptation Index (Sharp & Li, 1987): the geometric mean of relative
/// adaptiveness over the coding sequence's codons, excluding Met, Trp, and stop
/// codons (the standard convention). Returns `1.0` if no codon is counted.
///
/// `dna` length must be a multiple of three and every codon must be valid with
/// a weight present in `usage`.
pub fn cai(dna: &str, usage: &CodonUsage) -> Result<f64, CodonError> {
    if dna.is_empty() {
        return Err(CodonError::Empty { what: "dna" });
    }
    if dna.len() % 3 != 0 {
        return Err(CodonError::NotMultipleOfThree { len: dna.len() });
    }
    let mut sum_ln = 0.0;
    let mut counted = 0usize;
    for (index, chunk) in dna.as_bytes().chunks(3).enumerate() {
        let codon = std::str::from_utf8(chunk).unwrap_or("");
        let aa = translate_codon(codon).ok_or_else(|| CodonError::InvalidCodon {
            codon: codon.to_string(),
            index,
        })?;
        if aa == '*' || is_single_codon(aa) {
            continue;
        }
        let w = usage
            .weight(codon)
            .ok_or_else(|| CodonError::MissingWeight {
                codon: codon.to_string(),
            })?;
        sum_ln += w.ln();
        counted += 1;
    }
    if counted == 0 {
        return Ok(1.0);
    }
    Ok((sum_ln / counted as f64).exp())
}

/// The GC fraction of a DNA string, in `[0, 1]`. Rejects an empty string or any
/// non-`ACGT(U)` character.
pub fn gc_content(dna: &str) -> Result<f64, CodonError> {
    if dna.is_empty() {
        return Err(CodonError::Empty { what: "dna" });
    }
    let mut gc = 0usize;
    for (pos, b) in dna.bytes().enumerate() {
        match b.to_ascii_uppercase() {
            b'G' | b'C' => gc += 1,
            b'A' | b'T' | b'U' => {}
            other => {
                return Err(CodonError::InvalidCodon {
                    codon: (other as char).to_string(),
                    index: pos,
                })
            }
        }
    }
    Ok(gc as f64 / dna.len() as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::usage::illustrative_weights;

    #[test]
    fn reverse_translate_round_trips_to_protein() {
        let u = illustrative_weights();
        let protein = "MFLGKW";
        let dna = reverse_translate_optimal(protein, &u).unwrap();
        assert_eq!(dna.len(), protein.len() * 3);
        let back: String = dna
            .as_bytes()
            .chunks(3)
            .map(|c| translate_codon(std::str::from_utf8(c).unwrap()).unwrap())
            .collect();
        assert_eq!(back, protein);
        // M and W are single-codon
        assert!(dna.starts_with("ATG"));
        assert!(dna.ends_with("TGG"));
    }

    #[test]
    fn optimal_sequence_has_cai_one() {
        let u = illustrative_weights();
        // F, L, G are multi-codon; optimal picks weight-1.0 codons -> CAI 1.0
        let dna = reverse_translate_optimal("FLG", &u).unwrap();
        assert!((cai(&dna, &u).unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn suboptimal_codon_lowers_cai() {
        let u = illustrative_weights();
        // pick a deliberately non-optimal Phe codon (TTC, weight 0.5)
        assert!((u.weight("TTC").unwrap() - 0.5).abs() < 1e-12);
        let dna = "TTCTTC"; // two sub-optimal Phe codons
        let c = cai(dna, &u).unwrap();
        assert!((c - 0.5).abs() < 1e-12, "CAI = {c}");
    }

    #[test]
    fn gc_content_known_values() {
        assert!((gc_content("GGCC").unwrap() - 1.0).abs() < 1e-12);
        assert!(gc_content("ATAT").unwrap().abs() < 1e-12);
        assert!((gc_content("ATGC").unwrap() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn all_single_codon_sequence_has_cai_one() {
        let u = illustrative_weights();
        // Met+Trp only -> nothing counted -> CAI defined as 1.0
        assert!((cai("ATGTGG", &u).unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn rejects_bad_input() {
        let u = illustrative_weights();
        assert_eq!(
            reverse_translate_optimal("", &u).unwrap_err().code(),
            "empty"
        );
        assert_eq!(
            reverse_translate_optimal("MZ", &u).unwrap_err().code(),
            "invalid_residue"
        );
        assert_eq!(cai("ATGT", &u).unwrap_err().code(), "not_multiple_of_three");
        assert_eq!(gc_content("ATXG").unwrap_err().code(), "invalid_codon");
    }
}
