//! Molecular weight of DNA, RNA and protein sequences.
//!
//! - **Nucleic acids** — summed nucleotide-monophosphate residue
//!   masses minus the water lost per phosphodiester bond, plus one
//!   correction. The single-stranded weight is returned;
//!   [`molecular_weight_dna_double_stranded`] doubles + accounts for
//!   the complement.
//! - **Protein** — summed amino-acid residue masses plus one water
//!   (for the terminal `H` + `OH`). Both **average** and
//!   **monoisotopic** masses are provided.
//!
//! Masses are in daltons (g/mol).

use crate::error::{BioseqError, Result};
use crate::ops::revcomp::reverse_complement;
use crate::seq::{Seq, SeqKind};

/// Average residue mass of a DNA deoxynucleotide-monophosphate, Da.
/// These are the masses *in a chain* (the free dNMP minus one water).
fn dna_residue_mass(b: u8) -> Option<f64> {
    Some(match b {
        b'A' => 313.21,
        b'C' => 289.18,
        b'G' => 329.21,
        b'T' => 304.20,
        _ => return None,
    })
}

/// Average residue mass of an RNA ribonucleotide-monophosphate, Da.
fn rna_residue_mass(b: u8) -> Option<f64> {
    Some(match b {
        b'A' => 329.21,
        b'C' => 305.18,
        b'G' => 345.21,
        b'U' => 306.17,
        _ => return None,
    })
}

/// Average and monoisotopic residue masses of the 20 standard amino
/// acids (the mass contributed *in a peptide chain*, i.e. the free
/// amino acid minus one water). Returns `(average, monoisotopic)`.
fn aa_residue_mass(b: u8) -> Option<(f64, f64)> {
    // Average masses: IUPAC; monoisotopic: Unimod / standard tables.
    Some(match b {
        b'A' => (71.0788, 71.03711),
        b'R' => (156.1875, 156.10111),
        b'N' => (114.1038, 114.04293),
        b'D' => (115.0886, 115.02694),
        b'C' => (103.1388, 103.00919),
        b'E' => (129.1155, 129.04259),
        b'Q' => (128.1307, 128.05858),
        b'G' => (57.0519, 57.02146),
        b'H' => (137.1411, 137.05891),
        b'I' => (113.1594, 113.08406),
        b'L' => (113.1594, 113.08406),
        b'K' => (128.1741, 128.09496),
        b'M' => (131.1926, 131.04049),
        b'F' => (147.1766, 147.06841),
        b'P' => (97.1167, 97.05276),
        b'S' => (87.0782, 87.03203),
        b'T' => (101.1051, 101.04768),
        b'W' => (186.2132, 186.07931),
        b'Y' => (163.1760, 163.06333),
        b'V' => (99.1326, 99.06841),
        _ => return None,
    })
}

/// Mass of one water molecule, Da.
const WATER_AVG: f64 = 18.0153;
/// Monoisotopic mass of one water molecule, Da.
const WATER_MONO: f64 = 18.010565;

/// Single-stranded molecular weight of a nucleic-acid sequence, Da.
///
/// Uses average residue masses; the standard nucleic-acid formula
/// `Σ(residue masses) + 18.02` (one extra water for the free 5′/3′
/// termini, the OligoCalc convention). Returns
/// [`BioseqError::Invalid`] for a protein input or any non-canonical
/// base.
pub fn molecular_weight_nucleic(seq: &Seq) -> Result<f64> {
    let mass_fn: fn(u8) -> Option<f64> = match seq.kind() {
        SeqKind::Dna => dna_residue_mass,
        SeqKind::Rna => rna_residue_mass,
        SeqKind::Protein => {
            return Err(BioseqError::invalid(
                "kind",
                "use molecular_weight_protein for proteins",
            ))
        }
    };
    let mut sum = 0.0;
    for b in seq.iter() {
        match mass_fn(b) {
            Some(m) => sum += m,
            None => {
                return Err(BioseqError::invalid(
                    "sequence",
                    format!("non-canonical base `{}` has no tabulated mass", b as char),
                ))
            }
        }
    }
    if seq.is_empty() {
        return Ok(0.0);
    }
    Ok(sum + WATER_AVG)
}

/// Molecular weight of a double-stranded DNA molecule, Da — the
/// single-stranded weight of `seq` plus that of its reverse
/// complement.
pub fn molecular_weight_dna_double_stranded(seq: &Seq) -> Result<f64> {
    if seq.kind() != SeqKind::Dna {
        return Err(BioseqError::invalid("kind", "double-stranded weight needs DNA"));
    }
    let top = molecular_weight_nucleic(seq)?;
    let bottom = molecular_weight_nucleic(&reverse_complement(seq)?)?;
    Ok(top + bottom)
}

/// Average molecular weight of a protein, Da.
///
/// `Σ(average residue masses) + 18.0153` (one water for the chain
/// termini). Returns [`BioseqError::Invalid`] for a non-protein input
/// or any residue outside the 20 standard amino acids.
pub fn molecular_weight_protein(seq: &Seq) -> Result<f64> {
    protein_mass(seq, false)
}

/// Monoisotopic molecular weight of a protein, Da — like
/// [`molecular_weight_protein`] but using monoisotopic residue masses.
pub fn molecular_weight_protein_monoisotopic(seq: &Seq) -> Result<f64> {
    protein_mass(seq, true)
}

fn protein_mass(seq: &Seq, monoisotopic: bool) -> Result<f64> {
    if seq.kind() != SeqKind::Protein {
        return Err(BioseqError::invalid("kind", "expected a protein sequence"));
    }
    let mut sum = 0.0;
    for b in seq.iter() {
        if b == b'*' {
            continue; // a stop codon contributes nothing
        }
        match aa_residue_mass(b) {
            Some((avg, mono)) => sum += if monoisotopic { mono } else { avg },
            None => {
                return Err(BioseqError::invalid(
                    "sequence",
                    format!("residue `{}` is not a standard amino acid", b as char),
                ))
            }
        }
    }
    if sum == 0.0 {
        return Ok(0.0);
    }
    Ok(sum + if monoisotopic { WATER_MONO } else { WATER_AVG })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dna_weight_is_positive_and_scales() {
        let short = Seq::new(SeqKind::Dna, "ACGT").unwrap();
        let long = Seq::new(SeqKind::Dna, "ACGTACGT").unwrap();
        let w_short = molecular_weight_nucleic(&short).unwrap();
        let w_long = molecular_weight_nucleic(&long).unwrap();
        assert!(w_short > 0.0);
        assert!(w_long > w_short);
    }

    #[test]
    fn dna_weight_known_value() {
        // ACGT: 313.21+289.18+329.21+304.20 + 18.02 = 1253.82.
        let s = Seq::new(SeqKind::Dna, "ACGT").unwrap();
        let w = molecular_weight_nucleic(&s).unwrap();
        assert!((w - 1253.82).abs() < 0.1, "got {w}");
    }

    #[test]
    fn double_stranded_doubles_a_palindrome() {
        // For a palindrome, top == bottom, so ds == 2*ss.
        let s = Seq::new(SeqKind::Dna, "GAATTC").unwrap();
        let ss = molecular_weight_nucleic(&s).unwrap();
        let ds = molecular_weight_dna_double_stranded(&s).unwrap();
        assert!((ds - 2.0 * ss).abs() < 1e-6);
    }

    #[test]
    fn rna_weight_differs_from_dna() {
        let dna = Seq::new(SeqKind::Dna, "ACGA").unwrap();
        let rna = Seq::new(SeqKind::Rna, "ACGA").unwrap();
        // RNA residues are heavier (extra 2'-OH).
        assert!(
            molecular_weight_nucleic(&rna).unwrap()
                > molecular_weight_nucleic(&dna).unwrap()
        );
    }

    #[test]
    fn protein_weight_known_value() {
        // Glycine residue 57.0519 + water 18.0153 = 75.0672 (free Gly).
        let g = Seq::new(SeqKind::Protein, "G").unwrap();
        let w = molecular_weight_protein(&g).unwrap();
        assert!((w - 75.0672).abs() < 0.01, "got {w}");
    }

    #[test]
    fn monoisotopic_lighter_than_average() {
        let p = Seq::new(SeqKind::Protein, "MKVLAAGGWWYY").unwrap();
        let avg = molecular_weight_protein(&p).unwrap();
        let mono = molecular_weight_protein_monoisotopic(&p).unwrap();
        assert!(mono < avg, "mono {mono} should be < avg {avg}");
    }

    #[test]
    fn stop_codon_ignored_in_protein_mass() {
        let with = Seq::new(SeqKind::Protein, "MKV*").unwrap();
        let without = Seq::new(SeqKind::Protein, "MKV").unwrap();
        assert!(
            (molecular_weight_protein(&with).unwrap()
                - molecular_weight_protein(&without).unwrap())
            .abs()
                < 1e-9
        );
    }

    #[test]
    fn wrong_kind_rejected() {
        let p = Seq::new(SeqKind::Protein, "MK").unwrap();
        assert!(molecular_weight_nucleic(&p).is_err());
        let d = Seq::new(SeqKind::Dna, "ACGT").unwrap();
        assert!(molecular_weight_protein(&d).is_err());
    }

    #[test]
    fn empty_sequences_weigh_zero() {
        let d = Seq::new(SeqKind::Dna, "").unwrap();
        assert_eq!(molecular_weight_nucleic(&d).unwrap(), 0.0);
        let p = Seq::new(SeqKind::Protein, "").unwrap();
        assert_eq!(molecular_weight_protein(&p).unwrap(), 0.0);
    }
}
