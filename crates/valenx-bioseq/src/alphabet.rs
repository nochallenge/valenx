//! IUPAC nucleotide and protein alphabets.
//!
//! This module is the foundation every other module builds on. It
//! defines the three sequence kinds Valenx handles ([`SeqKind`]),
//! validates residues against the IUPAC code tables (including the
//! ambiguity codes `N R Y S W K M B D H V`), and provides the
//! nucleotide complement table used by reverse-complement.
//!
//! References: IUPAC-IUB Joint Commission on Biochemical Nomenclature
//! (1984), "Nomenclature for incompletely specified bases in nucleic
//! acid sequences."

use crate::error::{BioseqError, Result};
use serde::{Deserialize, Serialize};

/// The kind of biological sequence.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub enum SeqKind {
    /// Deoxyribonucleic acid — alphabet `ACGT` + ambiguity codes.
    Dna,
    /// Ribonucleic acid — alphabet `ACGU` + ambiguity codes.
    Rna,
    /// Protein — the 20 standard amino acids + `X B Z J U O *`.
    Protein,
}

impl SeqKind {
    /// Human-readable name (`"DNA"`, `"RNA"`, `"protein"`).
    pub fn name(self) -> &'static str {
        match self {
            SeqKind::Dna => "DNA",
            SeqKind::Rna => "RNA",
            SeqKind::Protein => "protein",
        }
    }

    /// `true` for [`SeqKind::Dna`] / [`SeqKind::Rna`].
    pub fn is_nucleotide(self) -> bool {
        matches!(self, SeqKind::Dna | SeqKind::Rna)
    }

    /// The unambiguous "canonical" residues for this kind.
    pub fn canonical_residues(self) -> &'static [u8] {
        match self {
            SeqKind::Dna => b"ACGT",
            SeqKind::Rna => b"ACGU",
            SeqKind::Protein => b"ACDEFGHIKLMNPQRSTVWY",
        }
    }
}

/// Returns `true` if `b` is a valid (case-insensitive) residue for
/// `kind`, counting IUPAC ambiguity codes and the gap character `-`.
pub fn is_valid_residue(kind: SeqKind, b: u8) -> bool {
    let u = b.to_ascii_uppercase();
    match kind {
        SeqKind::Dna => matches!(
            u,
            b'A' | b'C' | b'G' | b'T'
                | b'N' | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M'
                | b'B' | b'D' | b'H' | b'V'
                | b'-'
        ),
        SeqKind::Rna => matches!(
            u,
            b'A' | b'C' | b'G' | b'U'
                | b'N' | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M'
                | b'B' | b'D' | b'H' | b'V'
                | b'-'
        ),
        // The 20 standard amino acids plus the IUPAC extras (B Z J X)
        // and the rare-residue codes (U = selenocysteine, O =
        // pyrrolysine) cover every letter A–Z, so accept any ASCII
        // letter plus the stop `*` and gap `-`.
        SeqKind::Protein => u.is_ascii_uppercase() || u == b'*' || u == b'-',
    }
}

/// Validates an entire sequence against `kind`. Returns `Ok(())` or the
/// first offending residue as a [`BioseqError::Alphabet`].
pub fn validate(kind: SeqKind, residues: &[u8]) -> Result<()> {
    for &b in residues {
        if !is_valid_residue(kind, b) {
            return Err(BioseqError::alphabet(b as char, kind.name()));
        }
    }
    Ok(())
}

/// Complement of a single DNA base (IUPAC-aware), preserving case.
///
/// Ambiguity codes complement to their reverse-pairing code: `R` (A/G)
/// → `Y` (C/T), `B` (C/G/T) → `V` (A/C/G), etc. `N` and `-` are their
/// own complement. Returns `None` for residues that are not DNA codes.
pub fn complement_dna(b: u8) -> Option<u8> {
    let upper = b.is_ascii_uppercase() || !b.is_ascii_alphabetic();
    let c = match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'U' => b'A',
        b'N' => b'N',
        b'R' => b'Y', // A/G  -> T/C
        b'Y' => b'R', // C/T  -> G/A
        b'S' => b'S', // C/G  -> G/C
        b'W' => b'W', // A/T  -> T/A
        b'K' => b'M', // G/T  -> C/A
        b'M' => b'K', // A/C  -> T/G
        b'B' => b'V', // C/G/T -> G/C/A
        b'V' => b'B', // A/C/G -> T/G/C
        b'D' => b'H', // A/G/T -> T/C/A
        b'H' => b'D', // A/C/T -> T/G/A
        b'-' => b'-',
        _ => return None,
    };
    Some(if upper { c } else { c.to_ascii_lowercase() })
}

/// Complement of a single RNA base — like [`complement_dna`] but the
/// `A` partner is `U` instead of `T`.
pub fn complement_rna(b: u8) -> Option<u8> {
    let upper = b.is_ascii_uppercase() || !b.is_ascii_alphabetic();
    let c = match b.to_ascii_uppercase() {
        b'A' => b'U',
        b'U' => b'A',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'N' => b'N',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'-' => b'-',
        _ => return None,
    };
    Some(if upper { c } else { c.to_ascii_lowercase() })
}

/// Expands an IUPAC nucleotide ambiguity code into the set of
/// canonical bases it represents (always uppercase `ACGT`; `U` is
/// folded to `T`). Returns `None` for non-nucleotide input.
///
/// Used by restriction-site matching and primer binding so an `N` in a
/// recognition site matches any base.
pub fn expand_iupac(b: u8) -> Option<&'static [u8]> {
    Some(match b.to_ascii_uppercase() {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' | b'U' => b"T",
        b'R' => b"AG",
        b'Y' => b"CT",
        b'S' => b"CG",
        b'W' => b"AT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'B' => b"CGT",
        b'D' => b"AGT",
        b'H' => b"ACT",
        b'V' => b"ACG",
        b'N' => b"ACGT",
        _ => return None,
    })
}

/// Returns `true` if a concrete base `base` is one of the bases the
/// IUPAC code `code` stands for (both nucleotide; case-insensitive;
/// `U`/`T` interchangeable).
pub fn iupac_matches(code: u8, base: u8) -> bool {
    match expand_iupac(code) {
        Some(set) => set.contains(&base.to_ascii_uppercase().pipe_t()),
        None => false,
    }
}

/// Tiny extension to fold `U` to `T` so `iupac_matches` works for RNA
/// concrete bases too.
trait PipeT {
    fn pipe_t(self) -> u8;
}
impl PipeT for u8 {
    fn pipe_t(self) -> u8 {
        if self == b'U' {
            b'T'
        } else {
            self
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_validation() {
        assert!(validate(SeqKind::Dna, b"ACGTACGT").is_ok());
        assert!(validate(SeqKind::Rna, b"ACGUACGU").is_ok());
        assert!(validate(SeqKind::Protein, b"MKVLAAA").is_ok());
    }

    #[test]
    fn ambiguity_codes_accepted() {
        assert!(validate(SeqKind::Dna, b"ACGTNRYSWKMBDHV").is_ok());
        assert!(validate(SeqKind::Rna, b"ACGUNRYSWKMBDHV").is_ok());
    }

    #[test]
    fn rejects_wrong_alphabet() {
        let e = validate(SeqKind::Dna, b"ACGU").unwrap_err();
        assert!(matches!(e, BioseqError::Alphabet { .. }));
        assert!(validate(SeqKind::Rna, b"ACGT").is_err());
        assert!(validate(SeqKind::Protein, b"MKVL123").is_err());
    }

    #[test]
    fn dna_complement_table() {
        assert_eq!(complement_dna(b'A'), Some(b'T'));
        assert_eq!(complement_dna(b'g'), Some(b'c'));
        assert_eq!(complement_dna(b'N'), Some(b'N'));
        assert_eq!(complement_dna(b'R'), Some(b'Y'));
        assert_eq!(complement_dna(b'B'), Some(b'V'));
        assert_eq!(complement_dna(b'Z'), None);
    }

    #[test]
    fn rna_complement_table() {
        assert_eq!(complement_rna(b'A'), Some(b'U'));
        assert_eq!(complement_rna(b'U'), Some(b'A'));
        assert_eq!(complement_rna(b'g'), Some(b'c'));
    }

    #[test]
    fn iupac_expansion_and_match() {
        assert_eq!(expand_iupac(b'N'), Some(&b"ACGT"[..]));
        assert_eq!(expand_iupac(b'R'), Some(&b"AG"[..]));
        assert!(iupac_matches(b'N', b'G'));
        assert!(iupac_matches(b'R', b'A'));
        assert!(!iupac_matches(b'R', b'C'));
        assert!(iupac_matches(b'W', b'T'));
        assert!(iupac_matches(b'A', b'a'));
    }

    #[test]
    fn kind_helpers() {
        assert!(SeqKind::Dna.is_nucleotide());
        assert!(!SeqKind::Protein.is_nucleotide());
        assert_eq!(SeqKind::Rna.name(), "RNA");
        assert_eq!(SeqKind::Dna.canonical_residues(), b"ACGT");
    }
}
