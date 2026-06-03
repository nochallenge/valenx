//! Reverse complement and complement for DNA and RNA.
//!
//! Both honor the IUPAC ambiguity codes via [`crate::alphabet`]'s
//! complement tables. Protein sequences have no complement and yield
//! [`BioseqError::Invalid`].

use crate::alphabet::{self, SeqKind};
use crate::error::{BioseqError, Result};
use crate::seq::Seq;

/// Returns the complement of a nucleotide sequence (same 5′→3′ order,
/// every base replaced by its IUPAC complement).
pub fn complement(seq: &Seq) -> Result<Seq> {
    let comp_fn = match seq.kind() {
        SeqKind::Dna => alphabet::complement_dna,
        SeqKind::Rna => alphabet::complement_rna,
        SeqKind::Protein => {
            return Err(BioseqError::invalid(
                "kind",
                "protein sequences have no complement",
            ))
        }
    };
    let out: Vec<u8> = seq
        .as_bytes()
        .iter()
        .map(|&b| comp_fn(b).unwrap_or(b))
        .collect();
    Ok(Seq::new_unchecked(seq.kind(), out, seq.topology()))
}

/// Returns the reverse complement — the complement read 3′→5′, i.e.
/// what a sequencer would read off the opposite strand.
pub fn reverse_complement(seq: &Seq) -> Result<Seq> {
    let comp_fn = match seq.kind() {
        SeqKind::Dna => alphabet::complement_dna,
        SeqKind::Rna => alphabet::complement_rna,
        SeqKind::Protein => {
            return Err(BioseqError::invalid(
                "kind",
                "protein sequences have no reverse complement",
            ))
        }
    };
    let out: Vec<u8> = seq
        .as_bytes()
        .iter()
        .rev()
        .map(|&b| comp_fn(b).unwrap_or(b))
        .collect();
    Ok(Seq::new_unchecked(seq.kind(), out, seq.topology()))
}

/// Reverse complement of raw bytes assumed to be uppercase DNA. A
/// thin helper for the hot paths (restriction digest, PCR) that work
/// on `&[u8]` rather than [`Seq`].
pub fn reverse_complement_dna_bytes(bytes: &[u8]) -> Vec<u8> {
    bytes
        .iter()
        .rev()
        .map(|&b| alphabet::complement_dna(b).unwrap_or(b))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dna_reverse_complement() {
        let s = Seq::new(SeqKind::Dna, "ATGC").unwrap();
        assert_eq!(complement(&s).unwrap().as_str(), "TACG");
        assert_eq!(reverse_complement(&s).unwrap().as_str(), "GCAT");
    }

    #[test]
    fn rna_reverse_complement() {
        let s = Seq::new(SeqKind::Rna, "AUGC").unwrap();
        assert_eq!(reverse_complement(&s).unwrap().as_str(), "GCAU");
    }

    #[test]
    fn ambiguity_codes_complement() {
        let s = Seq::new(SeqKind::Dna, "RYSWKMN").unwrap();
        // reverse: N M K W S Y R  -> complement each: N K M W S R Y
        assert_eq!(reverse_complement(&s).unwrap().as_str(), "NKMWSRY");
    }

    #[test]
    fn protein_has_no_complement() {
        let s = Seq::new(SeqKind::Protein, "MKVL").unwrap();
        assert!(complement(&s).is_err());
        assert!(reverse_complement(&s).is_err());
    }

    #[test]
    fn palindrome_is_its_own_revcomp() {
        let s = Seq::new(SeqKind::Dna, "GAATTC").unwrap(); // EcoRI site
        assert_eq!(reverse_complement(&s).unwrap().as_str(), "GAATTC");
    }

    #[test]
    fn byte_helper_matches() {
        assert_eq!(reverse_complement_dna_bytes(b"ATGC"), b"GCAT");
    }
}
