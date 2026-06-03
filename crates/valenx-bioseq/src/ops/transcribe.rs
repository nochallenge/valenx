//! Transcription DNA→RNA and back-transcription RNA→DNA.
//!
//! Following the Biopython convention, [`transcribe`] takes the
//! *coding* strand of DNA and produces the mRNA — the operation is
//! simply `T → U`. [`back_transcribe`] is the inverse, `U → T`.
//!
//! (Real transcription reads the template strand; the coding strand is
//! identical to the mRNA except for `T`/`U`, so the byte-level
//! operation is the same — this is the standard library convention.)

use crate::alphabet::SeqKind;
use crate::error::{BioseqError, Result};
use crate::seq::Seq;

/// Transcribes a coding-strand DNA sequence into mRNA (`T → U`).
///
/// Returns [`BioseqError::Invalid`] if `seq` is not DNA.
pub fn transcribe(seq: &Seq) -> Result<Seq> {
    if seq.kind() != SeqKind::Dna {
        return Err(BioseqError::invalid(
            "kind",
            format!("transcribe expects DNA, got {}", seq.kind().name()),
        ));
    }
    let out: Vec<u8> = seq
        .as_bytes()
        .iter()
        .map(|&b| match b {
            b'T' => b'U',
            b't' => b'u',
            other => other,
        })
        .collect();
    Ok(Seq::new_unchecked(SeqKind::Rna, out, seq.topology()))
}

/// Back-transcribes mRNA into coding-strand DNA (`U → T`).
///
/// Returns [`BioseqError::Invalid`] if `seq` is not RNA.
pub fn back_transcribe(seq: &Seq) -> Result<Seq> {
    if seq.kind() != SeqKind::Rna {
        return Err(BioseqError::invalid(
            "kind",
            format!("back_transcribe expects RNA, got {}", seq.kind().name()),
        ));
    }
    let out: Vec<u8> = seq
        .as_bytes()
        .iter()
        .map(|&b| match b {
            b'U' => b'T',
            b'u' => b't',
            other => other,
        })
        .collect();
    Ok(Seq::new_unchecked(SeqKind::Dna, out, seq.topology()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dna_to_rna() {
        let dna = Seq::new(SeqKind::Dna, "ATGGCCTAA").unwrap();
        let rna = transcribe(&dna).unwrap();
        assert_eq!(rna.kind(), SeqKind::Rna);
        assert_eq!(rna.as_str(), "AUGGCCUAA");
    }

    #[test]
    fn rna_to_dna_roundtrip() {
        let dna = Seq::new(SeqKind::Dna, "ATGCGTACGT").unwrap();
        let rna = transcribe(&dna).unwrap();
        let back = back_transcribe(&rna).unwrap();
        assert_eq!(back.as_str(), dna.as_str());
        assert_eq!(back.kind(), SeqKind::Dna);
    }

    #[test]
    fn wrong_kinds_rejected() {
        let rna = Seq::new(SeqKind::Rna, "AUGC").unwrap();
        assert!(transcribe(&rna).is_err());
        let dna = Seq::new(SeqKind::Dna, "ATGC").unwrap();
        assert!(back_transcribe(&dna).is_err());
    }

    #[test]
    fn ambiguity_codes_pass_through() {
        let dna = Seq::new(SeqKind::Dna, "ATGNRY").unwrap();
        assert_eq!(transcribe(&dna).unwrap().as_str(), "AUGNRY");
    }
}
