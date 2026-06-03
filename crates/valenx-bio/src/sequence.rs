//! Canonical biological sequence: a name, an alphabet, and a
//! validated character buffer. All bytes are normalised to uppercase
//! on construction so equality compares cleanly across sources that
//! mix case (FASTA files, NCBI dumps, hand-typed input).

use serde::{Deserialize, Serialize};
use thiserror::Error;

pub use crate::alphabet::Alphabet;

/// Construction errors from [`Sequence::new`].
#[derive(Debug, Error, PartialEq, Eq)]
pub enum SequenceError {
    /// A byte of the input is not valid in the requested alphabet (or
    /// is outside the ASCII range).
    #[error("invalid byte 0x{byte:02x} ({char:?}) for {alphabet} at position {position}")]
    InvalidByte {
        /// The offending byte (after uppercasing).
        byte: u8,
        /// The original character.
        char: char,
        /// The alphabet's identifier ([`Alphabet::id`]).
        alphabet: &'static str,
        /// Zero-based character position in the input.
        position: usize,
    },
}

/// A name + alphabet + validated character buffer.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Sequence {
    /// Display name from the source (e.g. FASTA `>` line, NCBI accession).
    pub name: String,
    /// The sequence's residue alphabet — fixes the validation rules.
    pub alphabet: Alphabet,
    /// Uppercase-normalised characters. Always `alphabet.is_valid(b)` for
    /// every byte b; the constructor enforces this invariant.
    bytes: String,
}

impl Sequence {
    /// Build a new sequence, validating every byte against the alphabet.
    /// Returns `SequenceError::InvalidByte` at the first invalid position.
    pub fn new(
        name: impl Into<String>,
        alphabet: Alphabet,
        text: &str,
    ) -> Result<Self, SequenceError> {
        let mut bytes = String::with_capacity(text.len());
        for (i, c) in text.chars().enumerate() {
            let upper = c.to_ascii_uppercase();
            let b = upper as u32;
            if b > 0x7F {
                return Err(SequenceError::InvalidByte {
                    byte: 0,
                    char: c,
                    alphabet: alphabet.id(),
                    position: i,
                });
            }
            let byte = upper as u8;
            if !alphabet.is_valid(byte) {
                return Err(SequenceError::InvalidByte {
                    byte,
                    char: c,
                    alphabet: alphabet.id(),
                    position: i,
                });
            }
            bytes.push(upper);
        }
        Ok(Self {
            name: name.into(),
            alphabet,
            bytes,
        })
    }

    /// Length in characters (equivalent to bytes here — only ASCII is
    /// stored).
    pub fn len(&self) -> usize {
        self.bytes.len()
    }

    /// `true` if the sequence has zero characters.
    pub fn is_empty(&self) -> bool {
        self.bytes.is_empty()
    }

    /// Borrow the uppercase, validated character buffer.
    pub fn as_str(&self) -> &str {
        &self.bytes
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_validates_input() {
        // Valid DNA accepts.
        let s = Sequence::new("p53", Alphabet::Dna, "ACGTACGTN").unwrap();
        assert_eq!(s.len(), 9);
        // Invalid DNA rejects with the offending position. `E`
        // (glutamate) is a canonical amino acid that is NOT a valid
        // IUPAC nucleotide code, so it's a clean "protein-only"
        // residue for this assertion. (Don't use `M` here — it's a
        // valid IUPAC ambiguity code for `aMino` = A or C.)
        let err = Sequence::new("bad", Alphabet::Dna, "ACGTE").unwrap_err();
        assert!(matches!(
            err,
            SequenceError::InvalidByte { position: 4, .. }
        ));
    }

    #[test]
    fn empty_sequence_is_valid() {
        let s = Sequence::new("empty", Alphabet::Protein, "").unwrap();
        assert_eq!(s.len(), 0);
        assert!(s.is_empty());
    }

    #[test]
    fn case_insensitive_storage_normalises_to_upper() {
        let s = Sequence::new("p53", Alphabet::Dna, "acgt").unwrap();
        assert_eq!(s.as_str(), "ACGT");
    }

    #[test]
    fn serde_round_trip() {
        let s = Sequence::new("ubiq", Alphabet::Protein, "MKLI").unwrap();
        let json = serde_json::to_string(&s).unwrap();
        let back: Sequence = serde_json::from_str(&json).unwrap();
        assert_eq!(back, s);
    }
}
