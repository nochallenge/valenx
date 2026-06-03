//! [`RnaSeq`] — the validated RNA sequence every folder consumes.
//!
//! The folding, ensemble and interaction algorithms all need a
//! sequence over exactly the four bases `A C G U` encoded `0..4`. This
//! module wraps that representation and provides one set of
//! constructors:
//!
//! - [`RnaSeq::parse`] — from an ASCII string (`A C G U`, also accepts
//!   `T`, folding it to `U`, and lowercase).
//! - [`RnaSeq::from_seq`] — from a [`valenx_bioseq::Seq`], the
//!   Round 6 Block 1 sequence type. A DNA `Seq` is transcribed; an RNA
//!   `Seq` is used directly; a protein `Seq` is rejected.
//!
//! Ambiguity codes and gaps are rejected: a thermodynamic folder
//! cannot assign an energy to an `N`.

use crate::error::{Result, RnaStructError};
use crate::fold::energy::encode_seq;
use valenx_bioseq::{Seq, SeqKind};

/// A validated RNA sequence over `A C G U`.
///
/// Stores both the original uppercased ASCII bytes (for display and
/// I/O) and the internal `0..4` codes (for the energy model).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RnaSeq {
    ascii: Vec<u8>,
    codes: Vec<u8>,
}

impl RnaSeq {
    /// Parses an ASCII RNA sequence. Whitespace is stripped; `T`/`t`
    /// are folded to `U`; the result is uppercased.
    ///
    /// # Errors
    /// [`RnaStructError::Sequence`] if the input is empty or contains
    /// any character outside `A C G U T` (case-insensitive).
    pub fn parse(s: impl AsRef<[u8]>) -> Result<Self> {
        let ascii: Vec<u8> = s
            .as_ref()
            .iter()
            .filter(|b| !b.is_ascii_whitespace())
            .map(|b| {
                let u = b.to_ascii_uppercase();
                if u == b'T' {
                    b'U'
                } else {
                    u
                }
            })
            .collect();
        if ascii.is_empty() {
            return Err(RnaStructError::sequence("sequence is empty"));
        }
        let codes = encode_seq(&ascii)
            .map_err(|b| RnaStructError::sequence(format!("illegal base `{}`", b as char)))?;
        Ok(RnaSeq { ascii, codes })
    }

    /// Builds an [`RnaSeq`] from a [`valenx_bioseq::Seq`].
    ///
    /// A [`SeqKind::Rna`] sequence is taken as-is; a [`SeqKind::Dna`]
    /// sequence is transcribed (`T` → `U`); a protein sequence is
    /// rejected. Ambiguity codes still cause a [`RnaStructError`].
    ///
    /// # Errors
    /// [`RnaStructError::Sequence`] for a protein input or an
    /// ambiguous / empty nucleotide sequence.
    pub fn from_seq(seq: &Seq) -> Result<Self> {
        match seq.kind() {
            SeqKind::Protein => Err(RnaStructError::sequence(
                "cannot fold a protein sequence",
            )),
            SeqKind::Dna | SeqKind::Rna => Self::parse(seq.as_bytes()),
        }
    }

    /// Number of bases.
    pub fn len(&self) -> usize {
        self.codes.len()
    }

    /// `true` if the sequence has no bases.
    pub fn is_empty(&self) -> bool {
        self.codes.is_empty()
    }

    /// The internal `0..4` base codes (`A C G U`).
    pub fn codes(&self) -> &[u8] {
        &self.codes
    }

    /// The uppercased ASCII bytes.
    pub fn as_bytes(&self) -> &[u8] {
        &self.ascii
    }

    /// The sequence as a `&str`.
    pub fn as_str(&self) -> &str {
        // SAFETY: bytes are validated A/C/G/U on construction.
        std::str::from_utf8(&self.ascii).expect("ACGU is ASCII")
    }

    /// Concatenates two RNA sequences (used by cofolding to fold a
    /// duplex as one composite strand).
    pub fn concat(&self, other: &RnaSeq) -> RnaSeq {
        let mut ascii = self.ascii.clone();
        ascii.extend_from_slice(&other.ascii);
        let mut codes = self.codes.clone();
        codes.extend_from_slice(&other.codes);
        RnaSeq { ascii, codes }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_basic() {
        let r = RnaSeq::parse("GGGAAACCC").unwrap();
        assert_eq!(r.len(), 9);
        assert_eq!(r.as_str(), "GGGAAACCC");
    }

    #[test]
    fn parse_folds_t_and_lowercase() {
        let r = RnaSeq::parse("gattaca").unwrap();
        assert_eq!(r.as_str(), "GAUUACA");
    }

    #[test]
    fn parse_rejects_garbage() {
        assert!(RnaSeq::parse("").is_err());
        assert!(RnaSeq::parse("ACGN").is_err());
        assert!(RnaSeq::parse("ACG-").is_err());
    }

    #[test]
    fn from_bioseq_seq() {
        let dna = Seq::new(SeqKind::Dna, "ATGC").unwrap();
        assert_eq!(RnaSeq::from_seq(&dna).unwrap().as_str(), "AUGC");
        let rna = Seq::new(SeqKind::Rna, "ACGU").unwrap();
        assert_eq!(RnaSeq::from_seq(&rna).unwrap().as_str(), "ACGU");
        let prot = Seq::new(SeqKind::Protein, "MK").unwrap();
        assert!(RnaSeq::from_seq(&prot).is_err());
    }

    #[test]
    fn concat_joins() {
        let a = RnaSeq::parse("GGG").unwrap();
        let b = RnaSeq::parse("CCC").unwrap();
        assert_eq!(a.concat(&b).as_str(), "GGGCCC");
    }
}
