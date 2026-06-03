//! The [`Seq`] type — a validated biological sequence.
//!
//! A `Seq` carries its [`SeqKind`] (DNA / RNA / protein), its
//! [`Topology`] (linear or circular — circular matters for plasmids
//! and bacterial chromosomes), and its residues as raw ASCII bytes
//! (uppercased on construction). It is the Biopython `Seq` analogue.

use crate::alphabet;
use crate::error::{BioseqError, Result};
use serde::{Deserialize, Serialize};
use std::fmt;

// Re-export so the canonical `crate::seq::SeqKind` path resolves for
// the modules that pair `SeqKind` with `Seq` / `Topology`.
pub use crate::alphabet::SeqKind;

/// Linear vs. circular sequence topology.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize, Default)]
pub enum Topology {
    /// A linear molecule with distinct 5′ and 3′ ends.
    #[default]
    Linear,
    /// A covalently closed circular molecule (plasmid, bacterial
    /// chromosome). Slicing and feature coordinates may wrap.
    Circular,
}

/// A validated biological sequence.
///
/// Residues are stored uppercased; construction validates every
/// residue against the IUPAC table for the declared [`SeqKind`].
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Seq {
    kind: SeqKind,
    topology: Topology,
    residues: Vec<u8>,
}

impl Seq {
    /// Builds a linear sequence, validating every residue.
    ///
    /// Whitespace in the input is stripped; remaining residues are
    /// uppercased. Returns [`BioseqError::Alphabet`] on the first
    /// illegal residue.
    pub fn new(kind: SeqKind, residues: impl AsRef<[u8]>) -> Result<Self> {
        Self::with_topology(kind, residues, Topology::Linear)
    }

    /// Like [`Seq::new`] but with an explicit [`Topology`].
    pub fn with_topology(
        kind: SeqKind,
        residues: impl AsRef<[u8]>,
        topology: Topology,
    ) -> Result<Self> {
        let cleaned: Vec<u8> = residues
            .as_ref()
            .iter()
            .filter(|b| !b.is_ascii_whitespace())
            .map(|b| b.to_ascii_uppercase())
            .collect();
        alphabet::validate(kind, &cleaned)?;
        Ok(Seq {
            kind,
            topology,
            residues: cleaned,
        })
    }

    /// Builds a sequence without validation. The caller asserts that
    /// every byte is a legal residue for `kind`. Used internally by
    /// translation (output is protein, already guaranteed valid).
    pub(crate) fn new_unchecked(kind: SeqKind, residues: Vec<u8>, topology: Topology) -> Self {
        Seq {
            kind,
            topology,
            residues,
        }
    }

    /// The sequence kind.
    pub fn kind(&self) -> SeqKind {
        self.kind
    }

    /// The topology.
    pub fn topology(&self) -> Topology {
        self.topology
    }

    /// `true` if the sequence is circular.
    pub fn is_circular(&self) -> bool {
        self.topology == Topology::Circular
    }

    /// Returns a copy with the topology set to `topology`.
    pub fn with_topology_set(&self, topology: Topology) -> Self {
        Seq {
            topology,
            ..self.clone()
        }
    }

    /// Number of residues.
    pub fn len(&self) -> usize {
        self.residues.len()
    }

    /// `true` if the sequence has no residues.
    pub fn is_empty(&self) -> bool {
        self.residues.is_empty()
    }

    /// The residues as a byte slice (uppercased).
    pub fn as_bytes(&self) -> &[u8] {
        &self.residues
    }

    /// The residues as a `&str` (always valid UTF-8 — IUPAC codes are
    /// ASCII).
    pub fn as_str(&self) -> &str {
        // SAFETY: residues are validated ASCII on construction.
        std::str::from_utf8(&self.residues).expect("IUPAC residues are ASCII")
    }

    /// The residue at index `i`, or `None` if out of range. For a
    /// circular sequence, indices wrap modulo the length.
    pub fn get(&self, i: usize) -> Option<u8> {
        if self.residues.is_empty() {
            return None;
        }
        match self.topology {
            Topology::Linear => self.residues.get(i).copied(),
            Topology::Circular => Some(self.residues[i % self.residues.len()]),
        }
    }

    /// A half-open `[start, end)` slice as a new linear `Seq`.
    ///
    /// For a circular sequence, `end <= len` is required but a slice
    /// that wraps past the origin should use [`Seq::slice_circular`].
    /// Returns [`BioseqError::Invalid`] if the range is out of bounds
    /// or inverted.
    pub fn slice(&self, start: usize, end: usize) -> Result<Seq> {
        if start > end {
            return Err(BioseqError::invalid(
                "range",
                format!("start {start} > end {end}"),
            ));
        }
        if end > self.residues.len() {
            return Err(BioseqError::invalid(
                "range",
                format!("end {end} exceeds length {}", self.residues.len()),
            ));
        }
        Ok(Seq {
            kind: self.kind,
            topology: Topology::Linear,
            residues: self.residues[start..end].to_vec(),
        })
    }

    /// A slice of a circular sequence that may wrap past the origin.
    ///
    /// `start` and `end` are interpreted modulo the length; the slice
    /// is `length` residues long walking forward from `start`. The
    /// result is linear. Returns [`BioseqError::Invalid`] for an empty
    /// sequence or a `length` exceeding the molecule size.
    pub fn slice_circular(&self, start: usize, length: usize) -> Result<Seq> {
        let n = self.residues.len();
        if n == 0 {
            return Err(BioseqError::invalid("sequence", "empty"));
        }
        if length > n {
            return Err(BioseqError::invalid(
                "length",
                format!("slice length {length} exceeds molecule size {n}"),
            ));
        }
        let mut out = Vec::with_capacity(length);
        for k in 0..length {
            out.push(self.residues[(start + k) % n]);
        }
        Ok(Seq {
            kind: self.kind,
            topology: Topology::Linear,
            residues: out,
        })
    }

    /// Concatenates two sequences of the same kind into a new linear
    /// `Seq`. Returns [`BioseqError::Invalid`] on a kind mismatch.
    pub fn concat(&self, other: &Seq) -> Result<Seq> {
        if self.kind != other.kind {
            return Err(BioseqError::invalid(
                "kind",
                format!(
                    "cannot concat {} with {}",
                    self.kind.name(),
                    other.kind.name()
                ),
            ));
        }
        let mut residues = self.residues.clone();
        residues.extend_from_slice(&other.residues);
        Ok(Seq {
            kind: self.kind,
            topology: Topology::Linear,
            residues,
        })
    }

    /// Returns a copy rotated so that index `origin` becomes position
    /// 0. Only meaningful for circular sequences but works on linear
    /// ones too (treats them as if circular for the rotation). The
    /// result keeps the original topology.
    pub fn rotate(&self, origin: usize) -> Result<Seq> {
        let n = self.residues.len();
        if n == 0 {
            return Ok(self.clone());
        }
        let o = origin % n;
        let mut residues = Vec::with_capacity(n);
        residues.extend_from_slice(&self.residues[o..]);
        residues.extend_from_slice(&self.residues[..o]);
        Ok(Seq {
            kind: self.kind,
            topology: self.topology,
            residues,
        })
    }

    /// Counts occurrences of residue `b` (case-insensitive).
    pub fn count(&self, b: u8) -> usize {
        let u = b.to_ascii_uppercase();
        self.residues.iter().filter(|&&r| r == u).count()
    }

    /// An iterator over the residues.
    pub fn iter(&self) -> impl Iterator<Item = u8> + '_ {
        self.residues.iter().copied()
    }
}

impl fmt::Display for Seq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn construct_and_validate() {
        let s = Seq::new(SeqKind::Dna, "acgtACGT").unwrap();
        assert_eq!(s.len(), 8);
        assert_eq!(s.as_str(), "ACGTACGT");
        assert!(Seq::new(SeqKind::Dna, "ACGU").is_err());
    }

    #[test]
    fn whitespace_is_stripped() {
        let s = Seq::new(SeqKind::Dna, "AC GT\nAC\tGT").unwrap();
        assert_eq!(s.as_str(), "ACGTACGT");
    }

    #[test]
    fn slice_and_concat() {
        let a = Seq::new(SeqKind::Dna, "AAATTT").unwrap();
        assert_eq!(a.slice(1, 4).unwrap().as_str(), "AAT");
        assert!(a.slice(4, 2).is_err());
        assert!(a.slice(0, 99).is_err());
        let b = Seq::new(SeqKind::Dna, "GGG").unwrap();
        assert_eq!(a.concat(&b).unwrap().as_str(), "AAATTTGGG");
        let p = Seq::new(SeqKind::Protein, "MK").unwrap();
        assert!(a.concat(&p).is_err());
    }

    #[test]
    fn circular_indexing_and_rotation() {
        let c = Seq::with_topology(SeqKind::Dna, "ATGCAA", Topology::Circular).unwrap();
        assert!(c.is_circular());
        assert_eq!(c.get(6), Some(b'A')); // wraps to index 0
        assert_eq!(c.get(7), Some(b'T'));
        let r = c.rotate(2).unwrap();
        assert_eq!(r.as_str(), "GCAAAT");
        assert!(r.is_circular());
    }

    #[test]
    fn circular_slice_wraps() {
        let c = Seq::with_topology(SeqKind::Dna, "ATGCAA", Topology::Circular).unwrap();
        // start at index 4, take 4 residues: A A A T -> wraps the origin
        assert_eq!(c.slice_circular(4, 4).unwrap().as_str(), "AAAT");
        assert!(c.slice_circular(0, 99).is_err());
    }

    #[test]
    fn count_residues() {
        let s = Seq::new(SeqKind::Dna, "AAGGCC").unwrap();
        assert_eq!(s.count(b'a'), 2);
        assert_eq!(s.count(b'G'), 2);
        assert_eq!(s.count(b'T'), 0);
    }
}
