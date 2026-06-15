//! The position-specific scoring matrix and peptide scoring.

use serde::{Deserialize, Serialize};

use crate::aa::{aa_index, N_AA};
use crate::error::ImmunoError;

/// A position-specific scoring matrix (PSSM): one real weight for every
/// (position, residue) pair over a fixed window length.
///
/// `score(peptide) = sum over positions i of weight[i][residue_i]`. Higher is a
/// stronger predicted binder. Build one with [`Pssm::new`] (validated) or use a
/// ready-made illustrative matrix from [`crate::library`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Pssm {
    name: String,
    /// `weights[position][residue_column]`; one row per window position.
    weights: Vec<[f64; N_AA]>,
}

impl Pssm {
    /// Build a PSSM from per-position residue-weight rows.
    ///
    /// Each row holds one weight per residue, indexed by
    /// [`aa_index`]. Rejects an empty matrix
    /// ([`ImmunoError::Empty`]) and any non-finite weight
    /// ([`ImmunoError::NonFiniteWeight`]).
    pub fn new(name: impl Into<String>, weights: Vec<[f64; N_AA]>) -> Result<Self, ImmunoError> {
        if weights.is_empty() {
            return Err(ImmunoError::Empty { what: "matrix" });
        }
        for (pos, row) in weights.iter().enumerate() {
            for (aa, &w) in row.iter().enumerate() {
                if !w.is_finite() {
                    return Err(ImmunoError::NonFiniteWeight { pos, aa });
                }
            }
        }
        Ok(Self {
            name: name.into(),
            weights,
        })
    }

    /// The window length: the number of residues a scored peptide must have.
    pub fn length(&self) -> usize {
        self.weights.len()
    }

    /// The matrix's name / allele label.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// The stored weight for `residue` at `pos`, or `None` if `pos` is out of
    /// range or `residue` is not a standard amino acid.
    pub fn weight(&self, pos: usize, residue: u8) -> Option<f64> {
        let row = self.weights.get(pos)?;
        Some(row[aa_index(residue)?])
    }

    /// Score a peptide by summing its residues' per-position weights.
    ///
    /// The peptide must be exactly [`Pssm::length`] residues long
    /// ([`ImmunoError::LengthMismatch`]) and contain only standard amino acids
    /// ([`ImmunoError::InvalidResidue`]).
    pub fn score(&self, peptide: &str) -> Result<f64, ImmunoError> {
        let bytes = peptide.as_bytes();
        if bytes.is_empty() {
            return Err(ImmunoError::Empty { what: "peptide" });
        }
        if bytes.len() != self.weights.len() {
            return Err(ImmunoError::LengthMismatch {
                got: bytes.len(),
                expected: self.weights.len(),
            });
        }
        let mut total = 0.0;
        for (pos, &b) in bytes.iter().enumerate() {
            let idx = aa_index(b).ok_or(ImmunoError::InvalidResidue {
                residue: char::from(b),
                pos,
            })?;
            total += self.weights[pos][idx];
        }
        Ok(total)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A 2-position matrix: position 0 rewards `A` by 1.0, position 1 rewards
    /// `C` by 2.0; everything else is zero.
    fn two_pos() -> Pssm {
        let mut rows = vec![[0.0; N_AA]; 2];
        rows[0][aa_index(b'A').unwrap()] = 1.0;
        rows[1][aa_index(b'C').unwrap()] = 2.0;
        Pssm::new("test", rows).unwrap()
    }

    #[test]
    fn score_is_sum_of_position_weights() {
        let m = two_pos();
        assert!((m.score("AC").unwrap() - 3.0).abs() < 1e-12);
        assert!((m.score("AA").unwrap() - 1.0).abs() < 1e-12);
        assert!((m.score("GG").unwrap() - 0.0).abs() < 1e-12);
    }

    #[test]
    fn length_and_name() {
        let m = two_pos();
        assert_eq!(m.length(), 2);
        assert_eq!(m.name(), "test");
    }

    #[test]
    fn weight_lookup() {
        let m = two_pos();
        assert!((m.weight(0, b'A').unwrap() - 1.0).abs() < 1e-12);
        assert!((m.weight(1, b'C').unwrap() - 2.0).abs() < 1e-12);
        assert!((m.weight(0, b'C').unwrap() - 0.0).abs() < 1e-12);
        assert_eq!(m.weight(5, b'A'), None); // position out of range
        assert_eq!(m.weight(0, b'X'), None); // not a standard residue
    }

    #[test]
    fn new_rejects_empty_matrix() {
        let err = Pssm::new("x", vec![]).unwrap_err();
        assert_eq!(err.code(), "empty");
    }

    #[test]
    fn new_rejects_non_finite_weight() {
        let mut rows = vec![[0.0; N_AA]; 1];
        rows[0][3] = f64::NAN;
        let err = Pssm::new("x", rows).unwrap_err();
        assert_eq!(err.code(), "non_finite_weight");
    }

    #[test]
    fn score_rejects_length_mismatch() {
        let err = two_pos().score("ACA").unwrap_err();
        assert_eq!(err.code(), "length_mismatch");
    }

    #[test]
    fn score_rejects_invalid_residue() {
        let err = two_pos().score("AX").unwrap_err();
        assert_eq!(err.code(), "invalid_residue");
    }

    #[test]
    fn score_accepts_lowercase() {
        let m = two_pos();
        assert!((m.score("ac").unwrap() - 3.0).abs() < 1e-12);
    }

    #[test]
    fn serde_round_trips() {
        let m = two_pos();
        let json = serde_json::to_string(&m).unwrap();
        let back: Pssm = serde_json::from_str(&json).unwrap();
        assert_eq!(m, back);
    }
}
