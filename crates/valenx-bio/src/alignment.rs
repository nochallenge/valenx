//! Multiple-sequence alignment (MSA) canonical type.
//!
//! Holds a list of [`Sequence`]s of identical length over a single
//! [`Alphabet`]. Gap characters are conventionally `-` (dash) — the
//! [`Alphabet`] matchset accepts `-` for all three alphabets so an
//! aligned `Sequence` round-trips through [`Sequence::new`] unchanged.

use crate::sequence::{Alphabet, Sequence};
use serde::{Deserialize, Serialize};
use thiserror::Error;

/// A multiple-sequence alignment.
///
/// Invariants enforced by [`Alignment::new`]:
/// - At least one row.
/// - All rows share the same [`Alphabet`].
/// - All rows share the same length.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Alignment {
    rows: Vec<Sequence>,
}

/// Construction errors from [`Alignment::new`].
#[derive(Debug, Error, PartialEq, Eq)]
pub enum AlignmentError {
    /// The caller passed an empty row list — an alignment must hold at
    /// least one row.
    #[error("alignment must have at least one row")]
    Empty,
    /// One row's length disagrees with the first row's length. All rows
    /// must be padded to the same column count.
    #[error("row `{name}` has length {got}, expected {expected}")]
    UnequalRows {
        /// Identifier of the offending row.
        name: String,
        /// Length of the offending row.
        got: usize,
        /// Required length (the first row's length).
        expected: usize,
    },
    /// One row's alphabet disagrees with the first row's alphabet. All
    /// rows must share the same residue alphabet.
    #[error("row `{name}` has alphabet {got:?}, expected {expected:?}")]
    AlphabetMismatch {
        /// Identifier of the offending row.
        name: String,
        /// Alphabet of the offending row.
        got: Alphabet,
        /// Required alphabet (the first row's alphabet).
        expected: Alphabet,
    },
}

impl Alignment {
    /// Build an [`Alignment`] from `rows`, enforcing the type's invariants
    /// (non-empty, equal lengths, shared alphabet).
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError::Empty`], [`AlignmentError::UnequalRows`]
    /// or [`AlignmentError::AlphabetMismatch`] when the corresponding
    /// invariant is violated.
    pub fn new(rows: Vec<Sequence>) -> Result<Self, AlignmentError> {
        let first = rows.first().ok_or(AlignmentError::Empty)?;
        let expected_len = first.len();
        let expected_alphabet = first.alphabet;
        for row in &rows[1..] {
            if row.alphabet != expected_alphabet {
                return Err(AlignmentError::AlphabetMismatch {
                    name: row.name.clone(),
                    got: row.alphabet,
                    expected: expected_alphabet,
                });
            }
            if row.len() != expected_len {
                return Err(AlignmentError::UnequalRows {
                    name: row.name.clone(),
                    got: row.len(),
                    expected: expected_len,
                });
            }
        }
        Ok(Self { rows })
    }

    /// Borrow the rows in declaration order.
    pub fn rows(&self) -> &[Sequence] {
        &self.rows
    }

    /// Number of rows (sequences) in the alignment.
    pub fn row_count(&self) -> usize {
        self.rows.len()
    }

    /// Number of columns (residue + gap positions). Equal to the length
    /// of any row by the alignment's invariant.
    pub fn column_count(&self) -> usize {
        self.rows.first().map(|r| r.len()).unwrap_or(0)
    }

    /// Shared alphabet of every row.
    pub fn alphabet(&self) -> Alphabet {
        // First-row alphabet — invariant guarantees all rows match.
        self.rows[0].alphabet
    }

    /// Total `'-'` count across all rows.
    pub fn gap_count(&self) -> usize {
        self.rows
            .iter()
            .map(|r| r.as_str().bytes().filter(|&b| b == b'-').count())
            .sum()
    }
}
