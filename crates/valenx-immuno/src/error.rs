//! Error taxonomy for immunogenicity screening.
//!
//! Every fallible entry point returns [`ImmunoError`]. Inputs are validated up
//! front — matrices reject non-finite weights and empty position sets, peptides
//! must match the matrix length and contain only standard residues, and a
//! protein must be at least one window long.

use thiserror::Error;

/// Errors raised while building a matrix or screening a sequence.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum ImmunoError {
    /// A required collection (matrix rows, peptide, protein) was empty.
    #[error("empty {what}")]
    Empty {
        /// What was empty (e.g. `"matrix"`, `"peptide"`, `"protein"`).
        what: &'static str,
    },

    /// A peptide's length did not match the matrix's window length.
    #[error("peptide length {got} does not match matrix length {expected}")]
    LengthMismatch {
        /// The peptide length supplied.
        got: usize,
        /// The matrix (window) length required.
        expected: usize,
    },

    /// A residue was not one of the 20 standard amino acids.
    #[error(
        "invalid residue {residue:?} at position {pos} (expected one of ACDEFGHIKLMNPQRSTVWY)"
    )]
    InvalidResidue {
        /// The offending character.
        residue: char,
        /// Zero-based position within the peptide or protein.
        pos: usize,
    },

    /// A matrix weight was `NaN` or infinite.
    #[error("non-finite weight at position {pos}, residue index {aa}")]
    NonFiniteWeight {
        /// Zero-based matrix position (row).
        pos: usize,
        /// Zero-based residue column (`0..20`).
        aa: usize,
    },

    /// A score threshold was `NaN` or infinite.
    #[error("threshold must be finite, got {0}")]
    NonFiniteThreshold(f64),

    /// The protein was shorter than one matrix window.
    #[error("protein length {protein} is shorter than the window length {window}")]
    ProteinTooShort {
        /// The protein length supplied.
        protein: usize,
        /// The window (matrix) length.
        window: usize,
    },
}

impl ImmunoError {
    /// A short, stable machine-readable code for this error (handy for logs,
    /// tests, and UI without matching on the human message).
    pub fn code(&self) -> &'static str {
        match self {
            ImmunoError::Empty { .. } => "empty",
            ImmunoError::LengthMismatch { .. } => "length_mismatch",
            ImmunoError::InvalidResidue { .. } => "invalid_residue",
            ImmunoError::NonFiniteWeight { .. } => "non_finite_weight",
            ImmunoError::NonFiniteThreshold(_) => "non_finite_threshold",
            ImmunoError::ProteinTooShort { .. } => "protein_too_short",
        }
    }
}
