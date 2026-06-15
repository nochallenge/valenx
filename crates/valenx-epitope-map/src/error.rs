//! Error taxonomy for epitope mapping.

use thiserror::Error;

/// Errors raised while building a scale or mapping epitopes.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum EpitopeError {
    /// The sequence was empty.
    #[error("empty sequence")]
    EmptySequence,

    /// A residue was not one of the 20 standard amino acids.
    #[error("invalid residue {residue:?} at position {pos}")]
    InvalidResidue {
        /// The offending character.
        residue: char,
        /// Zero-based position.
        pos: usize,
    },

    /// A window length was zero or larger than the sequence.
    #[error("window {window} invalid for sequence length {len}")]
    BadWindow {
        /// Requested window length.
        window: usize,
        /// Sequence length.
        len: usize,
    },

    /// A custom scale was missing a value for a standard amino acid.
    #[error("scale missing a value for residue {residue:?}")]
    ScaleIncomplete {
        /// The residue with no scale value.
        residue: char,
    },

    /// A threshold or scale value was not finite.
    #[error("non-finite {what}")]
    NonFinite {
        /// What was non-finite.
        what: &'static str,
    },
}

impl EpitopeError {
    /// A short, stable machine-readable code for this error.
    pub fn code(&self) -> &'static str {
        match self {
            EpitopeError::EmptySequence => "empty_sequence",
            EpitopeError::InvalidResidue { .. } => "invalid_residue",
            EpitopeError::BadWindow { .. } => "bad_window",
            EpitopeError::ScaleIncomplete { .. } => "scale_incomplete",
            EpitopeError::NonFinite { .. } => "non_finite",
        }
    }
}
