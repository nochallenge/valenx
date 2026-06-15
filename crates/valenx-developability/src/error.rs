//! Error taxonomy for developability assessment.

use thiserror::Error;

/// Errors raised while assessing a sequence.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum DevelopabilityError {
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

    /// A pH or threshold value was not finite.
    #[error("non-finite {what}")]
    NonFinite {
        /// What was non-finite.
        what: &'static str,
    },
}

impl DevelopabilityError {
    /// A short, stable machine-readable code for this error.
    pub fn code(&self) -> &'static str {
        match self {
            DevelopabilityError::EmptySequence => "empty_sequence",
            DevelopabilityError::InvalidResidue { .. } => "invalid_residue",
            DevelopabilityError::BadWindow { .. } => "bad_window",
            DevelopabilityError::NonFinite { .. } => "non_finite",
        }
    }
}
