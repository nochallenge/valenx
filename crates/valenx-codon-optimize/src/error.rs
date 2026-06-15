//! Error taxonomy for codon optimization.

use thiserror::Error;

/// Errors raised while reverse-translating or scoring a sequence.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum CodonError {
    /// A required input (protein, DNA, usage table) was empty.
    #[error("empty {what}")]
    Empty {
        /// What was empty.
        what: &'static str,
    },

    /// A protein residue was not a codable amino acid.
    #[error("invalid amino acid {residue:?} at position {pos}")]
    InvalidResidue {
        /// The offending character.
        residue: char,
        /// Zero-based position.
        pos: usize,
    },

    /// A DNA codon was not three valid bases (T/C/A/G).
    #[error("invalid codon {codon:?} at codon index {index}")]
    InvalidCodon {
        /// The offending codon text.
        codon: String,
        /// Zero-based codon index.
        index: usize,
    },

    /// A DNA length was not a multiple of three.
    #[error("DNA length {len} is not a multiple of 3")]
    NotMultipleOfThree {
        /// The DNA length.
        len: usize,
    },

    /// A codon-usage weight was outside `(0, 1]`.
    #[error("weight {value} for codon {codon:?} is outside (0, 1]")]
    WeightOutOfRange {
        /// The offending codon.
        codon: String,
        /// The offending weight.
        value: f64,
    },

    /// No usage weight was available for a codon needed by the calculation.
    #[error("no usage weight for codon {codon:?}")]
    MissingWeight {
        /// The codon with no weight.
        codon: String,
    },
}

impl CodonError {
    /// A short, stable machine-readable code for this error.
    pub fn code(&self) -> &'static str {
        match self {
            CodonError::Empty { .. } => "empty",
            CodonError::InvalidResidue { .. } => "invalid_residue",
            CodonError::InvalidCodon { .. } => "invalid_codon",
            CodonError::NotMultipleOfThree { .. } => "not_multiple_of_three",
            CodonError::WeightOutOfRange { .. } => "weight_out_of_range",
            CodonError::MissingWeight { .. } => "missing_weight",
        }
    }
}
