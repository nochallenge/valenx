//! Error taxonomy for off-target screening.

use thiserror::Error;

/// Errors raised while screening a candidate against a reference set.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum OffTargetError {
    /// A required sequence or the reference set was empty.
    #[error("empty {what}")]
    Empty {
        /// What was empty (e.g. `"candidate"`, `"reference set"`).
        what: &'static str,
    },

    /// A sequence contained a non-standard amino-acid residue.
    #[error("invalid residue {residue:?} at position {pos} in the {which} sequence")]
    InvalidResidue {
        /// Which sequence (`"candidate"` or a reference id).
        which: String,
        /// Zero-based position of the offending residue.
        pos: usize,
        /// The offending character.
        residue: char,
    },

    /// The k-mer length `k` was zero.
    #[error("k-mer length k must be >= 1")]
    ZeroK,

    /// A similarity threshold was `NaN` or infinite.
    #[error("threshold must be finite, got {0}")]
    NonFiniteThreshold(f64),
}

impl OffTargetError {
    /// A short, stable machine-readable code for this error.
    pub fn code(&self) -> &'static str {
        match self {
            OffTargetError::Empty { .. } => "empty",
            OffTargetError::InvalidResidue { .. } => "invalid_residue",
            OffTargetError::ZeroK => "zero_k",
            OffTargetError::NonFiniteThreshold(_) => "non_finite_threshold",
        }
    }
}
