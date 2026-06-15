//! Error taxonomy for the selection funnel.

use thiserror::Error;

/// Errors raised during consensus ranking, diversity selection, or the funnel.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum SelectError {
    /// A required collection was empty.
    #[error("empty {what}")]
    Empty {
        /// What was empty (e.g. `"methods"`, `"candidates"`, `"features"`).
        what: &'static str,
    },

    /// Paired or parallel collections had inconsistent lengths/dimensions.
    #[error("inconsistent {what}")]
    Inconsistent {
        /// What was inconsistent (e.g. `"method length"`, `"feature dimension"`).
        what: &'static str,
    },

    /// A value that must be finite was `NaN` or infinite.
    #[error("non-finite {what}")]
    NonFinite {
        /// What was non-finite.
        what: &'static str,
    },

    /// The requested selection size `n` was zero.
    #[error("n must be >= 1")]
    ZeroN,

    /// A diversity radius was not strictly positive.
    #[error("radius {value} must be > 0")]
    NonPositiveRadius {
        /// The offending radius.
        value: f64,
    },

    /// A start index was out of range.
    #[error("start index {index} out of range for {len} items")]
    StartOutOfRange {
        /// The offending index.
        index: usize,
        /// The number of items.
        len: usize,
    },
}

impl SelectError {
    /// A short, stable machine-readable code for this error.
    pub fn code(&self) -> &'static str {
        match self {
            SelectError::Empty { .. } => "empty",
            SelectError::Inconsistent { .. } => "inconsistent",
            SelectError::NonFinite { .. } => "non_finite",
            SelectError::ZeroN => "zero_n",
            SelectError::NonPositiveRadius { .. } => "non_positive_radius",
            SelectError::StartOutOfRange { .. } => "start_out_of_range",
        }
    }
}
