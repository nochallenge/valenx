//! Error taxonomy for PPI inference.
//!
//! Mirrors the fail-loud contract of `valenx-binder-score` and the
//! coarse, stable-code style of `valenx-align`'s `AlignError`: a
//! wrong interaction call is worse than a refusal, so every
//! unsupported input returns one of these rather than a
//! plausible-but-wrong number.

use thiserror::Error;

/// Errors raised while predicting interface contacts or scoring a PPI.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum PpiError {
    /// A supplied alignment had no rows, or fewer than the minimum
    /// number of sequences mutual information needs to be meaningful.
    #[error("paired MSA has too few sequences: {got} (need >= {need})")]
    TooFewSequences {
        /// Sequences actually supplied.
        got: usize,
        /// Minimum required.
        need: usize,
    },

    /// A supplied alignment had zero columns (empty rows).
    #[error("paired MSA has zero alignment columns")]
    EmptyAlignment,

    /// The two halves of a paired MSA disagree in depth — they must
    /// have one row per *paired* organism, in the same order.
    #[error("paired MSA halves differ in depth: chain A has {a}, chain B has {b}")]
    DepthMismatch {
        /// Depth of the chain-A half.
        a: usize,
        /// Depth of the chain-B half.
        b: usize,
    },

    /// Alignment rows were not all the same length (not a valid MSA).
    #[error("alignment rows differ in length: {got} vs expected {expected}")]
    RaggedRows {
        /// The offending row length.
        got: usize,
        /// The width established by the first row.
        expected: usize,
    },

    /// The complementarity term was requested but a structure was
    /// missing or empty — the geometric term cannot be computed without
    /// coordinates for both chains.
    #[error("complementarity requested but {what} structure/chain is missing or empty")]
    MissingStructure {
        /// Which side was missing (`"chain_a"` / `"chain_b"`).
        what: &'static str,
    },

    /// A value that must be finite was `NaN` or infinite.
    #[error("non-finite {what}")]
    NonFinite {
        /// What was non-finite.
        what: &'static str,
    },

    /// A weight was negative or non-finite.
    #[error("weight {value} for {what} must be finite and >= 0")]
    BadWeight {
        /// Which weight.
        what: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A requested rank / count parameter was out of range (e.g. `L/5`
    /// precision on an interface with no columns).
    #[error("invalid `{what}`: {reason}")]
    Invalid {
        /// Logical parameter name.
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },
}

impl PpiError {
    /// A short, stable machine-readable code for this error.
    pub fn code(&self) -> &'static str {
        match self {
            PpiError::TooFewSequences { .. } => "too_few_sequences",
            PpiError::EmptyAlignment => "empty_alignment",
            PpiError::DepthMismatch { .. } => "depth_mismatch",
            PpiError::RaggedRows { .. } => "ragged_rows",
            PpiError::MissingStructure { .. } => "missing_structure",
            PpiError::NonFinite { .. } => "non_finite",
            PpiError::BadWeight { .. } => "bad_weight",
            PpiError::Invalid { .. } => "invalid",
        }
    }

    /// Convenience constructor for [`PpiError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        PpiError::Invalid {
            what,
            reason: reason.into(),
        }
    }
}
