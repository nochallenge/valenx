//! Error type for reproducibility-bundle construction.

use thiserror::Error;

/// Something was wrong while assembling a [`crate::ReproBundle`].
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum ReproError {
    /// A required text field was empty (or only whitespace).
    #[error("{field} must not be empty")]
    Empty {
        /// Which field was empty.
        field: &'static str,
    },
    /// Two provenance steps share the same ordinal.
    #[error("duplicate provenance step ordinal {0}")]
    DuplicateStepOrdinal(u32),
}
