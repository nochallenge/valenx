//! Error type for reproducibility-bundle construction.

use thiserror::Error;

/// Something was wrong while assembling a [`crate::ReproBundle`] or verifying a
/// [`crate::provenance::Manifest`].
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
    /// Two items in the same role (input or output) of a
    /// [`crate::provenance::Manifest`] share a name.
    #[error("duplicate {role} name {name:?}")]
    DuplicateName {
        /// Which role the clash was in (`"input"` or `"output"`).
        role: &'static str,
        /// The repeated name.
        name: String,
    },
    /// A [`crate::provenance::Manifest`]'s stored `params_hash` did not match the
    /// hash recomputed from its parameters — the parameters were tampered with.
    #[error("manifest params_hash mismatch: expected {expected}, found {found}")]
    ParamsHashMismatch {
        /// The hash recomputed from the recorded parameters.
        expected: String,
        /// The (stale / forged) hash stored on the manifest.
        found: String,
    },
    /// A [`crate::provenance::Manifest`]'s stored `digest` did not match the
    /// digest recomputed from its fields — the lineage was tampered with.
    #[error("manifest digest mismatch: expected {expected}, found {found}")]
    DigestMismatch {
        /// The digest recomputed from the recorded fields.
        expected: String,
        /// The (stale / forged) digest stored on the manifest.
        found: String,
    },
}

impl ReproError {
    /// A short, stable, machine-readable code for this error.
    pub fn code(&self) -> &'static str {
        match self {
            ReproError::Empty { .. } => "empty",
            ReproError::DuplicateStepOrdinal(_) => "duplicate-step-ordinal",
            ReproError::DuplicateName { .. } => "duplicate-name",
            ReproError::ParamsHashMismatch { .. } => "params-hash-mismatch",
            ReproError::DigestMismatch { .. } => "digest-mismatch",
        }
    }
}
