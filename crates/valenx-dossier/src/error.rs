//! Error taxonomy for dossier assembly.

use thiserror::Error;

/// Errors raised while building a [`crate::RunDossier`].
#[derive(Debug, Error)]
pub enum DossierError {
    /// A required text field (goal, id, flag) was empty.
    #[error("empty {what}")]
    Empty {
        /// What was empty.
        what: &'static str,
    },

    /// A score or component value was `NaN` or infinite.
    #[error("non-finite {what}")]
    NonFinite {
        /// What was non-finite.
        what: &'static str,
    },

    /// A calibrated confidence fell outside `[0, 1]`.
    #[error("confidence {value} is outside [0, 1]")]
    ConfidenceOutOfRange {
        /// The offending confidence.
        value: f64,
    },

    /// The underlying reproducibility bundle could not be built.
    #[error("reproducibility bundle error: {0}")]
    Repro(#[from] valenx_repro::ReproError),
}

impl DossierError {
    /// A short, stable machine-readable code for this error.
    pub fn code(&self) -> &'static str {
        match self {
            DossierError::Empty { .. } => "empty",
            DossierError::NonFinite { .. } => "non_finite",
            DossierError::ConfidenceOutOfRange { .. } => "confidence_out_of_range",
            DossierError::Repro(_) => "repro",
        }
    }
}
