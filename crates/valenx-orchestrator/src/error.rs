//! Error taxonomy for the design funnel.

use thiserror::Error;

/// Errors raised while running [`crate::run_funnel`].
///
/// Each stage's underlying error is wrapped (`#[from]`) so the failing stage is
/// always identifiable; [`OrchestratorError::code`] gives a short, stable label.
#[derive(Debug, Error)]
pub enum OrchestratorError {
    /// A required input was empty (e.g. no candidates were supplied).
    #[error("empty {what}")]
    Empty {
        /// What was empty.
        what: &'static str,
    },

    /// An invariant the orchestrator relies on was violated — for example a
    /// shortlist entry referenced an id absent from the input candidates. This
    /// signals a bug, not bad input.
    #[error("internal invariant violated: {what}")]
    Internal {
        /// The invariant that failed.
        what: &'static str,
    },

    /// The selection stage ([`valenx_select`]) failed.
    #[error("selection stage failed: {0}")]
    Select(#[from] valenx_select::SelectError),

    /// The safety-consolidation stage ([`valenx_safety`]) failed.
    #[error("safety stage failed: {0}")]
    Safety(#[from] valenx_safety::SafetyError),

    /// The dossier-assembly stage ([`valenx_dossier`]) failed.
    #[error("dossier stage failed: {0}")]
    Dossier(#[from] valenx_dossier::DossierError),
}

impl OrchestratorError {
    /// A short, stable machine-readable code for this error.
    pub fn code(&self) -> &'static str {
        match self {
            OrchestratorError::Empty { .. } => "empty",
            OrchestratorError::Internal { .. } => "internal",
            OrchestratorError::Select(_) => "select",
            OrchestratorError::Safety(_) => "safety",
            OrchestratorError::Dossier(_) => "dossier",
        }
    }
}
