//! Error taxonomy for `valenx-neuro`.
//!
//! Every fallible public function returns [`Result<_, NeuroError>`]. The
//! hand-rolled `Display` / `Error` impls mirror the rest of the workspace
//! (no `thiserror` dependency).

use std::fmt;

/// Errors produced by `valenx-neuro`.
#[derive(Debug, Clone, PartialEq)]
pub enum NeuroError {
    /// The simulation input is invalid: an empty scene, a non-positive
    /// timestep, a non-positive mesh resolution, or a degenerate geometry.
    Invalid {
        /// Human-readable reason, surfaced in the UI.
        reason: String,
    },
    /// The underlying FEM field / bioheat solver failed. The string is the
    /// underlying solver error, verbatim.
    Solver(String),
}

impl fmt::Display for NeuroError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NeuroError::Invalid { reason } => {
                write!(f, "invalid simulation input: {reason}")
            }
            NeuroError::Solver(e) => write!(f, "field/bioheat solver failed: {e}"),
        }
    }
}

impl std::error::Error for NeuroError {}

/// Convenience alias for `Result<T, NeuroError>`.
pub type Result<T> = std::result::Result<T, NeuroError>;
