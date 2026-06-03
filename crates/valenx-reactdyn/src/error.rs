//! Error taxonomy for `valenx-reactdyn`.
//!
//! Every fallible public function returns [`Result<_, ReactDynError>`].
//! The variants are deliberately coarse — a caller usually only needs to
//! know whether the electronic-structure step failed, whether the input
//! was nonsensical, or whether the run was refused for being too large.
//! The hand-rolled `Display` / `Error` impls mirror `valenx-qchem`'s
//! `QchemError` (no `thiserror` dependency).

use std::fmt;

/// Errors produced by `valenx-reactdyn`.
#[derive(Debug, Clone, PartialEq)]
pub enum ReactDynError {
    /// The electronic-structure (qchem) step failed — SCF did not
    /// converge, a basis was missing for an element, or the geometry was
    /// rejected. The string is the underlying qchem error, verbatim.
    Qchem(String),
    /// The simulation input is invalid: an empty system, a non-positive
    /// timestep, zero steps, or a non-positive finite-difference delta.
    Invalid {
        /// Human-readable reason, surfaced in the UI.
        reason: String,
    },
    /// The requested run exceeds the safety guard on `atoms × steps`,
    /// which bounds the (expensive) numerical-gradient compute cost.
    GuardExceeded {
        /// Number of atoms in the system.
        atoms: usize,
        /// Number of requested steps.
        steps: usize,
        /// The configured cap on `atoms × steps`.
        cap: usize,
    },
}

impl fmt::Display for ReactDynError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ReactDynError::Qchem(e) => {
                write!(f, "electronic-structure step failed: {e}")
            }
            ReactDynError::Invalid { reason } => {
                write!(f, "invalid simulation input: {reason}")
            }
            ReactDynError::GuardExceeded { atoms, steps, cap } => write!(
                f,
                "run too large: {atoms} atoms × {steps} steps exceeds the cost guard \
                 ({cap}); reduce the system size or the number of steps"
            ),
        }
    }
}

impl std::error::Error for ReactDynError {}

/// Convenience alias for `Result<T, ReactDynError>`.
pub type Result<T> = std::result::Result<T, ReactDynError>;
