//! Typed errors for the libigl algorithms port.

use thiserror::Error;

/// Errors raised by `valenx-libigl-port`.
#[derive(Debug, Error)]
pub enum LibiglError {
    /// Bad parameter.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Need at least N input vertices / triangles.
    #[error("need >= {needed}, got {given}")]
    NotEnough {
        /// What's missing (`"vertices"`, `"triangles"`, …).
        what: &'static str,
        /// Required count.
        needed: usize,
        /// Caller's count.
        given: usize,
    },

    /// Algorithm did not converge inside the configured iteration limit.
    #[error("did not converge after {iters} iters (residual {residual})")]
    DidNotConverge {
        /// Iterations performed.
        iters: usize,
        /// Final residual.
        residual: f64,
    },

    /// Linear-system solve failed (singular matrix).
    #[error("singular system in {algo}")]
    Singular {
        /// Algorithm label.
        algo: &'static str,
    },
}

/// Coarse category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input.
    Input,
    /// Numerics / algorithm.
    Algorithm,
}

impl LibiglError {
    /// Stable kebab code.
    pub fn code(&self) -> &'static str {
        match self {
            LibiglError::BadParameter { .. } => "libigl.bad_parameter",
            LibiglError::NotEnough { .. } => "libigl.not_enough",
            LibiglError::DidNotConverge { .. } => "libigl.did_not_converge",
            LibiglError::Singular { .. } => "libigl.singular",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            LibiglError::BadParameter { .. } | LibiglError::NotEnough { .. } => {
                ErrorCategory::Input
            }
            LibiglError::DidNotConverge { .. } | LibiglError::Singular { .. } => {
                ErrorCategory::Algorithm
            }
        }
    }
}
