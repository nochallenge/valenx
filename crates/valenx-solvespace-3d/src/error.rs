//! Typed error taxonomy for the 3D constraint solver.

use thiserror::Error;

/// Errors raised by `valenx-solvespace-3d`.
#[derive(Debug, Error)]
pub enum Solve3DError {
    /// Bad parameter (out-of-range count, NaN, etc.).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Entity reference doesn't resolve in the sketch's entity table.
    #[error("unknown entity {0:?}")]
    UnknownEntity(crate::entity::EntityId),

    /// Constraint references an entity of the wrong kind (e.g. asking
    /// for the X coordinate of a `Plane3`).
    #[error("entity kind mismatch on {id:?}: expected {expected}")]
    KindMismatch {
        /// Offending entity id.
        id: crate::entity::EntityId,
        /// Human-readable expected kind.
        expected: &'static str,
    },

    /// Jacobian factorisation gave up — usually because the system is
    /// singular (parallel constraints, degenerate workplane, etc.).
    #[error("singular system at iteration {0}")]
    Singular(usize),

    /// Solver could not converge within the configured iteration budget.
    #[error("did not converge after {0} iterations (residual {1})")]
    DidNotConverge(usize, f64),
}

/// Coarse error category, mirrors `valenx-sketch`.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User input (parameters, entity ids).
    Input,
    /// Solver / numerics.
    Algorithm,
}

impl Solve3DError {
    /// Stable kebab-cased code.
    pub fn code(&self) -> &'static str {
        match self {
            Solve3DError::BadParameter { .. } => "solvespace3d.bad_parameter",
            Solve3DError::UnknownEntity(_) => "solvespace3d.unknown_entity",
            Solve3DError::KindMismatch { .. } => "solvespace3d.kind_mismatch",
            Solve3DError::Singular(_) => "solvespace3d.singular",
            Solve3DError::DidNotConverge(_, _) => "solvespace3d.did_not_converge",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            Solve3DError::BadParameter { .. }
            | Solve3DError::UnknownEntity(_)
            | Solve3DError::KindMismatch { .. } => ErrorCategory::Input,
            Solve3DError::Singular(_) | Solve3DError::DidNotConverge(_, _) => {
                ErrorCategory::Algorithm
            }
        }
    }
}
