//! Error type for the CFD solver.
//!
//! The solver itself is total — [`crate::solve_simple`] always returns
//! a [`crate::FlowSolution`], with a `converged` flag rather than an
//! `Err` for an iteration-cap miss, because a non-converged field is
//! still useful (it can be inspected, or the run resumed with more
//! iterations). This type covers the genuine *input* failure modes a
//! caller building a job from external data must handle.

use thiserror::Error;

/// Errors raised when setting up a CFD case from external input.
#[derive(Debug, Error)]
pub enum CfdError {
    /// A grid or fluid parameter was out of range — a zero cell count,
    /// a non-positive domain length, a non-physical density / viscosity.
    #[error("invalid CFD parameter `{name}`: {reason}")]
    BadParameter {
        /// Name of the offending parameter.
        name: &'static str,
        /// Human-readable explanation.
        reason: String,
    },

    /// The boundary specification cannot drive a well-posed flow — for
    /// example, an inlet on every side with no outlet (mass cannot
    /// leave), or an all-outlet domain (mass cannot enter).
    #[error("ill-posed boundary conditions: {0}")]
    IllPosedBoundary(String),
}

impl CfdError {
    /// A stable, kebab-cased identifier for the error.
    pub fn code(&self) -> &'static str {
        match self {
            CfdError::BadParameter { .. } => "cfd.bad_parameter",
            CfdError::IllPosedBoundary(_) => "cfd.ill_posed_boundary",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_are_stable() {
        assert_eq!(
            CfdError::BadParameter {
                name: "nx",
                reason: "zero".into()
            }
            .code(),
            "cfd.bad_parameter"
        );
        assert_eq!(
            CfdError::IllPosedBoundary("no outlet".into()).code(),
            "cfd.ill_posed_boundary"
        );
    }

    #[test]
    fn display_carries_context() {
        let e = CfdError::BadParameter {
            name: "viscosity",
            reason: "must be positive".into(),
        };
        assert!(e.to_string().contains("viscosity"));
        assert!(e.to_string().contains("must be positive"));
    }
}
