//! Typed error taxonomy for the CGAL algorithms port.

use thiserror::Error;

/// Errors raised by `valenx-cgal-port`.
#[derive(Debug, Error)]
pub enum CgalError {
    /// Bad parameter.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Reason.
        reason: String,
    },

    /// Need at least N points for this algorithm.
    #[error("need >= {needed} points, got {given}")]
    NotEnoughPoints {
        /// Required count.
        needed: usize,
        /// Caller's count.
        given: usize,
    },

    /// Collinear input — e.g. all points on a single line, which
    /// blows up Bowyer-Watson.
    #[error("input is collinear / coplanar; cannot triangulate")]
    Degenerate,

    /// Algorithm produced a non-finite or otherwise invalid result.
    #[error("numerical failure in {algo}: {reason}")]
    Numerical {
        /// Algorithm label.
        algo: &'static str,
        /// Reason.
        reason: String,
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

impl CgalError {
    /// Stable kebab code.
    pub fn code(&self) -> &'static str {
        match self {
            CgalError::BadParameter { .. } => "cgal.bad_parameter",
            CgalError::NotEnoughPoints { .. } => "cgal.not_enough_points",
            CgalError::Degenerate => "cgal.degenerate",
            CgalError::Numerical { .. } => "cgal.numerical",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            CgalError::BadParameter { .. } | CgalError::NotEnoughPoints { .. } => {
                ErrorCategory::Input
            }
            CgalError::Degenerate | CgalError::Numerical { .. } => ErrorCategory::Algorithm,
        }
    }
}
