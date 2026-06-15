//! Inductor / coil error taxonomy.
//!
//! Every fallible constructor in this crate funnels its rejection
//! through [`CoilError`]. The type carries a stable
//! [`code`](CoilError::code) and a coarse
//! [`category`](CoilError::category) so callers can branch on the kind
//! of failure (bad user input vs. a degenerate geometry) without string
//! matching on the human-readable message.

use thiserror::Error;

/// Errors raised when validating coil parameters or evaluating coil
/// models.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum CoilError {
    /// A scalar parameter that must be strictly positive was zero or
    /// negative (e.g. the length `l`, the cross-sectional area `A`, the
    /// frequency `f`, or the resistance `R`).
    #[error("parameter `{name}` must be > 0, got {value}")]
    NonPositive {
        /// The name of the offending parameter, as used in the formulas.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A parameter that must be finite (not NaN, not infinite) was not.
    #[error("parameter `{name}` must be finite, got {value}")]
    NotFinite {
        /// The name of the offending parameter.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// The number of turns `N` was negative. Zero turns is permitted
    /// (it yields zero inductance), but a negative turn count has no
    /// physical meaning.
    #[error("turn count `N` must be >= 0, got {value}")]
    NegativeTurns {
        /// The rejected turn count.
        value: f64,
    },

    /// The relative permeability `mu_r` was below 1. For the
    /// non-ferromagnetic (and idealised ferromagnetic) regime this model
    /// targets, `mu_r >= 1`; a value below 1 (strong diamagnetism) is
    /// outside scope.
    #[error("relative permeability `mu_r` must be >= 1, got {value}")]
    BadPermeability {
        /// The rejected relative permeability.
        value: f64,
    },
}

/// Coarse classification of a [`CoilError`], for telemetry and
/// caller-side branching.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// The caller supplied an invalid value (out of range, non-finite).
    Input,
    /// The requested configuration falls outside the model's validity
    /// domain.
    Domain,
}

impl CoilError {
    /// A stable, kebab-cased identifier for this error variant.
    ///
    /// Unlike the [`Display`](std::fmt::Display) message, the code is
    /// part of the crate's contract and is safe to match on in tests and
    /// telemetry pipelines.
    pub fn code(&self) -> &'static str {
        match self {
            CoilError::NonPositive { .. } => "coil.non-positive",
            CoilError::NotFinite { .. } => "coil.not-finite",
            CoilError::NegativeTurns { .. } => "coil.negative-turns",
            CoilError::BadPermeability { .. } => "coil.bad-permeability",
        }
    }

    /// The coarse [`ErrorCategory`] this error belongs to.
    pub fn category(&self) -> ErrorCategory {
        match self {
            CoilError::NonPositive { .. }
            | CoilError::NotFinite { .. }
            | CoilError::NegativeTurns { .. } => ErrorCategory::Input,
            CoilError::BadPermeability { .. } => ErrorCategory::Domain,
        }
    }
}

/// Reject a value that is not strictly positive (and not finite).
///
/// Used by the validated constructors so that the same rule produces the
/// same [`CoilError`] everywhere.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64, CoilError> {
    if !value.is_finite() {
        return Err(CoilError::NotFinite { name, value });
    }
    if value <= 0.0 {
        return Err(CoilError::NonPositive { name, value });
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn require_positive_accepts_positive() {
        assert_eq!(require_positive("x", 2.5).unwrap(), 2.5);
    }

    #[test]
    fn require_positive_rejects_zero() {
        let err = require_positive("l", 0.0).unwrap_err();
        assert_eq!(err.code(), "coil.non-positive");
        assert_eq!(err.category(), ErrorCategory::Input);
    }

    #[test]
    fn require_positive_rejects_negative() {
        let err = require_positive("A", -1.0).unwrap_err();
        assert!(matches!(err, CoilError::NonPositive { .. }));
    }

    #[test]
    fn require_positive_rejects_nan() {
        let err = require_positive("f", f64::NAN).unwrap_err();
        assert_eq!(err.code(), "coil.not-finite");
        assert_eq!(err.category(), ErrorCategory::Input);
    }

    #[test]
    fn require_positive_rejects_infinity() {
        let err = require_positive("R", f64::INFINITY).unwrap_err();
        assert!(matches!(err, CoilError::NotFinite { .. }));
    }

    #[test]
    fn categories_partition_variants() {
        assert_eq!(
            CoilError::NegativeTurns { value: -1.0 }.category(),
            ErrorCategory::Input
        );
        assert_eq!(
            CoilError::BadPermeability { value: 0.5 }.category(),
            ErrorCategory::Domain
        );
        assert_eq!(
            CoilError::BadPermeability { value: 0.5 }.code(),
            "coil.bad-permeability"
        );
    }
}
