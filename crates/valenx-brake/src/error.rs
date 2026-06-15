//! Brake-model error taxonomy.
//!
//! A single [`BrakeError`] enum covers every way an input can be
//! rejected. The `check_*` constructors centralise the validation so
//! the [`crate::disc`], [`crate::band`] and [`crate::energy`] modules
//! share one consistent set of guards (finiteness, sign, and the
//! physical domain of the friction coefficient).

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Errors raised while evaluating a brake model.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum BrakeError {
    /// A parameter was not a finite number (`NaN` or `±∞`).
    #[error("parameter `{name}` must be finite, got {value}")]
    NotFinite {
        /// Offending parameter name.
        name: &'static str,
        /// The non-finite value supplied.
        value: f64,
    },

    /// A parameter that must be strictly positive was zero or negative.
    #[error("parameter `{name}` must be > 0, got {value}")]
    NonPositive {
        /// Offending parameter name.
        name: &'static str,
        /// The non-positive value supplied.
        value: f64,
    },

    /// A parameter that must be non-negative was negative.
    #[error("parameter `{name}` must be >= 0, got {value}")]
    Negative {
        /// Offending parameter name.
        name: &'static str,
        /// The negative value supplied.
        value: f64,
    },

    /// The friction coefficient `mu` was outside the modelled range
    /// `0 < mu <= mu_max`.
    #[error("friction coefficient `mu` must be in (0, {max}], got {value}")]
    FrictionOutOfRange {
        /// The `mu` value supplied.
        value: f64,
        /// The inclusive upper bound this crate accepts.
        max: f64,
    },

    /// A count (e.g. the number of pad faces) was zero.
    #[error("count `{name}` must be >= 1, got 0")]
    ZeroCount {
        /// Offending count name.
        name: &'static str,
    },
}

/// Coarse category for a [`BrakeError`], for callers that group errors
/// without matching every variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum ErrorCategory {
    /// A user-supplied input was invalid.
    Input,
    /// A value fell outside the model's valid physical domain.
    Domain,
}

/// The widest friction coefficient this crate will accept.
///
/// Coulomb friction coefficients for engineering brake materials sit
/// well below 1; a hard upper bound of `2.0` rejects nonsense input
/// (negative, zero, or absurdly large `mu`) while still admitting every
/// realistic — and a generous margin of exotic — value.
pub const MU_MAX: f64 = 2.0;

impl BrakeError {
    /// Stable kebab-cased identifier for logging / matching.
    pub fn code(&self) -> &'static str {
        match self {
            BrakeError::NotFinite { .. } => "brake.not_finite",
            BrakeError::NonPositive { .. } => "brake.non_positive",
            BrakeError::Negative { .. } => "brake.negative",
            BrakeError::FrictionOutOfRange { .. } => "brake.friction_out_of_range",
            BrakeError::ZeroCount { .. } => "brake.zero_count",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            BrakeError::FrictionOutOfRange { .. } => ErrorCategory::Domain,
            _ => ErrorCategory::Input,
        }
    }
}

/// Reject a non-finite value, returning it unchanged when it is finite.
///
/// # Errors
/// [`BrakeError::NotFinite`] if `value` is `NaN` or infinite.
pub fn check_finite(name: &'static str, value: f64) -> Result<f64, BrakeError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(BrakeError::NotFinite { name, value })
    }
}

/// Require a finite, strictly-positive value.
///
/// # Errors
/// [`BrakeError::NotFinite`] if `value` is not finite, or
/// [`BrakeError::NonPositive`] if `value <= 0`.
pub fn check_positive(name: &'static str, value: f64) -> Result<f64, BrakeError> {
    let value = check_finite(name, value)?;
    if value > 0.0 {
        Ok(value)
    } else {
        Err(BrakeError::NonPositive { name, value })
    }
}

/// Require a finite, non-negative value.
///
/// # Errors
/// [`BrakeError::NotFinite`] if `value` is not finite, or
/// [`BrakeError::Negative`] if `value < 0`.
pub fn check_non_negative(name: &'static str, value: f64) -> Result<f64, BrakeError> {
    let value = check_finite(name, value)?;
    if value >= 0.0 {
        Ok(value)
    } else {
        Err(BrakeError::Negative { name, value })
    }
}

/// Require the friction coefficient `mu` to lie in `(0, MU_MAX]`.
///
/// # Errors
/// [`BrakeError::NotFinite`] if `mu` is not finite, or
/// [`BrakeError::FrictionOutOfRange`] if `mu <= 0` or `mu > MU_MAX`.
pub fn check_friction(mu: f64) -> Result<f64, BrakeError> {
    let mu = check_finite("mu", mu)?;
    if mu > 0.0 && mu <= MU_MAX {
        Ok(mu)
    } else {
        Err(BrakeError::FrictionOutOfRange {
            value: mu,
            max: MU_MAX,
        })
    }
}

/// Require a count to be at least one.
///
/// # Errors
/// [`BrakeError::ZeroCount`] if `count == 0`.
pub fn check_count(name: &'static str, count: u32) -> Result<u32, BrakeError> {
    if count >= 1 {
        Ok(count)
    } else {
        Err(BrakeError::ZeroCount { name })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finite_accepts_and_rejects() {
        assert!((check_finite("x", 1.5).unwrap() - 1.5).abs() < 1e-12);
        // NaN != NaN, so we cannot `assert_eq!` the whole error; compare
        // the stable code/category instead.
        let nan_err = check_finite("x", f64::NAN).unwrap_err();
        assert_eq!(nan_err.code(), "brake.not_finite");
        assert_eq!(nan_err.category(), ErrorCategory::Input);

        let inf_err = check_finite("x", f64::INFINITY).unwrap_err();
        assert_eq!(inf_err.code(), "brake.not_finite");
        assert_eq!(
            check_finite("x", f64::NEG_INFINITY).unwrap_err().code(),
            "brake.not_finite"
        );
    }

    #[test]
    fn positive_guard() {
        assert!((check_positive("x", 2.0).unwrap() - 2.0).abs() < 1e-12);
        assert_eq!(
            check_positive("x", 0.0).unwrap_err(),
            BrakeError::NonPositive {
                name: "x",
                value: 0.0
            }
        );
        assert_eq!(
            check_positive("x", -3.0).unwrap_err(),
            BrakeError::NonPositive {
                name: "x",
                value: -3.0
            }
        );
    }

    #[test]
    fn non_negative_guard() {
        assert!((check_non_negative("x", 0.0).unwrap()).abs() < 1e-12);
        assert!((check_non_negative("x", 4.0).unwrap() - 4.0).abs() < 1e-12);
        assert_eq!(
            check_non_negative("x", -1e-9).unwrap_err().code(),
            "brake.negative"
        );
    }

    #[test]
    fn friction_domain() {
        assert!((check_friction(0.4).unwrap() - 0.4).abs() < 1e-12);
        assert!((check_friction(MU_MAX).unwrap() - MU_MAX).abs() < 1e-12);
        assert_eq!(
            check_friction(0.0).unwrap_err(),
            BrakeError::FrictionOutOfRange {
                value: 0.0,
                max: MU_MAX
            }
        );
        assert_eq!(
            check_friction(MU_MAX + 0.1).unwrap_err().code(),
            "brake.friction_out_of_range"
        );
        // The friction-domain error is categorised as a Domain error.
        assert_eq!(
            check_friction(-0.1).unwrap_err().category(),
            ErrorCategory::Domain
        );
    }

    #[test]
    fn count_guard() {
        assert_eq!(check_count("n", 1).unwrap(), 1);
        assert_eq!(check_count("n", 4).unwrap(), 4);
        assert_eq!(
            check_count("n", 0).unwrap_err(),
            BrakeError::ZeroCount { name: "n" }
        );
    }

    #[test]
    fn categories_and_codes() {
        assert_eq!(
            BrakeError::NotFinite {
                name: "x",
                value: f64::NAN
            }
            .category(),
            ErrorCategory::Input
        );
        assert_eq!(
            BrakeError::FrictionOutOfRange {
                value: 9.0,
                max: MU_MAX
            }
            .category(),
            ErrorCategory::Domain
        );
        assert_eq!(
            BrakeError::ZeroCount { name: "n" }.code(),
            "brake.zero_count"
        );
    }
}
