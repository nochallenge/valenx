//! Error taxonomy for the three-phase power module.
//!
//! A single [`ThreePhaseError`] enum covers every failure mode of the
//! validated constructors and computations. Each variant carries the
//! offending parameter name and a human-readable reason so callers can
//! surface actionable diagnostics.

use thiserror::Error;

/// Errors raised when validating three-phase inputs or computing
/// derived quantities.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum ThreePhaseError {
    /// A magnitude that must be strictly positive (a voltage or
    /// current) was zero or negative.
    #[error("`{name}` must be > 0, got {value}")]
    NonPositive {
        /// Name of the offending parameter.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A magnitude that must be finite (not NaN or infinite) was not.
    #[error("`{name}` must be finite, got {value}")]
    NotFinite {
        /// Name of the offending parameter.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A power factor `cos(phi)` outside the closed interval
    /// `[-1.0, 1.0]` was supplied.
    #[error("power factor must lie in [-1, 1], got {value}")]
    PowerFactorOutOfRange {
        /// The rejected power-factor value.
        value: f64,
    },
}

/// Coarse classification of a [`ThreePhaseError`], useful for grouping
/// diagnostics or deciding retry behaviour at a call site.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The supplied value violated a domain constraint (sign,
    /// finiteness, or range).
    Domain,
}

impl ThreePhaseError {
    /// Stable, kebab-cased identifier for this error, suitable for logs
    /// or machine matching.
    pub fn code(&self) -> &'static str {
        match self {
            ThreePhaseError::NonPositive { .. } => "threephase.non-positive",
            ThreePhaseError::NotFinite { .. } => "threephase.not-finite",
            ThreePhaseError::PowerFactorOutOfRange { .. } => "threephase.power-factor-out-of-range",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            ThreePhaseError::NonPositive { .. }
            | ThreePhaseError::NotFinite { .. }
            | ThreePhaseError::PowerFactorOutOfRange { .. } => ErrorCategory::Domain,
        }
    }
}

/// Validate that `value` is finite and strictly positive, returning the
/// value on success.
///
/// Used by the validated constructors for voltage and current
/// magnitudes, which are always strictly positive in a physical
/// balanced system.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64, ThreePhaseError> {
    if !value.is_finite() {
        return Err(ThreePhaseError::NotFinite { name, value });
    }
    if value <= 0.0 {
        return Err(ThreePhaseError::NonPositive { name, value });
    }
    Ok(value)
}

/// Validate that a power factor `cos(phi)` is finite and within the
/// closed interval `[-1.0, 1.0]`, returning it on success.
pub(crate) fn require_power_factor(value: f64) -> Result<f64, ThreePhaseError> {
    if !value.is_finite() {
        return Err(ThreePhaseError::NotFinite {
            name: "power_factor",
            value,
        });
    }
    if !(-1.0..=1.0).contains(&value) {
        return Err(ThreePhaseError::PowerFactorOutOfRange { value });
    }
    Ok(value)
}
