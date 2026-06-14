//! Error taxonomy for `valenx-endurance`.
//!
//! Every fallible constructor and function in this crate returns
//! [`Result<_, EnduranceError>`]. The variants are deliberately coarse:
//! an exercise-physiology caller usually only cares whether
//!
//! 1. a supplied quantity was out of its physical domain — a negative
//!    partial pressure, a non-positive `p50`, a hemoglobin concentration
//!    below zero, an exercise intensity outside `0..=1`
//!    ([`EnduranceError::OutOfDomain`]); or
//! 2. a value that must be finite was `NaN` or infinite
//!    ([`EnduranceError::NotFinite`]).
//!
//! Use [`EnduranceError::code`] for stable log / telemetry tagging and
//! [`EnduranceError::category`] to bucket failures without matching every
//! variant. The pattern mirrors the other Valenx domain crates
//! (`valenx-springs`, `valenx-md`).

use thiserror::Error;

/// Errors produced by `valenx-endurance`.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum EnduranceError {
    /// A quantity fell outside its physically valid domain — for
    /// example a negative partial pressure, a non-positive `p50` or Hill
    /// coefficient, a negative hemoglobin concentration, or an exercise
    /// intensity outside the closed interval `0..=1`.
    #[error("`{name}` out of domain: {value} ({reason})")]
    OutOfDomain {
        /// Logical parameter name (e.g. `"po2"`, `"p50"`, `"hb_g_dl"`).
        name: &'static str,
        /// The offending value, formatted for the message.
        value: f64,
        /// Human-readable explanation of the valid domain.
        reason: &'static str,
    },

    /// A value that must be finite was `NaN` or infinite.
    #[error("`{name}` must be finite, got {value}")]
    NotFinite {
        /// Logical parameter name.
        name: &'static str,
        /// The offending value.
        value: f64,
    },
}

/// Coarse error category, for bucketing without matching every variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The caller supplied an invalid input value.
    Input,
}

impl EnduranceError {
    /// Stable kebab-cased identifier for log / telemetry tagging.
    ///
    /// The string is part of the crate's observable contract and will
    /// not change for a given variant across patch releases.
    #[must_use]
    pub fn code(&self) -> &'static str {
        match self {
            EnduranceError::OutOfDomain { .. } => "endurance.out_of_domain",
            EnduranceError::NotFinite { .. } => "endurance.not_finite",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    #[must_use]
    pub fn category(&self) -> ErrorCategory {
        match self {
            EnduranceError::OutOfDomain { .. } | EnduranceError::NotFinite { .. } => {
                ErrorCategory::Input
            }
        }
    }
}

/// Return `Err(NotFinite)` if `value` is `NaN` or infinite.
///
/// A small shared guard used by the validated constructors throughout
/// the crate so every public entry point rejects non-finite input
/// uniformly.
pub(crate) fn require_finite(name: &'static str, value: f64) -> Result<(), EnduranceError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(EnduranceError::NotFinite { name, value })
    }
}

/// Return `Err(OutOfDomain)` unless `value >= 0.0` (and finite).
pub(crate) fn require_non_negative(
    name: &'static str,
    value: f64,
    reason: &'static str,
) -> Result<(), EnduranceError> {
    require_finite(name, value)?;
    if value < 0.0 {
        return Err(EnduranceError::OutOfDomain {
            name,
            value,
            reason,
        });
    }
    Ok(())
}

/// Return `Err(OutOfDomain)` unless `value > 0.0` (and finite).
pub(crate) fn require_positive(
    name: &'static str,
    value: f64,
    reason: &'static str,
) -> Result<(), EnduranceError> {
    require_finite(name, value)?;
    if value <= 0.0 {
        return Err(EnduranceError::OutOfDomain {
            name,
            value,
            reason,
        });
    }
    Ok(())
}

/// Return `Err(OutOfDomain)` unless `value` lies in the closed interval
/// `[lo, hi]` (and is finite).
pub(crate) fn require_in_closed(
    name: &'static str,
    value: f64,
    lo: f64,
    hi: f64,
    reason: &'static str,
) -> Result<(), EnduranceError> {
    require_finite(name, value)?;
    if value < lo || value > hi {
        return Err(EnduranceError::OutOfDomain {
            name,
            value,
            reason,
        });
    }
    Ok(())
}
