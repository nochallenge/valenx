//! Error taxonomy for dimensionless-group construction.
//!
//! Every public constructor in this crate validates its inputs and
//! returns a [`DimensionlessError`] on bad data rather than producing a
//! `NaN`, an infinity, or a silently-wrong value. The two failure modes
//! are a non-finite input (`NaN` / `±inf`) and an out-of-domain input
//! (for example a non-positive density or a zero denominator).

use thiserror::Error;

/// Errors raised when building a dimensionless group.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum DimensionlessError {
    /// A named input was not finite (it was `NaN` or `±inf`).
    #[error("parameter `{name}` must be finite, got {value}")]
    NotFinite {
        /// Name of the offending parameter.
        name: &'static str,
        /// The non-finite value that was supplied.
        value: f64,
    },

    /// A named input fell outside its required domain — for example a
    /// density that must be strictly positive but was zero or negative,
    /// or a length that must be non-zero.
    #[error("parameter `{name}` is out of domain: {reason} (got {value})")]
    OutOfDomain {
        /// Name of the offending parameter.
        name: &'static str,
        /// Human-readable statement of the required domain.
        reason: &'static str,
        /// The value that violated the domain.
        value: f64,
    },
}

impl DimensionlessError {
    /// Stable, kebab-cased identifier for programmatic matching and
    /// logging. Distinct from the `Display` message, which is meant for
    /// humans and may change wording.
    pub fn code(&self) -> &'static str {
        match self {
            DimensionlessError::NotFinite { .. } => "dimensionless.not-finite",
            DimensionlessError::OutOfDomain { .. } => "dimensionless.out-of-domain",
        }
    }

    /// Name of the parameter that triggered the error.
    pub fn parameter(&self) -> &'static str {
        match self {
            DimensionlessError::NotFinite { name, .. } => name,
            DimensionlessError::OutOfDomain { name, .. } => name,
        }
    }
}

/// Require that `value` is finite, otherwise return
/// [`DimensionlessError::NotFinite`]. Internal helper shared by every
/// module's constructors.
pub(crate) fn require_finite(name: &'static str, value: f64) -> Result<f64, DimensionlessError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DimensionlessError::NotFinite { name, value })
    }
}

/// Require that `value` is finite **and** strictly greater than zero.
/// Used for quantities like density, viscosity, conductivity, length,
/// and speed of sound that are physically positive and appear in a
/// denominator or under a square root.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64, DimensionlessError> {
    let value = require_finite(name, value)?;
    if value > 0.0 {
        Ok(value)
    } else {
        Err(DimensionlessError::OutOfDomain {
            name,
            reason: "must be strictly positive",
            value,
        })
    }
}

/// Require that `value` is finite and not negative (zero is allowed).
/// Used for quantities such as velocity magnitude or a heat-transfer
/// coefficient that may legitimately be zero but never negative.
pub(crate) fn require_non_negative(
    name: &'static str,
    value: f64,
) -> Result<f64, DimensionlessError> {
    let value = require_finite(name, value)?;
    if value >= 0.0 {
        Ok(value)
    } else {
        Err(DimensionlessError::OutOfDomain {
            name,
            reason: "must be zero or positive",
            value,
        })
    }
}
