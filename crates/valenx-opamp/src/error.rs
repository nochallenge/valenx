//! Error taxonomy for op-amp circuit models.
//!
//! Every fallible constructor in this crate funnels its rejections
//! through [`OpAmpError`]. The variants distinguish *why* a parameter
//! was rejected (non-finite, non-positive, empty input set) so a caller
//! can surface a precise message or branch on [`OpAmpError::code`].

use thiserror::Error;

/// A convenience alias for results produced by this crate.
pub type Result<T> = core::result::Result<T, OpAmpError>;

/// Errors raised when constructing or evaluating an op-amp circuit.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum OpAmpError {
    /// A parameter was `NaN` or `±∞`.
    ///
    /// Resistances, voltages and frequencies must be finite real
    /// numbers; a non-finite value almost always signals an upstream
    /// computation error and is rejected at the boundary.
    #[error("parameter `{name}` must be finite, got {value}")]
    NotFinite {
        /// The offending parameter's name.
        name: &'static str,
        /// The non-finite value that was supplied.
        value: f64,
    },

    /// A parameter that must be strictly positive was `<= 0`.
    ///
    /// Resistor values and the gain-bandwidth product are physical
    /// magnitudes: zero or negative values have no meaning in the ideal
    /// model and would produce a division by zero or a sign flip.
    #[error("parameter `{name}` must be > 0, got {value}")]
    NonPositive {
        /// The offending parameter's name.
        name: &'static str,
        /// The non-positive value that was supplied.
        value: f64,
    },

    /// A summing amplifier was built with no input branches.
    ///
    /// `out = -Rf * Σ(Vᵢ/Rᵢ)` over an empty set would be a trivial zero
    /// and is far more likely to be a caller mistake than an intent, so
    /// it is rejected.
    #[error("summing amplifier requires at least one input branch")]
    NoInputs,
}

impl OpAmpError {
    /// A stable, kebab-cased identifier for this error.
    ///
    /// Useful for logging, telemetry, or matching in tests without
    /// depending on the human-readable [`Display`](std::fmt::Display)
    /// string.
    pub fn code(&self) -> &'static str {
        match self {
            OpAmpError::NotFinite { .. } => "opamp.not-finite",
            OpAmpError::NonPositive { .. } => "opamp.non-positive",
            OpAmpError::NoInputs => "opamp.no-inputs",
        }
    }
}

/// Validate that `value` is finite, returning [`OpAmpError::NotFinite`]
/// otherwise.
///
/// This is the shared front-half of every numeric check in the crate:
/// callers run it first, then layer a sign / magnitude check on top.
pub(crate) fn ensure_finite(name: &'static str, value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(OpAmpError::NotFinite { name, value })
    }
}

/// Validate that `value` is finite and strictly positive.
///
/// Returns the value unchanged on success, or the appropriate
/// [`OpAmpError`] variant ([`NotFinite`](OpAmpError::NotFinite) takes
/// precedence over [`NonPositive`](OpAmpError::NonPositive)).
pub(crate) fn ensure_positive(name: &'static str, value: f64) -> Result<f64> {
    let value = ensure_finite(name, value)?;
    if value > 0.0 {
        Ok(value)
    } else {
        Err(OpAmpError::NonPositive { name, value })
    }
}
