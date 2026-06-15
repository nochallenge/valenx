//! Error taxonomy for the transmission-line calculator.
//!
//! All public constructors funnel their domain checks through
//! [`TlError`], so callers get a single typed error to match on rather
//! than panics or silent `NaN` results.

use thiserror::Error;

/// Errors raised while constructing or evaluating transmission-line
/// quantities.
///
/// Every variant carries enough context (the offending parameter name
/// and the rejected value) to render an actionable message without the
/// caller having to reconstruct it.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum TlError {
    /// A parameter that must be strictly positive was zero or negative.
    ///
    /// Distributed line constants (`L`, `C`) and the characteristic
    /// impedance must be `> 0`; a non-positive value has no physical
    /// meaning for a passive lossless line.
    #[error("parameter `{name}` must be strictly positive, got {value}")]
    NonPositive {
        /// Name of the rejected parameter (e.g. `"inductance_per_m"`).
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A parameter that must be non-negative (a resistance / impedance
    /// magnitude) was negative.
    ///
    /// A load resistance of exactly `0 Ω` (a short) is physically valid
    /// and allowed; only a *negative* resistance is rejected here, as a
    /// purely resistive passive load cannot present negative resistance.
    #[error("parameter `{name}` must be non-negative, got {value}")]
    Negative {
        /// Name of the rejected parameter (e.g. `"load_ohms"`).
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A supplied value was not finite (`NaN` or `±∞`).
    ///
    /// Floating-point inputs are validated up front so that downstream
    /// formulas never silently propagate a non-finite value.
    #[error("parameter `{name}` must be finite, got {value}")]
    NotFinite {
        /// Name of the rejected parameter.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A reflection-coefficient magnitude outside the passive range
    /// `0 ..= 1` was supplied to a constructor that expects `|gamma|`.
    ///
    /// For a passive termination on a lossless line `|gamma| <= 1`; a
    /// magnitude above unity would imply a reflected wave carrying more
    /// power than the incident wave.
    #[error("reflection magnitude must lie in 0..=1, got {value}")]
    GammaOutOfRange {
        /// The rejected magnitude.
        value: f64,
    },
}

impl TlError {
    /// Stable, kebab-cased identifier for this error.
    ///
    /// Useful for logging, metrics, or mapping to UI strings without
    /// matching on the (non-exhaustive in spirit) `Display` text.
    pub fn code(&self) -> &'static str {
        match self {
            TlError::NonPositive { .. } => "transmissionline.non-positive",
            TlError::Negative { .. } => "transmissionline.negative",
            TlError::NotFinite { .. } => "transmissionline.not-finite",
            TlError::GammaOutOfRange { .. } => "transmissionline.gamma-out-of-range",
        }
    }
}

/// Validate that `value` is finite, returning [`TlError::NotFinite`]
/// otherwise.
///
/// Shared helper used by the public constructors so every numeric entry
/// point rejects `NaN` / `±∞` consistently.
pub(crate) fn ensure_finite(name: &'static str, value: f64) -> Result<f64, TlError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(TlError::NotFinite { name, value })
    }
}

/// Validate that `value` is finite and strictly positive.
pub(crate) fn ensure_positive(name: &'static str, value: f64) -> Result<f64, TlError> {
    let value = ensure_finite(name, value)?;
    if value > 0.0 {
        Ok(value)
    } else {
        Err(TlError::NonPositive { name, value })
    }
}

/// Validate that `value` is finite and non-negative.
pub(crate) fn ensure_non_negative(name: &'static str, value: f64) -> Result<f64, TlError> {
    let value = ensure_finite(name, value)?;
    if value >= 0.0 {
        Ok(value)
    } else {
        Err(TlError::Negative { name, value })
    }
}
