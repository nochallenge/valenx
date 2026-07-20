//! Error taxonomy for the RLC-resonance models.
//!
//! Every public computation validates its inputs through these
//! constructors before touching a `sqrt`/`ln`/division, so a caller can
//! never drive a silent `NaN`/`Inf` into a result. The physics here
//! demands strictly positive inductance, capacitance and frequency, and
//! a non-negative resistance, so the constructors reject anything
//! non-finite or out of that domain.

use thiserror::Error;

/// Shorthand for `Result<T, RlcError>`.
pub type Result<T> = core::result::Result<T, RlcError>;

/// Anything that can go wrong validating or evaluating an RLC model.
///
/// This enum is `#[non_exhaustive]`: future releases may add variants
/// without it being a breaking change, so downstream `match` arms must
/// include a wildcard.
#[derive(Debug, Error, Clone, PartialEq)]
#[non_exhaustive]
pub enum RlcError {
    /// A quantity that physics requires to be strictly positive was zero,
    /// negative, or non-finite (`NaN`/`±Inf`). Inductance `L`,
    /// capacitance `C`, and any frequency fall in this class.
    #[error("`{name}` must be finite and > 0, got {value}")]
    NonPositive {
        /// Which quantity was bad (`"inductance"`, `"capacitance"`, ...).
        name: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A quantity that may be zero but never negative or non-finite was
    /// out of range. Resistance `R` falls in this class: `R = 0` is the
    /// ideal lossless limit (`Q → ∞`), which the API reports explicitly
    /// rather than dividing by zero.
    #[error("`{name}` must be finite and >= 0, got {value}")]
    Negative {
        /// Which quantity was bad (`"resistance"`).
        name: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A model was evaluated at a singular point where its closed form
    /// blows up — e.g. the series bandwidth `BW = R/(2 pi L)` is finite,
    /// but the loaded quality factor `Q = (1/R) sqrt(L/C)` diverges as
    /// `R → 0` (a lossless tank has no defined `−3 dB` bandwidth ratio).
    /// Carries a short reason.
    #[error("singular model: {0}")]
    Singular(&'static str),
}

impl RlcError {
    /// Stable, machine-readable identifier for this error.
    ///
    /// Kebab/dotted form so logs and tests can match on a constant
    /// instead of the human-readable `Display` string.
    pub fn code(&self) -> &'static str {
        match self {
            RlcError::NonPositive { .. } => "rlc.non_positive",
            RlcError::Negative { .. } => "rlc.negative",
            RlcError::Singular(_) => "rlc.singular",
        }
    }
}

/// Validate that `value` is finite and strictly positive.
///
/// Used for inductance, capacitance and frequency, all of which must be
/// `> 0` for the resonance formulas to be well defined.
///
/// # Errors
///
/// Returns [`RlcError::NonPositive`] when `value` is `NaN`, `±Inf`, zero
/// or negative.
pub fn require_positive(name: &'static str, value: f64) -> Result<f64> {
    if !value.is_finite() || value <= 0.0 {
        Err(RlcError::NonPositive { name, value })
    } else {
        Ok(value)
    }
}

/// Validate that `value` is finite and non-negative.
///
/// Used for resistance, where `0` is the meaningful lossless limit but a
/// negative or non-finite value is unphysical.
///
/// # Errors
///
/// Returns [`RlcError::Negative`] when `value` is `NaN`, `±Inf` or `< 0`.
pub fn require_non_negative(name: &'static str, value: f64) -> Result<f64> {
    if !value.is_finite() || value < 0.0 {
        Err(RlcError::Negative { name, value })
    } else {
        Ok(value)
    }
}
