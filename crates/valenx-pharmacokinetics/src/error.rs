//! Error taxonomy for `valenx-pharmacokinetics`.
//!
//! Every fallible constructor returns [`Result<_, PkError>`]. The model
//! parameters are physical quantities with sign and finiteness
//! constraints (a volume of distribution and a clearance are strictly
//! positive, a dose is non-negative, a Hill coefficient is strictly
//! positive); the validated constructors reject anything that violates
//! those so that the downstream closed-form expressions stay well-defined
//! (no division by zero, no negative concentrations).

use thiserror::Error;

/// Errors raised when constructing a pharmacokinetic / pharmacodynamic
/// model from raw parameters.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum PkError {
    /// A parameter that must be strictly positive (`> 0`) was not — for
    /// example a zero or negative volume of distribution, clearance,
    /// EC50, or Hill coefficient. A non-finite value (NaN / infinity)
    /// also lands here.
    #[error("parameter `{name}` must be strictly positive and finite, got {value}")]
    NotPositive {
        /// The offending parameter's name (stable, kebab-free identifier).
        name: &'static str,
        /// The value that was supplied.
        value: f64,
    },

    /// A parameter that must be non-negative (`>= 0`) was negative — for
    /// example a negative dose or a negative time. A non-finite value
    /// (NaN / infinity) also lands here.
    #[error("parameter `{name}` must be non-negative and finite, got {value}")]
    Negative {
        /// The offending parameter's name.
        name: &'static str,
        /// The value that was supplied.
        value: f64,
    },
}

/// Convenience alias for `Result<T, PkError>`.
pub type Result<T> = std::result::Result<T, PkError>;

/// Validate that `value` is strictly positive and finite, returning it on
/// success or [`PkError::NotPositive`] otherwise.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(PkError::NotPositive { name, value })
    }
}

/// Validate that `value` is non-negative and finite, returning it on
/// success or [`PkError::Negative`] otherwise.
pub(crate) fn require_non_negative(name: &'static str, value: f64) -> Result<f64> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(PkError::Negative { name, value })
    }
}
