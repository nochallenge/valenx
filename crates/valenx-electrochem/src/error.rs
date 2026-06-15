//! Electrochemistry error taxonomy.
//!
//! Every fallible constructor in this crate funnels through
//! [`ElectrochemError`]. The variants carry enough context (the offending
//! parameter name and a human-readable reason) to render a useful message,
//! and [`ElectrochemError::code`] gives a stable kebab-cased identifier for
//! programmatic matching, logging, and tests.

use thiserror::Error;

/// Errors raised when validating electrochemical inputs.
#[derive(Debug, Error)]
pub enum ElectrochemError {
    /// A numeric parameter fell outside its physically meaningful domain.
    ///
    /// Examples: a non-positive number of electrons transferred, a
    /// non-positive absolute temperature, or a reaction quotient `<= 0`
    /// (the natural log in the Nernst equation is only defined for a
    /// strictly positive argument).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// The name of the offending parameter (e.g. `"n"`, `"temperature_k"`).
        name: &'static str,
        /// A human-readable explanation of why the value was rejected.
        reason: String,
    },

    /// A supplied value was required to be finite but was `NaN` or infinite.
    ///
    /// Guards the public constructors against `NaN` / `inf` leaking into a
    /// calculation, where it would silently poison every downstream result.
    #[error("non-finite value for `{name}`")]
    NonFinite {
        /// The name of the parameter that was not finite.
        name: &'static str,
    },
}

/// A coarse classification of an [`ElectrochemError`].
///
/// Useful when a caller wants to react to broad failure classes (for
/// instance, surfacing input errors to a user while logging algorithm
/// errors) without matching every individual variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// The caller supplied an invalid value (out of domain, wrong sign).
    Input,
    /// A value was not a finite real number.
    Numeric,
}

impl ElectrochemError {
    /// A stable, kebab-cased identifier for this error.
    ///
    /// The returned string is part of the crate's API contract: it is safe
    /// to match on in tests and downstream code, and will not change for a
    /// given variant across patch releases.
    pub fn code(&self) -> &'static str {
        match self {
            ElectrochemError::BadParameter { .. } => "electrochem.bad_parameter",
            ElectrochemError::NonFinite { .. } => "electrochem.non_finite",
        }
    }

    /// The coarse [`ErrorCategory`] this error belongs to.
    pub fn category(&self) -> ErrorCategory {
        match self {
            ElectrochemError::BadParameter { .. } => ErrorCategory::Input,
            ElectrochemError::NonFinite { .. } => ErrorCategory::Numeric,
        }
    }
}

/// Internal helper: reject non-finite inputs at the crate boundary.
///
/// Returns `Ok(value)` when `value` is finite, otherwise an
/// [`ElectrochemError::NonFinite`] tagged with `name`.
pub(crate) fn require_finite(value: f64, name: &'static str) -> Result<f64, ElectrochemError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ElectrochemError::NonFinite { name })
    }
}

/// Internal helper: require `value > 0` for an already-finite input.
///
/// `name` labels the parameter and `reason` describes the constraint in the
/// error message when the check fails.
pub(crate) fn require_positive(
    value: f64,
    name: &'static str,
    reason: &str,
) -> Result<f64, ElectrochemError> {
    let value = require_finite(value, name)?;
    if value > 0.0 {
        Ok(value)
    } else {
        Err(ElectrochemError::BadParameter {
            name,
            reason: reason.to_string(),
        })
    }
}

/// Internal helper: require `value >= 0` for an already-finite input.
///
/// Used for quantities (such as elapsed time or transferred charge) where
/// zero is physically meaningful but a negative value is not.
pub(crate) fn require_non_negative(
    value: f64,
    name: &'static str,
    reason: &str,
) -> Result<f64, ElectrochemError> {
    let value = require_finite(value, name)?;
    if value >= 0.0 {
        Ok(value)
    } else {
        Err(ElectrochemError::BadParameter {
            name,
            reason: reason.to_string(),
        })
    }
}
