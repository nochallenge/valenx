//! Creep / stress-rupture error taxonomy.
//!
//! Every fallible model in this crate validates its inputs through one
//! of the constructors here, so an invalid physical quantity (a
//! non-positive temperature, a negative rupture time, a non-finite
//! coefficient) is surfaced as a typed [`CreepError`] rather than
//! producing a silent `NaN` / `inf`.

use thiserror::Error;

/// Errors raised by the creep and stress-rupture models.
#[derive(Debug, Error)]
pub enum CreepError {
    /// A scalar input was outside its valid domain (non-positive,
    /// negative, or otherwise unphysical).
    #[error("invalid value for `{name}`: {value} ({reason})")]
    InvalidValue {
        /// Name of the offending parameter.
        name: &'static str,
        /// The value that was supplied.
        value: f64,
        /// Why it was rejected.
        reason: &'static str,
    },

    /// A scalar input was not a finite number (`NaN` or `+/-inf`).
    #[error("non-finite value for `{name}`: {value}")]
    NotFinite {
        /// Name of the offending parameter.
        name: &'static str,
        /// The value that was supplied.
        value: f64,
    },
}

/// Coarse category for an error, useful for callers that want to react
/// differently to bad user input versus a malformed material constant.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The supplied value is the wrong magnitude or sign.
    Domain,
    /// The supplied value is `NaN` or infinite.
    NotFinite,
}

impl CreepError {
    /// Stable, kebab-cased identifier for the error variant.
    ///
    /// Useful for logging or matching without depending on the human
    /// readable [`Display`](std::fmt::Display) text.
    pub fn code(&self) -> &'static str {
        match self {
            CreepError::InvalidValue { .. } => "creep.invalid_value",
            CreepError::NotFinite { .. } => "creep.not_finite",
        }
    }

    /// Coarse [`ErrorCategory`] for the variant.
    pub fn category(&self) -> ErrorCategory {
        match self {
            CreepError::InvalidValue { .. } => ErrorCategory::Domain,
            CreepError::NotFinite { .. } => ErrorCategory::NotFinite,
        }
    }
}

/// Require that `value` is a finite number, returning a typed error
/// otherwise. Shared by the validating constructors throughout the
/// crate so the finiteness check is written once.
pub(crate) fn require_finite(name: &'static str, value: f64) -> Result<f64, CreepError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(CreepError::NotFinite { name, value })
    }
}

/// Require that `value` is finite and strictly greater than zero.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64, CreepError> {
    let value = require_finite(name, value)?;
    if value > 0.0 {
        Ok(value)
    } else {
        Err(CreepError::InvalidValue {
            name,
            value,
            reason: "must be strictly positive",
        })
    }
}

/// Require that `value` is finite and not negative (zero is allowed).
pub(crate) fn require_non_negative(name: &'static str, value: f64) -> Result<f64, CreepError> {
    let value = require_finite(name, value)?;
    if value >= 0.0 {
        Ok(value)
    } else {
        Err(CreepError::InvalidValue {
            name,
            value,
            reason: "must not be negative",
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn positive_accepts_positive_rejects_zero_and_negative() {
        assert_eq!(require_positive("x", 1.5).unwrap(), 1.5);
        assert!(matches!(
            require_positive("x", 0.0),
            Err(CreepError::InvalidValue { .. })
        ));
        assert!(matches!(
            require_positive("x", -2.0),
            Err(CreepError::InvalidValue { .. })
        ));
    }

    #[test]
    fn non_negative_accepts_zero_rejects_negative() {
        assert_eq!(require_non_negative("x", 0.0).unwrap(), 0.0);
        assert!(matches!(
            require_non_negative("x", -0.1),
            Err(CreepError::InvalidValue { .. })
        ));
    }

    #[test]
    fn non_finite_rejected() {
        assert!(matches!(
            require_finite("x", f64::NAN),
            Err(CreepError::NotFinite { .. })
        ));
        assert!(matches!(
            require_positive("x", f64::INFINITY),
            Err(CreepError::NotFinite { .. })
        ));
    }

    #[test]
    fn code_and_category_are_stable() {
        let domain = CreepError::InvalidValue {
            name: "x",
            value: -1.0,
            reason: "r",
        };
        assert_eq!(domain.code(), "creep.invalid_value");
        assert_eq!(domain.category(), ErrorCategory::Domain);

        let nf = CreepError::NotFinite {
            name: "x",
            value: f64::NAN,
        };
        assert_eq!(nf.code(), "creep.not_finite");
        assert_eq!(nf.category(), ErrorCategory::NotFinite);
    }
}
