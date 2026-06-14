//! Thermoregulation error taxonomy.
//!
//! Every validated constructor in this crate funnels its rejection
//! through [`ThermoregError`]. The enum is deliberately small: the
//! inputs are physical quantities, so the only ways they go wrong are
//! "non-finite", "out of physical range" (e.g. a negative mass), or a
//! genuinely degenerate combination.

use thiserror::Error;

/// Errors raised when constructing or evaluating a thermoregulation
/// model from caller-supplied numbers.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum ThermoregError {
    /// A supplied value was `NaN` or infinite.
    ///
    /// Carries the parameter name so the caller can see which input
    /// was bad.
    #[error("parameter `{name}` must be finite, got {value}")]
    NotFinite {
        /// Name of the offending parameter.
        name: &'static str,
        /// The non-finite value that was supplied.
        value: f64,
    },

    /// A value fell outside its allowed physical range.
    #[error("parameter `{name}` = {value} is out of range ({reason})")]
    OutOfRange {
        /// Name of the offending parameter.
        name: &'static str,
        /// The value that was supplied.
        value: f64,
        /// Human-readable statement of the permitted range.
        reason: &'static str,
    },
}

/// Coarse classification of a [`ThermoregError`], useful for callers
/// (UI / MCP) that want to react to a *class* of failure rather than
/// match every variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The input was malformed (non-finite).
    Malformed,
    /// The input was well-formed but physically invalid.
    Domain,
}

impl ThermoregError {
    /// Stable, kebab-cased identifier for this error, suitable for
    /// logging or as a machine-readable error code.
    pub fn code(&self) -> &'static str {
        match self {
            ThermoregError::NotFinite { .. } => "thermoreg.not-finite",
            ThermoregError::OutOfRange { .. } => "thermoreg.out-of-range",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            ThermoregError::NotFinite { .. } => ErrorCategory::Malformed,
            ThermoregError::OutOfRange { .. } => ErrorCategory::Domain,
        }
    }
}

/// Reject a non-finite value.
///
/// Returns `Ok(value)` when `value.is_finite()`, otherwise a
/// [`ThermoregError::NotFinite`] tagged with `name`.
pub(crate) fn finite(name: &'static str, value: f64) -> Result<f64, ThermoregError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ThermoregError::NotFinite { name, value })
    }
}

/// Reject a value that is not strictly greater than zero.
///
/// First checks finiteness, then positivity. Used for quantities such
/// as mass and surface area that are meaningless at zero or below.
pub(crate) fn positive(name: &'static str, value: f64) -> Result<f64, ThermoregError> {
    let value = finite(name, value)?;
    if value > 0.0 {
        Ok(value)
    } else {
        Err(ThermoregError::OutOfRange {
            name,
            value,
            reason: "must be > 0",
        })
    }
}

/// Reject a value that is negative (zero is permitted).
///
/// Used for quantities such as a sweat rate or a convective
/// coefficient that may legitimately be zero but never negative.
pub(crate) fn non_negative(name: &'static str, value: f64) -> Result<f64, ThermoregError> {
    let value = finite(name, value)?;
    if value >= 0.0 {
        Ok(value)
    } else {
        Err(ThermoregError::OutOfRange {
            name,
            value,
            reason: "must be >= 0",
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finite_passes_through_real_value() {
        assert_eq!(finite("x", 3.5).unwrap(), 3.5);
    }

    #[test]
    fn finite_rejects_nan_and_inf() {
        let e = finite("x", f64::NAN).unwrap_err();
        assert_eq!(e.code(), "thermoreg.not-finite");
        assert_eq!(e.category(), ErrorCategory::Malformed);
        assert!(finite("x", f64::INFINITY).is_err());
    }

    #[test]
    fn positive_rejects_zero_and_negative() {
        assert!(positive("m", 0.0).is_err());
        assert!(positive("m", -1.0).is_err());
        assert_eq!(positive("m", 70.0).unwrap(), 70.0);
    }

    #[test]
    fn non_negative_allows_zero_but_not_negative() {
        assert_eq!(non_negative("s", 0.0).unwrap(), 0.0);
        let e = non_negative("s", -0.1).unwrap_err();
        assert_eq!(e.category(), ErrorCategory::Domain);
        assert_eq!(e.code(), "thermoreg.out-of-range");
    }

    #[test]
    fn error_is_displayable() {
        let e = ThermoregError::OutOfRange {
            name: "m",
            value: -1.0,
            reason: "must be > 0",
        };
        let s = format!("{e}");
        assert!(s.contains("out of range"));
        assert!(s.contains("must be > 0"));
    }
}
