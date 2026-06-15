//! Thermistor modelling error taxonomy.
//!
//! Every fallible constructor and conversion in the crate returns
//! [`Result<_, ThermistorError>`]. The error type carries a stable
//! kebab-cased [`code`](ThermistorError::code) and a coarse
//! [`category`](ThermistorError::category) so callers (CLIs, UIs, an
//! LLM/MCP surface) can branch on the failure class without string
//! matching on the human-readable message.

use thiserror::Error;

/// Errors raised while building thermistor models or converting between
/// resistance and temperature.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum ThermistorError {
    /// A scalar parameter was outside its valid domain.
    ///
    /// Used for non-positive resistances, non-positive absolute
    /// temperatures, and a non-positive or non-finite `beta`.
    #[error("bad parameter `{name}` = {value}: {reason}")]
    BadParameter {
        /// Name of the offending parameter (e.g. `"r0"`, `"t0_kelvin"`).
        name: &'static str,
        /// The rejected value, formatted for display.
        value: f64,
        /// Why the value was rejected.
        reason: &'static str,
    },

    /// Two calibration points that must differ were equal (or too
    /// close to distinguish), so the model cannot be solved.
    ///
    /// Raised by beta calibration when the two temperatures coincide,
    /// and by Steinhart-Hart fitting when two of the three sample
    /// resistances coincide.
    #[error("degenerate calibration: {0}")]
    Degenerate(&'static str),

    /// A computed quantity left the representable range — typically an
    /// `exp`/`ln` argument that overflowed to a non-finite value, or a
    /// Steinhart-Hart polynomial that yielded a non-positive `1/T`.
    #[error("non-finite result: {0}")]
    NonFinite(&'static str),
}

/// Coarse classification of a [`ThermistorError`].
///
/// Mirrors the convention used across the Valenx numerical crates: a
/// small enum that groups failures by who is responsible for them.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// Caller supplied an out-of-domain input value.
    Input,
    /// The model is unsolvable as posed (degenerate calibration data).
    Calibration,
    /// The arithmetic left the representable floating-point range.
    Numeric,
}

impl ThermistorError {
    /// Stable, kebab-cased identifier for programmatic matching.
    ///
    /// The string is part of the crate's public contract and will not
    /// change for a given variant across patch releases.
    pub fn code(&self) -> &'static str {
        match self {
            ThermistorError::BadParameter { .. } => "thermistor.bad_parameter",
            ThermistorError::Degenerate(_) => "thermistor.degenerate",
            ThermistorError::NonFinite(_) => "thermistor.non_finite",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            ThermistorError::BadParameter { .. } => ErrorCategory::Input,
            ThermistorError::Degenerate(_) => ErrorCategory::Calibration,
            ThermistorError::NonFinite(_) => ErrorCategory::Numeric,
        }
    }
}

/// Validate that a resistance is strictly positive and finite, in ohms.
///
/// Returns the value unchanged on success so it can be used inline.
///
/// # Errors
///
/// Returns [`ThermistorError::BadParameter`] if `r` is not strictly
/// positive or is not finite.
pub(crate) fn check_resistance(name: &'static str, r: f64) -> Result<f64, ThermistorError> {
    if !r.is_finite() {
        return Err(ThermistorError::BadParameter {
            name,
            value: r,
            reason: "resistance must be finite",
        });
    }
    if r <= 0.0 {
        return Err(ThermistorError::BadParameter {
            name,
            value: r,
            reason: "resistance must be strictly positive (ohms)",
        });
    }
    Ok(r)
}

/// Validate that an absolute temperature is strictly positive and
/// finite, in kelvin.
///
/// Returns the value unchanged on success.
///
/// # Errors
///
/// Returns [`ThermistorError::BadParameter`] if `t` is not strictly
/// positive (absolute temperatures cannot be zero or negative) or is
/// not finite.
pub(crate) fn check_temperature(name: &'static str, t: f64) -> Result<f64, ThermistorError> {
    if !t.is_finite() {
        return Err(ThermistorError::BadParameter {
            name,
            value: t,
            reason: "temperature must be finite",
        });
    }
    if t <= 0.0 {
        return Err(ThermistorError::BadParameter {
            name,
            value: t,
            reason: "absolute temperature must be strictly positive (kelvin)",
        });
    }
    Ok(t)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_non_positive_resistance() {
        let err = check_resistance("r", 0.0).unwrap_err();
        assert_eq!(err.code(), "thermistor.bad_parameter");
        assert_eq!(err.category(), ErrorCategory::Input);
    }

    #[test]
    fn rejects_negative_temperature() {
        let err = check_temperature("t", -5.0).unwrap_err();
        assert_eq!(err.code(), "thermistor.bad_parameter");
        assert_eq!(err.category(), ErrorCategory::Input);
    }

    #[test]
    fn rejects_non_finite() {
        assert!(check_resistance("r", f64::NAN).is_err());
        assert!(check_temperature("t", f64::INFINITY).is_err());
    }

    #[test]
    fn accepts_valid_values_unchanged() {
        let r = check_resistance("r", 10_000.0).unwrap();
        assert!((r - 10_000.0).abs() < 1e-12);
        let t = check_temperature("t", 298.15).unwrap();
        assert!((t - 298.15).abs() < 1e-12);
    }

    #[test]
    fn distinct_codes_and_categories() {
        let deg = ThermistorError::Degenerate("x");
        assert_eq!(deg.code(), "thermistor.degenerate");
        assert_eq!(deg.category(), ErrorCategory::Calibration);

        let nf = ThermistorError::NonFinite("x");
        assert_eq!(nf.code(), "thermistor.non_finite");
        assert_eq!(nf.category(), ErrorCategory::Numeric);
    }
}
