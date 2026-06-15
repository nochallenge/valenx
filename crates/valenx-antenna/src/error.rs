//! Error taxonomy for the antenna-link calculations.
//!
//! Every fallible entry point in this crate returns
//! [`Result<_, AntennaError>`]. The variants are intentionally coarse —
//! the only failure modes for these closed-form link relations are
//! out-of-domain inputs (a non-positive frequency, distance, aperture
//! area or wavelength, or a non-finite quantity). Stable
//! [`code`](AntennaError::code) and [`category`](AntennaError::category)
//! accessors are provided for telemetry / host surfaces.

use thiserror::Error;

/// Errors raised by the antenna-link models.
#[derive(Debug, Error, PartialEq)]
pub enum AntennaError {
    /// A physical input that must be strictly positive was zero or
    /// negative (frequency, distance, aperture area, wavelength,
    /// diameter, ...).
    #[error("non-positive {name}: expected a value > 0, got {value}")]
    NonPositive {
        /// Name of the offending parameter.
        name: &'static str,
        /// The supplied value.
        value: f64,
    },

    /// A supplied quantity was `NaN` or infinite. The closed-form link
    /// relations are only defined for finite inputs.
    #[error("non-finite {name}: value is NaN or infinite")]
    NonFinite {
        /// Name of the offending parameter.
        name: &'static str,
    },

    /// A linear power-ratio / gain was negative. Gains expressed as a
    /// linear ratio (not yet in decibels) must be `>= 0`.
    #[error("negative gain ratio for {name}: expected >= 0, got {value}")]
    NegativeGain {
        /// Name of the offending parameter.
        name: &'static str,
        /// The supplied value.
        value: f64,
    },
}

/// Coarse error category, useful for host-side grouping / metrics.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// A user-supplied input was out of the model's valid domain.
    Input,
}

impl AntennaError {
    /// Stable kebab-cased identifier for this error, suitable for logs
    /// and telemetry keys.
    pub fn code(&self) -> &'static str {
        match self {
            AntennaError::NonPositive { .. } => "antenna.non_positive",
            AntennaError::NonFinite { .. } => "antenna.non_finite",
            AntennaError::NegativeGain { .. } => "antenna.negative_gain",
        }
    }

    /// Coarse category for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            AntennaError::NonPositive { .. }
            | AntennaError::NonFinite { .. }
            | AntennaError::NegativeGain { .. } => ErrorCategory::Input,
        }
    }
}

/// Validate that `value` is finite and strictly positive, returning it
/// unchanged on success.
///
/// This is the shared gate used by every constructor / free function in
/// the crate that requires a positive physical quantity.
///
/// # Errors
///
/// Returns [`AntennaError::NonFinite`] if `value` is `NaN` or infinite,
/// or [`AntennaError::NonPositive`] if `value <= 0`.
pub fn require_positive(name: &'static str, value: f64) -> Result<f64, AntennaError> {
    if !value.is_finite() {
        return Err(AntennaError::NonFinite { name });
    }
    if value <= 0.0 {
        return Err(AntennaError::NonPositive { name, value });
    }
    Ok(value)
}

/// Validate that `value` is finite and non-negative (used for linear
/// gain ratios, which may legitimately be zero).
///
/// # Errors
///
/// Returns [`AntennaError::NonFinite`] if `value` is `NaN` or infinite,
/// or [`AntennaError::NegativeGain`] if `value < 0`.
pub fn require_non_negative_gain(name: &'static str, value: f64) -> Result<f64, AntennaError> {
    if !value.is_finite() {
        return Err(AntennaError::NonFinite { name });
    }
    if value < 0.0 {
        return Err(AntennaError::NegativeGain { name, value });
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn positive_passes_through() {
        let v = require_positive("d", 3.0).unwrap();
        assert!((v - 3.0).abs() < 1e-12);
    }

    #[test]
    fn zero_is_non_positive() {
        let e = require_positive("d", 0.0).unwrap_err();
        assert_eq!(
            e,
            AntennaError::NonPositive {
                name: "d",
                value: 0.0
            }
        );
        assert_eq!(e.code(), "antenna.non_positive");
        assert_eq!(e.category(), ErrorCategory::Input);
    }

    #[test]
    fn negative_is_non_positive() {
        let e = require_positive("f", -1.5).unwrap_err();
        assert_eq!(
            e,
            AntennaError::NonPositive {
                name: "f",
                value: -1.5
            }
        );
    }

    #[test]
    fn nan_is_non_finite() {
        let e = require_positive("d", f64::NAN).unwrap_err();
        assert_eq!(e, AntennaError::NonFinite { name: "d" });
        assert_eq!(e.code(), "antenna.non_finite");
    }

    #[test]
    fn infinite_is_non_finite() {
        let e = require_positive("d", f64::INFINITY).unwrap_err();
        assert_eq!(e, AntennaError::NonFinite { name: "d" });
    }

    #[test]
    fn gain_zero_allowed() {
        let v = require_non_negative_gain("g", 0.0).unwrap();
        assert!((v - 0.0).abs() < 1e-12);
    }

    #[test]
    fn gain_negative_rejected() {
        let e = require_non_negative_gain("g", -0.1).unwrap_err();
        assert_eq!(
            e,
            AntennaError::NegativeGain {
                name: "g",
                value: -0.1
            }
        );
        assert_eq!(e.code(), "antenna.negative_gain");
        assert_eq!(e.category(), ErrorCategory::Input);
    }

    #[test]
    fn gain_nan_rejected() {
        let e = require_non_negative_gain("g", f64::NAN).unwrap_err();
        assert_eq!(e, AntennaError::NonFinite { name: "g" });
    }

    #[test]
    fn error_display_is_informative() {
        let e = AntennaError::NonPositive {
            name: "freq_hz",
            value: -2.0,
        };
        let s = format!("{e}");
        assert!(s.contains("freq_hz"));
        assert!(s.contains("-2"));
    }
}
