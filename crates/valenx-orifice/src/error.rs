//! Error taxonomy for the differential-pressure flow-meter models.
//!
//! Every fallible constructor in this crate returns [`OrificeError`].
//! The variants are deliberately fine-grained so callers (and tests)
//! can match on the precise failure rather than a single opaque string.

use thiserror::Error;

/// Errors raised while building or evaluating a flow meter.
///
/// All numeric inputs are physical quantities that must be finite and,
/// for most of them, strictly positive. The variants below name the
/// offending parameter and the constraint it violated.
#[derive(Debug, Error)]
pub enum OrificeError {
    /// A scalar parameter was not a finite number (it was `NaN` or an
    /// infinity). Stored is the parameter name.
    #[error("parameter `{name}` must be finite, got a non-finite value")]
    NotFinite {
        /// Name of the offending parameter.
        name: &'static str,
    },

    /// A scalar parameter that must be strictly positive was zero or
    /// negative.
    #[error("parameter `{name}` must be strictly positive, got {value}")]
    NonPositive {
        /// Name of the offending parameter.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A scalar parameter that must be non-negative (zero allowed) was
    /// negative. Used for the pressure drop and the flow rate, where a
    /// value of exactly zero is physically meaningful.
    #[error("parameter `{name}` must be non-negative, got {value}")]
    Negative {
        /// Name of the offending parameter.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// The throat (bore) diameter `d` was not smaller than the pipe
    /// diameter `D`. The diameter ratio `beta = d / D` must lie in the
    /// open interval `(0, 1)`, so `d` must be strictly less than `D`.
    #[error("throat diameter d = {d} must be strictly less than pipe diameter D = {pipe}")]
    ThroatNotSmaller {
        /// The throat / bore diameter.
        d: f64,
        /// The upstream pipe diameter.
        pipe: f64,
    },

    /// The discharge coefficient `Cd` fell outside the physically
    /// sensible half-open range `(0, 1]`. A real differential-pressure
    /// meter always loses energy, so `Cd <= 1`; a value of zero or
    /// negative would describe a meter that passes no (or reversed)
    /// flow at a finite pressure drop.
    #[error("discharge coefficient Cd = {value} must lie in (0, 1]")]
    DischargeCoefficientOutOfRange {
        /// The rejected value.
        value: f64,
    },
}

/// Coarse classification of an [`OrificeError`], handy for routing or
/// for deciding whether a failure is the caller's fault (bad input) or
/// a domain-boundary violation in the model.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// The value supplied by the caller was malformed (non-finite or out
    /// of the allowed sign range).
    Input,
    /// The combination of values is geometrically or physically
    /// inconsistent (e.g. throat larger than the pipe).
    Domain,
}

impl OrificeError {
    /// A stable, machine-readable kebab/dotted identifier for the
    /// variant. Useful for logging and for assertions in tests that do
    /// not want to depend on the human-readable [`Display`] string.
    ///
    /// [`Display`]: std::fmt::Display
    #[must_use]
    pub fn code(&self) -> &'static str {
        match self {
            OrificeError::NotFinite { .. } => "orifice.not_finite",
            OrificeError::NonPositive { .. } => "orifice.non_positive",
            OrificeError::Negative { .. } => "orifice.negative",
            OrificeError::ThroatNotSmaller { .. } => "orifice.throat_not_smaller",
            OrificeError::DischargeCoefficientOutOfRange { .. } => {
                "orifice.discharge_coefficient_out_of_range"
            }
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    #[must_use]
    pub fn category(&self) -> ErrorCategory {
        match self {
            OrificeError::NotFinite { .. }
            | OrificeError::NonPositive { .. }
            | OrificeError::Negative { .. }
            | OrificeError::DischargeCoefficientOutOfRange { .. } => ErrorCategory::Input,
            OrificeError::ThroatNotSmaller { .. } => ErrorCategory::Domain,
        }
    }
}

/// Validate that `value` is finite, returning [`OrificeError::NotFinite`]
/// otherwise. Internal helper shared by the public constructors.
pub(crate) fn require_finite(name: &'static str, value: f64) -> Result<f64, OrificeError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(OrificeError::NotFinite { name })
    }
}

/// Validate that `value` is finite and strictly positive.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64, OrificeError> {
    let value = require_finite(name, value)?;
    if value > 0.0 {
        Ok(value)
    } else {
        Err(OrificeError::NonPositive { name, value })
    }
}

/// Validate that `value` is finite and non-negative (zero allowed).
pub(crate) fn require_non_negative(name: &'static str, value: f64) -> Result<f64, OrificeError> {
    let value = require_finite(name, value)?;
    if value >= 0.0 {
        Ok(value)
    } else {
        Err(OrificeError::Negative { name, value })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn require_finite_accepts_finite_and_rejects_nan_inf() {
        assert!((require_finite("x", 1.5).unwrap() - 1.5).abs() < 1e-15);
        assert!(matches!(
            require_finite("x", f64::NAN).unwrap_err(),
            OrificeError::NotFinite { name: "x" }
        ));
        assert!(matches!(
            require_finite("x", f64::NEG_INFINITY).unwrap_err(),
            OrificeError::NotFinite { name: "x" }
        ));
    }

    #[test]
    fn require_positive_rejects_zero_and_negative() {
        assert!((require_positive("x", 2.0).unwrap() - 2.0).abs() < 1e-15);
        assert!(matches!(
            require_positive("x", 0.0).unwrap_err(),
            OrificeError::NonPositive { name: "x", .. }
        ));
        assert!(matches!(
            require_positive("x", -3.0).unwrap_err(),
            OrificeError::NonPositive { name: "x", .. }
        ));
    }

    #[test]
    fn require_non_negative_allows_zero_but_not_negative() {
        assert!(require_non_negative("x", 0.0).unwrap().abs() < 1e-15);
        assert!((require_non_negative("x", 4.0).unwrap() - 4.0).abs() < 1e-15);
        assert!(matches!(
            require_non_negative("x", -0.5).unwrap_err(),
            OrificeError::Negative { name: "x", .. }
        ));
    }

    #[test]
    fn codes_are_distinct_and_stable() {
        let errs = [
            OrificeError::NotFinite { name: "a" },
            OrificeError::NonPositive {
                name: "a",
                value: -1.0,
            },
            OrificeError::Negative {
                name: "a",
                value: -1.0,
            },
            OrificeError::ThroatNotSmaller { d: 2.0, pipe: 1.0 },
            OrificeError::DischargeCoefficientOutOfRange { value: 2.0 },
        ];
        let codes: Vec<&str> = errs.iter().map(OrificeError::code).collect();
        // All five codes are unique.
        for (i, a) in codes.iter().enumerate() {
            for b in codes.iter().skip(i + 1) {
                assert_ne!(a, b, "duplicate code {a}");
            }
        }
        assert_eq!(codes[0], "orifice.not_finite");
        assert_eq!(codes[3], "orifice.throat_not_smaller");
    }

    #[test]
    fn categories_partition_input_vs_domain() {
        assert_eq!(
            OrificeError::NotFinite { name: "a" }.category(),
            ErrorCategory::Input
        );
        assert_eq!(
            OrificeError::DischargeCoefficientOutOfRange { value: 9.0 }.category(),
            ErrorCategory::Input
        );
        assert_eq!(
            OrificeError::ThroatNotSmaller { d: 2.0, pipe: 1.0 }.category(),
            ErrorCategory::Domain
        );
    }

    #[test]
    fn display_includes_offending_values() {
        let msg = OrificeError::ThroatNotSmaller { d: 0.2, pipe: 0.1 }.to_string();
        assert!(msg.contains("0.2"), "message names d: {msg}");
        assert!(msg.contains("0.1"), "message names D: {msg}");
    }
}
