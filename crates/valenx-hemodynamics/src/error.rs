//! Error taxonomy for the hemodynamics models.
//!
//! Every fallible computation in this crate validates its inputs up
//! front and returns a [`HemodynamicsError`]. The physiological flow
//! formulas are only meaningful for strictly-positive geometry and
//! material properties (a radius, a viscosity, a length), so a
//! non-positive input is reported as [`HemodynamicsError::NonPositive`]
//! rather than silently producing an infinity or a NaN. Inputs that may
//! legitimately be zero but never negative (an elapsed time, a heart
//! rate) are reported as [`HemodynamicsError::Negative`].

use thiserror::Error;

/// Errors raised by the hemodynamics calculations.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum HemodynamicsError {
    /// A parameter that must be strictly positive (`> 0`) was not.
    ///
    /// Used for quantities that appear in a denominator or as a power
    /// base where zero is physically meaningless: a vessel radius, a
    /// dynamic viscosity, a vessel length, or a time constant.
    #[error("parameter `{name}` must be > 0, got {value}")]
    NonPositive {
        /// Name of the offending parameter.
        name: &'static str,
        /// The value that was supplied.
        value: f64,
    },

    /// A parameter that must be non-negative (`>= 0`) was negative.
    ///
    /// Used for quantities where zero is a valid (if degenerate) value
    /// but a negative number is not: an elapsed time, a heart rate, a
    /// stroke volume, or a flow rate.
    #[error("parameter `{name}` must be >= 0, got {value}")]
    Negative {
        /// Name of the offending parameter.
        name: &'static str,
        /// The value that was supplied.
        value: f64,
    },

    /// A parameter was not a finite number (it was `NaN` or infinite).
    #[error("parameter `{name}` must be finite, got {value}")]
    NotFinite {
        /// Name of the offending parameter.
        name: &'static str,
        /// The value that was supplied.
        value: f64,
    },
}

/// Coarse classification of a [`HemodynamicsError`], handy for
/// telemetry and for callers that want to branch on error class
/// without matching every variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The supplied input was outside its valid domain.
    Domain,
}

impl HemodynamicsError {
    /// Stable, kebab-cased identifier for the error variant.
    ///
    /// The string is part of the public contract and is intended for
    /// logs and dashboards; it will not change for a given variant.
    pub fn code(&self) -> &'static str {
        match self {
            HemodynamicsError::NonPositive { .. } => "hemodynamics.non-positive",
            HemodynamicsError::Negative { .. } => "hemodynamics.negative",
            HemodynamicsError::NotFinite { .. } => "hemodynamics.not-finite",
        }
    }

    /// Coarse [`ErrorCategory`] for the variant.
    pub fn category(&self) -> ErrorCategory {
        match self {
            HemodynamicsError::NonPositive { .. }
            | HemodynamicsError::Negative { .. }
            | HemodynamicsError::NotFinite { .. } => ErrorCategory::Domain,
        }
    }
}

/// Validate that `value` is finite and strictly positive (`> 0`).
///
/// Returns the value unchanged on success so it can be used inline.
///
/// # Errors
///
/// Returns [`HemodynamicsError::NotFinite`] if `value` is `NaN` or
/// infinite, or [`HemodynamicsError::NonPositive`] if it is `<= 0`.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64, HemodynamicsError> {
    if !value.is_finite() {
        return Err(HemodynamicsError::NotFinite { name, value });
    }
    if value <= 0.0 {
        return Err(HemodynamicsError::NonPositive { name, value });
    }
    Ok(value)
}

/// Validate that `value` is finite and non-negative (`>= 0`).
///
/// Returns the value unchanged on success so it can be used inline.
///
/// # Errors
///
/// Returns [`HemodynamicsError::NotFinite`] if `value` is `NaN` or
/// infinite, or [`HemodynamicsError::Negative`] if it is `< 0`.
pub(crate) fn require_non_negative(
    name: &'static str,
    value: f64,
) -> Result<f64, HemodynamicsError> {
    if !value.is_finite() {
        return Err(HemodynamicsError::NotFinite { name, value });
    }
    if value < 0.0 {
        return Err(HemodynamicsError::Negative { name, value });
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn require_positive_accepts_positive() {
        let v = require_positive("x", 2.5).expect("2.5 is positive");
        assert!((v - 2.5).abs() < 1e-12);
    }

    #[test]
    fn require_positive_rejects_zero_and_negative() {
        assert_eq!(
            require_positive("r", 0.0),
            Err(HemodynamicsError::NonPositive {
                name: "r",
                value: 0.0
            })
        );
        assert_eq!(
            require_positive("r", -1.0),
            Err(HemodynamicsError::NonPositive {
                name: "r",
                value: -1.0
            })
        );
    }

    #[test]
    fn require_positive_rejects_non_finite() {
        // NaN != NaN, so the error value can only be matched structurally,
        // never compared for equality.
        match require_positive("r", f64::NAN) {
            Err(HemodynamicsError::NotFinite { name, value }) => {
                assert_eq!(name, "r");
                assert!(value.is_nan());
            }
            other => panic!("expected NotFinite(NaN), got {other:?}"),
        }
        assert!(matches!(
            require_positive("r", f64::INFINITY),
            Err(HemodynamicsError::NotFinite { .. })
        ));
    }

    #[test]
    fn require_non_negative_accepts_zero() {
        let v = require_non_negative("t", 0.0).expect("0 is non-negative");
        assert!(v.abs() < 1e-12);
    }

    #[test]
    fn require_non_negative_rejects_negative() {
        assert_eq!(
            require_non_negative("t", -0.5),
            Err(HemodynamicsError::Negative {
                name: "t",
                value: -0.5
            })
        );
    }

    #[test]
    fn code_and_category_are_stable() {
        let e = HemodynamicsError::NonPositive {
            name: "r",
            value: 0.0,
        };
        assert_eq!(e.code(), "hemodynamics.non-positive");
        assert_eq!(e.category(), ErrorCategory::Domain);

        let e = HemodynamicsError::Negative {
            name: "t",
            value: -1.0,
        };
        assert_eq!(e.code(), "hemodynamics.negative");

        let e = HemodynamicsError::NotFinite {
            name: "r",
            value: f64::NAN,
        };
        assert_eq!(e.code(), "hemodynamics.not-finite");
    }

    #[test]
    fn display_is_readable() {
        let e = HemodynamicsError::NonPositive {
            name: "radius_m",
            value: -1.0,
        };
        let msg = e.to_string();
        assert!(msg.contains("radius_m"));
        assert!(msg.contains("> 0"));
    }
}
