//! Shaft-design error taxonomy.
//!
//! Every fallible entry point returns [`ShaftError`]. Validated
//! constructors (e.g. [`crate::ShaftSection::new`]) reject physically
//! meaningless inputs up front so the stress formulae only ever run on
//! a positive diameter and finite loads.

use thiserror::Error;

/// Errors raised when building or analysing a shaft section.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum ShaftError {
    /// A diameter, length, or other geometric quantity that must be
    /// strictly positive was zero or negative.
    #[error("non-positive {name}: {value} (must be > 0)")]
    NonPositive {
        /// Name of the offending quantity, e.g. `"diameter_m"`.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A supplied value was not finite (NaN or infinite).
    #[error("non-finite {name}: {value}")]
    NotFinite {
        /// Name of the offending quantity.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },
}

impl ShaftError {
    /// Reject a quantity that must be strictly positive and finite.
    ///
    /// Returns the value unchanged on success so it can be used inline.
    ///
    /// # Errors
    ///
    /// [`ShaftError::NotFinite`] if `value` is NaN or infinite, or
    /// [`ShaftError::NonPositive`] if `value <= 0`.
    pub fn require_positive(name: &'static str, value: f64) -> Result<f64, ShaftError> {
        if !value.is_finite() {
            return Err(ShaftError::NotFinite { name, value });
        }
        if value <= 0.0 {
            return Err(ShaftError::NonPositive { name, value });
        }
        Ok(value)
    }

    /// Reject a quantity that must be finite but may take any sign
    /// (e.g. a bending moment or torque, which can act either way).
    ///
    /// Returns the value unchanged on success.
    ///
    /// # Errors
    ///
    /// [`ShaftError::NotFinite`] if `value` is NaN or infinite.
    pub fn require_finite(name: &'static str, value: f64) -> Result<f64, ShaftError> {
        if !value.is_finite() {
            return Err(ShaftError::NotFinite { name, value });
        }
        Ok(value)
    }

    /// Stable kebab-cased identifier for telemetry / logs.
    pub fn code(&self) -> &'static str {
        match self {
            ShaftError::NonPositive { .. } => "shaftdesign.non_positive",
            ShaftError::NotFinite { .. } => "shaftdesign.not_finite",
        }
    }

    /// Coarse error category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            ShaftError::NonPositive { .. } | ShaftError::NotFinite { .. } => ErrorCategory::Input,
        }
    }
}

/// Coarse classification of a [`ShaftError`].
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// Caller-supplied input was invalid.
    Input,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn require_positive_passes_finite_positive() {
        assert_eq!(ShaftError::require_positive("d", 3.5).unwrap(), 3.5);
    }

    #[test]
    fn require_positive_rejects_zero_and_negative() {
        assert!(matches!(
            ShaftError::require_positive("d", 0.0),
            Err(ShaftError::NonPositive { .. })
        ));
        assert!(matches!(
            ShaftError::require_positive("d", -1.0),
            Err(ShaftError::NonPositive { .. })
        ));
    }

    #[test]
    fn require_positive_rejects_non_finite_before_sign() {
        assert!(matches!(
            ShaftError::require_positive("d", f64::NAN),
            Err(ShaftError::NotFinite { .. })
        ));
        assert!(matches!(
            ShaftError::require_positive("d", f64::INFINITY),
            Err(ShaftError::NotFinite { .. })
        ));
    }

    #[test]
    fn require_finite_allows_negative_and_zero() {
        assert_eq!(ShaftError::require_finite("m", -42.0).unwrap(), -42.0);
        assert_eq!(ShaftError::require_finite("m", 0.0).unwrap(), 0.0);
    }

    #[test]
    fn require_finite_rejects_nan_and_inf() {
        assert!(ShaftError::require_finite("m", f64::NAN).is_err());
        assert!(ShaftError::require_finite("m", f64::NEG_INFINITY).is_err());
    }

    #[test]
    fn codes_and_categories_are_stable() {
        let np = ShaftError::NonPositive {
            name: "d",
            value: -1.0,
        };
        let nf = ShaftError::NotFinite {
            name: "m",
            value: f64::NAN,
        };
        assert_eq!(np.code(), "shaftdesign.non_positive");
        assert_eq!(nf.code(), "shaftdesign.not_finite");
        assert_eq!(np.category(), ErrorCategory::Input);
        assert_eq!(nf.category(), ErrorCategory::Input);
    }

    #[test]
    fn display_renders_name_and_value() {
        let e = ShaftError::NonPositive {
            name: "diameter",
            value: -2.0,
        };
        let msg = e.to_string();
        assert!(msg.contains("diameter"));
        assert!(msg.contains("-2"));
    }
}
