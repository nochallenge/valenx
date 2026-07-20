//! Error taxonomy for `valenx-fin`.
//!
//! Every fallible public function returns [`Result<T, FinError>`].
//! The constructors here are the single choke point through which all
//! physical inputs pass: they reject non-finite values (`NaN`, infinities)
//! and values outside the physical domain (the fin model demands strictly
//! positive geometry and transport properties, and a non-negative tip
//! length `m * L`). Centralising the checks means the formula modules can
//! assume their arguments are already sane.

use thiserror::Error;

/// Errors raised by the extended-surface (fin) heat-transfer model.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum FinError {
    /// A supplied value was not finite (`NaN`, `+inf`, or `-inf`).
    #[error("parameter `{name}` must be finite, got {value}")]
    NotFinite {
        /// Name of the offending parameter.
        name: &'static str,
        /// The non-finite value received.
        value: f64,
    },

    /// A value the physics requires to be strictly positive was not.
    ///
    /// Examples: convection coefficient `h`, perimeter `P`, thermal
    /// conductivity `k`, cross-sectional area `A`, and fin length `L`.
    #[error("parameter `{name}` must be strictly positive, got {value}")]
    NonPositive {
        /// Name of the offending parameter.
        name: &'static str,
        /// The non-positive value received.
        value: f64,
    },

    /// A value the physics requires to be non-negative was negative.
    ///
    /// Used for the dimensionless fin length `mL`, which may legitimately
    /// be zero (the conduction-dominated limit) but never negative.
    #[error("parameter `{name}` must be non-negative, got {value}")]
    Negative {
        /// Name of the offending parameter.
        name: &'static str,
        /// The negative value received.
        value: f64,
    },
}

impl FinError {
    /// Stable, kebab-cased identifier for programmatic matching and logs.
    ///
    /// The string is part of the public contract and will not change for a
    /// given variant.
    pub fn code(&self) -> &'static str {
        match self {
            FinError::NotFinite { .. } => "fin.not-finite",
            FinError::NonPositive { .. } => "fin.non-positive",
            FinError::Negative { .. } => "fin.negative",
        }
    }

    /// Validate that `value` is finite, returning it unchanged on success.
    ///
    /// # Errors
    ///
    /// Returns [`FinError::NotFinite`] when `value` is `NaN` or infinite.
    pub fn finite(name: &'static str, value: f64) -> std::result::Result<f64, FinError> {
        if value.is_finite() {
            Ok(value)
        } else {
            Err(FinError::NotFinite { name, value })
        }
    }

    /// Validate that `value` is finite and strictly positive.
    ///
    /// # Errors
    ///
    /// Returns [`FinError::NotFinite`] when `value` is not finite, or
    /// [`FinError::NonPositive`] when `value <= 0`.
    pub fn positive(name: &'static str, value: f64) -> std::result::Result<f64, FinError> {
        let v = FinError::finite(name, value)?;
        if v > 0.0 {
            Ok(v)
        } else {
            Err(FinError::NonPositive { name, value: v })
        }
    }

    /// Validate that `value` is finite and non-negative (`>= 0`).
    ///
    /// # Errors
    ///
    /// Returns [`FinError::NotFinite`] when `value` is not finite, or
    /// [`FinError::Negative`] when `value < 0`.
    pub fn non_negative(name: &'static str, value: f64) -> std::result::Result<f64, FinError> {
        let v = FinError::finite(name, value)?;
        if v >= 0.0 {
            Ok(v)
        } else {
            Err(FinError::Negative { name, value: v })
        }
    }
}

/// Convenience alias for `Result<T, FinError>`.
pub type Result<T> = std::result::Result<T, FinError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finite_accepts_normal_value() {
        assert_eq!(FinError::finite("x", 1.5).unwrap(), 1.5);
    }

    #[test]
    fn finite_rejects_nan_and_inf() {
        assert_eq!(
            FinError::finite("x", f64::NAN).unwrap_err().code(),
            "fin.not-finite"
        );
        assert!(matches!(
            FinError::finite("x", f64::INFINITY),
            Err(FinError::NotFinite { name: "x", .. })
        ));
    }

    #[test]
    fn positive_accepts_positive_and_rejects_zero_and_negative() {
        assert_eq!(FinError::positive("h", 12.0).unwrap(), 12.0);
        assert_eq!(
            FinError::positive("h", 0.0).unwrap_err().code(),
            "fin.non-positive"
        );
        assert!(matches!(
            FinError::positive("h", -1.0),
            Err(FinError::NonPositive { .. })
        ));
    }

    #[test]
    fn non_negative_accepts_zero_but_rejects_negative() {
        assert_eq!(FinError::non_negative("mL", 0.0).unwrap(), 0.0);
        assert!(matches!(
            FinError::non_negative("mL", -0.001),
            Err(FinError::Negative { .. })
        ));
    }

    #[test]
    fn positive_propagates_not_finite() {
        // A non-finite input is reported as NotFinite, not NonPositive.
        assert!(matches!(
            FinError::positive("k", f64::NEG_INFINITY),
            Err(FinError::NotFinite { .. })
        ));
    }
}
