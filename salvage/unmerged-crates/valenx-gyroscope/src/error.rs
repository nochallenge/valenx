//! Error taxonomy for the gyroscope workbench.
//!
//! Every public computation routes its inputs through the validated
//! constructors here, so a non-finite (`NaN`/`Inf`) or out-of-domain
//! value is rejected up front rather than silently propagating into a
//! `NaN` couple or a divide-by-zero precession rate.

use thiserror::Error;

/// Shorthand for `Result<T, GyroError>`.
pub type Result<T> = core::result::Result<T, GyroError>;

/// Anything that can go wrong validating an input or evaluating a
/// gyroscopic relation.
///
/// This enum is `#[non_exhaustive]`: new variants may be added without it
/// being a breaking change, so downstream `match` arms must include a
/// wildcard.
#[derive(Debug, Error, Clone, PartialEq)]
#[non_exhaustive]
pub enum GyroError {
    /// A quantity that must be a real number was `NaN` or infinite.
    #[error("`{name}` must be finite, got {value}")]
    NonFinite {
        /// Which quantity was non-finite.
        name: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A quantity that physics requires to be strictly positive
    /// (a moment of inertia, a spin rate that sits in a denominator)
    /// was zero or negative.
    #[error("`{name}` must be > 0, got {value}")]
    NonPositive {
        /// Which quantity was out of domain.
        name: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A quantity that physics requires to be non-negative (a couple
    /// magnitude, a precession rate used as a magnitude) was negative.
    #[error("`{name}` must be >= 0, got {value}")]
    Negative {
        /// Which quantity was out of domain.
        name: &'static str,
        /// The offending value.
        value: f64,
    },
}

impl GyroError {
    /// Stable, dot-separated identifier for routing / logging. The string
    /// is part of the crate's public contract and will not change for an
    /// existing variant.
    pub fn code(&self) -> &'static str {
        match self {
            GyroError::NonFinite { .. } => "gyroscope.non_finite",
            GyroError::NonPositive { .. } => "gyroscope.non_positive",
            GyroError::Negative { .. } => "gyroscope.negative",
        }
    }
}

/// Reject a non-finite value for `name`.
#[inline]
pub(crate) fn finite(name: &'static str, value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(GyroError::NonFinite { name, value })
    }
}

/// Require `name` to be finite and strictly positive.
#[inline]
pub(crate) fn positive(name: &'static str, value: f64) -> Result<f64> {
    let value = finite(name, value)?;
    if value > 0.0 {
        Ok(value)
    } else {
        Err(GyroError::NonPositive { name, value })
    }
}

/// Require `name` to be finite and non-negative.
#[inline]
pub(crate) fn non_negative(name: &'static str, value: f64) -> Result<f64> {
    let value = finite(name, value)?;
    if value >= 0.0 {
        Ok(value)
    } else {
        Err(GyroError::Negative { name, value })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finite_accepts_real_numbers_and_rejects_nan_inf() {
        assert_eq!(finite("x", 3.5).unwrap(), 3.5);
        assert_eq!(finite("x", 0.0).unwrap(), 0.0);
        assert_eq!(finite("x", -7.0).unwrap(), -7.0);
        assert!(matches!(
            finite("x", f64::NAN),
            Err(GyroError::NonFinite { .. })
        ));
        assert!(matches!(
            finite("x", f64::INFINITY),
            Err(GyroError::NonFinite { .. })
        ));
        assert!(matches!(
            finite("x", f64::NEG_INFINITY),
            Err(GyroError::NonFinite { .. })
        ));
    }

    #[test]
    fn positive_rejects_zero_and_negative() {
        assert_eq!(positive("x", 1.0).unwrap(), 1.0);
        assert!(matches!(
            positive("x", 0.0),
            Err(GyroError::NonPositive { .. })
        ));
        assert!(matches!(
            positive("x", -1.0),
            Err(GyroError::NonPositive { .. })
        ));
    }

    #[test]
    fn non_negative_accepts_zero_rejects_negative() {
        assert_eq!(non_negative("x", 0.0).unwrap(), 0.0);
        assert_eq!(non_negative("x", 2.0).unwrap(), 2.0);
        assert!(matches!(
            non_negative("x", -0.001),
            Err(GyroError::Negative { .. })
        ));
    }

    #[test]
    fn codes_are_stable_and_distinct() {
        let nf = GyroError::NonFinite {
            name: "x",
            value: f64::NAN,
        };
        let np = GyroError::NonPositive {
            name: "x",
            value: 0.0,
        };
        let ng = GyroError::Negative {
            name: "x",
            value: -1.0,
        };
        assert_eq!(nf.code(), "gyroscope.non_finite");
        assert_eq!(np.code(), "gyroscope.non_positive");
        assert_eq!(ng.code(), "gyroscope.negative");
        // Distinct codes per variant.
        assert_ne!(nf.code(), np.code());
        assert_ne!(np.code(), ng.code());
    }
}
