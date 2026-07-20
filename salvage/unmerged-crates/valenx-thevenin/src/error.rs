//! Error taxonomy for the Thevenin / Norton equivalent models.
//!
//! Every public constructor in this crate funnels invalid input through
//! [`TheveninError`]. The validation rules encode the physics domain:
//! resistances must be strictly positive (a real two-terminal source has
//! a finite, non-zero internal resistance for the Norton / max-power
//! transforms to be defined), and every numeric input must be finite
//! (no `NaN`, no infinities).

use thiserror::Error;

/// Errors raised when constructing or transforming an equivalent network.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum TheveninError {
    /// A numeric input was `NaN` or infinite.
    #[error("non-finite value for `{name}`: {value}")]
    NonFinite {
        /// The offending parameter name.
        name: &'static str,
        /// The offending value (already known to be `NaN` or infinite).
        value: f64,
    },

    /// A resistance that physics requires to be strictly positive was
    /// zero or negative.
    #[error("non-positive resistance `{name}` = {value} ohm (must be > 0)")]
    NonPositiveResistance {
        /// The offending parameter name.
        name: &'static str,
        /// The offending value (zero or negative).
        value: f64,
    },

    /// A quantity that must not be negative (e.g. a load resistance, which
    /// may legitimately be zero for a short circuit) was negative.
    #[error("negative value `{name}` = {value} (must be >= 0)")]
    Negative {
        /// The offending parameter name.
        name: &'static str,
        /// The offending value (negative).
        value: f64,
    },
}

impl TheveninError {
    /// Stable kebab-cased identifier, suitable for logs and machine
    /// dispatch. The string is part of the crate's public contract and
    /// will not change for a given variant.
    pub fn code(&self) -> &'static str {
        match self {
            TheveninError::NonFinite { .. } => "thevenin.non-finite",
            TheveninError::NonPositiveResistance { .. } => "thevenin.non-positive-resistance",
            TheveninError::Negative { .. } => "thevenin.negative",
        }
    }
}

/// Reject `NaN` / infinite inputs. Returns the value unchanged when finite.
pub(crate) fn finite(name: &'static str, value: f64) -> Result<f64, TheveninError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(TheveninError::NonFinite { name, value })
    }
}

/// Require a finite, strictly positive resistance.
pub(crate) fn positive_resistance(name: &'static str, value: f64) -> Result<f64, TheveninError> {
    let value = finite(name, value)?;
    if value > 0.0 {
        Ok(value)
    } else {
        Err(TheveninError::NonPositiveResistance { name, value })
    }
}

/// Require a finite, non-negative value (zero allowed).
pub(crate) fn non_negative(name: &'static str, value: f64) -> Result<f64, TheveninError> {
    let value = finite(name, value)?;
    if value >= 0.0 {
        Ok(value)
    } else {
        Err(TheveninError::Negative { name, value })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finite_accepts_normal_and_rejects_nan_inf() {
        assert_eq!(finite("x", 3.5).unwrap(), 3.5);
        assert_eq!(finite("x", 0.0).unwrap(), 0.0);
        assert_eq!(finite("x", -2.0).unwrap(), -2.0);
        assert!(matches!(
            finite("x", f64::NAN),
            Err(TheveninError::NonFinite { .. })
        ));
        assert!(matches!(
            finite("x", f64::INFINITY),
            Err(TheveninError::NonFinite { .. })
        ));
        assert!(matches!(
            finite("x", f64::NEG_INFINITY),
            Err(TheveninError::NonFinite { .. })
        ));
    }

    #[test]
    fn positive_resistance_rejects_zero_and_negative() {
        assert_eq!(positive_resistance("r", 10.0).unwrap(), 10.0);
        assert!(matches!(
            positive_resistance("r", 0.0),
            Err(TheveninError::NonPositiveResistance { .. })
        ));
        assert!(matches!(
            positive_resistance("r", -1.0),
            Err(TheveninError::NonPositiveResistance { .. })
        ));
        // Non-finite is caught first, as a NonFinite error.
        assert!(matches!(
            positive_resistance("r", f64::NAN),
            Err(TheveninError::NonFinite { .. })
        ));
    }

    #[test]
    fn non_negative_allows_zero_rejects_negative() {
        assert_eq!(non_negative("r", 0.0).unwrap(), 0.0);
        assert_eq!(non_negative("r", 4.0).unwrap(), 4.0);
        assert!(matches!(
            non_negative("r", -0.001),
            Err(TheveninError::Negative { .. })
        ));
    }

    #[test]
    fn codes_are_stable_and_distinct() {
        let a = TheveninError::NonFinite {
            name: "x",
            value: f64::NAN,
        };
        let b = TheveninError::NonPositiveResistance {
            name: "r",
            value: 0.0,
        };
        let c = TheveninError::Negative {
            name: "r",
            value: -1.0,
        };
        assert_eq!(a.code(), "thevenin.non-finite");
        assert_eq!(b.code(), "thevenin.non-positive-resistance");
        assert_eq!(c.code(), "thevenin.negative");
        // All distinct.
        assert_ne!(a.code(), b.code());
        assert_ne!(b.code(), c.code());
        assert_ne!(a.code(), c.code());
    }
}
