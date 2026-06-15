//! Error taxonomy for the thermal-expansion models.
//!
//! Every fallible entry point in this crate validates its inputs through
//! the [`ThermalError`] constructors below, so a constructed value is
//! always physically admissible (finite, and positive where a physical
//! quantity must be positive). The asserts in the topic modules then
//! exercise the *numerics*, not the guards.

use thiserror::Error;

/// Errors raised when constructing or evaluating a thermal-expansion model.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum ThermalError {
    /// A quantity that must be strictly positive was zero or negative.
    ///
    /// Examples: a reference length, a coefficient of thermal expansion,
    /// or a Young's modulus supplied as `<= 0`.
    #[error("non-positive `{name}`: got {value}, expected a value > 0")]
    NonPositive {
        /// Name of the offending parameter.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A quantity that must be finite was `NaN` or an infinity.
    #[error("non-finite `{name}`: got {value}, expected a finite value")]
    NonFinite {
        /// Name of the offending parameter.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A named material was not present in the built-in library.
    #[error("unknown material `{name}`")]
    UnknownMaterial {
        /// The requested material key.
        name: String,
    },
}

/// Coarse classification of a [`ThermalError`], useful for callers that
/// want to branch on the *kind* of failure without matching every variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// The caller supplied an invalid numeric input.
    Input,
    /// The caller referenced something the library does not contain.
    Lookup,
}

impl ThermalError {
    /// Return a stable, kebab-cased identifier for this error.
    ///
    /// The string is suitable for logs and machine matching and will not
    /// change for a given variant across patch releases.
    pub fn code(&self) -> &'static str {
        match self {
            ThermalError::NonPositive { .. } => "thermalexpansion.non-positive",
            ThermalError::NonFinite { .. } => "thermalexpansion.non-finite",
            ThermalError::UnknownMaterial { .. } => "thermalexpansion.unknown-material",
        }
    }

    /// Return the coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            ThermalError::NonPositive { .. } | ThermalError::NonFinite { .. } => {
                ErrorCategory::Input
            }
            ThermalError::UnknownMaterial { .. } => ErrorCategory::Lookup,
        }
    }
}

/// Validate that `value` is finite, returning [`ThermalError::NonFinite`]
/// (tagged with `name`) otherwise.
///
/// This is the shared guard behind temperature changes, which may legally
/// be negative (cooling) but must always be finite.
pub(crate) fn require_finite(name: &'static str, value: f64) -> Result<f64, ThermalError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ThermalError::NonFinite { name, value })
    }
}

/// Validate that `value` is finite and strictly positive, returning the
/// appropriate [`ThermalError`] (tagged with `name`) otherwise.
///
/// Used for quantities — lengths, areas, volumes, coefficients of thermal
/// expansion, elastic moduli — that have no physical meaning at or below
/// zero in these models.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64, ThermalError> {
    let value = require_finite(name, value)?;
    if value > 0.0 {
        Ok(value)
    } else {
        Err(ThermalError::NonPositive { name, value })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn require_positive_accepts_positive() {
        assert_eq!(require_positive("len", 2.5).unwrap(), 2.5);
    }

    #[test]
    fn require_positive_rejects_zero_and_negative() {
        assert_eq!(
            require_positive("len", 0.0).unwrap_err(),
            ThermalError::NonPositive {
                name: "len",
                value: 0.0
            }
        );
        assert_eq!(
            require_positive("len", -1.0).unwrap_err(),
            ThermalError::NonPositive {
                name: "len",
                value: -1.0
            }
        );
    }

    #[test]
    fn require_finite_allows_negative_but_rejects_nan_and_inf() {
        assert_eq!(require_finite("dT", -40.0).unwrap(), -40.0);
        assert!(matches!(
            require_finite("dT", f64::NAN).unwrap_err(),
            ThermalError::NonFinite { name: "dT", .. }
        ));
        assert!(matches!(
            require_finite("dT", f64::INFINITY).unwrap_err(),
            ThermalError::NonFinite { .. }
        ));
    }

    #[test]
    fn codes_and_categories_are_stable() {
        let np = ThermalError::NonPositive {
            name: "x",
            value: -1.0,
        };
        let nf = ThermalError::NonFinite {
            name: "x",
            value: f64::NAN,
        };
        let um = ThermalError::UnknownMaterial {
            name: "unobtainium".to_string(),
        };
        assert_eq!(np.code(), "thermalexpansion.non-positive");
        assert_eq!(nf.code(), "thermalexpansion.non-finite");
        assert_eq!(um.code(), "thermalexpansion.unknown-material");
        assert_eq!(np.category(), ErrorCategory::Input);
        assert_eq!(nf.category(), ErrorCategory::Input);
        assert_eq!(um.category(), ErrorCategory::Lookup);
    }
}
