//! Thermocouple error taxonomy.

use thiserror::Error;

/// Errors raised when constructing or evaluating a [`Thermocouple`].
///
/// [`Thermocouple`]: crate::thermocouple::Thermocouple
#[derive(Debug, Error)]
pub enum ThermocoupleError {
    /// A parameter was non-finite (`NaN` or infinite).
    #[error("parameter `{name}` must be finite, got {value}")]
    NonFinite {
        /// Parameter name.
        name: &'static str,
        /// Offending value.
        value: f64,
    },

    /// The Seebeck sensitivity was zero or negative.
    ///
    /// A non-positive sensitivity makes the EMF map degenerate (and
    /// non-invertible at zero), so it is rejected at construction time.
    #[error("Seebeck sensitivity must be strictly positive, got {0} V/C")]
    NonPositiveSensitivity(f64),
}

impl ThermocoupleError {
    /// Stable kebab-cased identifier, handy for logging and tests.
    ///
    /// The string is part of the crate's public contract and will not
    /// change for a given variant.
    pub fn code(&self) -> &'static str {
        match self {
            ThermocoupleError::NonFinite { .. } => "thermocouple.non_finite",
            ThermocoupleError::NonPositiveSensitivity(_) => "thermocouple.non_positive_sensitivity",
        }
    }

    /// Coarse category for grouping errors in a UI or report.
    pub fn category(&self) -> ErrorCategory {
        match self {
            ThermocoupleError::NonFinite { .. } => ErrorCategory::Input,
            ThermocoupleError::NonPositiveSensitivity(_) => ErrorCategory::Config,
        }
    }
}

/// Coarse error category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// Bad caller-supplied measurement input (non-finite temperature or
    /// voltage).
    Input,
    /// Bad device configuration (an invalid Seebeck sensitivity).
    Config,
}
