//! Error type for the analog-filter models.

use thiserror::Error;

/// Shorthand for `Result<T, FilterError>`.
pub type Result<T> = core::result::Result<T, FilterError>;

/// Anything that can go wrong building or evaluating a filter model.
///
/// Every public constructor in this crate validates its inputs and
/// returns one of these variants rather than producing a silent
/// `NaN`/`Inf` (e.g. from a zero `R`, `L`, or `C` feeding a division or
/// a `sqrt`).
///
/// This enum is `#[non_exhaustive]`: new variants may be added in future
/// releases without it being a breaking change, so downstream `match`
/// arms must include a wildcard.
#[derive(Debug, Error, Clone, PartialEq)]
#[non_exhaustive]
pub enum FilterError {
    /// A component value (resistance, inductance, or capacitance) was
    /// non-physical — that is, not a strictly-positive finite number.
    ///
    /// Zero or negative values are rejected because they would drive a
    /// division by zero (the cutoff / resonance formulas all divide by
    /// `R`, `L`, or `C`) or take the square root of a non-positive
    /// number; `NaN` / `Inf` are rejected for the same reason.
    #[error("invalid component: {field} = {value} (must be a finite value > 0)")]
    InvalidComponent {
        /// Which component was bad (`"R"`, `"L"`, or `"C"`).
        field: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A frequency argument was non-physical — negative or non-finite.
    ///
    /// Zero is permitted (the DC point), but a negative or `NaN` / `Inf`
    /// frequency is rejected.
    #[error("invalid frequency: {value} Hz (must be a finite value >= 0)")]
    InvalidFrequency {
        /// The offending value, in hertz.
        value: f64,
    },
}

impl FilterError {
    /// Construct an [`FilterError::InvalidComponent`] for `field`.
    #[must_use]
    pub(crate) fn component(field: &'static str, value: f64) -> Self {
        Self::InvalidComponent { field, value }
    }

    /// Construct an [`FilterError::InvalidFrequency`].
    #[must_use]
    pub(crate) fn frequency(value: f64) -> Self {
        Self::InvalidFrequency { value }
    }
}

/// Validate that `value` is a strictly-positive finite component value,
/// returning it on success or an [`FilterError::InvalidComponent`].
pub(crate) fn check_component(field: &'static str, value: f64) -> Result<f64> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(FilterError::component(field, value))
    }
}

/// Validate that `value` is a finite, non-negative frequency in hertz,
/// returning it on success or an [`FilterError::InvalidFrequency`].
pub(crate) fn check_frequency(value: f64) -> Result<f64> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(FilterError::frequency(value))
    }
}
