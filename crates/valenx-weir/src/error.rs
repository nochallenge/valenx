//! Error taxonomy for weir-flow hydraulics.
//!
//! Every fallible constructor in this crate validates its inputs and
//! returns a [`WeirError`] describing exactly which physical quantity
//! was out of range, rather than panicking or silently producing a
//! nonsensical discharge.

use thiserror::Error;

/// Errors raised while constructing or evaluating a weir.
///
/// All variants are produced by the validated constructors in the
/// [`crate::rectangular`] and [`crate::vnotch`] modules. They carry the
/// offending quantity by name together with the rejected value so a
/// caller can surface a precise diagnostic.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum WeirError {
    /// A quantity that must be strictly positive was zero or negative.
    ///
    /// Crest length, head, discharge coefficient and gravitational
    /// acceleration must all be `> 0` for the weir formulae to be
    /// physically meaningful.
    #[error("non-positive `{name}`: got {value}, must be > 0")]
    NonPositive {
        /// Name of the offending quantity (e.g. `"head"`).
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A quantity required to be finite was `NaN` or infinite.
    #[error("non-finite `{name}`: got {value}")]
    NotFinite {
        /// Name of the offending quantity.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A V-notch vertex angle that was outside the open interval
    /// `(0, π)` radians — a triangular notch must subtend a strictly
    /// positive angle that is less than a straight line.
    #[error("V-notch angle out of range: got {radians} rad, must be in (0, π)")]
    NotchAngleOutOfRange {
        /// The rejected full vertex angle, in radians.
        radians: f64,
    },
}

/// Coarse category for a [`WeirError`], useful for routing or metrics.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A user-supplied input value was invalid.
    Input,
    /// A configuration / tuning constant was invalid.
    Config,
}

impl WeirError {
    /// Stable, kebab-cased identifier for this error.
    ///
    /// The string is part of the crate's public contract: callers may
    /// match on it for logging or localization without depending on the
    /// human-readable [`Display`](std::fmt::Display) text.
    pub fn code(&self) -> &'static str {
        match self {
            WeirError::NonPositive { .. } => "weir.non-positive",
            WeirError::NotFinite { .. } => "weir.not-finite",
            WeirError::NotchAngleOutOfRange { .. } => "weir.notch-angle-out-of-range",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    ///
    /// The discharge coefficient and gravitational acceleration are
    /// treated as configuration; geometry and head are user input.
    pub fn category(&self) -> ErrorCategory {
        match self {
            WeirError::NonPositive { name, .. } => match *name {
                "discharge_coefficient" | "gravity" => ErrorCategory::Config,
                _ => ErrorCategory::Input,
            },
            WeirError::NotFinite { name, .. } => match *name {
                "discharge_coefficient" | "gravity" => ErrorCategory::Config,
                _ => ErrorCategory::Input,
            },
            WeirError::NotchAngleOutOfRange { .. } => ErrorCategory::Input,
        }
    }
}

/// Validate that `value` is finite and strictly positive.
///
/// Returns `value` unchanged on success, otherwise the appropriate
/// [`WeirError`] tagged with `name`. Used by the public constructors so
/// the validation rules live in exactly one place.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64, WeirError> {
    if !value.is_finite() {
        return Err(WeirError::NotFinite { name, value });
    }
    if value <= 0.0 {
        return Err(WeirError::NonPositive { name, value });
    }
    Ok(value)
}
