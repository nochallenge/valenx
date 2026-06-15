//! Error taxonomy for the soil bearing-capacity crate.
//!
//! Every fallible entry point returns [`SoilBearingError`]. Construct
//! validated inputs through the checking constructors on the input
//! types ([`crate::soil::SoilProperties::new`],
//! [`crate::footing::Footing::new`]) rather than building the structs
//! by hand, so out-of-domain values are rejected at the boundary.

use thiserror::Error;

/// Errors raised while validating geotechnical inputs or evaluating the
/// Terzaghi bearing-capacity equation.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum SoilBearingError {
    /// A scalar parameter was outside its admissible domain.
    ///
    /// `name` is a stable, machine-readable identifier for the
    /// offending field; `value` is the rejected number; `reason`
    /// explains the constraint that was violated.
    #[error("parameter `{name}` = {value} is invalid: {reason}")]
    InvalidParameter {
        /// Stable identifier of the offending parameter (e.g. `"phi_deg"`).
        name: &'static str,
        /// The rejected value.
        value: f64,
        /// Human-readable description of the violated constraint.
        reason: &'static str,
    },

    /// A non-finite value (NaN or infinity) was supplied where a finite
    /// number is required.
    #[error("parameter `{name}` must be finite, got {value}")]
    NotFinite {
        /// Stable identifier of the offending parameter.
        name: &'static str,
        /// The rejected (non-finite) value.
        value: f64,
    },
}

impl SoilBearingError {
    /// Stable, kebab/dotted identifier for the error variant.
    ///
    /// Useful for logging, metrics, or matching in callers without
    /// depending on the human-readable [`std::fmt::Display`] text.
    pub fn code(&self) -> &'static str {
        match self {
            SoilBearingError::InvalidParameter { .. } => "soilbearing.invalid_parameter",
            SoilBearingError::NotFinite { .. } => "soilbearing.not_finite",
        }
    }

    /// Construct an [`SoilBearingError::InvalidParameter`].
    ///
    /// Internal helper used by the validated constructors.
    pub(crate) fn invalid(name: &'static str, value: f64, reason: &'static str) -> Self {
        SoilBearingError::InvalidParameter {
            name,
            value,
            reason,
        }
    }

    /// Reject `value` if it is not finite, returning
    /// [`SoilBearingError::NotFinite`].
    ///
    /// Internal helper used by the validated constructors.
    pub(crate) fn require_finite(name: &'static str, value: f64) -> Result<f64, Self> {
        if value.is_finite() {
            Ok(value)
        } else {
            Err(SoilBearingError::NotFinite { name, value })
        }
    }
}
