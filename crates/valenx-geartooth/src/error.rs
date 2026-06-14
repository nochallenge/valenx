//! Gear-tooth strength error taxonomy.
//!
//! A single [`GearToothError`] enum covers every failure mode of the
//! crate. Constructors elsewhere validate their inputs and return
//! [`GearToothError::BadParameter`] with the offending parameter name
//! and a human-readable reason; queries that fall outside a model's
//! valid domain return [`GearToothError::OutOfDomain`].

use thiserror::Error;

/// Errors raised by gear-tooth strength calculations.
#[derive(Debug, Error)]
pub enum GearToothError {
    /// A supplied parameter was non-finite, non-positive, or otherwise
    /// outside the physically meaningful range for its quantity.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Name of the offending parameter (stable, `snake_case`).
        name: &'static str,
        /// Human-readable explanation of why it was rejected.
        reason: String,
    },

    /// A query was made outside the model's valid domain — for example
    /// a Lewis form factor for a tooth count below the practical
    /// undercut-free minimum.
    #[error("out of domain: {0}")]
    OutOfDomain(String),
}

impl GearToothError {
    /// Build a [`GearToothError::BadParameter`] from a static name and an
    /// owned reason string.
    ///
    /// Small helper so call sites read as one line.
    pub fn bad_parameter(name: &'static str, reason: impl Into<String>) -> Self {
        Self::BadParameter {
            name,
            reason: reason.into(),
        }
    }

    /// Stable kebab-cased identifier for this error, suitable for logs
    /// and machine matching. Never changes for a given variant.
    pub fn code(&self) -> &'static str {
        match self {
            Self::BadParameter { .. } => "geartooth.bad_parameter",
            Self::OutOfDomain(_) => "geartooth.out_of_domain",
        }
    }

    /// Coarse [`ErrorCategory`] for grouping in a UI or report.
    pub fn category(&self) -> ErrorCategory {
        match self {
            Self::BadParameter { .. } => ErrorCategory::Input,
            Self::OutOfDomain(_) => ErrorCategory::Domain,
        }
    }
}

/// Coarse classification of a [`GearToothError`].
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The caller supplied an invalid input value.
    Input,
    /// The requested operation lies outside the model's valid domain.
    Domain,
}

/// Validate that `value` is finite and strictly greater than zero,
/// returning it unchanged on success.
///
/// Used by the public constructors so the rejection message is uniform.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64, GearToothError> {
    if !value.is_finite() {
        return Err(GearToothError::bad_parameter(
            name,
            format!("must be finite, got {value}"),
        ));
    }
    if value <= 0.0 {
        return Err(GearToothError::bad_parameter(
            name,
            format!("must be strictly positive, got {value}"),
        ));
    }
    Ok(value)
}

/// Validate that `value` is finite and lies within the inclusive range
/// `[lo, hi]`, returning it unchanged on success.
pub(crate) fn require_in_range(
    name: &'static str,
    value: f64,
    lo: f64,
    hi: f64,
) -> Result<f64, GearToothError> {
    if !value.is_finite() {
        return Err(GearToothError::bad_parameter(
            name,
            format!("must be finite, got {value}"),
        ));
    }
    if value < lo || value > hi {
        return Err(GearToothError::bad_parameter(
            name,
            format!("must lie in [{lo}, {hi}], got {value}"),
        ));
    }
    Ok(value)
}
