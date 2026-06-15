//! Flywheel sizing error taxonomy.
//!
//! Every fallible constructor and computation in this crate returns a
//! [`FlywheelError`]. The variants carry enough context (the offending
//! parameter name and a human-readable reason) to surface a useful
//! message, and [`FlywheelError::code`] gives a stable kebab-cased
//! identifier suitable for logs or UI lookups.

use thiserror::Error;

/// Errors raised while constructing or evaluating flywheel models.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum FlywheelError {
    /// A scalar input was outside its allowed range (for example a
    /// non-positive mass, radius, density, or angular speed where a
    /// strictly positive value is required).
    #[error("invalid parameter `{name}`: {reason} (got {value})")]
    InvalidParameter {
        /// The parameter that failed validation.
        name: &'static str,
        /// Why the value was rejected.
        reason: &'static str,
        /// The offending value, for diagnostics.
        value: f64,
    },

    /// A pair of inputs is individually fine but jointly inconsistent —
    /// for example an annulus whose inner radius is not strictly less
    /// than its outer radius, or a speed band whose minimum exceeds its
    /// maximum.
    #[error("inconsistent parameters: {0}")]
    Inconsistent(&'static str),
}

impl FlywheelError {
    /// Build an [`FlywheelError::InvalidParameter`] error.
    ///
    /// Centralising construction keeps the validation call sites terse
    /// and the message format consistent across modules.
    #[must_use]
    pub fn invalid(name: &'static str, reason: &'static str, value: f64) -> Self {
        Self::InvalidParameter {
            name,
            reason,
            value,
        }
    }

    /// Validate that `value` is finite and strictly positive, returning
    /// it on success or an [`FlywheelError::InvalidParameter`] otherwise.
    ///
    /// This is the workhorse guard used by the rotor / flywheel
    /// constructors for masses, radii, densities, and speeds.
    pub fn require_positive(name: &'static str, value: f64) -> Result<f64, Self> {
        if !value.is_finite() {
            return Err(Self::invalid(name, "must be finite", value));
        }
        if value <= 0.0 {
            return Err(Self::invalid(name, "must be > 0", value));
        }
        Ok(value)
    }

    /// Validate that `value` is finite and non-negative (zero allowed),
    /// returning it on success.
    ///
    /// Used for quantities such as an angular speed that may legitimately
    /// be zero (a flywheel at rest stores no energy but is still valid).
    pub fn require_non_negative(name: &'static str, value: f64) -> Result<f64, Self> {
        if !value.is_finite() {
            return Err(Self::invalid(name, "must be finite", value));
        }
        if value < 0.0 {
            return Err(Self::invalid(name, "must be >= 0", value));
        }
        Ok(value)
    }

    /// A stable kebab-cased identifier for the error variant, suitable
    /// for log filtering or localisation keys.
    #[must_use]
    pub fn code(&self) -> &'static str {
        match self {
            Self::InvalidParameter { .. } => "flywheel.invalid_parameter",
            Self::Inconsistent(_) => "flywheel.inconsistent",
        }
    }

    /// The coarse category this error belongs to.
    #[must_use]
    pub fn category(&self) -> ErrorCategory {
        match self {
            Self::InvalidParameter { .. } => ErrorCategory::Input,
            Self::Inconsistent(_) => ErrorCategory::Input,
        }
    }
}

/// A coarse grouping of [`FlywheelError`] variants for triage.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The caller supplied a bad input value or an inconsistent
    /// combination of inputs.
    Input,
    /// An internal algorithmic precondition was violated.
    Algorithm,
}
