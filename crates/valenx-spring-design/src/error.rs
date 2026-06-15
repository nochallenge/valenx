//! Error taxonomy for helical-spring design calculations.
//!
//! Every fallible constructor in this crate returns
//! [`SpringError`]. The variants distinguish *input* faults (a
//! non-finite or non-positive dimension) from *domain* faults (a
//! physically meaningless geometry, e.g. a wire thicker than the coil
//! it is wound into).

use thiserror::Error;

/// Errors raised when validating spring inputs or evaluating models.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum SpringError {
    /// A supplied parameter was not a finite, strictly-positive number.
    ///
    /// All physical spring dimensions (wire diameter, coil diameter,
    /// coil count, modulus, force) must be finite and `> 0`. `NaN`,
    /// infinities, zero, and negatives are rejected here.
    #[error("non-positive or non-finite parameter `{name}` = {value}: {reason}")]
    NonPositive {
        /// Name of the offending parameter, e.g. `"wire_diameter_mm"`.
        name: &'static str,
        /// The numeric value that failed validation.
        value: f64,
        /// Human-readable explanation of the requirement.
        reason: &'static str,
    },

    /// The geometry is physically impossible.
    ///
    /// Raised, for example, when the mean coil diameter is not strictly
    /// greater than the wire diameter (the spring index `C = D/d` would
    /// be `<= 1`, meaning the wire could not close into a coil).
    #[error("degenerate spring geometry: {0}")]
    Degenerate(String),
}

/// Coarse classification of a [`SpringError`], handy for UI grouping
/// and for deciding whether the fault is the user's input or an
/// algorithmic-domain limit.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The caller supplied a bad value (out of range / non-finite).
    Input,
    /// The combination of (individually valid) values is geometrically
    /// inconsistent.
    Geometry,
}

impl SpringError {
    /// Stable, kebab-cased identifier suitable for logs and tests.
    ///
    /// The string is part of the crate's public contract and will not
    /// change for a given variant across patch releases.
    pub fn code(&self) -> &'static str {
        match self {
            SpringError::NonPositive { .. } => "spring.non-positive",
            SpringError::Degenerate(_) => "spring.degenerate",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            SpringError::NonPositive { .. } => ErrorCategory::Input,
            SpringError::Degenerate(_) => ErrorCategory::Geometry,
        }
    }
}

/// Validate that `value` is finite and strictly positive.
///
/// Returns `value` unchanged on success, or a
/// [`SpringError::NonPositive`] carrying `name` and `reason` on
/// failure. Used by every public constructor in this crate so the
/// rejection message is uniform.
///
/// # Errors
///
/// Returns [`SpringError::NonPositive`] if `value` is `NaN`, infinite,
/// zero, or negative.
pub(crate) fn require_positive(
    value: f64,
    name: &'static str,
    reason: &'static str,
) -> Result<f64, SpringError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(SpringError::NonPositive {
            name,
            value,
            reason,
        })
    }
}
