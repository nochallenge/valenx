//! Capacitor calculator error taxonomy.
//!
//! Every fallible function in this crate returns
//! [`Result<_, CapacitorError>`](Result). The error type carries stable
//! [`code`](CapacitorError::code) and [`category`](CapacitorError::category)
//! accessors so callers can branch on failures without string-matching the
//! human-readable [`Display`](core::fmt::Display) message.

use thiserror::Error;

/// Convenience alias for results produced by this crate.
pub type Result<T> = core::result::Result<T, CapacitorError>;

/// Errors raised while validating inputs or evaluating capacitor models.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum CapacitorError {
    /// A scalar parameter fell outside its physically meaningful domain
    /// (for example a non-positive plate area, gap, frequency or
    /// capacitance, or a negative time / resistance).
    #[error("invalid parameter `{name}` = {value}: {reason}")]
    InvalidParameter {
        /// Name of the offending parameter.
        name: &'static str,
        /// The supplied value, echoed back for diagnostics.
        value: f64,
        /// Why the value is rejected.
        reason: &'static str,
    },

    /// A network combinator was handed an empty list of capacitors. A
    /// series or parallel combination needs at least one branch.
    #[error("empty capacitor network: {0}")]
    EmptyNetwork(&'static str),
}

impl CapacitorError {
    /// Construct an [`InvalidParameter`](CapacitorError::InvalidParameter)
    /// error if `value` is not strictly greater than zero.
    ///
    /// Returns `Ok(value)` when the value is positive and finite, so call
    /// sites can validate inline:
    ///
    /// ```
    /// use valenx_capacitor::error::CapacitorError;
    ///
    /// let a = CapacitorError::require_positive("area_m2", 1.0e-4, "plate area")?;
    /// assert_eq!(a, 1.0e-4);
    /// assert!(CapacitorError::require_positive("gap_m", 0.0, "gap").is_err());
    /// # Ok::<(), CapacitorError>(())
    /// ```
    pub fn require_positive(name: &'static str, value: f64, reason: &'static str) -> Result<f64> {
        if value.is_finite() && value > 0.0 {
            Ok(value)
        } else {
            Err(CapacitorError::InvalidParameter {
                name,
                value,
                reason,
            })
        }
    }

    /// Construct an [`InvalidParameter`](CapacitorError::InvalidParameter)
    /// error if `value` is negative (or non-finite). Zero is accepted —
    /// useful for quantities such as elapsed time `t = 0` that are
    /// physically valid at the lower bound.
    ///
    /// Returns `Ok(value)` when the value is non-negative and finite.
    ///
    /// ```
    /// use valenx_capacitor::error::CapacitorError;
    ///
    /// assert!(CapacitorError::require_non_negative("t_s", 0.0, "time").is_ok());
    /// assert!(CapacitorError::require_non_negative("t_s", -1.0, "time").is_err());
    /// ```
    pub fn require_non_negative(
        name: &'static str,
        value: f64,
        reason: &'static str,
    ) -> Result<f64> {
        if value.is_finite() && value >= 0.0 {
            Ok(value)
        } else {
            Err(CapacitorError::InvalidParameter {
                name,
                value,
                reason,
            })
        }
    }

    /// Stable kebab-cased identifier for this error, suitable for logs and
    /// telemetry. Unlike the [`Display`](core::fmt::Display) text, this string is part of the
    /// crate's API contract and will not change for cosmetic reasons.
    pub fn code(&self) -> &'static str {
        match self {
            CapacitorError::InvalidParameter { .. } => "capacitor.invalid_parameter",
            CapacitorError::EmptyNetwork(_) => "capacitor.empty_network",
        }
    }

    /// Coarse category for grouping failures in a UI or dashboard.
    pub fn category(&self) -> ErrorCategory {
        match self {
            CapacitorError::InvalidParameter { .. } => ErrorCategory::Input,
            CapacitorError::EmptyNetwork(_) => ErrorCategory::Input,
        }
    }
}

/// Coarse classification of [`CapacitorError`] variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// The caller supplied an out-of-domain or otherwise invalid input.
    Input,
}
