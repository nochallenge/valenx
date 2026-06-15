//! Heat-transfer error taxonomy.
//!
//! Every fallible constructor in this crate validates its physical
//! inputs (lengths, areas, conductivities and coefficients must be
//! strictly positive; temperatures and counts must be finite) and
//! returns a [`HeatTransferError`] on violation. The error type carries
//! stable [`code`](HeatTransferError::code) and
//! [`category`](HeatTransferError::category) accessors for telemetry.

use thiserror::Error;

/// Errors raised by the heat-transfer crate.
#[derive(Debug, Error)]
pub enum HeatTransferError {
    /// A physical parameter was outside its valid domain (for example a
    /// non-positive area, length, conductivity, or convection
    /// coefficient).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name (stable, machine-readable).
        name: &'static str,
        /// Human-readable reason the value was rejected.
        reason: String,
    },

    /// A resistance network was empty where at least one element is
    /// required (a series or parallel combination of zero resistors is
    /// undefined).
    #[error("empty resistance network: {0}")]
    EmptyNetwork(&'static str),
}

/// Coarse error category, useful for routing / metrics.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The caller supplied an invalid value.
    Input,
    /// The requested operation is not defined for the given structure.
    Domain,
}

impl HeatTransferError {
    /// Stable, kebab-cased identifier for this error variant.
    ///
    /// Suitable as a log/telemetry key; never changes for a given
    /// variant.
    pub fn code(&self) -> &'static str {
        match self {
            HeatTransferError::BadParameter { .. } => "heat-transfer.bad_parameter",
            HeatTransferError::EmptyNetwork(_) => "heat-transfer.empty_network",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            HeatTransferError::BadParameter { .. } => ErrorCategory::Input,
            HeatTransferError::EmptyNetwork(_) => ErrorCategory::Domain,
        }
    }
}

/// Convenience result alias used throughout the crate.
pub type Result<T> = core::result::Result<T, HeatTransferError>;

/// Internal helper: reject a non-positive (or non-finite) physical
/// quantity with a uniform [`HeatTransferError::BadParameter`].
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64> {
    if !value.is_finite() || value <= 0.0 {
        return Err(HeatTransferError::BadParameter {
            name,
            reason: format!("must be finite and > 0, got {value}"),
        });
    }
    Ok(value)
}

/// Internal helper: reject a non-finite quantity (used for values that
/// may legitimately be negative or zero, such as a temperature).
pub(crate) fn require_finite(name: &'static str, value: f64) -> Result<f64> {
    if !value.is_finite() {
        return Err(HeatTransferError::BadParameter {
            name,
            reason: format!("must be finite, got {value}"),
        });
    }
    Ok(value)
}
