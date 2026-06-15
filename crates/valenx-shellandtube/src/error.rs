//! Shell-and-tube sizing error taxonomy.
//!
//! Every fallible constructor in the crate funnels invalid input through
//! [`HxError`]. Each variant carries the offending parameter name and a
//! human-readable reason, and maps to a stable kebab-cased [`HxError::code`]
//! plus a coarse [`ErrorCategory`] for programmatic handling.

use thiserror::Error;

/// Errors raised while validating inputs or sizing a shell-and-tube
/// heat exchanger.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum HxError {
    /// A scalar parameter was outside its physically valid domain
    /// (for example a non-positive duty, area, diameter or length).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Name of the offending parameter.
        name: &'static str,
        /// Why the value was rejected.
        reason: String,
    },

    /// The LMTD correction factor `F` was outside the open-closed
    /// interval `(0, 1]` required by the LMTD method.
    #[error("correction factor F = {value} out of range: {reason}")]
    CorrectionFactorOutOfRange {
        /// Bit pattern of the rejected `F` value, rendered for context.
        value: f64,
        /// Why the value was rejected.
        reason: String,
    },

    /// The supplied terminal temperature differences are physically
    /// inconsistent — both must be strictly positive for a feasible
    /// exchanger (a sign change implies a temperature cross / pinch).
    #[error("infeasible temperature profile: {0}")]
    InfeasibleTemperatureProfile(String),
}

/// Coarse category for an [`HxError`], so callers can branch on the
/// *kind* of failure without matching every concrete variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// Caller-supplied numeric input was invalid.
    Input,
    /// The requested operating point is thermodynamically infeasible.
    Infeasible,
}

impl HxError {
    /// Stable, kebab-cased identifier for this error, suitable for logs,
    /// metrics labels or machine-readable diagnostics. The string is part
    /// of the crate's API contract and will not change for a given variant.
    pub fn code(&self) -> &'static str {
        match self {
            HxError::BadParameter { .. } => "shellandtube.bad_parameter",
            HxError::CorrectionFactorOutOfRange { .. } => {
                "shellandtube.correction_factor_out_of_range"
            }
            HxError::InfeasibleTemperatureProfile(_) => {
                "shellandtube.infeasible_temperature_profile"
            }
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            HxError::BadParameter { .. } | HxError::CorrectionFactorOutOfRange { .. } => {
                ErrorCategory::Input
            }
            HxError::InfeasibleTemperatureProfile(_) => ErrorCategory::Infeasible,
        }
    }

    /// Build a [`HxError::BadParameter`] from a static name and any
    /// displayable reason. Internal helper used by the sizing code.
    pub(crate) fn bad(name: &'static str, reason: impl Into<String>) -> Self {
        HxError::BadParameter {
            name,
            reason: reason.into(),
        }
    }

    /// Validate that a parameter is finite and strictly positive,
    /// returning it unchanged or a [`HxError::BadParameter`].
    ///
    /// # Errors
    ///
    /// Returns [`HxError::BadParameter`] when `value` is NaN, infinite,
    /// zero or negative.
    pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64, HxError> {
        if !value.is_finite() {
            return Err(HxError::bad(name, format!("must be finite, got {value}")));
        }
        if value <= 0.0 {
            return Err(HxError::bad(name, format!("must be > 0, got {value}")));
        }
        Ok(value)
    }
}
