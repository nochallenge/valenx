//! Combustion error taxonomy.

use thiserror::Error;

/// Errors raised by combustion stoichiometry and flame-temperature
/// calculations.
#[derive(Debug, Error)]
pub enum CombustionError {
    /// A fuel formula `CxHy` had a non-physical atom count.
    ///
    /// At least one carbon atom is required (`x >= 1`); hydrogen may be
    /// zero only for pure carbon, which is rejected here because the
    /// closed-form hydrocarbon model needs `y >= 1` to stay meaningful.
    #[error("invalid fuel formula C{carbon}H{hydrogen}: {reason}")]
    InvalidFuel {
        /// Number of carbon atoms supplied.
        carbon: u32,
        /// Number of hydrogen atoms supplied.
        hydrogen: u32,
        /// Why the formula was rejected.
        reason: &'static str,
    },

    /// A parameter that must be strictly positive was zero or negative.
    #[error("bad parameter `{name}` = {value}: {reason}")]
    BadParameter {
        /// Parameter name.
        name: &'static str,
        /// Offending value.
        value: f64,
        /// Reason it was rejected.
        reason: &'static str,
    },
}

/// Coarse error category, mirroring the sibling-crate convention.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// Caller-supplied input was invalid.
    Input,
    /// A tunable knob (heating value, cp, temperature) was invalid.
    Config,
}

impl CombustionError {
    /// Validate that `value` is finite and strictly positive, returning a
    /// [`CombustionError::BadParameter`] otherwise.
    ///
    /// Used by the stoichiometry and flame-temperature routines so the
    /// "must be > 0" guard is written once and rejects NaN, infinities,
    /// zero, and negatives uniformly. The positive comparison keeps the
    /// check readable (no negated partial-ordering operator).
    pub(crate) fn require_positive(
        name: &'static str,
        value: f64,
        reason: &'static str,
    ) -> Result<(), CombustionError> {
        if value.is_finite() && value > 0.0 {
            Ok(())
        } else {
            Err(CombustionError::BadParameter {
                name,
                value,
                reason,
            })
        }
    }

    /// Stable kebab-cased identifier for logs and tests.
    pub fn code(&self) -> &'static str {
        match self {
            CombustionError::InvalidFuel { .. } => "combustion.invalid_fuel",
            CombustionError::BadParameter { .. } => "combustion.bad_parameter",
        }
    }

    /// Coarse category for the error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            CombustionError::InvalidFuel { .. } => ErrorCategory::Input,
            CombustionError::BadParameter { .. } => ErrorCategory::Config,
        }
    }
}
