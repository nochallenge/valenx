//! Error taxonomy for the PID-tuning calculator.
//!
//! Every fallible entry point returns [`PidTuningError`]. The variants
//! are deliberately narrow: the only way a Ziegler-Nichols ultimate-gain
//! computation can fail is a non-physical input (a non-positive,
//! infinite, or NaN ultimate gain or oscillation period). The validated
//! constructors on [`crate::ultimate::UltimateMeasurement`] are the sole
//! producers of these errors, so callers always get a checked value or a
//! precise reason for rejection.

use thiserror::Error;

/// Errors raised while validating tuning inputs.
#[derive(Debug, Error)]
pub enum PidTuningError {
    /// A required parameter was not a finite, strictly positive number.
    ///
    /// The ultimate gain `Ku` and the ultimate period `Tu` must both be
    /// real, finite, and greater than zero: a gain of zero never
    /// sustains oscillation and a period of zero is unphysical. NaN and
    /// the infinities are rejected here as well.
    #[error("parameter `{name}` must be finite and > 0, got {value}")]
    NonPositive {
        /// Stable parameter name (`"Ku"` or `"Tu"`).
        name: &'static str,
        /// The offending value as supplied by the caller.
        value: f64,
    },
}

/// Coarse classification of a [`PidTuningError`].
///
/// Useful for UI / logging layers that want to bucket failures without
/// matching every concrete variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// The caller supplied an invalid measurement value.
    Input,
}

impl PidTuningError {
    /// Stable, kebab-cased identifier for this error.
    ///
    /// Intended for machine consumption (telemetry keys, test
    /// assertions) where the human-readable [`Display`](std::fmt::Display)
    /// string is unsuitable.
    pub fn code(&self) -> &'static str {
        match self {
            PidTuningError::NonPositive { .. } => "pidtuning.non-positive",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            PidTuningError::NonPositive { .. } => ErrorCategory::Input,
        }
    }
}
