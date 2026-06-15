//! Vibration analysis error taxonomy.
//!
//! Every fallible constructor in this crate returns
//! [`Result<_, VibrationError>`]. The error carries stable
//! [`code`](VibrationError::code) and [`category`](VibrationError::category)
//! accessors so callers (telemetry, UI) can branch on the failure kind
//! without string-matching the human-readable message.

use thiserror::Error;

/// Errors raised while building or evaluating a vibration model.
#[derive(Debug, Error)]
pub enum VibrationError {
    /// A physical parameter was outside its valid domain.
    ///
    /// Mass `m`, stiffness `k` and the damping coefficient `c` must all
    /// be physically meaningful: `m` and `k` must be strictly positive
    /// (a zero or negative mass / stiffness has no natural frequency)
    /// and `c` must be non-negative. The offending `name` and a short
    /// `reason` are reported.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name, e.g. `"mass_kg"`.
        name: &'static str,
        /// Why the value was rejected.
        reason: String,
    },

    /// The closed-form model does not apply for the requested regime.
    ///
    /// For example, the damped natural frequency
    /// `wd = wn*sqrt(1 - zeta^2)` is only real for an *underdamped*
    /// system (`zeta < 1`); asking for `wd` on a critically- or
    /// over-damped system raises this.
    #[error("not applicable: {0}")]
    NotApplicable(String),

    /// Two successive-peak amplitudes given to the logarithmic-decrement
    /// estimator were not a valid decaying pair (both must be strictly
    /// positive and the later peak must not exceed the earlier one).
    #[error("invalid decay data: {0}")]
    InvalidDecay(String),
}

/// Coarse category for a [`VibrationError`], for grouping in telemetry
/// or UI without inspecting the specific variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// Caller supplied an invalid input value.
    Input,
    /// The requested quantity does not exist for this model's regime.
    Domain,
}

impl VibrationError {
    /// A stable, kebab-cased identifier for this error.
    ///
    /// Unlike the [`Display`](std::fmt::Display) message, this string is
    /// part of the crate's contract and safe to match on.
    pub fn code(&self) -> &'static str {
        match self {
            VibrationError::BadParameter { .. } => "vibration.bad_parameter",
            VibrationError::NotApplicable(_) => "vibration.not_applicable",
            VibrationError::InvalidDecay(_) => "vibration.invalid_decay",
        }
    }

    /// The coarse [`ErrorCategory`] this error belongs to.
    pub fn category(&self) -> ErrorCategory {
        match self {
            VibrationError::BadParameter { .. } => ErrorCategory::Input,
            VibrationError::InvalidDecay(_) => ErrorCategory::Input,
            VibrationError::NotApplicable(_) => ErrorCategory::Domain,
        }
    }
}
