//! The closed-loop ultimate-gain measurement that drives every rule.
//!
//! The Ziegler-Nichols ultimate-gain ("ultimate sensitivity") method
//! starts from a single experiment: with integral and derivative action
//! disabled, the proportional gain of a closed loop is raised until the
//! output exhibits sustained, constant-amplitude oscillation. The gain
//! at that point is the *ultimate gain* `Ku` and the period of the
//! oscillation is the *ultimate period* `Tu`. Those two numbers are the
//! complete input to the tuning table.
//!
//! [`UltimateMeasurement`] is a validated newtype-style wrapper around
//! that `(Ku, Tu)` pair: it can only be constructed through
//! [`UltimateMeasurement::new`], which rejects any non-finite or
//! non-positive value, so downstream tuning math never has to re-check
//! its inputs.

use serde::{Deserialize, Serialize};

use crate::error::PidTuningError;

/// A validated closed-loop ultimate-gain measurement `(Ku, Tu)`.
///
/// Both fields are guaranteed finite and strictly positive once the
/// value exists, because the only constructor is the validated
/// [`UltimateMeasurement::new`].
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct UltimateMeasurement {
    ku: f64,
    tu: f64,
}

impl UltimateMeasurement {
    /// Build a measurement from an ultimate gain and ultimate period.
    ///
    /// `ku` is the dimensionless proportional gain at which the loop
    /// first sustains oscillation; `tu` is the oscillation period in
    /// seconds. Both must be finite and strictly greater than zero.
    ///
    /// # Errors
    ///
    /// Returns [`PidTuningError::NonPositive`] if `ku` or `tu` is NaN,
    /// infinite, or not strictly positive. The error names the offending
    /// parameter (`"Ku"` or `"Tu"`).
    pub fn new(ku: f64, tu: f64) -> Result<Self, PidTuningError> {
        if !ku.is_finite() || ku <= 0.0 {
            return Err(PidTuningError::NonPositive {
                name: "Ku",
                value: ku,
            });
        }
        if !tu.is_finite() || tu <= 0.0 {
            return Err(PidTuningError::NonPositive {
                name: "Tu",
                value: tu,
            });
        }
        Ok(Self { ku, tu })
    }

    /// The ultimate gain `Ku` (dimensionless).
    pub fn ultimate_gain(&self) -> f64 {
        self.ku
    }

    /// The ultimate period `Tu`, in seconds.
    pub fn ultimate_period(&self) -> f64 {
        self.tu
    }
}
