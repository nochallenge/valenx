//! Error taxonomy for the lead-screw workbench.
//!
//! Every fallible constructor in this crate returns
//! [`LeadScrewError`]. The variants are deliberately specific so a
//! caller (CLI, GUI field validator, or upstream solver) can map a
//! failure back to the exact offending input.

use thiserror::Error;

/// Errors raised while validating lead-screw inputs.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum LeadScrewError {
    /// A quantity that must be strictly positive was zero or negative.
    ///
    /// Lead, diameter, torque, RPM and the microstep count all feed
    /// into divisions or physically-meaningful ratios where a
    /// non-positive value is nonsensical (and, for `lead`, would divide
    /// by zero in the thrust / resolution formulas).
    #[error("`{name}` must be > 0, got {value}")]
    NotPositive {
        /// Name of the offending parameter.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// An efficiency was outside the open-closed interval `(0, 1]`.
    ///
    /// Mechanical efficiency `eta` is a dimensionless fraction of the
    /// applied screw torque that survives friction to become useful
    /// axial work. A value `<= 0` means no thrust is ever produced
    /// (degenerate), and a value `> 1` would manufacture energy.
    #[error("efficiency `{name}` must be in (0, 1], got {value}")]
    EfficiencyOutOfRange {
        /// Name of the offending efficiency parameter.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A friction coefficient was negative.
    ///
    /// The coefficient of friction `mu` between nut and screw thread is
    /// a non-negative ratio; `0` models the idealized frictionless
    /// (always back-drivable) limit.
    #[error("friction coefficient `{name}` must be >= 0, got {value}")]
    NegativeFriction {
        /// Name of the offending parameter.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A non-finite (NaN or infinite) value was supplied.
    ///
    /// Caught explicitly so that a stray `NaN` never silently
    /// propagates through the closed-form formulas (where, e.g.,
    /// `NaN > 0.0` is `false` and would slip past a naive range check).
    #[error("`{name}` must be finite, got {value}")]
    NotFinite {
        /// Name of the offending parameter.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// The microstep count rounded to zero steps per revolution.
    ///
    /// Distinct from [`LeadScrewError::NotPositive`]: the raw input was
    /// a valid positive integer-like count but is reported separately so
    /// a driver-configuration UI can flag "microstepping" specifically.
    #[error("microsteps per revolution must be >= 1, got {0}")]
    ZeroMicrosteps(u32),
}

impl LeadScrewError {
    /// Stable, kebab-cased identifier for this error.
    ///
    /// Useful for logging / telemetry where the human-readable
    /// [`Display`](std::fmt::Display) string is not a stable key.
    pub fn code(&self) -> &'static str {
        match self {
            LeadScrewError::NotPositive { .. } => "leadscrew.not-positive",
            LeadScrewError::EfficiencyOutOfRange { .. } => "leadscrew.efficiency-out-of-range",
            LeadScrewError::NegativeFriction { .. } => "leadscrew.negative-friction",
            LeadScrewError::NotFinite { .. } => "leadscrew.not-finite",
            LeadScrewError::ZeroMicrosteps(_) => "leadscrew.zero-microsteps",
        }
    }
}

/// Internal helper: reject non-finite and non-positive values in one
/// place so every constructor enforces the same contract.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64, LeadScrewError> {
    if !value.is_finite() {
        return Err(LeadScrewError::NotFinite { name, value });
    }
    if value <= 0.0 {
        return Err(LeadScrewError::NotPositive { name, value });
    }
    Ok(value)
}

/// Internal helper: reject non-finite and out-of-range efficiencies.
pub(crate) fn require_efficiency(name: &'static str, value: f64) -> Result<f64, LeadScrewError> {
    if !value.is_finite() {
        return Err(LeadScrewError::NotFinite { name, value });
    }
    if value <= 0.0 || value > 1.0 {
        return Err(LeadScrewError::EfficiencyOutOfRange { name, value });
    }
    Ok(value)
}

/// Internal helper: reject non-finite and negative friction values.
pub(crate) fn require_non_negative_friction(
    name: &'static str,
    value: f64,
) -> Result<f64, LeadScrewError> {
    if !value.is_finite() {
        return Err(LeadScrewError::NotFinite { name, value });
    }
    if value < 0.0 {
        return Err(LeadScrewError::NegativeFriction { name, value });
    }
    Ok(value)
}
