//! Error type for the projectile-ballistics crate.

use thiserror::Error;

/// Shorthand for `Result<T, ProjectileError>`.
pub type Result<T> = core::result::Result<T, ProjectileError>;

/// Anything that can go wrong validating inputs or running a trajectory.
///
/// This enum is `#[non_exhaustive]`: new variants may be added in future
/// releases without it being a breaking change, so downstream `match`
/// arms must include a wildcard.
#[derive(Debug, Error, Clone, PartialEq)]
#[non_exhaustive]
pub enum ProjectileError {
    /// A scalar parameter was non-physical in a way that would otherwise
    /// feed a silent `NaN`/`Inf` into the equations — e.g. a non-finite,
    /// negative, or zero value where a strictly positive one is required
    /// (speed, gravity, mass, …). Carries the field name and the value.
    #[error("invalid parameter `{field}` = {value} ({reason})")]
    InvalidParameter {
        /// Which quantity was out of range (e.g. `"speed"`, `"gravity"`).
        field: &'static str,
        /// The offending value.
        value: f64,
        /// Short human-readable reason (e.g. `"must be > 0"`).
        reason: &'static str,
    },

    /// A launch angle (radians) was outside the supported `[0, π/2]`
    /// quadrant for which the elementary range/apex formulae and the
    /// drag integrator are defined.
    #[error("launch angle {value} rad is outside the supported [0, π/2] range")]
    AngleOutOfRange {
        /// The offending angle in radians.
        value: f64,
    },

    /// The numerical integrator hit its hard step budget before the
    /// projectile returned to (or fell below) the launch height. Carries
    /// the budget that was exhausted.
    #[error("integration exceeded the maximum step count ({0})")]
    StepBudgetExceeded(u64),

    /// The golden-section optimal-angle search was given an invalid
    /// bracket (lower bound not strictly below the upper bound, or a
    /// bound outside `[0, π/2]`). Carries a short reason.
    #[error("invalid optimisation bracket: {0}")]
    InvalidBracket(&'static str),
}

impl ProjectileError {
    /// Validate that `value` is finite and strictly positive, returning it
    /// unchanged or an [`InvalidParameter`](ProjectileError::InvalidParameter).
    ///
    /// Used by the public constructors so every entry point rejects
    /// non-physical magnitudes up front rather than silently producing a
    /// `NaN` result downstream.
    pub fn require_positive(field: &'static str, value: f64) -> Result<f64> {
        if !value.is_finite() {
            Err(ProjectileError::InvalidParameter {
                field,
                value,
                reason: "must be finite",
            })
        } else if value <= 0.0 {
            Err(ProjectileError::InvalidParameter {
                field,
                value,
                reason: "must be > 0",
            })
        } else {
            Ok(value)
        }
    }

    /// Validate that `value` is finite and non-negative (`>= 0`).
    ///
    /// Used for quantities such as the drag coefficient where exactly zero
    /// is physically meaningful (the vacuum limit) but a negative or
    /// non-finite value is not.
    pub fn require_non_negative(field: &'static str, value: f64) -> Result<f64> {
        if !value.is_finite() {
            Err(ProjectileError::InvalidParameter {
                field,
                value,
                reason: "must be finite",
            })
        } else if value < 0.0 {
            Err(ProjectileError::InvalidParameter {
                field,
                value,
                reason: "must be >= 0",
            })
        } else {
            Ok(value)
        }
    }

    /// Validate that `angle_rad` lies in the closed `[0, π/2]` quadrant.
    ///
    /// A non-finite angle (`NaN`/`±∞`) also fails this check, since it is
    /// never contained in the range.
    pub fn require_angle(angle_rad: f64) -> Result<f64> {
        if (0.0..=core::f64::consts::FRAC_PI_2).contains(&angle_rad) {
            Ok(angle_rad)
        } else {
            Err(ProjectileError::AngleOutOfRange { value: angle_rad })
        }
    }
}
