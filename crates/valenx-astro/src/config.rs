//! Run configuration for an ascent simulation.

use serde::{Deserialize, Serialize};

use crate::error::AstroError;
use crate::guidance::GuidanceProgram;
use crate::sim::{MAX_SAMPLES, MAX_SIM_STEPS};
use crate::wind::WindModel;

/// Smallest accepted integration step (s). A 0.1 ms step is already far
/// finer than any ascent dynamics require; anything below this is treated
/// as an input error (it only serves to explode the step count).
pub const MIN_TIME_STEP: f64 = 1e-4;

/// Largest accepted simulated duration (s) — about four months. No
/// launch ascent or insertion runs anywhere near this; a larger value is
/// rejected as an input error.
pub const MAX_SIM_TIME: f64 = 1e7;

/// How the vehicle is steered and when its engines cut off.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum GuidanceMode {
    /// Fly the open-loop gravity turn and burn every stage to
    /// depletion. Simple and self-stabilising, but an overpowered
    /// vehicle ends on an eccentric orbit (no insertion targeting).
    OpenLoopGravityTurn,

    /// Closed-loop orbital insertion to a target circular altitude.
    ///
    /// The vehicle flies the gravity turn until its osculating apoapsis
    /// reaches `target_altitude_m`, then **cuts off and coasts** to
    /// apoapsis, then **reignites** and burns prograde-horizontal to
    /// raise periapsis up to the target — circularising the orbit. The
    /// leftover propellant from the early cutoff is what funds the
    /// circularisation burn.
    ClosedLoopInsertion {
        /// Target circular-orbit altitude above the equatorial radius (m).
        target_altitude_m: f64,
    },
}

impl Default for GuidanceMode {
    fn default() -> Self {
        Self::OpenLoopGravityTurn
    }
}

/// Everything that parameterises a run other than the vehicle itself.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct AscentConfig {
    /// Launch-site elevation above the equatorial radius (m).
    pub launch_altitude_m: f64,
    /// Open-loop pitch program.
    pub guidance: GuidanceProgram,
    /// Fixed integration step (s).
    pub time_step: f64,
    /// Hard cap on simulated time (s); the run stops at the latest here
    /// even if no other termination fires.
    pub max_time: f64,
    /// Minimum spacing between recorded trajectory samples (s). Keeps
    /// the output series small for long flights.
    pub sample_interval: f64,
    /// Steering / cutoff strategy.
    #[serde(default)]
    pub mode: GuidanceMode,
    /// Horizontal wind model applied to the drag calculation.
    #[serde(default)]
    pub wind: WindModel,
}

impl Default for AscentConfig {
    fn default() -> Self {
        Self {
            launch_altitude_m: 0.0,
            guidance: GuidanceProgram::default(),
            time_step: 0.1,
            max_time: 2_000.0,
            sample_interval: 1.0,
            mode: GuidanceMode::default(),
            wind: WindModel::default(),
        }
    }
}

impl AscentConfig {
    /// Validate the run settings.
    pub fn validate(&self) -> Result<(), AstroError> {
        if !self.time_step.is_finite() || self.time_step <= 0.0 {
            return Err(AstroError::InvalidIntegration("time_step must be > 0"));
        }
        if !self.max_time.is_finite() || self.max_time <= 0.0 {
            return Err(AstroError::InvalidIntegration("max_time must be > 0"));
        }
        if !self.sample_interval.is_finite() || self.sample_interval <= 0.0 {
            return Err(AstroError::InvalidIntegration(
                "sample_interval must be > 0",
            ));
        }
        // Absolute bounds so a finite-but-absurd config cannot drive the
        // fixed-step loop or the sample Vec unbounded.
        if self.time_step < MIN_TIME_STEP {
            return Err(AstroError::InvalidIntegration(
                "time_step below the minimum (1e-4 s)",
            ));
        }
        if self.max_time > MAX_SIM_TIME {
            return Err(AstroError::InvalidIntegration(
                "max_time above the maximum (1e7 s)",
            ));
        }
        // Cap the derived step and sample counts. Compute in f64 and
        // compare against the cap before any `as u64` cast so an absurd
        // ratio can't silently wrap. (`as u64` saturates, but we want a
        // clean OutOfRange rather than a saturated budget.)
        let step_count = (self.max_time / self.time_step).ceil();
        if step_count > MAX_SIM_STEPS as f64 {
            return Err(AstroError::OutOfRange {
                what: "steps",
                value: step_count.min(u64::MAX as f64) as u64,
                max: MAX_SIM_STEPS,
            });
        }
        let sample_count = (self.max_time / self.sample_interval).ceil();
        if sample_count > MAX_SAMPLES as f64 {
            return Err(AstroError::OutOfRange {
                what: "samples",
                value: sample_count.min(u64::MAX as f64) as u64,
                max: MAX_SAMPLES,
            });
        }
        if !self.launch_altitude_m.is_finite() {
            return Err(AstroError::InvalidIntegration(
                "launch_altitude_m not finite",
            ));
        }
        if let GuidanceMode::ClosedLoopInsertion { target_altitude_m } = self.mode {
            if !target_altitude_m.is_finite() || target_altitude_m <= 0.0 {
                return Err(AstroError::InvalidGuidance("target_altitude_m must be > 0"));
            }
        }
        self.guidance.validate()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::{MAX_SAMPLES, MAX_SIM_STEPS};

    #[test]
    fn default_and_presets_validate() {
        assert!(AscentConfig::default().validate().is_ok());
        assert!(crate::presets::leo_ascent_config().validate().is_ok());
        assert!(crate::presets::leo_insertion_config().validate().is_ok());
    }

    #[test]
    fn rejects_absurdly_small_time_step() {
        // A 1 ns step is physically absurd for an ascent and would, with
        // any sane max_time, blow past the step budget.
        let cfg = AscentConfig {
            time_step: 1e-9,
            ..AscentConfig::default()
        };
        assert!(cfg.validate().is_err(), "1 ns step must be rejected");
    }

    #[test]
    fn rejects_enormous_max_time() {
        let cfg = AscentConfig {
            max_time: 1e15,
            ..AscentConfig::default()
        };
        assert!(cfg.validate().is_err(), "1e15 s max_time must be rejected");
    }

    #[test]
    fn rejects_step_count_above_max() {
        // dt=1e-9, max_time=1e15 -> ~1e24 steps: must be a clean Err,
        // not a hang. (Each individual bound also trips, but the step
        // count is the OOM/hang guard the finding is about.)
        let cfg = AscentConfig {
            time_step: 1e-3,
            max_time: 1e7, // both individually in-range...
            sample_interval: 1.0,
            ..AscentConfig::default()
        };
        // 1e7 / 1e-3 = 1e10 steps > MAX_SIM_STEPS (1e8) -> reject.
        match cfg.validate() {
            Err(AstroError::OutOfRange { what, max, .. }) => {
                assert_eq!(what, "steps");
                assert_eq!(max, MAX_SIM_STEPS);
            }
            other => panic!("expected OutOfRange(steps), got {other:?}"),
        }
    }

    #[test]
    fn rejects_sample_count_above_max() {
        // Keep the step count in range but make the sample count explode.
        let cfg = AscentConfig {
            time_step: 1.0,
            max_time: 1e7,        // 1e7 steps, in-range
            sample_interval: 1e-3, // 1e10 samples > MAX_SAMPLES -> reject
            ..AscentConfig::default()
        };
        match cfg.validate() {
            Err(AstroError::OutOfRange { what, max, .. }) => {
                assert_eq!(what, "samples");
                assert_eq!(max, MAX_SAMPLES);
            }
            other => panic!("expected OutOfRange(samples), got {other:?}"),
        }
    }
}
