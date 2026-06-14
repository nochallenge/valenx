//! Discrete-time PID controller.
//!
//! Implements the textbook proportional-integral-derivative control law
//! in its **positional** form
//!
//! ```text
//! u(t) = Kp * e + Ki * I + Kd * D
//! ```
//!
//! where `e = setpoint - measurement` is the error, `I` is the running
//! integral of the error (accumulated by the rectangle rule
//! `I += e * dt`) and `D` is the discrete derivative of the error
//! (`D = (e - e_prev) / dt`). The controller is stepped once per sample
//! at a fixed sample time `dt`.
//!
//! ## Model
//!
//! Two extras make the controller usable rather than a toy:
//!
//! - **Output clamping** — the command is saturated to an optional
//!   `[out_min, out_max]` range.
//! - **Anti-windup (integral clamping)** — when the raw command
//!   saturates, the integral term is held / clamped so the accumulator
//!   does not "wind up" and cause a large overshoot when the error
//!   finally reverses. We use the simplest robust scheme: clamp the
//!   *integral contribution* `Ki * I` to the output range so it can never
//!   alone push the command past the limits.
//!
//! ## Honest scope
//!
//! This is a standard fixed-step digital PID — the form found in any
//! controls textbook and in countless embedded loops. It does **not**
//! model derivative filtering, setpoint weighting / two-degree-of-freedom
//! structures, bumpless auto/manual transfer, or gain scheduling. It is a
//! research/educational building block, not a certified industrial
//! controller.

use serde::{Deserialize, Serialize};

use crate::error::{ControlsError, Result};

/// Gains and limits for a discrete PID controller.
///
/// Construct with [`PidConfig::new`], which validates that every gain is
/// finite and the sample time is finite and strictly positive.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PidConfig {
    /// Proportional gain `Kp`.
    pub kp: f64,
    /// Integral gain `Ki`.
    pub ki: f64,
    /// Derivative gain `Kd`.
    pub kd: f64,
    /// Sample time `dt` (s), strictly positive.
    pub dt: f64,
    /// Optional lower output clamp. `None` means unbounded below.
    pub out_min: Option<f64>,
    /// Optional upper output clamp. `None` means unbounded above.
    pub out_max: Option<f64>,
}

impl PidConfig {
    /// Construct a PID configuration with the given gains and sample time
    /// and **no** output limits.
    ///
    /// # Errors
    ///
    /// Returns [`ControlsError::InvalidParameter`] if any gain is
    /// non-finite, or if `dt` is non-finite or not strictly positive.
    pub fn new(kp: f64, ki: f64, kd: f64, dt: f64) -> Result<Self> {
        for (name, g) in [("kp", kp), ("ki", ki), ("kd", kd)] {
            if !g.is_finite() {
                return Err(ControlsError::invalid(name, "gain must be finite"));
            }
        }
        if !dt.is_finite() || dt <= 0.0 {
            return Err(ControlsError::invalid(
                "dt",
                "sample time must be finite and > 0",
            ));
        }
        Ok(Self {
            kp,
            ki,
            kd,
            dt,
            out_min: None,
            out_max: None,
        })
    }

    /// Set symmetric / asymmetric output clamps, returning the updated
    /// config (builder style).
    ///
    /// # Errors
    ///
    /// Returns [`ControlsError::InvalidParameter`] if either bound is
    /// non-finite, or if `min > max`.
    pub fn with_limits(mut self, min: f64, max: f64) -> Result<Self> {
        if !min.is_finite() || !max.is_finite() {
            return Err(ControlsError::invalid(
                "out_limits",
                "output limits must be finite",
            ));
        }
        if min > max {
            return Err(ControlsError::invalid(
                "out_limits",
                "out_min must not exceed out_max",
            ));
        }
        self.out_min = Some(min);
        self.out_max = Some(max);
        Ok(self)
    }
}

/// A stateful discrete PID controller.
///
/// Create one from a [`PidConfig`] with [`Pid::new`], then drive it with
/// repeated [`Pid::update`] calls — one per sample. State (the running
/// integral and the previous error) lives in the struct; call
/// [`Pid::reset`] to clear it.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Pid {
    cfg: PidConfig,
    integral: f64,
    prev_error: f64,
    initialized: bool,
}

impl Pid {
    /// Build a controller from a validated [`PidConfig`].
    pub fn new(cfg: PidConfig) -> Self {
        Self {
            cfg,
            integral: 0.0,
            prev_error: 0.0,
            initialized: false,
        }
    }

    /// The configuration this controller was built with.
    pub fn config(&self) -> &PidConfig {
        &self.cfg
    }

    /// Current value of the running error integral `I`.
    pub fn integral(&self) -> f64 {
        self.integral
    }

    /// Clear the integral accumulator and the stored previous error.
    ///
    /// Use this on a setpoint change or when re-engaging the loop so the
    /// derivative term does not see a spurious one-sample jump and the
    /// integral starts fresh.
    pub fn reset(&mut self) {
        self.integral = 0.0;
        self.prev_error = 0.0;
        self.initialized = false;
    }

    /// Advance the controller one sample and return the command `u`.
    ///
    /// `setpoint` is the reference and `measurement` the latest plant
    /// output; the error is `e = setpoint - measurement`. The integral is
    /// updated by the rectangle rule and the derivative by a backward
    /// difference. On the very first call after construction or
    /// [`reset`](Self::reset) the derivative term is taken as zero (no
    /// previous error exists), avoiding a derivative kick.
    ///
    /// The returned command is clamped to the configured output limits;
    /// the integral is held (anti-windup) so its contribution `Ki * I`
    /// can never alone exceed those limits.
    pub fn update(&mut self, setpoint: f64, measurement: f64) -> f64 {
        let error = setpoint - measurement;
        let dt = self.cfg.dt;

        // Derivative: backward difference; zero on the first sample.
        let derivative = if self.initialized {
            (error - self.prev_error) / dt
        } else {
            0.0
        };

        // Integrate with the rectangle rule.
        self.integral += error * dt;

        // Anti-windup: clamp the *integral contribution* Ki*I to the
        // output range so the accumulator cannot wind up beyond what the
        // actuator can deliver. Only meaningful when Ki != 0.
        if self.cfg.ki != 0.0 {
            let i_term = self.cfg.ki * self.integral;
            let clamped = clamp_opt(i_term, self.cfg.out_min, self.cfg.out_max);
            if clamped != i_term {
                self.integral = clamped / self.cfg.ki;
            }
        }

        let raw = self.cfg.kp * error + self.cfg.ki * self.integral + self.cfg.kd * derivative;

        self.prev_error = error;
        self.initialized = true;

        clamp_opt(raw, self.cfg.out_min, self.cfg.out_max)
    }
}

/// Clamp `value` to the optional `[min, max]` interval. A `None` bound
/// means unbounded on that side.
fn clamp_opt(value: f64, min: Option<f64>, max: Option<f64>) -> f64 {
    let mut v = value;
    if let Some(lo) = min {
        if v < lo {
            v = lo;
        }
    }
    if let Some(hi) = max {
        if v > hi {
            v = hi;
        }
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn rejects_bad_config() {
        assert!(PidConfig::new(f64::NAN, 0.0, 0.0, 0.1).is_err());
        assert!(PidConfig::new(1.0, f64::INFINITY, 0.0, 0.1).is_err());
        assert!(PidConfig::new(1.0, 0.0, 0.0, 0.0).is_err());
        assert!(PidConfig::new(1.0, 0.0, 0.0, -0.1).is_err());
        assert!(PidConfig::new(1.0, 0.0, 0.0, f64::NAN).is_err());
        assert!(PidConfig::new(1.0, 0.5, 0.05, 0.1).is_ok());
    }

    #[test]
    fn limits_must_be_ordered_and_finite() {
        let cfg = PidConfig::new(1.0, 0.0, 0.0, 0.1).unwrap();
        assert!(cfg.with_limits(5.0, -5.0).is_err());
        assert!(cfg.with_limits(f64::NAN, 1.0).is_err());
        let ok = cfg.with_limits(-10.0, 10.0).unwrap();
        assert_eq!(ok.out_min, Some(-10.0));
        assert_eq!(ok.out_max, Some(10.0));
    }

    #[test]
    fn pure_proportional_matches_kp_times_error() {
        // With Ki = Kd = 0, u = Kp * e exactly.
        let cfg = PidConfig::new(2.0, 0.0, 0.0, 0.1).unwrap();
        let mut pid = Pid::new(cfg);
        let u = pid.update(10.0, 4.0); // e = 6
        assert!((u - 12.0).abs() < EPS, "u = {u}");
    }

    #[test]
    fn integral_accumulates_error_over_time() {
        // Constant error of 1.0 each step, Ki = 1, dt = 0.5:
        // after n steps I = n * e * dt = n * 0.5, u = Ki*I (Kp=Kd=0).
        let cfg = PidConfig::new(0.0, 1.0, 0.0, 0.5).unwrap();
        let mut pid = Pid::new(cfg);
        let u1 = pid.update(1.0, 0.0);
        let u2 = pid.update(1.0, 0.0);
        let u3 = pid.update(1.0, 0.0);
        assert!((u1 - 0.5).abs() < EPS, "u1 = {u1}");
        assert!((u2 - 1.0).abs() < EPS, "u2 = {u2}");
        assert!((u3 - 1.5).abs() < EPS, "u3 = {u3}");
        assert!((pid.integral() - 1.5).abs() < EPS);
    }

    #[test]
    fn derivative_is_zero_on_first_sample_then_backward_difference() {
        // Kd = 2, dt = 0.1. First sample: D = 0 -> u = 0 (Kp=Ki=0).
        // Second sample: error jumps 0 -> 5, D = (5-0)/0.1 = 50,
        // u = Kd*D = 100. (Measurement chosen so error is exactly known.)
        let cfg = PidConfig::new(0.0, 0.0, 2.0, 0.1).unwrap();
        let mut pid = Pid::new(cfg);
        let u0 = pid.update(0.0, 0.0); // e = 0, first sample D = 0
        assert!(u0.abs() < EPS, "u0 = {u0}");
        let u1 = pid.update(5.0, 0.0); // e = 5, D = (5-0)/0.1 = 50
        assert!((u1 - 100.0).abs() < EPS, "u1 = {u1}");
    }

    #[test]
    fn output_is_clamped_to_limits() {
        let cfg = PidConfig::new(10.0, 0.0, 0.0, 0.1)
            .unwrap()
            .with_limits(-5.0, 5.0)
            .unwrap();
        let mut pid = Pid::new(cfg);
        // e = 100 -> raw = 1000, clamped to +5.
        let hi = pid.update(100.0, 0.0);
        assert!((hi - 5.0).abs() < EPS, "hi = {hi}");
        pid.reset();
        // e = -100 -> raw = -1000, clamped to -5.
        let lo = pid.update(-100.0, 0.0);
        assert!((lo - (-5.0)).abs() < EPS, "lo = {lo}");
    }

    #[test]
    fn anti_windup_bounds_the_integral_contribution() {
        // Drive a persistent error into a saturating integrator; the
        // integral term must not wind up past the output limit.
        let cfg = PidConfig::new(0.0, 1.0, 0.0, 1.0)
            .unwrap()
            .with_limits(-2.0, 2.0)
            .unwrap();
        let mut pid = Pid::new(cfg);
        for _ in 0..50 {
            pid.update(10.0, 0.0); // large constant error
        }
        // Ki * I is clamped to the +2 output limit, so I itself is <= 2.
        assert!(
            pid.integral() <= 2.0 + EPS,
            "integral wound up to {}",
            pid.integral()
        );
        // And the command stays clamped.
        let u = pid.update(10.0, 0.0);
        assert!((u - 2.0).abs() < EPS, "u = {u}");
    }

    #[test]
    fn reset_clears_state() {
        let cfg = PidConfig::new(1.0, 1.0, 0.0, 0.1).unwrap();
        let mut pid = Pid::new(cfg);
        pid.update(5.0, 0.0);
        assert!(pid.integral().abs() > 0.0);
        pid.reset();
        assert!(pid.integral().abs() < EPS);
    }

    #[test]
    fn full_pid_equals_sum_of_three_terms() {
        // Hand-compute u = Kp*e + Ki*I + Kd*D on the second sample.
        // Kp=1, Ki=2, Kd=3, dt=0.5.
        // Sample 1: e1 = 4. I = 4*0.5 = 2. D = 0 (first). u1 = 1*4 + 2*2 + 0 = 8.
        // Sample 2: e2 = 6. I = 2 + 6*0.5 = 5. D = (6-4)/0.5 = 4.
        //           u2 = 1*6 + 2*5 + 3*4 = 6 + 10 + 12 = 28.
        let cfg = PidConfig::new(1.0, 2.0, 3.0, 0.5).unwrap();
        let mut pid = Pid::new(cfg);
        let u1 = pid.update(4.0, 0.0);
        assert!((u1 - 8.0).abs() < EPS, "u1 = {u1}");
        let u2 = pid.update(6.0, 0.0);
        assert!((u2 - 28.0).abs() < EPS, "u2 = {u2}");
    }

    #[test]
    fn config_round_trips_through_json() {
        let cfg = PidConfig::new(1.0, 0.5, 0.05, 0.02)
            .unwrap()
            .with_limits(-1.0, 1.0)
            .unwrap();
        let json = serde_json::to_string(&cfg).unwrap();
        let back: PidConfig = serde_json::from_str(&json).unwrap();
        assert_eq!(cfg, back);
    }
}
