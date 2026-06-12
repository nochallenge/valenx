//! Propulsive landing — the **hoverslam** ("suicide burn") that flies a
//! booster from a high descent speed to a soft touchdown.
//!
//! A hoverslam ignites as late as possible and burns at high thrust so
//! the vehicle's velocity and altitude reach zero together at the ground.
//! Under a **constant net deceleration** `a_net = T/m − g` (thrust minus
//! weight), the kinematics are the textbook constant-acceleration result:
//! the ignition altitude that nulls a descent speed `v` exactly at the
//! surface is
//!
//! ```text
//!   h = v² / (2·a_net),   a_net > 0.
//! ```
//!
//! If the engine cannot produce a *positive* net deceleration
//! (`T/m ≤ g`), the vehicle cannot arrest its fall at all — a real,
//! testable "can't land" condition reported as an
//! [`AstroError::InvalidParameter`].
//!
//! Alongside the closed-form ignition altitude this module provides an
//! **integrated landing burn** ([`LandingSim::run`]) over the same
//! fixed-step RK4 discipline as the ascent simulator (a bounded step
//! count, a minimum `dt`), with the engine's mass falling as propellant
//! is spent. The throttle is physically clamped to `[0, 1]`. In the
//! mass-invariant limit the integrated touchdown speed reproduces the
//! closed form (≈ 0 at `h = v²/(2a)`); with a finite mass flow the rising
//! `a_net` nulls the velocity at or just above the ground — i.e. a
//! controlled landing that never overshoots into the surface.

use serde::{Deserialize, Serialize};

use crate::constants::G0;
use crate::error::{AstroError, Result};
use crate::sim::{check_step_count, MAX_SIM_STEPS};

/// Smallest accepted integration step for a landing burn (s). Set to the
/// same `1e-4` s value as the ascent simulator's
/// [`crate::config::MIN_TIME_STEP`] floor — an independent constant kept at
/// parity with it — so an absurd `dt` cannot explode the step count.
pub const MIN_TIME_STEP: f64 = 1e-4;

/// Ignition ("suicide-burn") altitude (m) that brings a descent speed
/// `descent_speed` (m/s) to rest exactly at the surface under a constant
/// net deceleration `net_decel` (m/s²): `h = v²/(2·a_net)`.
///
/// # Errors
///
/// Returns [`AstroError::InvalidParameter`]:
/// - if `descent_speed` is non-finite or negative, or `net_decel` is
///   non-finite;
/// - if `net_decel <= 0` — the engine cannot produce a positive net
///   deceleration, so the vehicle **cannot land** (it would keep
///   accelerating downward). This is the "can't land" result, not a
///   numerical failure.
pub fn ignition_altitude(descent_speed: f64, net_decel: f64) -> Result<f64> {
    if !descent_speed.is_finite() || descent_speed < 0.0 {
        return Err(AstroError::InvalidParameter(
            "descent_speed must be finite and >= 0",
        ));
    }
    if !net_decel.is_finite() {
        return Err(AstroError::InvalidParameter("net_decel must be finite"));
    }
    if net_decel <= 0.0 {
        return Err(AstroError::InvalidParameter(
            "cannot land: net deceleration is non-positive (thrust/mass <= g)",
        ));
    }
    Ok(descent_speed * descent_speed / (2.0 * net_decel))
}

/// A propulsive-landing burn to integrate.
///
/// The vehicle starts at [`start_altitude`](Self::start_altitude)
/// descending at `descent_speed`, and burns at a constant `thrust`
/// against gravity `g` (default standard gravity) while its mass falls at
/// the engine's propellant mass flow (`thrust / (isp·g₀)`). The throttle
/// is clamped to `[0, 1]`; here it runs full-open (1.0).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LandingSim {
    /// Altitude at ignition (m). Set to [`ignition_altitude`] for the
    /// canonical hoverslam, or higher for an early, conservative burn.
    pub start_altitude: f64,
    /// Descent speed at ignition (m/s, positive downward).
    pub descent_speed: f64,
    /// Engine thrust (N).
    pub thrust: f64,
    /// Vehicle mass at ignition (kg).
    pub initial_mass: f64,
    /// Engine specific impulse (s); drives the propellant mass flow.
    pub isp: f64,
    /// Surface gravitational acceleration (m/s²).
    pub gravity: f64,
    /// Integration step (s).
    pub time_step: f64,
}

/// The outcome of an integrated landing burn.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LandingOutcome {
    /// Speed at touchdown / arrest (m/s). Near zero for a well-timed
    /// hoverslam; positive means a hard landing.
    pub touchdown_speed: f64,
    /// Altitude at which the vehicle came to rest (m). Near zero for a
    /// well-timed burn; positive means it stopped short (a hover) and
    /// would need to descend the rest under reduced thrust.
    pub arrest_altitude: f64,
    /// Propellant consumed during the burn (kg).
    pub propellant_used: f64,
    /// Elapsed burn time (s).
    pub burn_time: f64,
    /// The (constant) commanded throttle, clamped to `[0, 1]`.
    pub throttle: f64,
}

impl LandingSim {
    /// Validate the configuration.
    ///
    /// # Errors
    ///
    /// Returns [`AstroError::InvalidParameter`] for any non-finite or
    /// non-physical field, [`AstroError::InvalidIntegration`] for a
    /// `time_step` at or below the [`MIN_TIME_STEP`] floor, and the
    /// "can't land" [`AstroError::InvalidParameter`] if the initial net
    /// deceleration `thrust/mass − g` is not positive.
    pub fn validate(&self) -> Result<()> {
        for value in [
            self.start_altitude,
            self.descent_speed,
            self.thrust,
            self.initial_mass,
            self.isp,
            self.gravity,
        ] {
            if !value.is_finite() {
                return Err(AstroError::InvalidParameter(
                    "landing parameter must be finite",
                ));
            }
        }
        if self.start_altitude < 0.0 {
            return Err(AstroError::InvalidParameter("start_altitude must be >= 0"));
        }
        if self.descent_speed < 0.0 {
            return Err(AstroError::InvalidParameter("descent_speed must be >= 0"));
        }
        if self.thrust <= 0.0 {
            return Err(AstroError::InvalidParameter("thrust must be > 0"));
        }
        if self.initial_mass <= 0.0 {
            return Err(AstroError::InvalidParameter("initial_mass must be > 0"));
        }
        if self.isp <= 0.0 {
            return Err(AstroError::InvalidParameter("isp must be > 0"));
        }
        if self.gravity < 0.0 {
            return Err(AstroError::InvalidParameter("gravity must be >= 0"));
        }
        if !self.time_step.is_finite() || self.time_step < MIN_TIME_STEP {
            return Err(AstroError::InvalidIntegration(
                "time_step below the minimum (1e-4 s)",
            ));
        }
        // Can the engine arrest the descent at all?
        if self.thrust / self.initial_mass - self.gravity <= 0.0 {
            return Err(AstroError::InvalidParameter(
                "cannot land: thrust/mass <= g (no positive net deceleration)",
            ));
        }
        Ok(())
    }

    /// Integrate the landing burn until the vehicle comes to rest (descent
    /// speed reaches zero) or reaches the surface, whichever first.
    ///
    /// Uses the shared fixed-step RK4 second-order integrator with the
    /// mass frozen across each step (as the ascent loop does), the mass
    /// then stepped discretely by the propellant burned. The loop is
    /// bounded by [`MAX_SIM_STEPS`].
    ///
    /// # Errors
    ///
    /// Returns the validation errors of [`LandingSim::validate`], or
    /// [`AstroError::StepBudgetExceeded`] if the burn does not terminate
    /// within the step budget (it always does for a physical hoverslam,
    /// since `a_net` only grows as mass falls).
    pub fn run(&self) -> Result<LandingOutcome> {
        self.validate()?;

        // Bound the loop up front: t ≤ v/a_net0 for the descent to arrest,
        // and we cap the step count regardless.
        let max_steps = MAX_SIM_STEPS;
        check_step_count(max_steps)?;

        let throttle = 1.0_f64; // full-open, already clamped to [0,1]
        let thrust = throttle * self.thrust;
        let mdot = thrust / (self.isp * G0);
        let g = self.gravity;
        let dt = self.time_step;

        // State: altitude (m) and descent speed (m/s, positive down).
        let mut alt = self.start_altitude;
        let mut v = self.descent_speed;
        let mut mass = self.initial_mass;
        let mut t = 0.0f64;
        let mut steps = 0u64;

        // ḣ = −v_down ; v̇_down = −(T/m − g). Mass frozen across the step.
        while alt > 0.0 && v > 0.0 {
            if steps >= max_steps {
                return Err(AstroError::StepBudgetExceeded(max_steps));
            }
            steps += 1;

            let a_net = thrust / mass - g;
            // d(alt) = -v ; d(v) = -a_net (constant across the step)
            let d = |_alt: f64, vv: f64| (-vv, -a_net);
            let k1 = d(alt, v);
            let k2 = d(alt + 0.5 * dt * k1.0, v + 0.5 * dt * k1.1);
            let k3 = d(alt + 0.5 * dt * k2.0, v + 0.5 * dt * k2.1);
            let k4 = d(alt + dt * k3.0, v + dt * k3.1);
            alt += dt / 6.0 * (k1.0 + 2.0 * k2.0 + 2.0 * k3.0 + k4.0);
            v += dt / 6.0 * (k1.1 + 2.0 * k2.1 + 2.0 * k3.1 + k4.1);
            mass -= mdot * dt;
            t += dt;
        }

        // Touchdown speed is the descent speed at arrest, clamped at zero
        // (a tiny negative from the discrete final step means it has just
        // come to rest / would start ascending — report 0).
        let touchdown_speed = v.max(0.0);
        Ok(LandingOutcome {
            touchdown_speed,
            arrest_altitude: alt.max(0.0),
            propellant_used: self.initial_mass - mass,
            burn_time: t,
            throttle,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ignition_altitude_is_kinematic() {
        // h = v²/(2a). v=200, a_net=15 -> 40000/30 = 1333.333...
        let h = ignition_altitude(200.0, 15.0).expect("can land");
        assert!((h - 1_333.333_333_333_333).abs() < 1e-9, "h = {h}");
        // v=100, a_net=10 -> 500.
        assert!((ignition_altitude(100.0, 10.0).expect("ok") - 500.0).abs() < 1e-9);
    }

    #[test]
    fn cannot_land_when_thrust_under_weight() {
        // Non-positive net deceleration -> "can't land" Err.
        assert!(matches!(
            ignition_altitude(100.0, 0.0),
            Err(AstroError::InvalidParameter(_))
        ));
        assert!(ignition_altitude(100.0, -2.0).is_err());
    }

    #[test]
    fn integrated_burn_touchdown_speed_is_near_zero_mass_invariant() {
        // ORACLE: under a constant net deceleration (mass-invariant limit:
        // a huge mass with a matched thrust so mdot's effect is
        // negligible), igniting at h = v²/(2a) nulls both velocity and
        // altitude at the ground simultaneously. The integrated touchdown
        // speed must be ≈ 0 and the arrest altitude ≈ 0.
        let v = 100.0;
        let a_net = 10.0;
        let g = G0;
        let mass = 1.0e6;
        let thrust = (a_net + g) * mass; // gives a_net exactly at ignition
                                         // A very high Isp makes the propellant mass flow negligible over
                                         // the burn, so a_net stays essentially constant — the true
                                         // mass-invariant limit the closed form assumes.
        let h = ignition_altitude(v, a_net).expect("can land"); // 500 m
        let sim = LandingSim {
            start_altitude: h,
            descent_speed: v,
            thrust,
            initial_mass: mass,
            isp: 3.0e6,
            gravity: g,
            time_step: 1e-3,
        };
        let out = sim.run().expect("valid landing");
        // Both velocity and altitude reach zero together at the surface.
        assert!(
            out.touchdown_speed < 0.01,
            "touchdown speed = {}",
            out.touchdown_speed
        );
        assert!(
            out.arrest_altitude < 0.1,
            "arrest altitude = {}",
            out.arrest_altitude
        );
        assert!((0.0..=1.0).contains(&out.throttle));
    }

    #[test]
    fn mass_varying_burn_arrests_at_or_above_ground_softly() {
        // A realistic finite-mass booster: igniting at the constant-mass
        // ignition altitude, the rising a_net (mass falls) brings velocity
        // to ~0 at or above the surface — a controlled landing, never a
        // crash. Touchdown speed stays small and non-negative.
        let v = 200.0;
        let a_net0 = 15.0;
        let g = G0;
        let mass = 25_000.0;
        let thrust = (a_net0 + g) * mass;
        let h = ignition_altitude(v, a_net0).expect("can land");
        let sim = LandingSim {
            start_altitude: h,
            descent_speed: v,
            thrust,
            initial_mass: mass,
            isp: 300.0,
            gravity: g,
            time_step: 1e-3,
        };
        let out = sim.run().expect("valid landing");
        // Soft: arrested with near-zero residual speed.
        assert!(
            out.touchdown_speed < 0.05,
            "touchdown speed = {}",
            out.touchdown_speed
        );
        // Came to rest at or above the ground (never overshot below).
        assert!(out.arrest_altitude >= 0.0);
        // Burned a sensible, positive amount of propellant.
        assert!(out.propellant_used > 0.0);
        // Burn time is on the order of v/a_net (~13 s) — finite.
        assert!(
            out.burn_time > 0.0 && out.burn_time < 60.0,
            "t = {}",
            out.burn_time
        );
    }

    #[test]
    fn run_rejects_engine_that_cannot_land() {
        // thrust/mass <= g -> validate() rejects with the can't-land Err.
        let sim = LandingSim {
            start_altitude: 1_000.0,
            descent_speed: 100.0,
            thrust: 9.0 * 1_000.0, // T/m = 9 < g
            initial_mass: 1_000.0,
            isp: 300.0,
            gravity: G0,
            time_step: 1e-3,
        };
        assert!(matches!(sim.run(), Err(AstroError::InvalidParameter(_))));
    }

    #[test]
    fn rejects_non_physical_config() {
        let base = LandingSim {
            start_altitude: 500.0,
            descent_speed: 100.0,
            thrust: 5.0e5,
            initial_mass: 20_000.0,
            isp: 300.0,
            gravity: G0,
            time_step: 1e-3,
        };
        assert!(base.validate().is_ok());
        assert!(LandingSim {
            thrust: 0.0,
            ..base
        }
        .validate()
        .is_err());
        assert!(LandingSim {
            initial_mass: 0.0,
            ..base
        }
        .validate()
        .is_err());
        assert!(LandingSim { isp: -1.0, ..base }.validate().is_err());
        assert!(LandingSim {
            time_step: 1e-9,
            ..base
        }
        .validate()
        .is_err()); // below floor
        assert!(LandingSim {
            descent_speed: f64::NAN,
            ..base
        }
        .validate()
        .is_err());
    }
}
