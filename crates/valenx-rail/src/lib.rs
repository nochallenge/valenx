//! # valenx-rail
//!
//! Rail-vehicle (train) **running resistance and tractive performance** by
//! the classical Davis equation.
//!
//! ## What
//!
//! Given a train's mass, its three Davis resistance coefficients and a
//! constant tractive effort, this crate computes:
//!
//! - the running (rolling + flange + aerodynamic) resistance
//!   `R(v) = A + B*v + C*v^2` ([`Train::running_resistance`]);
//! - the grade resistance `m*g*grade` ([`Train::grade_resistance`]) and the
//!   total resistance ([`Train::total_resistance`]);
//! - the net tractive force `TE - R - grade` ([`Train::net_tractive_force`])
//!   and the resulting acceleration ([`Train::acceleration`]);
//! - the drawbar power needed to hold a speed on a grade
//!   ([`Train::drawbar_power`]);
//! - the steepest grade the train can start on from rest
//!   ([`Train::max_startable_grade`]); and
//! - the balancing (steady) speed on a grade ([`Train::balancing_speed`]),
//!   where tractive effort exactly equals the total resistance.
//!
//! ## Model
//!
//! Running resistance uses the Davis (1926) quadratic
//!
//! `R(v) = A + B*v + C*v^2`
//!
//! with `A` a speed-independent term (bearing / rolling / flange), `B*v` a
//! term linear in speed and `C*v^2` the aerodynamic term. Grade resistance
//! is the along-track component of weight, `R_g = m*g*grade`, using the
//! small-angle convention where `grade` is the rise-over-run fraction
//! (positive uphill, negative downhill). The net accelerating force is
//!
//! `F_net = TE - R(v) - R_g`,
//!
//! and the acceleration is `F_net / m`. At rest the train can just start on
//! a grade when `TE = A + m*g*grade`, giving the maximum startable grade
//! `(TE - A) / (m*g)`. The balancing speed solves `TE = R(v) + R_g`, i.e.
//! the positive root of `C*v^2 + B*v + (A + R_g - TE) = 0`.
//!
//! Units are SI: kg for mass, newtons for `A` and `TE`, N*s/m for `B`,
//! N*s^2/m^2 for `C`, m/s for speed; results are newtons, m/s^2 and watts.
//!
//! ## Honest scope
//!
//! Research/educational grade. Resistance is the lumped Davis quadratic and
//! tractive effort is treated as **constant** — a real locomotive is
//! power-limited, so its tractive effort falls roughly as `P/v` above its
//! base speed. The constant-TE [`Train::balancing_speed`] is therefore only
//! physical in the low-speed, roughly-constant-TE region; at high speed a
//! real train's balancing speed is set by its power, not by this formula.
//! It also omits curve resistance, tunnel resistance, rotational-inertia
//! mass factors, adhesion (wheel-slip) limits, and any signalling, braking
//! or timetable consideration. It is not a substitute for a railway
//! traction simulator or certified rolling-stock engineering.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use serde::{Deserialize, Serialize};

/// Standard gravity (m/s^2).
pub const GRAVITY: f64 = 9.80665;

/// An out-of-domain train input.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum RailError {
    /// A quantity that must be finite and strictly positive was not.
    #[error("{quantity} must be finite and positive, got {value}")]
    NonPositive {
        /// The offending quantity's name.
        quantity: &'static str,
        /// The offending value.
        value: f64,
    },
    /// A quantity that must be finite and non-negative was not.
    #[error("{quantity} must be finite and non-negative, got {value}")]
    Negative {
        /// The offending quantity's name.
        quantity: &'static str,
        /// The offending value.
        value: f64,
    },
}

/// Return `value` when finite and strictly positive, else an error.
fn require_positive(quantity: &'static str, value: f64) -> Result<f64, RailError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(RailError::NonPositive { quantity, value })
    }
}

/// Return `value` when finite and non-negative, else an error.
fn require_non_negative(quantity: &'static str, value: f64) -> Result<f64, RailError> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(RailError::Negative { quantity, value })
    }
}

/// A train, validated on construction.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Train {
    /// Total mass `m` (kg).
    pub mass_kg: f64,
    /// Davis constant term `A` (N), the speed-independent resistance.
    pub davis_a_n: f64,
    /// Davis linear coefficient `B` (N*s/m).
    pub davis_b_n_s_per_m: f64,
    /// Davis quadratic (aerodynamic) coefficient `C` (N*s^2/m^2).
    pub davis_c_n_s2_per_m2: f64,
    /// Constant tractive effort `TE` (N).
    pub max_tractive_effort_n: f64,
}

impl Train {
    /// Build a validated train. Mass, the aerodynamic coefficient `C` and
    /// the tractive effort must be finite and positive; `A` and `B` must be
    /// finite and non-negative.
    ///
    /// # Errors
    ///
    /// Returns [`RailError`] when any input is out of its physical domain.
    pub fn new(
        mass_kg: f64,
        davis_a_n: f64,
        davis_b_n_s_per_m: f64,
        davis_c_n_s2_per_m2: f64,
        max_tractive_effort_n: f64,
    ) -> Result<Self, RailError> {
        let mass_kg = require_positive("mass", mass_kg)?;
        let davis_a_n = require_non_negative("Davis A", davis_a_n)?;
        let davis_b_n_s_per_m = require_non_negative("Davis B", davis_b_n_s_per_m)?;
        let davis_c_n_s2_per_m2 = require_positive("Davis C", davis_c_n_s2_per_m2)?;
        let max_tractive_effort_n = require_positive("tractive effort", max_tractive_effort_n)?;
        Ok(Self {
            mass_kg,
            davis_a_n,
            davis_b_n_s_per_m,
            davis_c_n_s2_per_m2,
            max_tractive_effort_n,
        })
    }

    /// Weight `m*g` (N).
    pub fn weight(&self) -> f64 {
        self.mass_kg * GRAVITY
    }

    /// Davis running resistance `R(v) = A + B*v + C*v^2` (N) at speed `v`
    /// (m/s, clamped to non-negative).
    pub fn running_resistance(&self, speed_m_s: f64) -> f64 {
        let v = speed_m_s.max(0.0);
        self.davis_a_n + self.davis_b_n_s_per_m * v + self.davis_c_n_s2_per_m2 * v * v
    }

    /// Grade resistance `m*g*grade` (N), where `grade` is the rise/run
    /// fraction (positive uphill). Negative on a downgrade.
    pub fn grade_resistance(&self, grade: f64) -> f64 {
        self.weight() * grade
    }

    /// Total resistance opposing motion, `R(v) + m*g*grade` (N).
    pub fn total_resistance(&self, speed_m_s: f64, grade: f64) -> f64 {
        self.running_resistance(speed_m_s) + self.grade_resistance(grade)
    }

    /// Net tractive (accelerating) force `TE - R(v) - m*g*grade` (N); may be
    /// negative when resistance exceeds the tractive effort.
    pub fn net_tractive_force(&self, speed_m_s: f64, grade: f64) -> f64 {
        self.max_tractive_effort_n - self.total_resistance(speed_m_s, grade)
    }

    /// Acceleration `F_net / m` (m/s^2) at a speed and grade.
    pub fn acceleration(&self, speed_m_s: f64, grade: f64) -> f64 {
        self.net_tractive_force(speed_m_s, grade) / self.mass_kg
    }

    /// Drawbar power to hold a steady speed on a grade,
    /// `total_resistance(v, grade) * v` (W).
    pub fn drawbar_power(&self, speed_m_s: f64, grade: f64) -> f64 {
        self.total_resistance(speed_m_s, grade) * speed_m_s.max(0.0)
    }

    /// Steepest grade (rise/run fraction) the train can start on from rest,
    /// where standing tractive effort balances `A` plus grade resistance:
    /// `(TE - A) / (m*g)`.
    pub fn max_startable_grade(&self) -> f64 {
        (self.max_tractive_effort_n - self.davis_a_n) / self.weight()
    }

    /// Balancing (steady) speed (m/s) on a grade — the positive root of
    /// `C*v^2 + B*v + (A + m*g*grade - TE) = 0`. `None` when the train
    /// cannot reach a positive steady speed (tractive effort never exceeds
    /// the resistance). Constant-TE model — see the crate honest scope.
    pub fn balancing_speed(&self, grade: f64) -> Option<f64> {
        let c = self.davis_c_n_s2_per_m2;
        let b = self.davis_b_n_s_per_m;
        let c0 = self.davis_a_n + self.grade_resistance(grade) - self.max_tractive_effort_n;
        let disc = b * b - 4.0 * c * c0;
        if disc < 0.0 {
            return None;
        }
        let v = (-b + disc.sqrt()) / (2.0 * c);
        if v > 0.0 && v.is_finite() {
            Some(v)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-6 * b.abs().max(1.0)
    }

    /// A heavy freight train (about 4000 t) with strong tractive effort.
    fn freight() -> Train {
        Train::new(4_000_000.0_f64, 30000.0, 350.0, 12.0, 600000.0).unwrap()
    }

    #[test]
    fn resistance_at_rest_is_the_constant_term_and_grows_with_speed() {
        let t = freight();
        assert!(close(t.running_resistance(0.0), 30000.0));
        assert!(t.running_resistance(30.0) > t.running_resistance(10.0));
        // R(10) = 30000 + 350*10 + 12*100 = 34700.
        assert!(close(t.running_resistance(10.0), 34700.0));
    }

    #[test]
    fn balancing_speed_zeroes_the_net_force() {
        let t = freight();
        for &grade in &[0.0_f64, 0.005, 0.01] {
            let v = t
                .balancing_speed(grade)
                .expect("freight reaches a balancing speed");
            assert!(v > 0.0);
            assert!(
                t.net_tractive_force(v, grade).abs() < 1.0,
                "net force at balancing speed should be ~0, grade {grade}"
            );
        }
    }

    #[test]
    fn max_startable_grade_zeroes_the_standing_net_force() {
        let t = freight();
        let g = t.max_startable_grade();
        // At rest on exactly that grade, net force is zero.
        assert!(t.net_tractive_force(0.0, g).abs() < 1.0);
        // A steeper grade cannot be started (negative net force at rest).
        assert!(t.net_tractive_force(0.0, g + 0.005) < 0.0);
    }

    #[test]
    fn a_grade_too_steep_to_climb_has_no_balancing_speed() {
        let t = freight();
        let g = t.max_startable_grade() + 0.01;
        assert!(t.balancing_speed(g).is_none());
    }

    #[test]
    fn downgrade_raises_the_balancing_speed_over_level() {
        let t = freight();
        let level = t.balancing_speed(0.0).unwrap();
        let down = t.balancing_speed(-0.005).unwrap();
        assert!(down > level);
    }

    #[test]
    fn acceleration_is_net_force_over_mass() {
        let t = freight();
        let expected = t.net_tractive_force(5.0, 0.01) / t.mass_kg;
        assert!(close(t.acceleration(5.0, 0.01), expected));
    }

    #[test]
    fn drawbar_power_is_resistance_times_speed() {
        let t = freight();
        let v = 20.0;
        assert!(close(
            t.drawbar_power(v, 0.0),
            t.total_resistance(v, 0.0) * v
        ));
    }

    #[test]
    fn rejects_out_of_domain_inputs() {
        assert!(Train::new(0.0, 30000.0, 350.0, 12.0, 600000.0).is_err());
        assert!(Train::new(4_000_000.0, -1.0, 350.0, 12.0, 600000.0).is_err());
        assert!(Train::new(4_000_000.0, 30000.0, 350.0, 0.0, 600000.0).is_err());
        assert!(Train::new(4_000_000.0, 30000.0, 350.0, 12.0, 0.0).is_err());
    }

    #[test]
    fn train_is_serde_round_trippable() {
        let t = freight();
        let json = serde_json::to_string(&t).unwrap();
        let back: Train = serde_json::from_str(&json).unwrap();
        // Tolerance compare (serde_json float parse is ~1 ULP approximate).
        assert!(close(back.mass_kg, t.mass_kg));
        assert!(close(back.davis_c_n_s2_per_m2, t.davis_c_n_s2_per_m2));
        assert!(close(back.max_tractive_effort_n, t.max_tractive_effort_n));
    }
}
