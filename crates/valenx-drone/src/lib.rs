//! # valenx-drone
//!
//! Multirotor (drone / quadcopter) **hover performance** by momentum
//! (actuator-disk) theory.
//!
//! ## What
//!
//! Given a multirotor's rotor count, rotor radius, all-up mass, rotor
//! figure of merit and the air density, this crate computes:
//!
//! - the total rotor disk area `A = n * pi * R^2` ([`Multirotor::disk_area`])
//!   and the disk loading `T / A` ([`Multirotor::disk_loading`]);
//! - the hover thrust per rotor `T / n` ([`Multirotor::hover_thrust_per_rotor`]);
//! - the hover induced velocity `vi = sqrt(T / (2 * rho * A))`
//!   ([`Multirotor::induced_velocity_hover`]);
//! - the ideal hover power `P = T * vi = T^1.5 / sqrt(2 * rho * A)`
//!   ([`Multirotor::ideal_hover_power`]) and the actual shaft power
//!   `P / FM` ([`Multirotor::actual_hover_power`]); and
//! - the hover endurance from a battery's usable energy
//!   ([`Multirotor::hover_endurance_minutes`]) and the thrust-to-weight
//!   ratio for a given installed thrust ([`Multirotor::thrust_to_weight`]).
//!
//! ## Model
//!
//! In a steady hover the rotors must produce a total thrust equal to the
//! weight, `T = m * g`. Actuator-disk (momentum) theory models the rotor
//! system as an ideal disk of total area `A = n * pi * R^2` accelerating
//! air downward. Conservation of momentum gives the induced (downwash)
//! velocity at the disk
//!
//! `vi = sqrt(T / (2 * rho * A))`
//!
//! and the ideal induced power
//!
//! `P_ideal = T * vi = T^1.5 / sqrt(2 * rho * A)`.
//!
//! Real rotors need more than this; the **figure of merit** `FM` in
//! `(0, 1]` is the ratio of ideal to actual hover power, so the shaft power
//! is `P_actual = P_ideal / FM` (good model rotors are around 0.7). Using
//! the total thrust with the total disk area gives the correct system
//! power because the per-rotor powers sum back to `T^1.5 / sqrt(2 rho A)`.
//!
//! Hover endurance from a battery of usable energy `E` (watt-hours) is
//! `t = E / P_actual`, i.e. `60 * E / P_actual` minutes with `E` in Wh and
//! `P_actual` in watts.
//!
//! Units are SI: metres for `R`, kg for mass, kg/m^3 for `rho`. Results are
//! m^2 for area, N/m^2 (Pa) for disk loading, m/s for velocity and watts
//! for power.
//!
//! ## Honest scope
//!
//! Research/educational grade. This is ideal hover momentum theory for an
//! axial (vertical) hover only. It deliberately omits everything a real
//! flight-performance estimate needs, including but not limited to: forward
//! flight and climb, blade-element aerodynamics and tip losses beyond the
//! lumped figure of merit, propeller / motor / ESC efficiency curves and
//! their variation with RPM, battery voltage sag and usable-capacity
//! derating, rotor-to-rotor and ground-effect interference, control and
//! payload power, and any airframe structural or regulatory limit. It is
//! not flight-grade and is not a substitute for measured performance,
//! detailed analysis, or airworthiness certification.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

/// Standard gravity (m/s^2).
pub const GRAVITY: f64 = 9.80665;
/// Nominal sea-level air density (kg/m^3, ISA 15 C).
pub const SEA_LEVEL_AIR_DENSITY: f64 = 1.225;

/// An out-of-domain multirotor input.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum DroneError {
    /// A quantity that must be finite and strictly positive was not.
    #[error("{quantity} must be finite and positive, got {value}")]
    NonPositive {
        /// The offending quantity's name.
        quantity: &'static str,
        /// The offending value.
        value: f64,
    },
    /// The figure of merit was outside the physical range `(0, 1]`.
    #[error("figure of merit must be in (0, 1], got {value}")]
    FigureOfMerit {
        /// The offending value.
        value: f64,
    },
    /// The rotor count was zero.
    #[error("rotor count must be at least 1")]
    NoRotors,
}

/// Return `value` when it is finite and strictly positive, else an error.
fn require_positive(quantity: &'static str, value: f64) -> Result<f64, DroneError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(DroneError::NonPositive { quantity, value })
    }
}

/// A multirotor aircraft, validated on construction.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Multirotor {
    /// Number of lifting rotors `n` (>= 1).
    pub rotor_count: u32,
    /// Rotor radius `R` (m).
    pub rotor_radius_m: f64,
    /// All-up mass `m` (kg).
    pub mass_kg: f64,
    /// Rotor figure of merit `FM` in `(0, 1]` (ideal / actual hover power).
    pub figure_of_merit: f64,
    /// Air density `rho` (kg/m^3).
    pub air_density: f64,
}

impl Multirotor {
    /// Build a validated multirotor. The rotor count must be at least 1;
    /// rotor radius, mass and air density must be finite and positive; the
    /// figure of merit must lie in `(0, 1]`.
    ///
    /// # Errors
    ///
    /// Returns [`DroneError`] when any input is out of its physical domain.
    pub fn new(
        rotor_count: u32,
        rotor_radius_m: f64,
        mass_kg: f64,
        figure_of_merit: f64,
        air_density: f64,
    ) -> Result<Self, DroneError> {
        if rotor_count == 0 {
            return Err(DroneError::NoRotors);
        }
        let rotor_radius_m = require_positive("rotor radius", rotor_radius_m)?;
        let mass_kg = require_positive("mass", mass_kg)?;
        let air_density = require_positive("air density", air_density)?;
        if !(figure_of_merit.is_finite() && figure_of_merit > 0.0 && figure_of_merit <= 1.0) {
            return Err(DroneError::FigureOfMerit {
                value: figure_of_merit,
            });
        }
        Ok(Self {
            rotor_count,
            rotor_radius_m,
            mass_kg,
            figure_of_merit,
            air_density,
        })
    }

    /// All-up weight `W = m * g` (N) — the total thrust required to hover.
    pub fn weight(&self) -> f64 {
        self.mass_kg * GRAVITY
    }

    /// Total rotor disk area `A = n * pi * R^2` (m^2).
    pub fn disk_area(&self) -> f64 {
        f64::from(self.rotor_count) * PI * self.rotor_radius_m * self.rotor_radius_m
    }

    /// Disk loading `W / A` (N/m^2).
    pub fn disk_loading(&self) -> f64 {
        self.weight() / self.disk_area()
    }

    /// Hover thrust required per rotor `W / n` (N).
    pub fn hover_thrust_per_rotor(&self) -> f64 {
        self.weight() / f64::from(self.rotor_count)
    }

    /// Hover induced (downwash) velocity `vi = sqrt(W / (2 * rho * A))`
    /// (m/s).
    pub fn induced_velocity_hover(&self) -> f64 {
        (self.weight() / (2.0 * self.air_density * self.disk_area())).sqrt()
    }

    /// Ideal hover power `P = W * vi = W^1.5 / sqrt(2 * rho * A)` (W).
    pub fn ideal_hover_power(&self) -> f64 {
        self.weight() * self.induced_velocity_hover()
    }

    /// Actual hover shaft power `P_ideal / FM` (W).
    pub fn actual_hover_power(&self) -> f64 {
        self.ideal_hover_power() / self.figure_of_merit
    }

    /// Thrust-to-weight ratio for a given installed total thrust (N).
    pub fn thrust_to_weight(&self, installed_thrust_n: f64) -> f64 {
        installed_thrust_n / self.weight()
    }

    /// Hover endurance (minutes) from a battery of `usable_energy_wh`
    /// watt-hours: `60 * E / P_actual`.
    pub fn hover_endurance_minutes(&self, usable_energy_wh: f64) -> f64 {
        60.0 * usable_energy_wh / self.actual_hover_power()
    }

    /// Compute the full hover-performance report in one call.
    pub fn hover_performance(&self) -> HoverPerformance {
        HoverPerformance {
            weight_n: self.weight(),
            disk_area_m2: self.disk_area(),
            disk_loading_pa: self.disk_loading(),
            induced_velocity_m_s: self.induced_velocity_hover(),
            ideal_power_w: self.ideal_hover_power(),
            actual_power_w: self.actual_hover_power(),
        }
    }
}

/// A multirotor's computed hover performance.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct HoverPerformance {
    /// All-up weight `W` (N).
    pub weight_n: f64,
    /// Total rotor disk area `A` (m^2).
    pub disk_area_m2: f64,
    /// Disk loading `W / A` (N/m^2).
    pub disk_loading_pa: f64,
    /// Hover induced velocity `vi` (m/s).
    pub induced_velocity_m_s: f64,
    /// Ideal hover power (W).
    pub ideal_power_w: f64,
    /// Actual hover shaft power (W).
    pub actual_power_w: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-6 * b.abs().max(1.0)
    }

    /// A 1.5 kg quadcopter with 15 cm rotors at sea level.
    fn quad() -> Multirotor {
        Multirotor::new(4, 0.15, 1.5, 0.7, SEA_LEVEL_AIR_DENSITY).unwrap()
    }

    #[test]
    fn ideal_power_two_ways_agree() {
        // P = W*vi must equal the closed form W^1.5 / sqrt(2 rho A).
        let m = quad();
        let closed = m.weight().powf(1.5) / (2.0 * m.air_density * m.disk_area()).sqrt();
        assert!(close(m.ideal_hover_power(), closed));
    }

    #[test]
    fn hover_thrust_and_disk_quantities_match_hand_calc() {
        let m = quad();
        assert!(close(m.weight(), 1.5 * GRAVITY));
        // A = 4 * pi * 0.15^2.
        assert!(close(m.disk_area(), 4.0 * PI * 0.15 * 0.15));
        assert!(close(m.hover_thrust_per_rotor(), m.weight() / 4.0));
        assert!(close(m.disk_loading(), m.weight() / m.disk_area()));
    }

    #[test]
    fn induced_velocity_satisfies_momentum_thrust() {
        // T = 2 rho A vi^2 must recover the weight.
        let m = quad();
        let vi = m.induced_velocity_hover();
        let t = 2.0 * m.air_density * m.disk_area() * vi * vi;
        assert!(close(t, m.weight()));
    }

    #[test]
    fn figure_of_merit_scales_actual_power() {
        // FM = 1 -> actual == ideal; FM = 0.5 -> actual == 2 * ideal.
        let ideal = quad().ideal_hover_power();
        let perfect = Multirotor::new(4, 0.15, 1.5, 1.0, SEA_LEVEL_AIR_DENSITY).unwrap();
        let half = Multirotor::new(4, 0.15, 1.5, 0.5, SEA_LEVEL_AIR_DENSITY).unwrap();
        assert!(close(perfect.actual_hover_power(), ideal));
        assert!(close(half.actual_hover_power(), 2.0 * ideal));
    }

    #[test]
    fn bigger_rotors_lower_the_hover_power() {
        // Same mass, larger disk -> lower disk loading -> less induced power.
        let small = Multirotor::new(4, 0.12, 1.5, 0.7, SEA_LEVEL_AIR_DENSITY).unwrap();
        let big = Multirotor::new(4, 0.25, 1.5, 0.7, SEA_LEVEL_AIR_DENSITY).unwrap();
        assert!(big.actual_hover_power() < small.actual_hover_power());
    }

    #[test]
    fn endurance_is_energy_over_power() {
        let m = quad();
        let wh = 50.0;
        let expected = 60.0 * wh / m.actual_hover_power();
        assert!(close(m.hover_endurance_minutes(wh), expected));
        // Twice the battery, twice the endurance.
        assert!(close(
            m.hover_endurance_minutes(2.0 * wh),
            2.0 * m.hover_endurance_minutes(wh)
        ));
    }

    #[test]
    fn thrust_to_weight_ratio() {
        let m = quad();
        // Installed thrust equal to twice the weight gives T/W = 2.
        assert!((m.thrust_to_weight(2.0 * m.weight()) - 2.0).abs() < EPS);
    }

    #[test]
    fn rejects_out_of_domain_inputs() {
        assert!(Multirotor::new(0, 0.15, 1.5, 0.7, SEA_LEVEL_AIR_DENSITY).is_err());
        assert!(Multirotor::new(4, 0.0, 1.5, 0.7, SEA_LEVEL_AIR_DENSITY).is_err());
        assert!(Multirotor::new(4, 0.15, -1.0, 0.7, SEA_LEVEL_AIR_DENSITY).is_err());
        assert!(Multirotor::new(4, 0.15, 1.5, 1.5, SEA_LEVEL_AIR_DENSITY).is_err());
        assert!(Multirotor::new(4, 0.15, 1.5, 0.0, SEA_LEVEL_AIR_DENSITY).is_err());
    }

    #[test]
    fn hover_performance_is_serde_round_trippable() {
        let h = quad().hover_performance();
        let json = serde_json::to_string(&h).unwrap();
        let back: HoverPerformance = serde_json::from_str(&json).unwrap();
        // serde_json default float parse is ~1 ULP approximate -> tolerance
        // compare, not bit / string compare.
        assert!(close(back.weight_n, h.weight_n));
        assert!(close(back.disk_area_m2, h.disk_area_m2));
        assert!(close(back.induced_velocity_m_s, h.induced_velocity_m_s));
        assert!(close(back.ideal_power_w, h.ideal_power_w));
        assert!(close(back.actual_power_w, h.actual_power_w));
    }
}
