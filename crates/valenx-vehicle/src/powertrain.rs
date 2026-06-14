//! **Electric powertrain + energy/range** model — the EV "total simulation"
//! layer on top of the [`crate::Car`] chassis.
//!
//! A combustion car is power-limited everywhere; an EV is **torque-limited**
//! off the line (the motor delivers near-constant peak torque up to a base
//! speed) and **power-limited** above it. This module models that envelope, a
//! battery (usable energy + a discharge-power ceiling), a single-speed
//! reduction, regenerative braking, and the steady-cruise energy consumption
//! that sets **range**.
//!
//! It reuses the [`crate::Car`] for the rolling chassis: aerodynamic drag /
//! downforce, rolling resistance, and the friction-circle traction limit with
//! longitudinal weight transfer (a `Car` whose `peak_power` is set effectively
//! infinite returns exactly the traction limit, which we cap the motor force
//! against). Acceleration is then `(min(motor force, traction) − drag −
//! rolling) / m`.
//!
//! Validated against Tesla-class numbers (0–100 km/h ≈ 3–5 s, an electronically
//! limited top speed, and ~400–600 km highway range).
//!
//! Honest scope: a point-mass EV performance/energy model — research /
//! preliminary-design grade. No battery thermal/SOC-dependent power fade, no
//! motor-efficiency map, no drive-cycle (WLTP/EPA) integration, no inverter
//! detail. Those are the documented next steps.

use serde::{Deserialize, Serialize};

use crate::{Car, Sprint};

/// An electric drive unit: a torque-then-power envelope through a single-speed
/// reduction to the wheels.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ElectricPowertrain {
    /// Motor/inverter peak power (W).
    pub peak_power: f64,
    /// Motor peak torque (N·m) — held from zero up to the base speed.
    pub peak_torque: f64,
    /// Single-speed reduction ratio (motor rev per wheel rev).
    pub gear_ratio: f64,
    /// Driven-wheel radius (m).
    pub wheel_radius: f64,
    /// Driveline efficiency motor → wheels (0–1).
    pub efficiency: f64,
}

impl ElectricPowertrain {
    /// Peak wheel force from the torque limit (N): `T·ratio·η / r`.
    fn torque_limited_force(&self) -> f64 {
        self.peak_torque * self.gear_ratio * self.efficiency / self.wheel_radius
    }
}

/// A traction battery.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Battery {
    /// Usable energy (kWh).
    pub usable_kwh: f64,
    /// Maximum discharge power (W) — caps the powertrain.
    pub max_discharge_power: f64,
}

/// An electric vehicle: a rolling chassis ([`Car`]) plus an electric
/// powertrain, battery, accessory load and regenerative braking.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct EvCar {
    /// Rolling chassis (mass, aero, tires, weight transfer). Its `peak_power`
    /// is set effectively infinite so [`Car::tractive_force`] returns the pure
    /// friction-circle traction limit.
    pub chassis: Car,
    /// Electric drive unit.
    pub powertrain: ElectricPowertrain,
    /// Traction battery.
    pub battery: Battery,
    /// Accessory/HVAC load drawn while driving (W).
    pub aux_power: f64,
    /// Fraction of braking kinetic energy recovered by regen (0–1).
    pub regen_efficiency: f64,
    /// Electronically limited top speed (m/s), if any.
    pub top_speed_limit: Option<f64>,
}

impl EvCar {
    /// Effective power ceiling — the lesser of the motor and battery limits (W).
    pub fn available_power(&self) -> f64 {
        self.powertrain
            .peak_power
            .min(self.battery.max_discharge_power)
    }

    /// Wheel tractive force at speed `v` (N): the motor envelope (torque-limited
    /// then power-limited) capped by the chassis traction limit.
    pub fn tractive_force(&self, v: f64) -> f64 {
        let power_limited = if v > 0.1 {
            self.available_power() * self.powertrain.efficiency / v
        } else {
            f64::INFINITY
        };
        let motor = self.powertrain.torque_limited_force().min(power_limited);
        motor.min(self.chassis.tractive_force(v))
    }

    /// Net longitudinal acceleration available at speed `v` (m/s²).
    pub fn acceleration_at(&self, v: f64) -> f64 {
        (self.tractive_force(v) - self.chassis.drag_force(v) - self.chassis.rolling_force())
            / self.chassis.mass
    }

    /// Top speed the powertrain+aero balance allows, ignoring any limiter (m/s).
    pub fn unrestricted_top_speed(&self) -> f64 {
        let (mut lo, mut hi) = (0.0_f64, 150.0_f64);
        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            if self.acceleration_at(mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// Top speed (m/s) — the unrestricted top, capped by the electronic limiter.
    pub fn top_speed(&self) -> f64 {
        let phys = self.unrestricted_top_speed();
        self.top_speed_limit.map_or(phys, |l| phys.min(l))
    }

    /// Time + distance from rest to `target_speed` (m/s), by integration.
    pub fn accelerate_to(&self, target_speed: f64) -> Sprint {
        let dt = 1.0e-3;
        let target = target_speed.min(self.top_speed() * 0.999);
        let (mut v, mut t, mut d) = (0.0_f64, 0.0_f64, 0.0_f64);
        while v < target && t < 120.0 {
            let a = self.acceleration_at(v);
            if a <= 1e-4 {
                break;
            }
            v += a * dt;
            d += v * dt;
            t += dt;
        }
        Sprint {
            time: t,
            distance: d,
        }
    }

    /// Battery power drawn to cruise at constant speed `v` (W): the resistive
    /// power over the driveline, plus accessories.
    pub fn cruise_battery_power(&self, v: f64) -> f64 {
        let resist = self.chassis.drag_force(v) + self.chassis.rolling_force();
        resist * v / self.powertrain.efficiency + self.aux_power
    }

    /// Steady-cruise energy consumption at speed `v` (Wh/km).
    pub fn consumption_wh_per_km(&self, v: f64) -> f64 {
        self.cruise_battery_power(v) / (v * 3.6)
    }

    /// Range (km) cruising at constant speed `v`.
    pub fn range_km(&self, v: f64) -> f64 {
        self.battery.usable_kwh * 1000.0 / self.consumption_wh_per_km(v)
    }

    /// Kinetic energy regen can recover braking from `v0` to rest (kWh).
    pub fn regen_recovered_kwh(&self, v0: f64) -> f64 {
        self.regen_efficiency * 0.5 * self.chassis.mass * v0 * v0 / 3.6e6
    }
}

/// A Tesla-class dual-motor EV sedan.
pub fn tesla() -> EvCar {
    let chassis = Car {
        mass: 2_100.0,
        peak_power: 1.0e12, // effectively ∞ ⇒ chassis.tractive_force = traction limit
        drivetrain_efficiency: 0.90,
        drivetrain: crate::Drivetrain::AllWheel,
        tire_friction: 1.00,
        weight_dist_front: 0.48,
        wheelbase: 2.96,
        cg_height: 0.46,
        frontal_area: 2.22,
        drag_coefficient: 0.23,
        downforce_cla: 0.0,
        rolling_resistance: 0.012,
    };
    EvCar {
        chassis,
        powertrain: ElectricPowertrain {
            peak_power: 300_000.0, // ~400 hp combined
            peak_torque: 640.0,
            gear_ratio: 9.0,
            wheel_radius: 0.34,
            efficiency: 0.90,
        },
        battery: Battery {
            usable_kwh: 75.0,
            max_discharge_power: 320_000.0,
        },
        aux_power: 500.0,
        regen_efficiency: 0.60,
        top_speed_limit: Some(250.0 / 3.6), // electronically limited
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn kmh(v: f64) -> f64 {
        v / 3.6
    }

    #[test]
    fn tesla_class_ev_hits_known_numbers() {
        let ev = tesla();
        let sprint = ev.accelerate_to(kmh(100.0));
        assert!(
            (2.5..6.0).contains(&sprint.time),
            "0–100 km/h = {} s",
            sprint.time
        );
        let top_kmh = ev.top_speed() * 3.6;
        assert!((230.0..270.0).contains(&top_kmh), "top = {top_kmh} km/h");
        // The limiter bites: unrestricted physics top is higher than the cap.
        assert!(
            ev.unrestricted_top_speed() > ev.top_speed() + 0.5,
            "limiter should leave headroom"
        );
        let range = ev.range_km(kmh(100.0));
        assert!(
            (380.0..680.0).contains(&range),
            "highway range = {range} km"
        );
    }

    #[test]
    fn launch_is_torque_limited_then_power_limited() {
        let ev = tesla();
        // Off the line the motor torque envelope (not the battery power) sets the
        // force, so a low-speed force ≈ T·ratio·η/r.
        let f0 = ev.tractive_force(1.0);
        let expected = 640.0 * 9.0 * 0.9 / 0.34;
        assert!((f0 - expected).abs() / expected < 0.05, "launch force {f0}");
        // High up, force has dropped (power-limited).
        assert!(ev.tractive_force(60.0) < f0, "force should fall with speed");
    }

    #[test]
    fn bigger_battery_extends_range_proportionally() {
        let base = tesla();
        let mut big = tesla();
        big.battery.usable_kwh *= 2.0;
        let r = big.range_km(kmh(100.0)) / base.range_km(kmh(100.0));
        assert!(
            (r - 2.0).abs() < 1e-6,
            "double the kWh ⇒ ~double range, got {r}"
        );
    }

    #[test]
    fn regen_recovers_a_fraction_of_kinetic_energy() {
        let ev = tesla();
        let v0 = kmh(100.0);
        let ke_kwh = 0.5 * ev.chassis.mass * v0 * v0 / 3.6e6;
        let regen = ev.regen_recovered_kwh(v0);
        assert!(
            regen > 0.0 && regen < ke_kwh,
            "regen {regen} vs KE {ke_kwh}"
        );
        assert!((regen - 0.60 * ke_kwh).abs() < 1e-9);
    }

    #[test]
    fn is_deterministic() {
        let ev = tesla();
        assert_eq!(ev.accelerate_to(kmh(100.0)), ev.accelerate_to(kmh(100.0)));
        assert_eq!(ev.range_km(kmh(110.0)), ev.range_km(kmh(110.0)));
    }
}
