//! # valenx-vehicle
//!
//! Native **ground-vehicle (car) performance dynamics** — the "total
//! simulation" core of the Valenx car suite. Describe a car by its
//! performance-relevant parameters and read back the engineering answers a
//! vehicle dynamicist cares about: **0–100 km/h**, **top speed**, **braking
//! distance**, and **steady-state cornering** (max lateral g and skidpad
//! speed).
//!
//! ## The model
//!
//! A longitudinal/lateral **point-mass** model with the physics that actually
//! sets a car's numbers:
//!
//! - **Powertrain** as a peak-power limit: the tractive force the engine can
//!   deliver at speed `v` is `P·η / v` — so launch is power-rich and high
//!   speed is power-limited, exactly the shape of a real acceleration curve.
//! - **Tire friction circle**: traction is capped at `μ·N`, where `N` is the
//!   vertical load on the driven axle including **longitudinal weight
//!   transfer** (`m·a·h/L`, solved self-consistently with the acceleration it
//!   produces) and aerodynamic **downforce**. A rear-drive car squats and
//!   grips harder under power; a front-drive car unloads.
//! - **Aerodynamics**: drag `½ρv²·Cd·A` and downforce `½ρv²·ClA`, the latter
//!   pressing the tires down so grip rises with speed.
//! - **Rolling resistance** `Crr·m·g`.
//!
//! Top speed is the power-vs-drag balance; the sprint and braking numbers come
//! from integrating the available longitudinal acceleration; cornering comes
//! from the friction circle with speed-dependent downforce.
//!
//! ## Honest scope
//!
//! This is a **preliminary-design / research-grade** point-mass model. It
//! captures the dominant physics (power limit, grip limit, weight transfer,
//! aero) and lands on the right ballparks for real cars, but it is **not**:
//! a gear-by-gear powertrain with a real torque curve and shift logic; a
//! transient tire model (slip ratio/angle, Pacejka); a multi-body suspension
//! or chassis-flex model; or a commercial vehicle-dynamics suite (CarSim,
//! Adams/Car) or a driving game. Those are documented extensions on the way to
//! a fuller suite — the powertrain torque curve, a quasi-static lap simulator,
//! and the suspension/aero subsystems come next.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use serde::{Deserialize, Serialize};

/// Standard gravity (m/s²).
const G: f64 = 9.806_65;
/// Sea-level air density (kg/m³).
const RHO_AIR: f64 = 1.225;

/// Which axle(s) the engine drives — sets the share of vertical load available
/// for traction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Drivetrain {
    /// Front-wheel drive — the front axle drives (and unloads under accel).
    FrontWheel,
    /// Rear-wheel drive — the rear axle drives (and gains load under accel).
    RearWheel,
    /// All-wheel drive — all four tires put down power.
    AllWheel,
}

/// A car described by its performance-relevant parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Car {
    /// Total mass including driver (kg).
    pub mass: f64,
    /// Engine peak power (W).
    pub peak_power: f64,
    /// Driveline efficiency engine → wheels (0–1).
    pub drivetrain_efficiency: f64,
    /// Driven axle layout.
    pub drivetrain: Drivetrain,
    /// Tire–road friction coefficient `μ` (grip).
    pub tire_friction: f64,
    /// Static front weight fraction (0–1).
    pub weight_dist_front: f64,
    /// Wheelbase (m) — the lever arm for longitudinal weight transfer.
    pub wheelbase: f64,
    /// Centre-of-gravity height (m).
    pub cg_height: f64,
    /// Frontal area (m²).
    pub frontal_area: f64,
    /// Aerodynamic drag coefficient `Cd`.
    pub drag_coefficient: f64,
    /// Downforce coefficient–area product `Cl·A` (m², ≥ 0 = downforce).
    pub downforce_cla: f64,
    /// Rolling-resistance coefficient `Crr`.
    pub rolling_resistance: f64,
}

/// A 0-to-target acceleration result.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Sprint {
    /// Elapsed time (s).
    pub time: f64,
    /// Distance covered (m).
    pub distance: f64,
}

impl Car {
    /// Aerodynamic drag force at speed `v` (N).
    pub fn drag_force(&self, v: f64) -> f64 {
        0.5 * RHO_AIR * v * v * self.drag_coefficient * self.frontal_area
    }

    /// Aerodynamic downforce at speed `v` (N, ≥ 0).
    pub fn downforce(&self, v: f64) -> f64 {
        0.5 * RHO_AIR * v * v * self.downforce_cla
    }

    /// Rolling-resistance force (N) — first-order constant `Crr·m·g`.
    pub fn rolling_force(&self) -> f64 {
        self.rolling_resistance * self.mass * G
    }

    /// Vertical load on the driven axle (N) at longitudinal acceleration
    /// `a_long` and speed `v`: static distribution + dynamic weight transfer
    /// (`m·a·h/L`) + aerodynamic downforce.
    fn driven_axle_load(&self, a_long: f64, v: f64) -> f64 {
        let total = self.mass * G + self.downforce(v);
        let transfer = self.mass * a_long * self.cg_height / self.wheelbase.max(1e-6);
        let load = match self.drivetrain {
            Drivetrain::RearWheel => total * (1.0 - self.weight_dist_front) + transfer,
            Drivetrain::FrontWheel => total * self.weight_dist_front - transfer,
            Drivetrain::AllWheel => total,
        };
        load.max(0.0)
    }

    /// Maximum tractive force at the wheels at speed `v` (N) — the lesser of
    /// the power limit `P·η/v` and the traction limit `μ·N`, with the traction
    /// limit's weight transfer solved self-consistently with the acceleration
    /// it produces.
    pub fn tractive_force(&self, v: f64) -> f64 {
        let power_limited = if v > 0.1 {
            self.peak_power * self.drivetrain_efficiency / v
        } else {
            f64::INFINITY
        };
        // Fixed-point iterate: traction depends on a_long via weight transfer,
        // and a_long depends on the (capped) tractive force.
        let mut f = self.tire_friction * self.driven_axle_load(0.0, v);
        for _ in 0..10 {
            let f_eff = f.min(power_limited);
            let a = (f_eff - self.drag_force(v) - self.rolling_force()) / self.mass;
            f = self.tire_friction * self.driven_axle_load(a, v);
        }
        f.min(power_limited)
    }

    /// Net longitudinal acceleration available at speed `v` (m/s²) — tractive
    /// force minus drag and rolling resistance, over mass.
    pub fn acceleration_at(&self, v: f64) -> f64 {
        (self.tractive_force(v) - self.drag_force(v) - self.rolling_force()) / self.mass
    }

    /// Top speed (m/s) — the speed at which the available longitudinal
    /// acceleration falls to zero (power balanced by drag + rolling).
    pub fn top_speed(&self) -> f64 {
        let (mut lo, mut hi) = (0.0_f64, 200.0_f64); // 0 … 720 km/h bracket
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

    /// Time + distance to accelerate from rest to `target_speed` (m/s), by
    /// explicit integration of the available acceleration. Capped just below
    /// top speed so an unreachable target still returns.
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

    /// Braking distance from `v0` (m/s) to a stop (m) — friction-circle decel
    /// `μ·(m·g + downforce)/m` plus aero drag, integrated.
    pub fn braking_distance(&self, v0: f64) -> f64 {
        let dt = 1.0e-3;
        let (mut v, mut d) = (v0, 0.0_f64);
        while v > 0.0 && d < 2_000.0 {
            let grip = self.tire_friction * (self.mass * G + self.downforce(v)) / self.mass;
            let decel = grip + self.drag_force(v) / self.mass;
            v -= decel * dt;
            d += v.max(0.0) * dt;
        }
        d
    }

    /// Maximum lateral acceleration at speed `v`, in g — the friction circle
    /// with speed-dependent downforce: `μ·(m·g + downforce)/(m·g)`.
    pub fn max_lateral_g(&self, v: f64) -> f64 {
        self.tire_friction * (self.mass * G + self.downforce(v)) / (self.mass * G)
    }

    /// Steady-state cornering (skidpad) speed (m/s) for a corner of radius
    /// `radius` (m), solving `v²/R = a_lat(v)` self-consistently (downforce
    /// raises grip with speed).
    pub fn skidpad_speed(&self, radius: f64) -> f64 {
        let mut v = 10.0_f64;
        for _ in 0..100 {
            let a_lat = self.max_lateral_g(v) * G;
            v = (a_lat * radius.max(1e-6)).sqrt();
        }
        v
    }
}

/// A representative ~400 hp rear-drive sports car (≈ 1500 kg).
pub fn sports_car() -> Car {
    Car {
        mass: 1_500.0,
        peak_power: 300_000.0, // ~400 hp
        drivetrain_efficiency: 0.90,
        drivetrain: Drivetrain::RearWheel,
        tire_friction: 1.05,
        weight_dist_front: 0.45,
        wheelbase: 2.6,
        cg_height: 0.50,
        frontal_area: 2.2,
        drag_coefficient: 0.32,
        downforce_cla: 0.6,
        rolling_resistance: 0.012,
    }
}

/// A representative all-wheel-drive hypercar (≈ 1600 kg, ~1000 hp).
pub fn hypercar() -> Car {
    Car {
        mass: 1_600.0,
        peak_power: 750_000.0, // ~1000 hp
        drivetrain_efficiency: 0.92,
        drivetrain: Drivetrain::AllWheel,
        tire_friction: 1.20,
        weight_dist_front: 0.43,
        wheelbase: 2.7,
        cg_height: 0.46,
        frontal_area: 2.1,
        drag_coefficient: 0.35,
        downforce_cla: 2.5,
        rolling_resistance: 0.012,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// km/h → m/s.
    fn kmh(v: f64) -> f64 {
        v / 3.6
    }

    #[test]
    fn sports_car_hits_known_performance_ballparks() {
        let c = sports_car();
        let sprint = c.accelerate_to(kmh(100.0));
        assert!(
            (3.0..6.5).contains(&sprint.time),
            "0–100 km/h = {} s",
            sprint.time
        );
        let top_kmh = c.top_speed() * 3.6;
        assert!(
            (230.0..340.0).contains(&top_kmh),
            "top speed = {top_kmh} km/h"
        );
        let brake = c.braking_distance(kmh(100.0));
        assert!((28.0..48.0).contains(&brake), "100–0 km/h = {brake} m");
        let lat = c.max_lateral_g(kmh(80.0));
        assert!((0.95..1.40).contains(&lat), "lateral = {lat} g");
    }

    #[test]
    fn top_speed_is_the_power_drag_balance() {
        let c = sports_car();
        let vmax = c.top_speed();
        // At top speed the net acceleration is essentially zero.
        assert!(c.acceleration_at(vmax).abs() < 0.05, "a(vmax) should ≈ 0");
        // Just below it the car can still pull; just above it cannot.
        assert!(c.acceleration_at(vmax * 0.9) > 0.0);
        assert!(c.acceleration_at(vmax * 1.1) < 0.0);
    }

    #[test]
    fn more_power_accelerates_harder() {
        let base = sports_car();
        let mut strong = base;
        strong.peak_power *= 1.6;
        assert!(
            strong.accelerate_to(kmh(100.0)).time < base.accelerate_to(kmh(100.0)).time,
            "more power should cut the 0–100 time"
        );
        assert!(strong.top_speed() > base.top_speed());
    }

    #[test]
    fn awd_launches_at_least_as_hard_as_rwd_for_a_powerful_car() {
        // A high-power car is traction-limited off the line; putting all four
        // tires down (AWD) should not be slower than driving one axle (RWD).
        let mut rwd = hypercar();
        rwd.drivetrain = Drivetrain::RearWheel;
        let awd = hypercar(); // AllWheel
        assert!(
            awd.accelerate_to(kmh(100.0)).time <= rwd.accelerate_to(kmh(100.0)).time + 1e-6,
            "AWD {} should be ≤ RWD {}",
            awd.accelerate_to(kmh(100.0)).time,
            rwd.accelerate_to(kmh(100.0)).time
        );
    }

    #[test]
    fn downforce_raises_cornering_speed() {
        let base = sports_car();
        let mut winged = base;
        winged.downforce_cla = 3.0; // big wing
        assert!(
            winged.skidpad_speed(60.0) > base.skidpad_speed(60.0),
            "downforce should raise skidpad speed"
        );
        // And lateral grip should climb with speed when there is downforce.
        assert!(winged.max_lateral_g(kmh(200.0)) > winged.max_lateral_g(kmh(50.0)));
    }

    #[test]
    fn is_deterministic_and_finite() {
        let c = hypercar();
        let a = c.accelerate_to(kmh(100.0));
        let b = c.accelerate_to(kmh(100.0));
        assert_eq!(a, b);
        assert!(a.time.is_finite() && c.top_speed().is_finite());
        assert!(c.skidpad_speed(50.0).is_finite() && c.braking_distance(kmh(100.0)).is_finite());
    }
}
