//! # valenx-conveyor
//!
//! Closed-form steady-state models for a belt conveyor handling bulk
//! material.
//!
//! ## What
//!
//! Given a belt's bulk-material load (density, cross-sectional load
//! area, speed) and its drive geometry (incline, effective tension),
//! this crate computes the quantities a conveyor sizing exercise needs:
//!
//! - mass flow rate `mdot = rho * A * v` and volumetric throughput
//!   `Q_v = A * v`,
//! - rated capacity in tonnes per hour (which scales linearly with belt
//!   speed and load area),
//! - lift power against an incline `P_lift = mdot * g * h`,
//! - total drive power from belt tension `P = F * v`, and a transparent
//!   friction-plus-lift power breakdown.
//!
//! ## Model
//!
//! The material on the belt is treated as a moving prism of constant
//! cross-section `A` translated at speed `v`; continuity then gives the
//! mass flow `rho * A * v`. Lift power is the rate of work raising that
//! mass over a vertical rise `h`. Drive power is the effective belt
//! tension times belt speed. For a pure (frictionless) lift the two
//! power statements coincide; with friction present, `F` carries the
//! extra resistance and `P = F * v` is the general form. Every formula
//! is exercised against hand-computed ground-truth values in the unit
//! tests.
//!
//! ## Honest scope
//!
//! Research/educational grade. The relations here are introductory
//! textbook algebra (CEMA / ISO 5048 style at the first-principles
//! level). They deliberately do not model the detailed resistance
//! build-up of a real installation — idler rolling resistance, belt
//! flexure, skirtboard and pulley drag, take-up, drive efficiency or
//! material surcharge geometry — all of which a certified conveyor-design
//! package folds into the effective tension. This crate is NOT a
//! clinical/medical/production engineering tool and is not a substitute
//! for a certified design package or a professional engineer's sign-off.
//!
//! ## Surface
//!
//! - [`belt::Belt`] — a transported stream (`density`, `load_area`,
//!   `speed`) with [`belt::Belt::mass_flow`],
//!   [`belt::Belt::volumetric_flow`] and [`belt::Belt::capacity_tph`].
//! - [`belt::mass_flow_rate`] — standalone `rho * A * v`.
//! - [`power::lift_power`] / [`power::lift_power_g`] — `mdot * g * h`.
//! - [`power::drive_power`] — `F * v`.
//! - [`power::incline_rise`] — `L * sin(angle)`.
//! - [`power::PowerBreakdown`] — friction-plus-lift power resolution.
//! - [`error::ConveyorError`] — validated-input error taxonomy.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod belt;
pub mod error;
pub mod power;

pub use belt::{mass_flow_rate, Belt, G};
pub use error::{ConveyorError, ErrorCategory};
pub use power::{drive_power, incline_rise, lift_power, lift_power_g, PowerBreakdown};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn belt_round_trips_through_json() {
        let belt = Belt::new(1500.0, 0.2, 2.0).expect("valid belt");
        let json = serde_json::to_string(&belt).expect("serialize");
        let back: Belt = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(belt, back);
        assert!((back.mass_flow() - 600.0).abs() < 1e-9);
    }

    #[test]
    fn power_breakdown_round_trips_through_json() {
        let bd = PowerBreakdown::solve(600.0, 2.0, 100.0, 0.0, 5000.0).expect("valid");
        let json = serde_json::to_string(&bd).expect("serialize");
        let back: PowerBreakdown = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(bd, back);
    }

    #[test]
    fn end_to_end_capacity_and_power() {
        // A 1.2 m wide flat belt, 25% filled, sand at 1600 kg/m^3,
        // running at 2.5 m/s up a 100 m / 30 deg incline.
        let belt = Belt::from_flat_fill(1600.0, 1.2, 0.25, 2.5).expect("valid belt");
        // A = 0.25 * 1.2^2 = 0.36 m^2; mdot = 1600 * 0.36 * 2.5 = 1440 kg/s.
        assert!((belt.mass_flow() - 1440.0).abs() < 1e-6);
        // Capacity = 1440 * 3.6 = 5184 t/h.
        assert!((belt.capacity_tph() - 5184.0).abs() < 1e-3);

        let h = incline_rise(100.0, std::f64::consts::FRAC_PI_6).expect("valid"); // 50 m
        let p_lift = lift_power(belt.mass_flow(), h).expect("valid");
        // 1440 * 9.80665 * 50 = 706078.8 W.
        assert!((p_lift - 1440.0 * G * 50.0).abs() < 1e-3);
    }
}
