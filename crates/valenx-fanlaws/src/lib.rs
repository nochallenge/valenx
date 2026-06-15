//! # valenx-fanlaws
//!
//! Closed-form fan and centrifugal-blower **affinity laws** ("fan
//! laws") plus the fundamental fan **power identity**.
//!
//! ## What
//!
//! Given one operating point of a single, geometrically fixed fan, this
//! crate computes how its volumetric flow, pressure rise, and shaft
//! power scale when you change the impeller speed and/or the gas
//! density, and how shaft power relates to flow, pressure, and total
//! efficiency. Concretely it provides:
//!
//! - [`scale_flow`], [`scale_pressure`], [`scale_power`] — scale a
//!   single quantity between two speeds.
//! - [`correct_pressure_for_density`], [`correct_power_for_density`] —
//!   density correction at fixed speed.
//! - [`OperatingPoint`] + [`scale_operating_point`] — apply the full
//!   speed-and-density transform to a complete operating point.
//! - [`air_power`], [`shaft_power`], [`implied_efficiency`] and the
//!   validated [`Efficiency`] type — the `P = Q * dP / eta` identity
//!   and its inverse.
//!
//! ## Model
//!
//! For a fixed fan, with `r = n2 / n1` the speed ratio and
//! `d = rho2 / rho1` the density ratio:
//!
//! ```text
//! flow:      Q2 = Q1 * r                 (linear in speed)
//! pressure:  P2 = P1 * r^2 * d           (square of speed, linear in density)
//! power:     W2 = W1 * r^3 * d           (cube of speed, linear in density)
//! air power: P_air  = Q * dP
//! shaft:     P_shaft = Q * dP / eta,  eta in (0, 1]
//! ```
//!
//! The affinity exponents (1, 2, 3) follow from dimensional similarity:
//! flow scales with tip speed, pressure with the dynamic head
//! (proportional to speed squared and to density), and power with their
//! product. Efficiency is constrained to the half-open interval
//! `(0, 1]` so the power identity is always finite and shaft power
//! never falls below air power.
//!
//! ## Honest scope
//!
//! Research/educational grade. These are textbook closed-form
//! similarity relations under idealising assumptions: the fan's total
//! efficiency is treated as constant across the speed change (valid
//! only over modest speed ratios where Reynolds-number, tip-clearance,
//! and compressibility effects stay small), flow is taken as
//! incompressible so density enters as a simple linear multiplier, and
//! the fan is assumed to ride a fixed quadratic system-resistance curve
//! so the operating point stays on one similarity ray. The crate is
//! dimension-agnostic and performs no unit checking — you must keep
//! units self-consistent across the two points. It is NOT a clinical,
//! medical, or production engineering tool: do not use it for
//! ventilator design, HVAC code-compliance, or fan-selection sign-off
//! without an independent, qualified engineering review.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod affinity;
pub mod error;
pub mod power;

pub use affinity::{
    correct_power_for_density, correct_pressure_for_density, scale_flow, scale_operating_point,
    scale_power, scale_pressure, OperatingPoint,
};
pub use error::{ErrorCategory, FanLawError};
pub use power::{air_power, implied_efficiency, shaft_power, Efficiency};

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for analytic float comparisons.
    const EPS: f64 = 1e-9;

    #[test]
    fn end_to_end_fan_speed_up_with_power_check() {
        // A fan delivers Q=8 m^3/s at dP=600 Pa, drawing W shaft at
        // 60% efficiency, running at 1000 rev/min in 1.2 kg/m^3 air.
        // Step it to 1500 rev/min (r = 1.5) at the same density.
        let eta = Efficiency::new(0.6).unwrap();
        let p0 = OperatingPoint::new(
            8.0,
            600.0,
            shaft_power(8.0, 600.0, eta).unwrap(),
            1000.0,
            1.2,
        )
        .unwrap();

        let p1 = scale_operating_point(&p0, 1500.0, 1.2).unwrap();

        // Flow up by 1.5, pressure by 2.25, power by 3.375.
        assert!((p1.flow - 12.0).abs() < EPS);
        assert!((p1.pressure - 1350.0).abs() < 1e-7);
        assert!((p1.power - p0.power * 1.5f64.powi(3)).abs() < 1e-6);

        // The new shaft power, recomputed from the scaled flow and
        // pressure at the SAME efficiency, must agree with the affinity
        // result — the two routes are consistent.
        let p1_from_identity = shaft_power(p1.flow, p1.pressure, eta).unwrap();
        assert!(
            (p1.power - p1_from_identity).abs() < 1e-6,
            "affinity power {affinity} vs identity {identity}",
            affinity = p1.power,
            identity = p1_from_identity
        );
    }

    #[test]
    fn public_reexports_are_wired() {
        // Smoke test that the flattened re-exports resolve.
        assert!((scale_flow(1.0, 1.0, 2.0).unwrap() - 2.0).abs() < EPS);
        let eta = Efficiency::from_percent(50.0).unwrap();
        assert!((shaft_power(1.0, 1.0, eta).unwrap() - 2.0).abs() < EPS);
        assert_eq!(
            FanLawError::EfficiencyOutOfRange { value: 2.0 }.category(),
            ErrorCategory::Config
        );
    }
}
