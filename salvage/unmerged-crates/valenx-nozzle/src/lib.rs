//! # valenx-nozzle
//!
//! Isentropic compressible flow through a converging-diverging (de Laval)
//! nozzle, in closed form.
//!
//! ## What
//!
//! Given a Mach number `M` and the ratio of specific heats `gamma`, this
//! crate evaluates the textbook one-dimensional isentropic-flow relations:
//! the stagnation property ratios `T0/T`, `p0/p` and `rho0/rho`; the
//! area-Mach ratio `A/A*`; the choked critical-pressure ratio `p*/p0`; and a
//! dimensionless mass-flow parameter. It also inverts the area-Mach relation
//! (sub- and supersonic branches) and the pressure ratio back to `M`.
//!
//! ## Model
//!
//! Steady, adiabatic, reversible (isentropic) flow of a calorically perfect
//! gas with constant `gamma`. Writing `g` for `gamma`, the governing
//! equations are:
//!
//! ```text
//! T0/T   = 1 + (g-1)/2 * M^2
//! p0/p   = (T0/T)^( g/(g-1) )
//! rho0/rho = (T0/T)^( 1/(g-1) )
//! A/A*   = (1/M) * [ (2/(g+1)) (1 + (g-1)/2 M^2) ] ^ ( (g+1)/(2(g-1)) )
//! p*/p0  = (2/(g+1)) ^ ( g/(g-1) )
//! ```
//!
//! `A/A*` has a global minimum of `1` at the sonic throat `M = 1`, so the
//! nozzle is choked there and every `A/A* > 1` has one subsonic and one
//! supersonic solution.
//!
//! ## Honest scope
//!
//! Research/educational grade. These are standard textbook closed-form
//! models, validated in the test suite against analytic ground truth
//! (`M = 1` gives `A/A* = 1`; the `gamma = 1.4` air table values for `T0/T`,
//! `p0/p`, `A/A*`). It is NOT a clinical/medical or production-certified
//! engineering tool: there is no modelling of friction, heat transfer,
//! boundary layers, shocks, real-gas effects, two-phase flow, nor any
//! component tolerances, safety factors, or thermal/structural limits. Do
//! not use it to certify hardware.
//!
//! ## Example
//!
//! ```
//! use valenx_nozzle::{area_ratio, pressure_ratio, temperature_ratio, mach_from_area_ratio, Branch};
//!
//! let gamma = 1.4; // air
//!
//! // At the throat the flow is sonic and A/A* = 1 exactly.
//! assert!((area_ratio(1.0, gamma).unwrap() - 1.0).abs() < 1e-12);
//!
//! // Standard air-table values at M = 2.
//! assert!((temperature_ratio(2.0, gamma).unwrap() - 1.8).abs() < 1e-12);
//! assert!((pressure_ratio(2.0, gamma).unwrap() - 7.824_449).abs() < 1e-5);
//! assert!((area_ratio(2.0, gamma).unwrap() - 1.687_500).abs() < 1e-6);
//!
//! // Recover the supersonic Mach from an expansion area ratio of 1.6875.
//! let m = mach_from_area_ratio(1.6875, gamma, Branch::Supersonic).unwrap();
//! assert!((m - 2.0).abs() < 1e-6);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod area;
pub mod error;
pub mod ratios;

pub use area::{area_ratio, mach_from_area_ratio, mass_flow_parameter, Branch};
pub use error::{NozzleError, Result};
pub use ratios::{
    critical_pressure_ratio, density_ratio, mach_from_pressure_ratio, pressure_ratio,
    temperature_ratio,
};

#[cfg(test)]
mod tests {
    use super::*;

    const GAMMA_AIR: f64 = 1.4;

    // ----- stagnation ratios vs hand-worked / table ground truth ----------

    #[test]
    fn temperature_ratio_at_mach_zero_is_one() {
        let r = temperature_ratio(0.0, GAMMA_AIR).unwrap();
        assert!((r - 1.0).abs() < 1e-15, "got {r}");
    }

    #[test]
    fn temperature_ratio_air_table_values() {
        // T0/T = 1 + 0.2 M^2 for air.
        // M = 1 -> 1.2; M = 2 -> 1.8; M = 3 -> 2.8.
        assert!((temperature_ratio(1.0, GAMMA_AIR).unwrap() - 1.2).abs() < 1e-12);
        assert!((temperature_ratio(2.0, GAMMA_AIR).unwrap() - 1.8).abs() < 1e-12);
        assert!((temperature_ratio(3.0, GAMMA_AIR).unwrap() - 2.8).abs() < 1e-12);
    }

    #[test]
    fn pressure_ratio_air_table_values() {
        // NACA 1135 / Anderson table: p0/p for air.
        // M = 1 -> 1.892929; M = 2 -> 7.824449; M = 0.5 -> 1.186212.
        assert!((pressure_ratio(0.5, GAMMA_AIR).unwrap() - 1.186_212).abs() < 1e-5);
        assert!((pressure_ratio(1.0, GAMMA_AIR).unwrap() - 1.892_929).abs() < 1e-5);
        assert!((pressure_ratio(2.0, GAMMA_AIR).unwrap() - 7.824_449).abs() < 1e-5);
    }

    #[test]
    fn density_ratio_air_table_values() {
        // rho0/rho for air: M = 1 -> 1.577439; M = 2 -> 4.346916.
        assert!((density_ratio(1.0, GAMMA_AIR).unwrap() - 1.577_439).abs() < 1e-5);
        assert!((density_ratio(2.0, GAMMA_AIR).unwrap() - 4.346_916).abs() < 1e-5);
    }

    #[test]
    fn ratios_satisfy_perfect_gas_identity() {
        // The state equation links the three ratios:
        //   (p0/p) = (rho0/rho) * (T0/T).
        for &m in &[0.3_f64, 0.8, 1.0, 1.5, 2.5] {
            let p = pressure_ratio(m, GAMMA_AIR).unwrap();
            let rho = density_ratio(m, GAMMA_AIR).unwrap();
            let t = temperature_ratio(m, GAMMA_AIR).unwrap();
            assert!((p - rho * t).abs() < 1e-9, "m = {m}: {p} vs {}", rho * t);
        }
    }

    #[test]
    fn critical_pressure_ratio_air() {
        // p*/p0 = 0.528282 for air; equals 1 / (p0/p at M = 1).
        let crit = critical_pressure_ratio(GAMMA_AIR).unwrap();
        assert!((crit - 0.528_282).abs() < 1e-6, "got {crit}");
        let inv = 1.0 / pressure_ratio(1.0, GAMMA_AIR).unwrap();
        assert!((crit - inv).abs() < 1e-12);
    }

    #[test]
    fn critical_pressure_ratio_monatomic() {
        // gamma = 5/3: (2/(8/3))^(2.5) = 0.75^2.5 = 0.487139.
        let crit = critical_pressure_ratio(5.0 / 3.0).unwrap();
        assert!((crit - 0.487_139).abs() < 1e-6, "got {crit}");
    }

    // ----- area-Mach relation vs ground truth -----------------------------

    #[test]
    fn area_ratio_unity_at_sonic_throat() {
        // The defining ground-truth case across all gammas.
        for &g in &[1.2_f64, 1.4, 5.0 / 3.0] {
            let r = area_ratio(1.0, g).unwrap();
            assert!((r - 1.0).abs() < 1e-12, "gamma = {g}: got {r}");
        }
    }

    #[test]
    fn area_ratio_air_table_values() {
        // A/A* for air: M = 2 -> 1.687500; M = 3 -> 4.234568;
        // M = 0.5 -> 1.339844.
        assert!((area_ratio(0.5, GAMMA_AIR).unwrap() - 1.339_844).abs() < 1e-6);
        assert!((area_ratio(2.0, GAMMA_AIR).unwrap() - 1.687_500).abs() < 1e-6);
        assert!((area_ratio(3.0, GAMMA_AIR).unwrap() - 4.234_568).abs() < 1e-6);
    }

    #[test]
    fn area_ratio_is_minimised_at_throat() {
        // A/A* > 1 just either side of the sonic point.
        let throat = area_ratio(1.0, GAMMA_AIR).unwrap();
        assert!(area_ratio(0.95, GAMMA_AIR).unwrap() > throat);
        assert!(area_ratio(1.05, GAMMA_AIR).unwrap() > throat);
    }

    // ----- inversion: round-trips and table values ------------------------

    #[test]
    fn invert_area_ratio_supersonic_branch() {
        // A/A* = 1.6875 (air) -> M = 2 on the supersonic branch.
        let m = mach_from_area_ratio(1.6875, GAMMA_AIR, Branch::Supersonic).unwrap();
        assert!((m - 2.0).abs() < 1e-6, "got {m}");
    }

    #[test]
    fn invert_area_ratio_subsonic_branch() {
        // A/A* = 1.339844 (air) -> M = 0.5 on the subsonic branch.
        let m = mach_from_area_ratio(1.339_844, GAMMA_AIR, Branch::Subsonic).unwrap();
        assert!((m - 0.5).abs() < 1e-5, "got {m}");
    }

    #[test]
    fn invert_area_ratio_round_trips_both_branches() {
        for &m in &[0.1_f64, 0.4, 0.7, 0.95] {
            let ar = area_ratio(m, GAMMA_AIR).unwrap();
            let back = mach_from_area_ratio(ar, GAMMA_AIR, Branch::Subsonic).unwrap();
            assert!((back - m).abs() < 1e-6, "subsonic m = {m}: back = {back}");
        }
        for &m in &[1.05_f64, 1.5, 2.0, 3.5, 5.0] {
            let ar = area_ratio(m, GAMMA_AIR).unwrap();
            let back = mach_from_area_ratio(ar, GAMMA_AIR, Branch::Supersonic).unwrap();
            assert!((back - m).abs() < 1e-6, "supersonic m = {m}: back = {back}");
        }
    }

    #[test]
    fn invert_area_ratio_unity_gives_sonic() {
        let sub = mach_from_area_ratio(1.0, GAMMA_AIR, Branch::Subsonic).unwrap();
        let sup = mach_from_area_ratio(1.0, GAMMA_AIR, Branch::Supersonic).unwrap();
        assert!((sub - 1.0).abs() < 1e-12);
        assert!((sup - 1.0).abs() < 1e-12);
    }

    #[test]
    fn invert_pressure_ratio_round_trip() {
        for &m in &[0.2_f64, 0.6, 1.0, 1.8, 3.0] {
            let p_over_p0 = 1.0 / pressure_ratio(m, GAMMA_AIR).unwrap();
            let back = mach_from_pressure_ratio(p_over_p0, GAMMA_AIR).unwrap();
            assert!((back - m).abs() < 1e-9, "m = {m}: back = {back}");
        }
    }

    #[test]
    fn invert_pressure_ratio_unity_is_stagnation() {
        // p/p0 = 1 means the gas is at rest: M = 0.
        let m = mach_from_pressure_ratio(1.0, GAMMA_AIR).unwrap();
        assert!(m.abs() < 1e-12, "got {m}");
    }

    // ----- mass-flow parameter --------------------------------------------

    #[test]
    fn mass_flow_parameter_peaks_at_choke() {
        let throat = mass_flow_parameter(1.0, GAMMA_AIR).unwrap();
        for &m in &[0.1_f64, 0.5, 0.9, 1.1, 2.0, 4.0] {
            let v = mass_flow_parameter(m, GAMMA_AIR).unwrap();
            assert!(v < throat, "m = {m}: {v} should be < choke {throat}");
        }
    }

    #[test]
    fn mass_flow_parameter_zero_at_rest() {
        let v = mass_flow_parameter(0.0, GAMMA_AIR).unwrap();
        assert!(v.abs() < 1e-15, "got {v}");
    }

    #[test]
    fn choked_mass_flow_parameter_air_value() {
        // MFP at M = 1 for air:
        //   sqrt(1.4) * (1.2)^(-3) = 1.183216 * 0.578704 = 0.684731.
        let v = mass_flow_parameter(1.0, GAMMA_AIR).unwrap();
        assert!((v - 0.684_731).abs() < 1e-6, "got {v}");
    }

    // ----- error / validation paths ---------------------------------------

    #[test]
    fn rejects_negative_mach() {
        let e = temperature_ratio(-0.5, GAMMA_AIR).unwrap_err();
        assert_eq!(e.code(), "negative_mach");
    }

    #[test]
    fn rejects_non_finite_inputs() {
        assert_eq!(
            area_ratio(f64::NAN, GAMMA_AIR).unwrap_err().code(),
            "not_finite"
        );
        assert_eq!(
            temperature_ratio(f64::INFINITY, GAMMA_AIR)
                .unwrap_err()
                .code(),
            "not_finite"
        );
        assert_eq!(
            pressure_ratio(1.0, f64::NAN).unwrap_err().code(),
            "not_finite"
        );
    }

    #[test]
    fn rejects_out_of_range_gamma() {
        assert_eq!(
            temperature_ratio(1.0, 1.0).unwrap_err().code(),
            "gamma_out_of_range"
        );
        assert_eq!(
            temperature_ratio(1.0, 0.9).unwrap_err().code(),
            "gamma_out_of_range"
        );
        assert_eq!(
            temperature_ratio(1.0, 2.0).unwrap_err().code(),
            "gamma_out_of_range"
        );
    }

    #[test]
    fn area_ratio_rejects_zero_mach() {
        // A/A* is infinite at M = 0; rejected rather than returning inf.
        assert!(area_ratio(0.0, GAMMA_AIR).is_err());
    }

    #[test]
    fn invert_area_ratio_rejects_below_choke() {
        let e = mach_from_area_ratio(0.5, GAMMA_AIR, Branch::Supersonic).unwrap_err();
        assert_eq!(e.code(), "area_ratio_below_choke");
    }

    #[test]
    fn invert_pressure_ratio_rejects_out_of_range() {
        assert_eq!(
            mach_from_pressure_ratio(0.0, GAMMA_AIR).unwrap_err().code(),
            "pressure_ratio_out_of_range"
        );
        assert_eq!(
            mach_from_pressure_ratio(1.5, GAMMA_AIR).unwrap_err().code(),
            "pressure_ratio_out_of_range"
        );
    }

    #[test]
    fn error_codes_are_stable_strings() {
        // Guard against an accidental rename of a public code.
        let e = NozzleError::NotConverged {
            residual: 1.0,
            iters: 7,
        };
        assert_eq!(e.code(), "not_converged");
    }
}
