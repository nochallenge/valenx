//! # valenx-gasdynamics — one-dimensional compressible-flow relations
//!
//! Closed-form gas-dynamics relations for a calorically-perfect ideal
//! gas: the isentropic stagnation ratios and area-Mach relation, and the
//! stationary normal-shock jumps. Pure `f64` algebra over two inputs — a
//! Mach number `M` and the specific-heat ratio `gamma` — with no
//! iteration, no external process, and no platform dependency.
//!
//! ## What
//!
//! - **Isentropic relations** ([`isentropic`]) — the stagnation-to-static
//!   ratios `T0/T` ([`isentropic::temperature_ratio`]), `p0/p`
//!   ([`isentropic::pressure_ratio`]) and `rho0/rho`
//!   ([`isentropic::density_ratio`]), bundled by [`stagnation_ratios`],
//!   plus the area-Mach relation `A/A*` ([`area_mach_ratio`]).
//! - **Normal shock** ([`mod@normal_shock`]) — the jump relations across
//!   a stationary normal shock as a function of the upstream Mach number
//!   `M1` and `gamma`: the downstream Mach `M2`
//!   ([`normal_shock::downstream_mach`]), the static ratios `p2/p1`,
//!   `T2/T1`, `rho2/rho1`, and the stagnation-pressure ratio `p02/p01`,
//!   bundled by [`normal_shock()`] into a [`NormalShock`].
//! - **Errors** ([`error`]) — a [`thiserror`]-backed [`GasError`] with
//!   validated constructors and stable [`GasError::code`] /
//!   [`GasError::category`] accessors.
//!
//! ## Model
//!
//! The gas is calorically perfect (constant specific heats, fixed
//! `gamma = cp / cv > 1`) and ideal (`p = rho R T`). Under those
//! assumptions the relations are the standard textbook results
//! (Anderson, *Modern Compressible Flow*; NACA Report 1135):
//!
//! ```text
//! Isentropic:  T0/T = 1 + (gamma-1)/2 * M^2
//!              p0/p = (T0/T)^( gamma/(gamma-1) )
//!              rho0/rho = (T0/T)^( 1/(gamma-1) )
//!              A/A* = (1/M) * [ (2/(gamma+1))(1 + (gamma-1)/2 M^2) ]
//!                              ^ ( (gamma+1)/(2(gamma-1)) )
//!
//! Normal shock (M1 >= 1):
//!              M2^2 = (1 + (gamma-1)/2 M1^2) / (gamma M1^2 - (gamma-1)/2)
//!              p2/p1 = 1 + 2 gamma/(gamma+1) (M1^2 - 1)
//!              rho2/rho1 = (gamma+1)M1^2 / ((gamma-1)M1^2 + 2)
//!              T2/T1 = (p2/p1)/(rho2/rho1)
//!              p02/p01 <= 1   (the entropy-rise signature)
//! ```
//!
//! Limiting behaviour the tests pin down: every stagnation ratio is one
//! at `M = 0` and rises monotonically; `A/A* = 1` exactly at `M = 1`;
//! across a shock `M1 > 1` gives `M2 < 1`, `p2 > p1`, and `p02/p01 < 1`;
//! and the Rankine-Hugoniot entropy change is non-negative for every
//! admissible shock.
//!
//! ## Honest scope
//!
//! Research / educational grade. This crate implements textbook
//! closed-form analytic models for a *single* perfect gas with constant
//! specific heats. It deliberately ignores real-gas effects (variable
//! `cp`, dissociation, vibrational excitation, condensation), viscosity
//! and heat conduction, finite shock thickness, oblique / curved shocks,
//! and heat or mass addition. It is intended for learning and first-order
//! estimates — it is **NOT a clinical/medical or production engineering
//! tool** and must not be used as the basis for certified design.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod fanno;
pub mod isentropic;
pub mod normal_shock;
pub mod oblique_shock;
pub mod prandtl_meyer;
pub mod rayleigh;

pub use error::{ErrorCategory, GasError, Result};
pub use fanno::{fanno_state, FannoState};
pub use isentropic::{
    area_mach_ratio, density_ratio, pressure_ratio, stagnation_ratios, temperature_ratio,
    StagnationRatios,
};
pub use normal_shock::{normal_shock, NormalShock};
pub use oblique_shock::{
    deflection_angle, max_deflection_angle, oblique_shock, shock_angle, ObliqueShock,
};
pub use prandtl_meyer::{
    mach_after_expansion, mach_from_prandtl_meyer, nu_max, prandtl_meyer_angle,
};
pub use rayleigh::{rayleigh_state, RayleighState};

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute tolerance for the cross-module integration checks.
    const EPS: f64 = 1e-9;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < EPS
    }

    /// Cross-module consistency: the static temperature ratio across a
    /// normal shock equals the ratio of the two isentropic temperature
    /// ratios evaluated at `M1` and `M2`, because the stagnation
    /// temperature is conserved across an adiabatic shock
    /// (`T0` is the same on both sides):
    ///
    /// ```text
    /// T2/T1 = (T0/T1) / (T0/T2) = (1 + k M1^2) / (1 + k M2^2),  k=(g-1)/2
    /// ```
    #[test]
    fn shock_temperature_ratio_consistent_with_isentropic() {
        let g = 1.4;
        let m1 = 2.5;
        let shock = normal_shock(m1, g).unwrap();
        let stag1 = temperature_ratio(m1, g).unwrap();
        let stag2 = temperature_ratio(shock.downstream_mach, g).unwrap();
        // T2/T1 = (T0/T at M1) / (T0/T at M2) since T0 is conserved.
        assert!(
            close(shock.temperature_ratio, stag1 / stag2),
            "shock T2/T1 = {}, isentropic = {}",
            shock.temperature_ratio,
            stag1 / stag2
        );
    }

    /// The bundled re-exports resolve and agree with the per-quantity
    /// functions for a representative supersonic state.
    #[test]
    fn reexports_resolve_and_agree() {
        // `density_ratio` re-exported at the crate root is the
        // *isentropic* stagnation ratio; it must agree with the
        // corresponding field of the bundled struct.
        let r: StagnationRatios = stagnation_ratios(1.8, 1.4).unwrap();
        assert!(close(r.t0_over_t, temperature_ratio(1.8, 1.4).unwrap()));
        assert!(close(r.p0_over_p, pressure_ratio(1.8, 1.4).unwrap()));
        assert!(close(r.rho0_over_rho, density_ratio(1.8, 1.4).unwrap()));
        assert!(area_mach_ratio(1.0, 1.4).unwrap() > 0.0);

        // The normal-shock bundle resolves through the crate root too.
        let s: NormalShock = normal_shock(1.8, 1.4).unwrap();
        assert!(s.downstream_mach < 1.0 && s.pressure_ratio > 1.0);

        // Error category surface is reachable through the crate root.
        assert_eq!(
            GasError::subsonic_shock(0.5).category(),
            ErrorCategory::Domain
        );
    }
}
