//! # valenx-refrigeration — vapor-compression refrigeration thermodynamics
//!
//! Closed-form coefficient-of-performance and single-stage cycle models
//! for vapor-compression refrigerators and heat pumps. Pure scalar
//! `f64` algorithms with no external processes, property libraries or
//! platform dependencies.
//!
//! ## What
//!
//! Three small modules cover the textbook treatment of the
//! vapor-compression cycle:
//!
//! - [`cop`] — the coefficient-of-performance definitions: the cooling
//!   COP [`cop::cop_cool`] as the refrigerating effect over the
//!   compressor work, the heating COP [`cop::cop_heat`] as the condenser
//!   heat rejection over the compressor work, the energy-balance
//!   identity [`cop::cop_heat_from_cool`]
//!   (`COP_heat = COP_cool + 1`), and the same numbers formed directly
//!   from the cycle-corner enthalpies via
//!   [`cop::cop_cool_from_enthalpies`] and
//!   [`cop::cop_heat_from_enthalpies`].
//! - [`carnot`] — the reversible upper limits
//!   [`carnot::carnot_cop_cool`] (`Tc / (Th - Tc)`) and
//!   [`carnot::carnot_cop_heat`] (`Th / (Th - Tc)`) between two
//!   reservoirs, plus the second-law (exergetic) efficiencies
//!   [`carnot::second_law_efficiency_cool`] and
//!   [`carnot::second_law_efficiency_heat`] of a real cycle against its
//!   Carnot bound.
//! - [`cycle`] — a four-corner [`cycle::Cycle`] fixed by three specific
//!   enthalpies (with the throttle pinned isenthalpically, `h4 = h3`)
//!   that derives the specific refrigerating effect, compressor work and
//!   condenser rejection, both COPs, and — scaled to a target
//!   refrigeration duty — the refrigerant mass-flow rate and the
//!   absolute power and heat-flow rates, all collected in a
//!   [`cycle::CycleReport`].
//!
//! Every fallible entry point returns [`error::Result`], and
//! [`error::RefrigError`] exposes stable [`code`](error::RefrigError::code)
//! and [`category`](error::RefrigError::category) accessors for telemetry.
//!
//! ## Model
//!
//! Numbering the cycle corners in the direction of refrigerant flow —
//! state 1 the evaporator outlet / compressor inlet, state 2 the
//! compressor outlet, state 3 the condenser outlet and state 4 the
//! throttle outlet — the governing relations are:
//!
//! ```text
//! COP_cool = Q_evap / W_comp = (h1 - h4) / (h2 - h1)
//! COP_heat = Q_cond / W_comp = (h2 - h3) / (h2 - h1)
//! Q_cond   = Q_evap + W_comp                  (steady-flow energy balance)
//! COP_heat = COP_cool + 1                      (follows from the balance)
//! h4       = h3                                (isenthalpic expansion valve)
//! ```
//!
//! and, for the reversible device between a cold reservoir at absolute
//! temperature `Tc` and a hot reservoir at `Th > Tc`:
//!
//! ```text
//! COP_cool,Carnot = Tc / (Th - Tc)
//! COP_heat,Carnot = Th / (Th - Tc)
//! ```
//!
//! Both Carnot limits fall as the temperature lift `Th - Tc` grows, so a
//! larger lift always lowers the achievable COP, and no real cycle can
//! exceed its Carnot bound (the second-law efficiency lies in `(0, 1]`).
//!
//! ## Honest scope
//!
//! Research/educational grade: these are textbook closed-form and
//! idealised numerical models. The compression is treated as ideal
//! isentropic or fixed-efficiency, the expansion as a perfectly
//! isenthalpic throttle, and the reservoirs as lumped at single
//! temperatures; pressure drops, heat leaks, transient behaviour and
//! refrigerant-specific property surfaces are out of scope. This crate
//! is **not** a clinical/medical or production engineering tool, and is
//! not a substitute for a validated refrigerant-property library such as
//! CoolProp or REFPROP, nor for HVAC equipment-selection software. The
//! enthalpies it consumes must come from such a source; this crate only
//! combines them into the standard performance figures.

#![forbid(unsafe_code)]
#![warn(missing_docs)]
#![warn(unused_imports)]

pub mod carnot;
pub mod cop;
pub mod cycle;
pub mod error;

pub use carnot::{
    carnot_cop_cool, carnot_cop_heat, second_law_efficiency_cool, second_law_efficiency_heat,
};
pub use cop::{
    cop_cool, cop_cool_from_enthalpies, cop_heat, cop_heat_from_cool, cop_heat_from_enthalpies,
};
pub use cycle::{Cycle, CycleReport};
pub use error::{ErrorCategory, RefrigError, Result};

#[cfg(test)]
mod tests {
    //! Cross-module integration checks tying the COP, Carnot and cycle
    //! layers together on a single consistent scenario.

    use super::*;

    /// Absolute tolerance for floating-point comparisons.
    const EPS: f64 = 1e-9;

    #[test]
    fn cycle_cop_is_bounded_by_carnot() {
        // Take the reference R-134a cycle (evaporating around -20 C,
        // condensing around 31 C for 0.8 MPa R-134a) and confirm its
        // cooling COP sits below the Carnot limit for those reservoir
        // temperatures, giving a physical second-law efficiency.
        let cycle = Cycle::new(239.16, 275.39, 95.47).unwrap();
        let cop = cycle.cop_cool().unwrap();

        let t_cold = 253.15; // -20 C
        let t_hot = 304.15; //  31 C
        let carnot = carnot_cop_cool(t_cold, t_hot).unwrap();

        assert!(
            cop < carnot,
            "cycle COP {cop} should be below Carnot {carnot}"
        );

        let eta = second_law_efficiency_cool(cop, t_cold, t_hot).unwrap();
        assert!(eta > 0.0 && eta < 1.0, "eta={eta}");
    }

    #[test]
    fn all_three_paths_to_heating_cop_agree() {
        // (1) from heats, (2) from cooling-COP + 1, (3) from enthalpies —
        // all three must coincide for one consistent cycle.
        let cycle = Cycle::new(239.16, 275.39, 95.47).unwrap();

        let q_evap = cycle.refrigerating_effect();
        let w = cycle.compressor_work();
        let q_cond = cycle.heat_rejected();

        let from_heats = cop_heat(q_cond, w).unwrap();
        let from_plus_one = cop_heat_from_cool(cop_cool(q_evap, w).unwrap()).unwrap();
        let from_enthalpies = cop_heat_from_enthalpies(cycle.h1, cycle.h2, cycle.h3).unwrap();

        assert!((from_heats - from_plus_one).abs() < EPS);
        assert!((from_heats - from_enthalpies).abs() < EPS);
    }
}
