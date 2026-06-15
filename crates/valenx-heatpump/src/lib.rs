//! # valenx-heatpump — Carnot heat-pump thermodynamics
//!
//! Closed-form reversible-cycle performance of a heat pump (or, run in
//! reverse, a chiller) plus a load-versus-capacity *balance-point*
//! solver for an air-source unit in winter.
//!
//! ## What
//!
//! The [`cop`] module — **coefficient of performance** — provides the
//! reversible (Carnot) heating and cooling COPs from two absolute
//! reservoir temperatures, the [`CarnotCop`] value type that bundles
//! them, the `COP_heat = COP_cool + 1` identity, and a second-law /
//! Carnot-fraction derating that scales the ideal limits down to
//! real-machine numbers.
//!
//! The [`balance`] module — **balance point** — crosses a building
//! [`LoadLine`] (heating load rising as it gets colder) against an
//! air-source [`CapacityLine`] (heating capacity falling as it gets
//! colder), with both a bracketed-bisection root-find
//! ([`solve_balance_point`]) and the affine closed-form cross-check
//! ([`balance_point_linear`]) for the outdoor temperature at which load
//! equals capacity.
//!
//! The [`error`] module supplies a [`HeatPumpError`] enum with validated
//! constructors and stable [`code`](HeatPumpError::code) /
//! [`category`](HeatPumpError::category) accessors for telemetry.
//!
//! Every fallible function returns [`Result<_, HeatPumpError>`].
//!
//! ## Model
//!
//! A heat pump moves heat `Q_c` out of a cold reservoir at absolute
//! temperature `T_c` and delivers `Q_h` into a hot reservoir at `T_h`,
//! consuming compressor work `W = Q_h - Q_c`. The reversible Carnot
//! cycle sets the efficiency ceiling. Writing the lift as
//! `dT = T_h - T_c > 0`:
//!
//! ```text
//! COP_heat = Q_h / W = T_h / (T_h - T_c)        (useful output = heating)
//! COP_cool = Q_c / W = T_c / (T_h - T_c)        (useful output = cooling)
//! ```
//!
//! Because `T_h = T_c + dT`, the two differ by exactly one:
//! `COP_heat - COP_cool = (T_h - T_c) / dT = 1`, i.e.
//! `COP_heat = COP_cool + 1` — the compressor work itself shows up as
//! heat in the hot reservoir. Both COPs grow without bound as the lift
//! shrinks and fall toward `1` (heating) / `0` (cooling) as the lift
//! grows, so a heat pump always beats resistive heating
//! (`COP_heat > 1`).
//!
//! The balance-point model is two straight lines on a (outdoor
//! temperature, kilowatts) plot. The building load is
//! `load(T) = UA * (T_balance - T)` for `T < T_balance`, rising as the
//! outdoor temperature `T` drops. The air-source capacity is
//! `capacity(T) = cap_ref + slope * (T - T_ref)` with a positive
//! `slope`, so capacity *falls* as `T` drops. The balance point is where
//! the lines cross; below it a backup heat source must cover the
//! deficit. The crossing is solved by bisection on the residual
//! `capacity - load` (a method that survives a swap to the non-linear
//! capacity tables real tools use) and cross-checked against the affine
//! closed form.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are textbook closed-form and
//! simple numerical models — reversible-cycle COP limits and linear
//! load / capacity lines — **NOT a clinical/medical or production
//! engineering tool**. The Carnot COP is a thermodynamic *upper bound*
//! that no real vapour-compression machine reaches; this crate does not
//! model refrigerant properties, compressor maps, superheat / subcool,
//! defrost cycles, part-load degradation, or the temperature glide of a
//! finite-capacity heat exchanger. The balance-point lines are a
//! deliberately simple affine idealisation; real load and capacity
//! curves bend. Nothing here sizes, rates, certifies, or commissions
//! actual HVAC equipment — do not use it for that. Treat the outputs as
//! illustrative thermodynamic limits, not design figures.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod balance;
pub mod cop;
pub mod error;

// --- Convenience re-exports of the most-used types --------------------

pub use balance::{
    balance_point_linear, solve_balance_point, BalancePoint, CapacityLine, LoadLine,
};
pub use cop::{carnot_cop_cool, carnot_cop_heat, CarnotCop, Derated};
pub use error::{ErrorCategory, HeatPumpError, Result};

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute tolerance for floating-point comparisons.
    const EPS: f64 = 1e-9;

    /// End-to-end: pick a winter design point, compute the Carnot
    /// heating COP, and locate the balance point for a matched
    /// building / heat-pump pair — confirming the headline relations all
    /// hold together.
    #[test]
    fn cop_and_balance_point_end_to_end() {
        // Heating from -5 °C outdoor (268.15 K) to a 45 °C supply
        // (318.15 K): lift = 50 K.
        let c = CarnotCop::new(268.15, 318.15).unwrap();
        assert!(c.cop_heat > 1.0);
        assert!((c.cop_heat - c.cop_cool - 1.0).abs() < EPS);
        assert!((c.cop_heat - 318.15 / 50.0).abs() < EPS);

        // A real unit at ~50% of Carnot.
        let real = c.derated_heat(0.5).unwrap();
        assert!(real > 1.0 && real < c.cop_heat);

        // Balance point of an 18 °C / 0.5 kW-per-K building against a
        // 10 kW (rated at 8.3 °C, -0.25 kW/K) air-source unit.
        let load = LoadLine::new(18.0, 0.5).unwrap();
        let cap = CapacityLine::new(8.3, 10.0, 0.25).unwrap();
        let bp = solve_balance_point(&load, &cap, -25.0, 18.0, 1e-12).unwrap();
        assert!((bp.load_kw - bp.capacity_kw).abs() < 1e-6);

        // The bisection root agrees with the analytic closed form.
        let t_closed = balance_point_linear(&load, &cap).unwrap();
        assert!((bp.t_balance_c - t_closed).abs() < 1e-6);
    }
}
