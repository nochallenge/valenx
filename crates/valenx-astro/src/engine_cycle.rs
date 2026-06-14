//! Full-flow staged-combustion (FFSC) engine **power-cycle balance** — the
//! defining capability of a Raptor-class methalox engine, and the answer to
//! "can this engine's turbopumps actually sustain the chamber pressure it is
//! designed for?"
//!
//! In a full-flow staged-combustion cycle (SpaceX's Raptor is the only flying
//! example) *all* the propellant passes through one of two preburners before
//! reaching the main chamber:
//!
//! - an **oxidizer-rich** preburner burns all the oxygen with a little fuel,
//!   and its hot gas drives the **oxidizer turbopump**;
//! - a **fuel-rich** preburner burns all the fuel with a little oxygen, and
//!   its gas drives the **fuel turbopump**.
//!
//! Each preburner runs at an extreme mixture ratio so its outlet (the
//! turbine-inlet) temperature stays low enough — roughly 700–900 K — to keep
//! the turbine blades intact. The turbine exhaust, still oxidizer-rich or
//! fuel-rich, is injected into the main chamber where final combustion at the
//! design mixture ratio happens. Because both turbines run on the full
//! propellant flow at modest temperature, the pumps can reach the very high
//! discharge pressures that make a ~300-bar chamber — and the high `Isp` that
//! follows — possible.
//!
//! ## The balance this models
//!
//! For each shaft (a pump and a turbine rigidly coupled), the pump's work per
//! kilogram of propellant is `Δp / (ρ · η_pump)`, and the turbine's work per
//! kilogram of gas is `c_p · T_in · (1 − Π^(−(γ−1)/γ)) · η_turbine` — an
//! isentropic expansion across the turbine pressure ratio `Π`. In a full-flow
//! cycle the gas mass through a turbine is essentially the propellant mass
//! through its own pump, so the shaft *closes* when the turbine's specific
//! work meets the pump's. That immediately yields the **maximum pump
//! discharge pressure** a given turbine-inlet temperature can sustain, and
//! hence the **maximum chamber pressure** the whole cycle can run — the real
//! limit a staged-combustion engine is designed against.
//!
//! ## Honest scope
//!
//! This is a **0-D steady power balance**, not a transient engine model. It
//! treats each shaft's turbine and pump mass flows as equal (the small
//! preburner cross-flows are a second-order correction), takes the
//! turbine-inlet temperature as the design parameter it really is (set by the
//! preburner mixture ratio), and uses representative pump/turbine efficiencies
//! plus a lumped injector-stiffness pressure budget. It predicts whether a
//! cycle closes and the chamber pressure it tops out at; it does **not** model
//! turbopump cavitation, bearing or seal losses, real turbomachinery maps,
//! start transients, or combustion stability. For Raptor-class methalox inputs
//! it lands on ~300-bar closure with ~30–40 MW per turbopump — the regime of
//! the real engine. The preburner gas properties can be informed by
//! [`crate::thermochem`] via [`ShaftInputs::gas_props_from`].

use serde::{Deserialize, Serialize};

use crate::thermochem::CombustionResult;

/// Universal gas constant (J/(mol·K)).
const R_UNIVERSAL: f64 = 8.314_462_618;

/// Liquid-oxygen density at storage conditions (kg/m³).
pub const RHO_LOX: f64 = 1_141.0;
/// Liquid-methane density at storage conditions (kg/m³).
pub const RHO_LCH4: f64 = 423.0;
/// RP-1 (kerosene) density at storage conditions (kg/m³).
pub const RHO_RP1: f64 = 810.0;
/// Liquid-hydrogen density at storage conditions (kg/m³).
pub const RHO_LH2: f64 = 71.0;

/// One turbopump shaft of the cycle — a pump and the turbine that drives it.
///
/// In a full-flow cycle there are two of these: the oxidizer shaft (pump
/// moving LOX, turbine fed by the oxidizer-rich preburner) and the fuel shaft
/// (pump moving fuel, turbine fed by the fuel-rich preburner).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ShaftInputs {
    /// Density of the fluid this pump moves (kg/m³).
    pub propellant_density: f64,
    /// Pump inlet (tank/feed) pressure (Pa).
    pub pump_inlet_pressure: f64,
    /// Pump isentropic efficiency (0–1).
    pub pump_efficiency: f64,
    /// Turbine-inlet (preburner outlet) temperature (K) — the design knob,
    /// held low by running the preburner very far from stoichiometric.
    pub turbine_inlet_temperature: f64,
    /// Preburner-gas specific heat at constant pressure (J/(kg·K)).
    pub turbine_gas_cp: f64,
    /// Preburner-gas ratio of specific heats (dimensionless).
    pub turbine_gas_gamma: f64,
    /// Turbine pressure ratio `Π = p_in / p_out` (> 1) — staged-combustion
    /// turbines run a modest expansion because the exhaust is reinjected at
    /// near-chamber pressure.
    pub turbine_pressure_ratio: f64,
    /// Turbine isentropic efficiency (0–1).
    pub turbine_efficiency: f64,
}

impl ShaftInputs {
    /// Turbine specific work per unit gas mass (J/kg) from an isentropic
    /// expansion across the turbine pressure ratio.
    pub fn turbine_specific_work(&self) -> f64 {
        let g = self.turbine_gas_gamma;
        let drop = 1.0 - self.turbine_pressure_ratio.powf(-(g - 1.0) / g);
        self.turbine_gas_cp * self.turbine_inlet_temperature * drop * self.turbine_efficiency
    }

    /// Pump specific work per unit propellant mass (J/kg) to reach a given
    /// discharge pressure.
    pub fn pump_specific_work(&self, discharge_pressure: f64) -> f64 {
        (discharge_pressure - self.pump_inlet_pressure)
            / (self.propellant_density * self.pump_efficiency)
    }

    /// Pump discharge pressure (Pa) required to land the propellant in the
    /// chamber at `chamber_pressure` after the turbine expansion and the main
    /// injector pressure drop: `p_c · stiffness · Π`.
    pub fn required_discharge(&self, chamber_pressure: f64, injector_stiffness: f64) -> f64 {
        chamber_pressure * injector_stiffness * self.turbine_pressure_ratio
    }

    /// Highest pump discharge pressure (Pa) this shaft can sustain — the point
    /// where the turbine's specific work exactly equals the pump's.
    pub fn max_discharge(&self) -> f64 {
        self.pump_inlet_pressure
            + self.propellant_density * self.pump_efficiency * self.turbine_specific_work()
    }

    /// Highest chamber pressure (Pa) this shaft alone can support, given the
    /// shared main-injector stiffness budget.
    pub fn max_chamber_pressure(&self, injector_stiffness: f64) -> f64 {
        self.max_discharge() / (injector_stiffness * self.turbine_pressure_ratio)
    }

    /// Derive `(c_p [J/(kg·K)], γ)` of a preburner gas from a
    /// [`CombustionResult`], so the turbine gas properties can come straight
    /// from the [`crate::thermochem`] equilibrium solver instead of being
    /// guessed. `c_p = γ/(γ−1) · R / M`.
    pub fn gas_props_from(result: &CombustionResult) -> (f64, f64) {
        let g = result.gamma;
        let molar_mass_kg = result.molar_mass / 1_000.0; // g/mol → kg/mol
        let cp = g / (g - 1.0) * R_UNIVERSAL / molar_mass_kg;
        (cp, g)
    }
}

/// Inputs to a full-flow staged-combustion power-balance solve.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CycleInputs {
    /// Main-chamber oxidizer/fuel mixture ratio (mass), used to split the
    /// total flow between the two shafts.
    pub mixture_ratio: f64,
    /// Target main-chamber stagnation pressure (Pa).
    pub chamber_pressure: f64,
    /// Total propellant mass flow through the engine (kg/s).
    pub total_mass_flow: f64,
    /// Main-injector stiffness — the ratio of the injector-inlet pressure to
    /// the chamber pressure (~1.2), the margin needed for combustion
    /// stability.
    pub main_injector_stiffness: f64,
    /// Oxidizer turbopump shaft.
    pub ox: ShaftInputs,
    /// Fuel turbopump shaft.
    pub fuel: ShaftInputs,
}

impl CycleInputs {
    /// A Raptor-class methalox full-flow staged-combustion design point: a
    /// ~300-bar chamber at MR 3.6 with ~700 kg/s of propellant, oxidizer- and
    /// fuel-rich preburners held to ~750 K, and representative pump/turbine
    /// efficiencies.
    pub fn raptor_methalox() -> Self {
        let ox = ShaftInputs {
            propellant_density: RHO_LOX,
            pump_inlet_pressure: 5.0e5,
            pump_efficiency: 0.75,
            turbine_inlet_temperature: 800.0,
            // Oxidizer-rich preburner gas: mostly O2 with some H2O/CO2 — a
            // comparatively low specific heat.
            turbine_gas_cp: 1_300.0,
            turbine_gas_gamma: 1.3,
            turbine_pressure_ratio: 1.5,
            turbine_efficiency: 0.78,
        };
        let fuel = ShaftInputs {
            propellant_density: RHO_LCH4,
            pump_inlet_pressure: 5.0e5,
            pump_efficiency: 0.75,
            turbine_inlet_temperature: 800.0,
            // Fuel-rich preburner gas: unburned CH4 plus H2/CO/H2O — a much
            // higher specific heat than the oxidizer-rich gas.
            turbine_gas_cp: 3_500.0,
            turbine_gas_gamma: 1.3,
            turbine_pressure_ratio: 1.5,
            turbine_efficiency: 0.78,
        };
        Self {
            mixture_ratio: 3.6,
            chamber_pressure: 300.0e5, // 300 bar
            total_mass_flow: 700.0,
            main_injector_stiffness: 1.2,
            ox,
            fuel,
        }
    }
}

/// Per-shaft outcome of the power balance.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ShaftResult {
    /// Propellant mass flow through this shaft's pump (kg/s).
    pub mass_flow: f64,
    /// Pump discharge pressure required for the target chamber pressure (Pa).
    pub required_discharge_pressure: f64,
    /// Highest discharge pressure the turbine can drive the pump to (Pa).
    pub max_discharge_pressure: f64,
    /// Shaft power demanded by the pump at the target chamber pressure (W).
    pub pump_power: f64,
    /// Shaft power the turbine produces (W).
    pub turbine_power: f64,
    /// Highest chamber pressure this shaft alone can support (Pa).
    pub max_chamber_pressure: f64,
    /// Whether this shaft closes — turbine power ≥ pump power.
    pub closes: bool,
}

/// Outcome of a full-flow staged-combustion power balance.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CycleResult {
    /// Oxidizer shaft result.
    pub ox: ShaftResult,
    /// Fuel shaft result.
    pub fuel: ShaftResult,
    /// Whether the *whole* cycle closes at the target chamber pressure — both
    /// shafts must close.
    pub closes: bool,
    /// Highest chamber pressure the cycle can sustain (Pa) — the lower of the
    /// two shaft limits.
    pub max_chamber_pressure: f64,
}

/// Solve one shaft's power balance for a given mass flow and target chamber
/// pressure.
fn solve_shaft(shaft: &ShaftInputs, mass_flow: f64, inputs: &CycleInputs) -> ShaftResult {
    let required =
        shaft.required_discharge(inputs.chamber_pressure, inputs.main_injector_stiffness);
    let max_discharge = shaft.max_discharge();
    let pump_power = mass_flow * shaft.pump_specific_work(required);
    let turbine_power = mass_flow * shaft.turbine_specific_work();
    let max_chamber_pressure = shaft.max_chamber_pressure(inputs.main_injector_stiffness);
    ShaftResult {
        mass_flow,
        required_discharge_pressure: required,
        max_discharge_pressure: max_discharge,
        pump_power,
        turbine_power,
        max_chamber_pressure,
        // Closes when the turbine can drive the pump to (at least) the
        // discharge the target chamber pressure needs.
        closes: max_discharge >= required,
    }
}

/// Run the full-flow staged-combustion power balance.
///
/// Splits the total flow between the two shafts by the mixture ratio, solves
/// each shaft, and reports whether the cycle closes at the target chamber
/// pressure together with the maximum chamber pressure it could run.
pub fn solve_cycle(inputs: &CycleInputs) -> CycleResult {
    let mr = inputs.mixture_ratio.max(1e-6);
    let m_ox = inputs.total_mass_flow * mr / (1.0 + mr);
    let m_fuel = inputs.total_mass_flow / (1.0 + mr);

    let ox = solve_shaft(&inputs.ox, m_ox, inputs);
    let fuel = solve_shaft(&inputs.fuel, m_fuel, inputs);

    CycleResult {
        closes: ox.closes && fuel.closes,
        max_chamber_pressure: ox.max_chamber_pressure.min(fuel.max_chamber_pressure),
        ox,
        fuel,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::thermochem::{combust, Propellant};

    #[test]
    fn raptor_methalox_cycle_closes_near_300_bar() {
        let r = solve_cycle(&CycleInputs::raptor_methalox());
        assert!(r.closes, "Raptor-class cycle should close at 300 bar");
        assert!(r.ox.closes && r.fuel.closes);
        // The cycle should top out close to the real Raptor regime
        // (~300–330 bar); a first-order balance lands a little under.
        let pc_bar = r.max_chamber_pressure / 1.0e5;
        assert!(
            (250.0..400.0).contains(&pc_bar),
            "max chamber pressure {pc_bar} bar"
        );
        // Turbopump powers land in the realistic tens-of-MW range.
        let ox_mw = r.ox.turbine_power / 1.0e6;
        let fuel_mw = r.fuel.turbine_power / 1.0e6;
        assert!((10.0..80.0).contains(&ox_mw), "ox turbine power {ox_mw} MW");
        assert!(
            (5.0..60.0).contains(&fuel_mw),
            "fuel turbine power {fuel_mw} MW"
        );
        // Each turbine must out-produce its pump for the shaft to close.
        assert!(r.ox.turbine_power >= r.ox.pump_power);
        assert!(r.fuel.turbine_power >= r.fuel.pump_power);
    }

    #[test]
    fn cycle_does_not_close_at_an_overambitious_chamber_pressure() {
        // The same hardware cannot reach 600 bar — the pumps would need a
        // discharge the turbines can't drive.
        let mut inputs = CycleInputs::raptor_methalox();
        inputs.chamber_pressure = 600.0e5;
        let r = solve_cycle(&inputs);
        assert!(
            !r.closes,
            "600 bar should not close on Raptor-class hardware"
        );
    }

    #[test]
    fn colder_turbines_support_a_lower_chamber_pressure() {
        // The turbine-inlet temperature is the cycle's master limit: cooler
        // preburners extract less work, so the sustainable chamber pressure
        // falls.
        let hot = solve_cycle(&CycleInputs::raptor_methalox());
        let mut cold_inputs = CycleInputs::raptor_methalox();
        cold_inputs.ox.turbine_inlet_temperature = 500.0;
        cold_inputs.fuel.turbine_inlet_temperature = 500.0;
        let cold = solve_cycle(&cold_inputs);
        assert!(
            cold.max_chamber_pressure < hot.max_chamber_pressure,
            "cold {} should be < hot {}",
            cold.max_chamber_pressure,
            hot.max_chamber_pressure
        );
    }

    #[test]
    fn lower_pump_efficiency_lowers_the_ceiling() {
        let base = solve_cycle(&CycleInputs::raptor_methalox());
        let mut worse = CycleInputs::raptor_methalox();
        worse.ox.pump_efficiency = 0.5;
        worse.fuel.pump_efficiency = 0.5;
        let r = solve_cycle(&worse);
        assert!(
            r.max_chamber_pressure < base.max_chamber_pressure,
            "less efficient pumps should lower the ceiling"
        );
    }

    #[test]
    fn thermochem_can_supply_the_preburner_gas_properties() {
        // The preburner gas c_p/γ can come straight from the equilibrium
        // combustion solver rather than being hand-set. Use a representative
        // methalox combustion to drive the conversion and confirm it yields a
        // physical specific heat and the same γ.
        let comb = combust(Propellant::Ch4Lox, 3.6, 300.0);
        let (cp, gamma) = ShaftInputs::gas_props_from(&comb);
        assert!(cp.is_finite() && cp > 0.0, "c_p = {cp}");
        assert_eq!(gamma, comb.gamma);
        // Feed those properties into a shaft and confirm it still produces a
        // finite, positive turbine work.
        let mut shaft = CycleInputs::raptor_methalox().fuel;
        shaft.turbine_gas_cp = cp;
        shaft.turbine_gas_gamma = gamma;
        assert!(shaft.turbine_specific_work() > 0.0);
    }

    #[test]
    fn is_deterministic_and_finite() {
        let a = solve_cycle(&CycleInputs::raptor_methalox());
        let b = solve_cycle(&CycleInputs::raptor_methalox());
        assert_eq!(a, b);
        assert!(a.max_chamber_pressure.is_finite());
        assert!(a.ox.turbine_power.is_finite() && a.fuel.pump_power.is_finite());
    }
}
