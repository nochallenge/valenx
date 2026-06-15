//! Single-stage vapor-compression cycle from four corner states.
//!
//! The standard textbook cycle is fixed by the specific enthalpy at four
//! corners, numbered in the direction of refrigerant flow:
//!
//! - **State 1** — evaporator outlet / compressor inlet (low-pressure vapor).
//! - **State 2** — compressor outlet / condenser inlet (high-pressure superheated vapor).
//! - **State 3** — condenser outlet (high-pressure liquid).
//! - **State 4** — expansion-valve (throttle) outlet / evaporator inlet.
//!
//! The expansion valve is modelled as an isenthalpic throttle, so the
//! state-4 enthalpy equals the state-3 enthalpy (`h4 = h3`); the cycle is
//! therefore fully specified by `h1`, `h2` and `h3`. From these three
//! the per-unit-mass effects follow directly:
//!
//! ```text
//! refrigerating effect   q_evap = h1 - h4 = h1 - h3
//! compressor work        w_comp = h2 - h1
//! condenser rejection    q_cond = h2 - h3
//! ```
//!
//! and the steady-flow energy balance closes exactly:
//! `q_cond = q_evap + w_comp`.
//!
//! [`Cycle`] holds the three independent enthalpies; [`CycleReport`]
//! holds the derived specific effects, both coefficients of performance,
//! and — when a refrigeration duty is supplied — the required refrigerant
//! mass-flow rate and the absolute heat / power flows.

use crate::cop::{cop_cool_from_enthalpies, cop_heat_from_enthalpies};
use crate::error::{RefrigError, Result};
use serde::{Deserialize, Serialize};

/// The three independent specific enthalpies (kJ/kg) that fix a
/// single-stage vapor-compression cycle.
///
/// State 4 is omitted because the isenthalpic throttle pins it to state
/// 3 (`h4 = h3`).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Cycle {
    /// Specific enthalpy at the evaporator outlet / compressor inlet (state 1), kJ/kg.
    pub h1: f64,
    /// Specific enthalpy at the compressor outlet / condenser inlet (state 2), kJ/kg.
    pub h2: f64,
    /// Specific enthalpy at the condenser outlet (state 3), kJ/kg.
    /// Also the throttle-outlet enthalpy (state 4), since `h4 = h3`.
    pub h3: f64,
}

impl Cycle {
    /// Construct a cycle from its three corner enthalpies, validating
    /// that they are finite and thermodynamically ordered for a
    /// vapor-compression refrigeration loop.
    ///
    /// The ordering required is `h2 > h1 > h3`: the compressor must add
    /// enthalpy (`h2 > h1`) and the refrigerant must be cooler/heavier
    /// after condensation than the vapor leaving the evaporator
    /// (`h1 > h3`, equivalently `h1 > h4`).
    ///
    /// # Errors
    ///
    /// Returns [`RefrigError::BadParameter`] if any enthalpy is
    /// non-finite, and [`RefrigError::InconsistentCycle`] if the
    /// ordering `h2 > h1 > h3` is violated.
    pub fn new(h1: f64, h2: f64, h3: f64) -> Result<Self> {
        let h1 = RefrigError::require_finite("h1", h1)?;
        let h2 = RefrigError::require_finite("h2", h2)?;
        let h3 = RefrigError::require_finite("h3", h3)?;

        if h2 <= h1 {
            return Err(RefrigError::InconsistentCycle(
                "compressor must raise enthalpy: require h2 > h1",
            ));
        }
        if h1 <= h3 {
            return Err(RefrigError::InconsistentCycle(
                "refrigerating effect must be positive: require h1 > h3 (= h4)",
            ));
        }
        Ok(Cycle { h1, h2, h3 })
    }

    /// The throttle-outlet enthalpy (state 4), equal to `h3` for the
    /// isenthalpic expansion modelled here.
    pub fn h4(&self) -> f64 {
        self.h3
    }

    /// Specific refrigerating effect `q_evap = h1 - h4 = h1 - h3` (kJ/kg).
    pub fn refrigerating_effect(&self) -> f64 {
        self.h1 - self.h3
    }

    /// Specific compressor work `w_comp = h2 - h1` (kJ/kg).
    pub fn compressor_work(&self) -> f64 {
        self.h2 - self.h1
    }

    /// Specific condenser heat rejection `q_cond = h2 - h3` (kJ/kg).
    pub fn heat_rejected(&self) -> f64 {
        self.h2 - self.h3
    }

    /// Cooling coefficient of performance `(h1 - h4) / (h2 - h1)`.
    ///
    /// # Errors
    ///
    /// Propagates [`cop_cool_from_enthalpies`]; for a [`Cycle`] built via
    /// [`Cycle::new`] the inputs are already valid so this never fails,
    /// but the `Result` is preserved for uniformity.
    pub fn cop_cool(&self) -> Result<f64> {
        cop_cool_from_enthalpies(self.h1, self.h2, self.h4())
    }

    /// Heating coefficient of performance `(h2 - h3) / (h2 - h1)`.
    ///
    /// # Errors
    ///
    /// Propagates [`cop_heat_from_enthalpies`].
    pub fn cop_heat(&self) -> Result<f64> {
        cop_heat_from_enthalpies(self.h1, self.h2, self.h3)
    }

    /// Build a [`CycleReport`] of the specific effects and both COPs,
    /// without any absolute-duty scaling (the flow fields are `None`).
    ///
    /// # Errors
    ///
    /// Propagates the COP computations.
    pub fn report(&self) -> Result<CycleReport> {
        Ok(CycleReport {
            refrigerating_effect: self.refrigerating_effect(),
            compressor_work: self.compressor_work(),
            heat_rejected: self.heat_rejected(),
            cop_cool: self.cop_cool()?,
            cop_heat: self.cop_heat()?,
            mass_flow: None,
            cooling_capacity: None,
            compressor_power: None,
            heat_rejection_rate: None,
        })
    }

    /// Build a [`CycleReport`] scaled to a target cooling capacity
    /// (refrigeration duty) `q_dot_evap`, in kilowatts.
    ///
    /// The required refrigerant mass-flow rate is
    /// `m_dot = q_dot_evap / q_evap` (kg/s, with `q_evap` in kJ/kg), and
    /// the absolute compressor power and condenser heat-rejection rate
    /// follow by multiplying the corresponding specific quantities by
    /// `m_dot`.
    ///
    /// # Errors
    ///
    /// Returns [`RefrigError::BadParameter`] if `q_dot_evap` is negative
    /// or non-finite, and propagates the COP computations.
    pub fn report_for_duty(&self, q_dot_evap: f64) -> Result<CycleReport> {
        let duty = RefrigError::require_non_negative("q_dot_evap", q_dot_evap)?;
        let q_evap = self.refrigerating_effect();
        // `q_evap > 0` is guaranteed by `Cycle::new`.
        let m_dot = duty / q_evap;

        let mut report = self.report()?;
        report.mass_flow = Some(m_dot);
        report.cooling_capacity = Some(duty);
        report.compressor_power = Some(m_dot * self.compressor_work());
        report.heat_rejection_rate = Some(m_dot * self.heat_rejected());
        Ok(report)
    }
}

/// Derived results for a single-stage vapor-compression cycle.
///
/// The first three fields are specific (per unit mass, kJ/kg); the COPs
/// are dimensionless. The four optional flow fields are populated only by
/// [`Cycle::report_for_duty`] and carry absolute rates: mass flow in
/// kg/s and the three energy rates in kilowatts.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CycleReport {
    /// Specific refrigerating effect `h1 - h4` (kJ/kg).
    pub refrigerating_effect: f64,
    /// Specific compressor work `h2 - h1` (kJ/kg).
    pub compressor_work: f64,
    /// Specific condenser heat rejection `h2 - h3` (kJ/kg).
    pub heat_rejected: f64,
    /// Cooling coefficient of performance (dimensionless).
    pub cop_cool: f64,
    /// Heating coefficient of performance (dimensionless).
    pub cop_heat: f64,
    /// Refrigerant mass-flow rate (kg/s), present only for a duty-scaled report.
    pub mass_flow: Option<f64>,
    /// Cooling capacity / refrigeration duty (kW), present only for a duty-scaled report.
    pub cooling_capacity: Option<f64>,
    /// Compressor shaft power (kW), present only for a duty-scaled report.
    pub compressor_power: Option<f64>,
    /// Condenser heat-rejection rate (kW), present only for a duty-scaled report.
    pub heat_rejection_rate: Option<f64>,
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute tolerance for floating-point comparisons in this module.
    const EPS: f64 = 1e-9;

    /// The R-134a worked example from `cop.rs`, reused here as the
    /// reference cycle: h1 = 239.16, h2 = 275.39, h3 = 95.47 kJ/kg.
    fn reference_cycle() -> Cycle {
        Cycle::new(239.16, 275.39, 95.47).unwrap()
    }

    #[test]
    fn specific_effects_match_hand_calc() {
        let c = reference_cycle();
        assert!((c.refrigerating_effect() - 143.69).abs() < 1e-9, "{c:?}");
        assert!((c.compressor_work() - 36.23).abs() < 1e-9, "{c:?}");
        assert!((c.heat_rejected() - 179.92).abs() < 1e-9, "{c:?}");
        // h4 is pinned to h3.
        assert!((c.h4() - c.h3).abs() < EPS);
    }

    #[test]
    fn energy_balance_closes_exactly() {
        // The defining check: q_cond = q_evap + w_comp.
        let c = reference_cycle();
        let lhs = c.heat_rejected();
        let rhs = c.refrigerating_effect() + c.compressor_work();
        assert!((lhs - rhs).abs() < EPS, "lhs={lhs} rhs={rhs}");
    }

    #[test]
    fn report_cops_obey_plus_one_identity() {
        let r = reference_cycle().report().unwrap();
        assert!((r.cop_heat - (r.cop_cool + 1.0)).abs() < EPS, "{r:?}");
        // Cooling COP equals the textbook ~3.966.
        assert!((r.cop_cool - 3.9660_f64).abs() < 1e-3, "{r:?}");
        // Un-scaled report has no flow data.
        assert!(r.mass_flow.is_none());
        assert!(r.compressor_power.is_none());
    }

    #[test]
    fn duty_scaling_gives_consistent_rates() {
        // Scale to a 10 kW cooling load.
        let c = reference_cycle();
        let r = c.report_for_duty(10.0).unwrap();

        let m_dot = r.mass_flow.unwrap();
        // m_dot = Q_evap / q_evap.
        assert!(
            (m_dot - 10.0 / c.refrigerating_effect()).abs() < EPS,
            "{r:?}"
        );

        // The cooling capacity field echoes the requested duty.
        assert!((r.cooling_capacity.unwrap() - 10.0).abs() < EPS);

        // Absolute power = m_dot * specific work, and the COP recovered
        // from the absolute rates must match the specific-basis COP.
        let power = r.compressor_power.unwrap();
        assert!((power - m_dot * c.compressor_work()).abs() < EPS, "{r:?}");
        let cop_from_rates = r.cooling_capacity.unwrap() / power;
        assert!((cop_from_rates - r.cop_cool).abs() < 1e-9, "{r:?}");

        // Rate-level energy balance: Q_cond_rate = Q_evap_rate + W_dot.
        let q_cond_rate = r.heat_rejection_rate.unwrap();
        assert!(
            (q_cond_rate - (r.cooling_capacity.unwrap() + power)).abs() < 1e-9,
            "{r:?}"
        );
    }

    #[test]
    fn zero_duty_is_allowed_and_gives_zero_flow() {
        let r = reference_cycle().report_for_duty(0.0).unwrap();
        assert!((r.mass_flow.unwrap() - 0.0).abs() < EPS);
        assert!((r.compressor_power.unwrap() - 0.0).abs() < EPS);
    }

    #[test]
    fn constructor_rejects_unordered_enthalpies() {
        // h2 <= h1.
        assert!(matches!(
            Cycle::new(240.0, 240.0, 95.0),
            Err(RefrigError::InconsistentCycle(_))
        ));
        // h1 <= h3.
        assert!(matches!(
            Cycle::new(95.0, 275.0, 95.0),
            Err(RefrigError::InconsistentCycle(_))
        ));
        // Non-finite.
        assert!(matches!(
            Cycle::new(f64::INFINITY, 275.0, 95.0),
            Err(RefrigError::BadParameter { .. })
        ));
    }

    #[test]
    fn report_for_duty_rejects_negative_duty() {
        assert!(matches!(
            reference_cycle().report_for_duty(-1.0),
            Err(RefrigError::BadParameter { .. })
        ));
    }
}
