//! First-order Thevenin equivalent-circuit cell model.
//!
//! ## Model
//!
//! A lithium-ion (or comparable) cell is approximated by the standard
//! **first-order (one-RC) Thevenin equivalent circuit**:
//!
//! ```text
//!            R0           R1
//!     o----[ R0 ]----+--[ R1 ]--+----o  (+) terminal
//!                    |          |
//!   OCV(SoC)        ===C1       |
//!     (source)       |          |
//!     o--------------+----------+----o  (-) terminal
//! ```
//!
//! - **`OCV(SoC)`** — the open-circuit (rest) voltage source, a monotone
//!   function of state of charge supplied by an [`OcvSocTable`].
//! - **`R0`** — the series *ohmic* resistance. It produces an
//!   *instantaneous* `I·R0` voltage drop the moment a load is applied.
//! - **`R1 ∥ C1`** — one resistor-capacitor pair modelling the cell's
//!   *diffusion / polarisation* over-potential. Under a constant current
//!   its voltage relaxes towards `I·R1` with time constant `τ = R1·C1`.
//!
//! ### Sign convention
//!
//! Current `I` is **positive on discharge** (conventional current
//! leaving the positive terminal). Consequently:
//!
//! - the terminal voltage **drops** below OCV under discharge, and
//! - the state of charge **falls** while discharging.
//!
//! A negative `I` is a charge current: the terminal voltage rises above
//! OCV and the SoC climbs.
//!
//! ### Governing equations
//!
//! Terminal voltage (algebraic, holds at every instant):
//!
//! ```text
//! V_t = OCV(SoC) - I·R0 - V_rc
//! ```
//!
//! State of charge by **coulomb counting** (`Q` is the usable capacity
//! in ampere-seconds, i.e. `capacity_Ah · 3600`):
//!
//! ```text
//! SoC(t) = SoC0 - (1/Q) ∫ I dt
//! ```
//!
//! RC-pair over-potential, integrated **exactly** across a step of
//! constant current `I` (an exponential map, unconditionally stable):
//!
//! ```text
//! V_rc(t+Δt) = V_rc(t)·e^(-Δt/τ) + I·R1·(1 - e^(-Δt/τ)),   τ = R1·C1
//! ```
//!
//! At `I = 0` this reduces to pure relaxation `V_rc·e^(-Δt/τ)`, so the
//! over-potential decays to `1/e` of its value after exactly one `τ`.

use crate::error::{EcmError, Result};
use crate::ocv::OcvSocTable;
use serde::{Deserialize, Serialize};

/// Seconds per hour — converts a capacity in ampere-hours to the
/// ampere-seconds (coulombs) used by the coulomb-counter.
pub const SECONDS_PER_HOUR: f64 = 3600.0;

/// Static parameters of a first-order Thevenin cell.
///
/// All resistances and the capacitance are strictly positive; the
/// capacity is in ampere-hours. Build with [`CellParams::new`], which
/// validates every field, then attach an [`OcvSocTable`] when
/// constructing a [`CellState`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CellParams {
    /// Series ohmic resistance `R0` (ohms). Drives the instantaneous
    /// `I·R0` step.
    pub r0: f64,
    /// Polarisation resistance `R1` (ohms) of the RC pair.
    pub r1: f64,
    /// Polarisation capacitance `C1` (farads) of the RC pair.
    pub c1: f64,
    /// Usable cell capacity (ampere-hours) for the coulomb counter.
    pub capacity_ah: f64,
}

impl CellParams {
    /// Build a validated parameter set.
    ///
    /// # Errors
    ///
    /// Returns [`EcmError::Invalid`] if any of `r0`, `r1`, `c1` or
    /// `capacity_ah` is non-finite or not strictly positive.
    pub fn new(r0: f64, r1: f64, c1: f64, capacity_ah: f64) -> Result<Self> {
        check_positive("r0", r0)?;
        check_positive("r1", r1)?;
        check_positive("c1", c1)?;
        check_positive("capacity_ah", capacity_ah)?;
        Ok(Self {
            r0,
            r1,
            c1,
            capacity_ah,
        })
    }

    /// RC time constant `τ = R1·C1` (seconds).
    pub fn tau(&self) -> f64 {
        self.r1 * self.c1
    }

    /// Usable capacity expressed in coulombs (ampere-seconds).
    pub fn capacity_coulombs(&self) -> f64 {
        self.capacity_ah * SECONDS_PER_HOUR
    }
}

/// A Thevenin cell together with its mutable state (SoC and RC
/// over-potential) and its OCV-SoC curve.
///
/// Construct with [`CellState::new`]; advance it with
/// [`CellState::step`]; read the loaded terminal voltage with
/// [`CellState::terminal_voltage`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CellState {
    params: CellParams,
    ocv: OcvSocTable,
    /// State of charge, a fraction in `[0, 1]`.
    soc: f64,
    /// Voltage across the RC pair (volts).
    v_rc: f64,
}

impl CellState {
    /// Build a cell at a given initial state of charge, fully relaxed
    /// (`V_rc = 0`).
    ///
    /// # Errors
    ///
    /// Returns [`EcmError::SocOutOfRange`] if `soc0` is outside
    /// `[0, 1]`, or [`EcmError::Invalid`] if it is non-finite.
    pub fn new(params: CellParams, ocv: OcvSocTable, soc0: f64) -> Result<Self> {
        Self::with_overpotential(params, ocv, soc0, 0.0)
    }

    /// Build a cell at a given SoC *and* a given initial RC
    /// over-potential (volts). Useful for resuming a paused simulation.
    ///
    /// # Errors
    ///
    /// Returns [`EcmError::SocOutOfRange`] if `soc0 ∉ [0, 1]`, or
    /// [`EcmError::Invalid`] if `soc0` or `v_rc0` is non-finite.
    pub fn with_overpotential(
        params: CellParams,
        ocv: OcvSocTable,
        soc0: f64,
        v_rc0: f64,
    ) -> Result<Self> {
        check_soc(soc0, "initial SoC")?;
        if !v_rc0.is_finite() {
            return Err(EcmError::invalid("v_rc0", "must be finite"));
        }
        Ok(Self {
            params,
            ocv,
            soc: soc0,
            v_rc: v_rc0,
        })
    }

    /// Borrow the static parameters.
    pub fn params(&self) -> &CellParams {
        &self.params
    }

    /// Borrow the OCV-SoC table.
    pub fn ocv_table(&self) -> &OcvSocTable {
        &self.ocv
    }

    /// Current state of charge (fraction in `[0, 1]`).
    pub fn soc(&self) -> f64 {
        self.soc
    }

    /// Current RC-pair over-potential (volts).
    pub fn v_rc(&self) -> f64 {
        self.v_rc
    }

    /// Open-circuit voltage at the present state of charge.
    pub fn ocv(&self) -> f64 {
        self.ocv.ocv_at(self.soc)
    }

    /// Terminal voltage **under a given load current** at the present
    /// state, without advancing time.
    ///
    /// `V_t = OCV(SoC) − I·R0 − V_rc`, with `I > 0` on discharge. This
    /// captures the *instantaneous* `I·R0` drop the moment a current is
    /// applied (the RC term `V_rc` is whatever the pair has built up so
    /// far — zero immediately after a step change from rest).
    pub fn terminal_voltage(&self, current: f64) -> f64 {
        self.ocv() - current * self.params.r0 - self.v_rc
    }

    /// Terminal voltage at the present state with **no external load**
    /// (`I = 0`): `OCV(SoC) − V_rc`. Equal to [`ocv`](Self::ocv) only
    /// when the cell is fully relaxed.
    pub fn open_terminal_voltage(&self) -> f64 {
        self.terminal_voltage(0.0)
    }

    /// Advance the cell by `dt` seconds under a **constant** current
    /// `current` (amperes, `> 0` on discharge), and return the terminal
    /// voltage *at the start of the step* (i.e. with the over-potential
    /// the cell carried *into* the step, reflecting the instantaneous
    /// `I·R0` drop).
    ///
    /// The update is:
    ///
    /// - **SoC**: `SoC -= I·Δt / Q` (coulomb counting), then clamped to
    ///   `[0, 1]`.
    /// - **`V_rc`**: the exact exponential map
    ///   `V_rc·e^(-Δt/τ) + I·R1·(1 − e^(-Δt/τ))`.
    ///
    /// # Errors
    ///
    /// Returns [`EcmError::Invalid`] if `dt` is negative or non-finite,
    /// or if `current` is non-finite.
    pub fn step(&mut self, current: f64, dt: f64) -> Result<f64> {
        if !dt.is_finite() || dt < 0.0 {
            return Err(EcmError::invalid("dt", "must be finite and non-negative"));
        }
        if !current.is_finite() {
            return Err(EcmError::invalid("current", "must be finite"));
        }

        // Terminal voltage seen at the instant the step begins.
        let v_terminal = self.terminal_voltage(current);

        // Coulomb-counted SoC update. I>0 (discharge) lowers SoC.
        let dq = current * dt; // coulombs drawn
        self.soc -= dq / self.params.capacity_coulombs();
        self.soc = self.soc.clamp(0.0, 1.0);

        // Exact exponential update of the RC over-potential.
        let tau = self.params.tau();
        let decay = (-dt / tau).exp();
        self.v_rc = self.v_rc * decay + current * self.params.r1 * (1.0 - decay);

        Ok(v_terminal)
    }
}

/// Validate that a named scalar is finite and strictly positive.
fn check_positive(what: &'static str, x: f64) -> Result<()> {
    if !x.is_finite() {
        return Err(EcmError::invalid(what, "must be finite"));
    }
    if x <= 0.0 {
        return Err(EcmError::invalid(what, "must be strictly positive"));
    }
    Ok(())
}

/// Validate that a state-of-charge value is finite and within `[0, 1]`.
fn check_soc(soc: f64, context: &'static str) -> Result<()> {
    if !soc.is_finite() {
        return Err(EcmError::invalid("soc", "must be finite"));
    }
    if !(0.0..=1.0).contains(&soc) {
        return Err(EcmError::SocOutOfRange { soc, context });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A representative small cell: 10 mΩ ohmic, a 20 mΩ / 1000 F RC pair
    /// (τ = 20 s), 2 Ah capacity, on a 3.0–4.2 V curve.
    fn cell(soc0: f64) -> CellState {
        let params = CellParams::new(0.01, 0.02, 1000.0, 2.0).unwrap();
        let ocv = OcvSocTable::new(
            vec![0.0, 0.25, 0.5, 0.75, 1.0],
            vec![3.0, 3.5, 3.7, 3.9, 4.2],
        )
        .unwrap();
        CellState::new(params, ocv, soc0).unwrap()
    }

    #[test]
    fn tau_and_capacity_helpers() {
        let p = CellParams::new(0.01, 0.02, 1000.0, 2.0).unwrap();
        assert!((p.tau() - 20.0).abs() < 1e-12);
        assert!((p.capacity_coulombs() - 7200.0).abs() < 1e-9);
    }

    #[test]
    fn terminal_equals_ocv_at_zero_current_when_relaxed() {
        // Validation: at I = 0 (and fully relaxed) terminal V == OCV.
        let c = cell(0.5);
        let ocv = c.ocv();
        assert!((c.terminal_voltage(0.0) - ocv).abs() < 1e-12);
        assert!((c.open_terminal_voltage() - ocv).abs() < 1e-12);
        // OCV at SoC=0.5 is exactly the tabulated 3.7 V.
        assert!((ocv - 3.7).abs() < 1e-12);
    }

    #[test]
    fn instantaneous_ir0_drop_under_discharge() {
        // Validation: V drops by exactly I*R0 the instant a discharge
        // current is applied (RC pair still at zero -> no V_rc yet).
        let c = cell(0.5);
        let i = 5.0; // 5 A discharge
        let v_loaded = c.terminal_voltage(i);
        let expected = c.ocv() - i * c.params().r0; // 3.7 - 5*0.01 = 3.65
        assert!((v_loaded - expected).abs() < 1e-12, "got {v_loaded}");
        assert!((v_loaded - 3.65).abs() < 1e-12, "got {v_loaded}");
        // And it is genuinely below OCV.
        assert!(v_loaded < c.ocv());
    }

    #[test]
    fn charge_current_raises_terminal_voltage() {
        // Negative current (charge) pushes terminal V above OCV by I*R0.
        let c = cell(0.5);
        let i = -5.0;
        let v_loaded = c.terminal_voltage(i);
        assert!(v_loaded > c.ocv());
        assert!((v_loaded - 3.75).abs() < 1e-12, "got {v_loaded}");
    }

    #[test]
    fn soc_falls_when_discharging() {
        // Validation: SoC decreases under a positive (discharge) current.
        let mut c = cell(0.8);
        let before = c.soc();
        c.step(2.0, 60.0).unwrap(); // 2 A for 60 s
        assert!(
            c.soc() < before,
            "soc did not fall: {} -> {}",
            before,
            c.soc()
        );
    }

    #[test]
    fn soc_rises_when_charging() {
        let mut c = cell(0.5);
        let before = c.soc();
        c.step(-2.0, 60.0).unwrap();
        assert!(c.soc() > before);
    }

    #[test]
    fn coulomb_counting_conserves_charge() {
        // Validation: ΔSoC matches the exact coulomb integral I·t / Q.
        // Draw 1 A for 1 hour from a 2 Ah cell -> exactly 0.5 SoC drop.
        let mut c = cell(1.0);
        let q = c.params().capacity_coulombs();
        let i = 1.0;
        let dt = SECONDS_PER_HOUR; // 1 hour
        c.step(i, dt).unwrap();
        let expected_drop = i * dt / q; // = 0.5
        assert!((expected_drop - 0.5).abs() < 1e-12);
        assert!(
            ((1.0 - c.soc()) - expected_drop).abs() < 1e-12,
            "soc now {}",
            c.soc()
        );
    }

    #[test]
    fn coulomb_counting_additive_over_substeps() {
        // Charge conservation: one big step == many small steps summing
        // to the same charge (the counter is exactly linear in I·dt).
        let mut one = cell(0.9);
        let mut many = cell(0.9);
        one.step(1.5, 100.0).unwrap();
        for _ in 0..100 {
            many.step(1.5, 1.0).unwrap();
        }
        assert!(
            (one.soc() - many.soc()).abs() < 1e-12,
            "{} vs {}",
            one.soc(),
            many.soc()
        );
    }

    #[test]
    fn soc_clamps_at_empty() {
        // Over-discharge cannot drive SoC below zero.
        let mut c = cell(0.05);
        c.step(10.0, SECONDS_PER_HOUR).unwrap(); // way more than available
        assert!((c.soc() - 0.0).abs() < 1e-12, "soc = {}", c.soc());
    }

    #[test]
    fn soc_clamps_at_full() {
        let mut c = cell(0.98);
        c.step(-10.0, SECONDS_PER_HOUR).unwrap();
        assert!((c.soc() - 1.0).abs() < 1e-12, "soc = {}", c.soc());
    }

    #[test]
    fn rc_relaxes_to_one_over_e_after_one_tau() {
        // Validation: with no current, V_rc decays to 1/e after τ.
        // Seed an over-potential, then relax for exactly one tau.
        let params = CellParams::new(0.01, 0.02, 1000.0, 2.0).unwrap();
        let ocv = OcvSocTable::new(vec![0.0, 1.0], vec![3.0, 4.2]).unwrap();
        let v0 = 0.1;
        let mut c = CellState::with_overpotential(params, ocv, 0.5, v0).unwrap();
        let tau = c.params().tau();
        c.step(0.0, tau).unwrap(); // relax one time constant at I=0
        let expected = v0 * (-1.0f64).exp(); // v0 / e
        assert!(
            (c.v_rc() - expected).abs() < 1e-12,
            "v_rc = {}, expected {}",
            c.v_rc(),
            expected
        );
    }

    #[test]
    fn rc_charges_toward_i_r1_steady_state() {
        // Under sustained constant current the RC over-potential climbs
        // monotonically toward its analytic steady state I*R1.
        let mut c = cell(0.5);
        let i = 3.0;
        let steady = i * c.params().r1; // 3 * 0.02 = 0.06 V
        let mut prev = c.v_rc();
        // 600 s at τ = 20 s is 30 time constants; the residual
        // e^(-30) ≈ 1e-13, so the over-potential reaches I·R1 to well
        // within the tolerance below.
        for _ in 0..600 {
            c.step(i, 1.0).unwrap();
            assert!(c.v_rc() >= prev - 1e-12, "v_rc not monotone");
            assert!(c.v_rc() <= steady + 1e-12, "v_rc overshot steady state");
            prev = c.v_rc();
        }
        // After many time constants it has essentially reached I*R1.
        assert!((c.v_rc() - steady).abs() < 1e-9, "v_rc = {}", c.v_rc());
    }

    #[test]
    fn rc_step_matches_closed_form_after_known_interval() {
        // Single step from rest: V_rc(dt) = I*R1*(1 - e^(-dt/tau)) exactly.
        let mut c = cell(0.5);
        let i = 4.0;
        let dt = 10.0;
        let tau = c.params().tau();
        c.step(i, dt).unwrap();
        let expected = i * c.params().r1 * (1.0 - (-dt / tau).exp());
        assert!((c.v_rc() - expected).abs() < 1e-12, "v_rc = {}", c.v_rc());
    }

    #[test]
    fn terminal_voltage_returned_by_step_is_start_of_step_value() {
        // The value `step` returns is the terminal V at the step's start,
        // i.e. uses the pre-step v_rc. From rest that is OCV - I*R0.
        let mut c = cell(0.5);
        let i = 5.0;
        let ocv0 = c.ocv();
        let v = c.step(i, 5.0).unwrap();
        assert!((v - (ocv0 - i * c.params().r0)).abs() < 1e-12, "v = {v}");
    }

    #[test]
    fn step_rejects_negative_dt() {
        let mut c = cell(0.5);
        let err = c.step(1.0, -1.0).unwrap_err();
        assert_eq!(err.code(), "battery_ecm.invalid");
    }

    #[test]
    fn step_rejects_non_finite_current() {
        let mut c = cell(0.5);
        let err = c.step(f64::NAN, 1.0).unwrap_err();
        assert_eq!(err.code(), "battery_ecm.invalid");
    }

    #[test]
    fn params_reject_non_positive() {
        assert!(CellParams::new(0.0, 0.02, 1000.0, 2.0).is_err());
        assert!(CellParams::new(0.01, -1.0, 1000.0, 2.0).is_err());
        assert!(CellParams::new(0.01, 0.02, 0.0, 2.0).is_err());
        assert!(CellParams::new(0.01, 0.02, 1000.0, -2.0).is_err());
        assert!(CellParams::new(f64::INFINITY, 0.02, 1000.0, 2.0).is_err());
    }

    #[test]
    fn state_rejects_soc_out_of_range() {
        let params = CellParams::new(0.01, 0.02, 1000.0, 2.0).unwrap();
        let ocv = OcvSocTable::new(vec![0.0, 1.0], vec![3.0, 4.2]).unwrap();
        let err = CellState::new(params, ocv, 1.5).unwrap_err();
        assert_eq!(err.code(), "battery_ecm.soc_out_of_range");
    }

    #[test]
    fn zero_current_step_only_relaxes_no_soc_change() {
        // I = 0: SoC is unchanged, V_rc decays. (Charge conserved.)
        let mut c = cell(0.6);
        let soc0 = c.soc();
        c.step(0.0, 5.0).unwrap();
        assert!((c.soc() - soc0).abs() < 1e-12);
    }

    #[test]
    fn zero_duration_step_is_identity() {
        // dt = 0: nothing changes, terminal V is the instantaneous value.
        let mut c = cell(0.7);
        let soc0 = c.soc();
        let vrc0 = c.v_rc();
        let v = c.step(3.0, 0.0).unwrap();
        assert!((c.soc() - soc0).abs() < 1e-12);
        assert!((c.v_rc() - vrc0).abs() < 1e-12);
        assert!((v - c.terminal_voltage(3.0)).abs() < 1e-12);
    }
}
