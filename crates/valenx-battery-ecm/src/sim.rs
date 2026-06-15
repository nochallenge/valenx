//! Constant-current discharge / charge driver and a worked example cell.
//!
//! ## Model
//!
//! [`simulate_constant_current`] repeatedly applies [`CellState::step`]
//! with a fixed current and a fixed time step, recording one
//! [`SamplePoint`] per step. It is a thin convenience over the cell's
//! own integrator — every physical decision (the `I·R0` drop, the
//! coulomb-counted SoC, the exact RC relaxation) lives in
//! [`crate::thevenin`]; this module only schedules the steps and stops
//! when the run finishes or the configured cut-off voltage is reached.
//!
//! The discharge stops early if the loaded terminal voltage falls to or
//! below `cutoff_v` — the textbook way a constant-current discharge ends
//! at a protection threshold rather than at empty.

use crate::error::{EcmError, Result};
use crate::ocv::OcvSocTable;
use crate::thevenin::{CellParams, CellState};
use serde::{Deserialize, Serialize};

/// One recorded instant of a simulated profile.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SamplePoint {
    /// Elapsed time at the *start* of this step (seconds).
    pub t: f64,
    /// State of charge at the start of this step (fraction in `[0, 1]`).
    pub soc: f64,
    /// Terminal voltage seen at the start of this step (volts).
    pub terminal_v: f64,
    /// RC over-potential at the start of this step (volts).
    pub v_rc: f64,
}

/// Settings for a constant-current run.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct DischargeConfig {
    /// Applied current (amperes, `> 0` discharges, `< 0` charges).
    pub current: f64,
    /// Fixed time step (seconds, `> 0`).
    pub dt: f64,
    /// Number of steps to attempt.
    pub steps: usize,
    /// Lower terminal-voltage cut-off (volts); the run stops once the
    /// loaded terminal voltage reaches or falls below this. Use
    /// `f64::NEG_INFINITY` to disable.
    pub cutoff_v: f64,
}

/// Run a constant-current profile and return the per-step samples.
///
/// The cell is advanced in place. Each [`SamplePoint`] holds the state
/// *entering* its step; the run ends after `cfg.steps` steps or as soon
/// as the loaded terminal voltage reaches `cfg.cutoff_v`, whichever comes
/// first. At least the initial sample is always returned.
///
/// # Errors
///
/// Returns [`EcmError::Invalid`] if `cfg.dt` is not strictly positive or
/// `cfg.current` is non-finite. Per-step errors from
/// [`CellState::step`] are propagated.
pub fn simulate_constant_current(
    cell: &mut CellState,
    cfg: &DischargeConfig,
) -> Result<Vec<SamplePoint>> {
    if !cfg.dt.is_finite() || cfg.dt <= 0.0 {
        return Err(EcmError::invalid("dt", "must be strictly positive"));
    }
    if !cfg.current.is_finite() {
        return Err(EcmError::invalid("current", "must be finite"));
    }

    let mut out = Vec::with_capacity(cfg.steps + 1);
    let mut t = 0.0;
    for _ in 0..cfg.steps {
        let terminal_v = cell.terminal_voltage(cfg.current);
        out.push(SamplePoint {
            t,
            soc: cell.soc(),
            terminal_v,
            v_rc: cell.v_rc(),
        });
        if terminal_v <= cfg.cutoff_v {
            return Ok(out);
        }
        cell.step(cfg.current, cfg.dt)?;
        t += cfg.dt;
    }
    // Final state after the last step, so callers see where it ended.
    out.push(SamplePoint {
        t,
        soc: cell.soc(),
        terminal_v: cell.terminal_voltage(cfg.current),
        v_rc: cell.v_rc(),
    });
    Ok(out)
}

/// A worked example: a generic ~2 Ah lithium-ion-like cell.
///
/// Parameters are representative round numbers for a small cylindrical
/// cell (10 mΩ ohmic resistance, a 15 mΩ / 2000 F polarisation pair,
/// τ = 30 s) on a smooth 3.0–4.2 V rest-voltage curve. Intended for
/// demos and tests, **not** as identified parameters of any real
/// product.
pub fn example_li_ion_cell(soc0: f64) -> Result<CellState> {
    let params = CellParams::new(0.010, 0.015, 2000.0, 2.0)?;
    let ocv = OcvSocTable::new(
        vec![0.00, 0.10, 0.30, 0.50, 0.70, 0.90, 1.00],
        vec![3.00, 3.45, 3.60, 3.70, 3.85, 4.05, 4.20],
    )?;
    CellState::new(params, ocv, soc0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn discharge_lowers_soc_and_voltage_over_time() {
        let mut c = example_li_ion_cell(1.0).unwrap();
        let cfg = DischargeConfig {
            current: 2.0,
            dt: 10.0,
            steps: 50,
            cutoff_v: f64::NEG_INFINITY,
        };
        let samples = simulate_constant_current(&mut c, &cfg).unwrap();
        assert_eq!(samples.len(), cfg.steps + 1);

        // SoC is monotonically non-increasing across a discharge.
        for w in samples.windows(2) {
            assert!(
                w[1].soc <= w[0].soc + 1e-12,
                "soc rose: {} -> {}",
                w[0].soc,
                w[1].soc
            );
        }
        // The cell ends below where it started, in both SoC and voltage.
        assert!(samples.last().unwrap().soc < samples[0].soc);
        assert!(samples.last().unwrap().terminal_v < samples[0].terminal_v);
    }

    #[test]
    fn first_sample_shows_instant_ir0_drop_from_rest() {
        // Entering the discharge from rest, the very first terminal
        // voltage is OCV - I*R0 (no RC over-potential yet).
        let mut c = example_li_ion_cell(0.7).unwrap();
        let ocv0 = c.ocv();
        let r0 = c.params().r0;
        let i = 3.0;
        let cfg = DischargeConfig {
            current: i,
            dt: 1.0,
            steps: 5,
            cutoff_v: f64::NEG_INFINITY,
        };
        let samples = simulate_constant_current(&mut c, &cfg).unwrap();
        assert!((samples[0].v_rc - 0.0).abs() < 1e-12);
        assert!(
            (samples[0].terminal_v - (ocv0 - i * r0)).abs() < 1e-12,
            "got {}",
            samples[0].terminal_v
        );
    }

    #[test]
    fn cutoff_stops_the_run_early() {
        // A high cut-off voltage trips immediately at the first sample.
        let mut c = example_li_ion_cell(0.5).unwrap();
        let v0 = c.terminal_voltage(2.0);
        let cfg = DischargeConfig {
            current: 2.0,
            dt: 10.0,
            steps: 100,
            cutoff_v: v0 + 0.01, // already below cut-off at step 0
        };
        let samples = simulate_constant_current(&mut c, &cfg).unwrap();
        assert_eq!(samples.len(), 1, "should stop at the first sample");
    }

    #[test]
    fn ah_throughput_matches_soc_drop() {
        // Charge conservation end-to-end: the ampere-hours pushed through
        // equal the SoC drop times the capacity.
        let mut c = example_li_ion_cell(1.0).unwrap();
        let cap_ah = c.params().capacity_ah;
        let i = 1.0;
        let dt = 36.0; // s
        let steps = 100usize; // total 3600 s = 1 h at 1 A -> 1 Ah
        let cfg = DischargeConfig {
            current: i,
            dt,
            steps,
            cutoff_v: f64::NEG_INFINITY,
        };
        let soc_before = c.soc();
        simulate_constant_current(&mut c, &cfg).unwrap();
        let ah_drawn = i * (dt * steps as f64) / 3600.0; // 1.0 Ah
        let soc_drop = soc_before - c.soc();
        assert!(
            (soc_drop - ah_drawn / cap_ah).abs() < 1e-12,
            "drop {soc_drop}"
        );
    }

    #[test]
    fn rejects_non_positive_dt() {
        let mut c = example_li_ion_cell(0.5).unwrap();
        let cfg = DischargeConfig {
            current: 1.0,
            dt: 0.0,
            steps: 10,
            cutoff_v: f64::NEG_INFINITY,
        };
        assert!(simulate_constant_current(&mut c, &cfg).is_err());
    }
}
