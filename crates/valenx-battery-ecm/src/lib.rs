//! # valenx-battery-ecm — battery equivalent-circuit model
//!
//! A native, dependency-light **first-order Thevenin** battery-cell
//! model: open-circuit voltage minus an ohmic drop minus one relaxing
//! RC-pair over-potential, with state of charge tracked by coulomb
//! counting and the rest voltage read from a monotone OCV-SoC table.
//!
//! ## What
//!
//! Describe a cell by its [`thevenin::CellParams`] (series resistance
//! `R0`, polarisation pair `R1`/`C1`, capacity in ampere-hours) and an
//! [`ocv::OcvSocTable`] (the rest-voltage-versus-charge curve), pick an
//! initial state of charge, and you get a [`thevenin::CellState`] you can
//! drive:
//!
//! - [`CellState::terminal_voltage`](thevenin::CellState::terminal_voltage)
//!   — the loaded terminal voltage at the present instant for a given
//!   current.
//! - [`CellState::step`](thevenin::CellState::step) — advance the cell by
//!   `Δt` under a constant current, updating SoC and the RC
//!   over-potential.
//! - [`simulate_constant_current`] — run a full constant-current
//!   discharge / charge profile (with an optional voltage cut-off) and
//!   read back a per-step trace.
//!
//! ```
//! use valenx_battery_ecm::{example_li_ion_cell, DischargeConfig, simulate_constant_current};
//!
//! // A worked ~2 Ah cell, starting full.
//! let mut cell = example_li_ion_cell(1.0).expect("valid cell");
//!
//! // 2 A constant-current discharge, 10 s steps, stop at 3.0 V.
//! let cfg = DischargeConfig { current: 2.0, dt: 10.0, steps: 500, cutoff_v: 3.0 };
//! let trace = simulate_constant_current(&mut cell, &cfg).expect("valid run");
//!
//! // Voltage and SoC both fall over the discharge.
//! assert!(trace.last().unwrap().terminal_v < trace[0].terminal_v);
//! assert!(trace.last().unwrap().soc < trace[0].soc);
//! ```
//!
//! ## Model
//!
//! The cell is the standard **one-RC Thevenin equivalent circuit**. With
//! current `I` taken **positive on discharge**:
//!
//! ```text
//! V_terminal = OCV(SoC) - I·R0 - V_rc                 (terminal voltage)
//! SoC(t)     = SoC0 - (1/Q) ∫ I dt                     (coulomb counting)
//! V_rc(t+Δt) = V_rc·e^(-Δt/τ) + I·R1·(1 - e^(-Δt/τ)),  τ = R1·C1   (RC pair)
//! ```
//!
//! `Q = capacity_Ah · 3600` is the usable capacity in coulombs. The
//! ohmic term `I·R0` is the *instantaneous* step the moment a load is
//! applied; the RC term `V_rc` is the slowly-relaxing diffusion
//! over-potential, integrated with the **exact exponential map** (so it
//! is unconditionally stable and decays to `1/e` after one `τ`). `OCV`
//! is read from a piecewise-linear, strictly-increasing-in-SoC,
//! non-decreasing-in-voltage [`ocv::OcvSocTable`].
//!
//! These behaviours are pinned by the per-module tests: terminal voltage
//! equals OCV at `I = 0`; it drops by exactly `I·R0` at the instant of a
//! step change; SoC falls (and conserves charge) on discharge; the RC
//! over-potential relaxes with `τ = R1·C1`; and the OCV interpolation is
//! monotone.
//!
//! ## Honest scope
//!
//! Research / educational grade. Every formula here is the textbook
//! closed-form / well-established lumped-parameter model — the canonical
//! first-order Thevenin ECM that appears in every battery-modelling text
//! and is the workhorse of state-of-charge estimators. It is **not** a
//! clinical / medical / production engineering tool, and in particular:
//!
//! - **One RC pair only.** Real high-fidelity ECMs use two or three RC
//!   pairs (fast + slow diffusion) or a full impedance / fractional-order
//!   model; this crate ships the single-pair case.
//! - **No temperature.** Every parameter is constant. There is no
//!   thermal model, no Arrhenius temperature dependence of `R0`/`R1`, and
//!   no self-heating — a real cell's resistances vary strongly with
//!   temperature.
//! - **No ageing, no hysteresis.** Capacity and resistances do not drift
//!   with cycle count, and the OCV curve has no charge / discharge
//!   hysteresis loop.
//! - **No parameter identification.** The parameters and the example
//!   OCV curve are representative round numbers, **not** values fitted to
//!   a measured cell. Coulomb counting also assumes a perfectly known
//!   current with no sensor bias or coulombic-efficiency loss.
//! - **No protection / safety logic.** The optional voltage cut-off in
//!   [`simulate_constant_current`] is a convenience stop, not a battery
//!   management system; nothing here models over-current, over-voltage or
//!   thermal-runaway protection.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod ocv;
pub mod sim;
pub mod thevenin;

pub use error::{EcmError, ErrorCategory, Result};
pub use ocv::OcvSocTable;
pub use sim::{example_li_ion_cell, simulate_constant_current, DischargeConfig, SamplePoint};
pub use thevenin::{CellParams, CellState, SECONDS_PER_HOUR};

#[cfg(test)]
mod tests {
    use super::*;

    /// End-to-end: a full constant-current discharge of the example cell
    /// hits its 3.0 V cut-off, and every validated behaviour holds along
    /// the way — the entry voltage is the instantaneous `OCV - I·R0`, SoC
    /// and voltage fall monotonically, and the recovered ampere-hours
    /// match the SoC drop.
    #[test]
    fn full_discharge_end_to_end() {
        let mut cell = example_li_ion_cell(1.0).unwrap();
        let ocv_full = cell.ocv();
        let r0 = cell.params().r0;
        let cap_ah = cell.params().capacity_ah;
        let i = 1.0;
        let dt = 30.0;

        let cfg = DischargeConfig {
            current: i,
            dt,
            steps: 100_000, // generously many; the cut-off ends it first
            cutoff_v: 3.0,
        };
        let trace = simulate_constant_current(&mut cell, &cfg).unwrap();

        // First sample: instantaneous I*R0 drop from rest, no RC yet.
        assert!((trace[0].v_rc - 0.0).abs() < 1e-12);
        assert!((trace[0].terminal_v - (ocv_full - i * r0)).abs() < 1e-12);

        // The run actually terminated on the cut-off, below full SoC.
        let last = *trace.last().unwrap();
        assert!(
            last.terminal_v <= 3.0 + 1e-9,
            "ended at {}",
            last.terminal_v
        );
        assert!(last.soc < 1.0);

        // Monotone discharge: SoC and terminal voltage never rise.
        for w in trace.windows(2) {
            assert!(w[1].soc <= w[0].soc + 1e-12);
        }

        // Charge bookkeeping: Ah drawn up to the final sample equals the
        // SoC drop times capacity (within one step's worth of charge).
        let ah_drawn = i * last.t / SECONDS_PER_HOUR;
        let soc_drop = 1.0 - last.soc;
        let step_ah = i * dt / SECONDS_PER_HOUR;
        assert!(
            (ah_drawn - soc_drop * cap_ah).abs() < step_ah + 1e-9,
            "ah_drawn={ah_drawn}, soc_drop*cap={}",
            soc_drop * cap_ah
        );
    }

    /// The whole cell state round-trips through JSON (serde derives on
    /// every public data type).
    #[test]
    fn cell_state_serde_round_trip() {
        let cell = example_li_ion_cell(0.6).unwrap();
        let json = serde_json::to_string(&cell).unwrap();
        let back: CellState = serde_json::from_str(&json).unwrap();
        assert_eq!(cell, back);
    }
}
