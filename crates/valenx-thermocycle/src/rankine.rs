//! The basic ideal Rankine cycle — the ideal steam power plant.
//!
//! Unlike the gas cycles in this crate, the Rankine cycle's working fluid
//! changes phase (water ⇌ steam), so its efficiency cannot be written as
//! a closed form in a single ratio — it depends on the actual fluid
//! enthalpies at the four corner states, which in practice are looked up
//! from steam tables. This module therefore works from **supplied
//! enthalpies** (kJ/kg) rather than computing them.
//!
//! The four states of the simple ideal cycle are:
//!
//! - **1** — saturated liquid leaving the condenser (pump inlet).
//! - **2** — compressed liquid leaving the pump (boiler inlet).
//! - **3** — superheated (or saturated) steam leaving the boiler
//!   (turbine inlet).
//! - **4** — wet steam leaving the turbine (condenser inlet).
//!
//! The energy balance over one cycle gives
//!
//! ```text
//! w_turbine = h3 - h4      (work out of the turbine)
//! w_pump    = h2 - h1      (work into the pump)
//! q_in      = h3 - h2      (heat into the boiler)
//! w_net     = w_turbine - w_pump
//! η_rankine = w_net / q_in
//! ```
//!
//! For an ideal pump the pump work can instead be approximated from the
//! incompressible relation `w_pump ≈ v₁ (P₂ - P₁)`, but here we take it
//! directly from the enthalpy difference `h2 - h1` so the calculation is
//! exact for whatever state data the caller provides.

use crate::error::{CycleError, Result};
use serde::{Deserialize, Serialize};

/// The four corner-state specific enthalpies (kJ/kg) of a basic ideal
/// Rankine cycle.
///
/// Fields are named for the conventional state numbering described in the
/// [module documentation](self).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Rankine {
    /// State 1 — saturated liquid at the condenser exit / pump inlet
    /// (kJ/kg).
    pub h1: f64,
    /// State 2 — compressed liquid at the pump exit / boiler inlet
    /// (kJ/kg).
    pub h2: f64,
    /// State 3 — steam at the boiler exit / turbine inlet (kJ/kg).
    pub h3: f64,
    /// State 4 — wet steam at the turbine exit / condenser inlet (kJ/kg).
    pub h4: f64,
}

impl Rankine {
    /// Build a Rankine cycle from its four state enthalpies in kJ/kg.
    ///
    /// # Errors
    ///
    /// - [`CycleError::NotFinite`] if any enthalpy is non-finite.
    /// - [`CycleError::NoHeatInput`] if the boiler heat input
    ///   `q_in = h3 - h2` is not strictly positive, which would make the
    ///   efficiency undefined (the cycle absorbs no heat) — for a real
    ///   cycle the turbine-inlet enthalpy `h3` always exceeds the
    ///   pump-outlet enthalpy `h2`.
    pub fn new(h1: f64, h2: f64, h3: f64, h4: f64) -> Result<Self> {
        let h1 = crate::error::finite("h1", h1)?;
        let h2 = crate::error::finite("h2", h2)?;
        let h3 = crate::error::finite("h3", h3)?;
        let h4 = crate::error::finite("h4", h4)?;
        // `h2` and `h3` are both finite (checked above), so `q_in` is
        // finite and this `<=` comparison can never see a NaN.
        let q_in = h3 - h2;
        if q_in <= 0.0 {
            return Err(CycleError::NoHeatInput { q_in });
        }
        Ok(Self { h1, h2, h3, h4 })
    }

    /// Specific work delivered by the turbine, `w_turbine = h3 - h4`
    /// (kJ/kg).
    pub fn turbine_work(self) -> f64 {
        self.h3 - self.h4
    }

    /// Specific work consumed by the feed pump, `w_pump = h2 - h1`
    /// (kJ/kg).
    pub fn pump_work(self) -> f64 {
        self.h2 - self.h1
    }

    /// Specific heat added in the boiler, `q_in = h3 - h2` (kJ/kg).
    ///
    /// Guaranteed strictly positive by [`Rankine::new`].
    pub fn heat_in(self) -> f64 {
        self.h3 - self.h2
    }

    /// Specific heat rejected in the condenser, `q_out = h4 - h1`
    /// (kJ/kg).
    pub fn heat_out(self) -> f64 {
        self.h4 - self.h1
    }

    /// Net specific work, `w_net = w_turbine - w_pump` (kJ/kg).
    ///
    /// By the first law this also equals `q_in - q_out`, a relation the
    /// crate's tests check directly.
    pub fn net_work(self) -> f64 {
        self.turbine_work() - self.pump_work()
    }

    /// The thermal efficiency of the cycle, `η = w_net / q_in`.
    ///
    /// The back-work ratio (pump work relative to turbine work) is small
    /// for steam cycles, so this is typically a little below the turbine-
    /// only ratio `(h3 - h4) / (h3 - h2)`.
    pub fn efficiency(self) -> f64 {
        self.net_work() / self.heat_in()
    }

    /// The **back-work ratio** `w_pump / w_turbine` — the fraction of the
    /// turbine's gross output that is recirculated to drive the feed
    /// pump. Small (a few percent) for steam Rankine cycles, which is one
    /// of their practical advantages over gas cycles.
    pub fn back_work_ratio(self) -> f64 {
        self.pump_work() / self.turbine_work()
    }
}

/// The ideal Rankine thermal efficiency from the four state enthalpies
/// (kJ/kg), `η = ((h3 - h4) - (h2 - h1)) / (h3 - h2)`.
///
/// Convenience wrapper around [`Rankine::new`] + [`Rankine::efficiency`].
///
/// # Errors
///
/// Propagates the validation errors of [`Rankine::new`].
pub fn rankine_efficiency(h1: f64, h2: f64, h3: f64, h4: f64) -> Result<f64> {
    Ok(Rankine::new(h1, h2, h3, h4)?.efficiency())
}
