//! # valenx-threephase
//!
//! Balanced three-phase AC power relations: wye (star) and delta (mesh)
//! line/phase conversions plus real-power computation.
//!
//! ## What
//!
//! Closed-form helpers for the textbook quantities of a *balanced*
//! three-phase system:
//!
//! - [`Connection`] — wye or delta wiring, with line/phase voltage and
//!   current conversions in both directions.
//! - [`BalancedLoad`] — a validated load (connection plus per-element
//!   voltage, current, and power factor) that derives line quantities,
//!   per-phase power, and total power.
//! - [`power_from_line`] — total real power directly from line-to-line
//!   voltage, line current, and power factor.
//! - [`SQRT_3`] — the `sqrt(3)` line-to-phase conversion constant.
//!
//! ## Model
//!
//! The system is assumed *balanced*: three identical loads, equal
//! phase magnitudes, 120 degrees apart, with all magnitudes given as
//! RMS values. Under that assumption the relations reduce to:
//!
//! - Wye (star): `V_line = sqrt(3) * V_phase`, `I_line = I_phase`.
//! - Delta (mesh): `V_line = V_phase`, `I_line = sqrt(3) * I_phase`.
//! - Real power: `P = sqrt(3) * V_line * I_line * cos(phi)`, which for a
//!   balanced load is identically three times the per-phase power
//!   `V_phase * I_phase * cos(phi)`.
//!
//! The angle `phi` is the displacement between phase voltage and phase
//! current; `cos(phi)` is the power factor, constrained to `[-1, 1]`.
//!
//! ## Honest scope
//!
//! Research/educational grade. These are steady-state, single-frequency
//! phasor models for *balanced* systems only: there is no treatment of
//! unbalance, harmonics, transients, neutral currents, mutual coupling,
//! transmission-line effects, or protection. This crate is NOT a
//! clinical/medical/production engineering tool and must not be relied
//! on for the design, sizing, or safety certification of real
//! electrical installations.

#![forbid(unsafe_code)]
#![warn(missing_docs)]
#![warn(unused_imports)]

pub mod error;
pub mod threephase;

pub use error::{ErrorCategory, ThreePhaseError};
pub use threephase::{power_from_line, BalancedLoad, Connection, SQRT_3};
