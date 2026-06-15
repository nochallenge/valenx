//! # valenx-resistor-network
//!
//! Closed-form helpers for linear resistor networks: series and
//! parallel equivalent resistance, voltage and current dividers, and
//! the Wheatstone-bridge balance condition.
//!
//! ## What
//!
//! A small, dependency-light toolbox for the textbook DC
//! resistor-network identities, with validated inputs and exact
//! ground-truth tests:
//!
//! - [`combination::series`] / [`combination::parallel`] /
//!   [`combination::parallel_pair`] — equivalent resistance of
//!   resistors in series or parallel.
//! - [`divider::voltage_divider`] — `Vout = Vin * R2 / (R1 + R2)`.
//! - [`divider::current_divider_i1`] /
//!   [`divider::current_divider_i2`] — the two-branch current split.
//! - [`bridge::is_balanced`] / [`bridge::detector_voltage`] /
//!   [`bridge::balancing_r4`] — Wheatstone-bridge balance and output.
//! - [`network::Combination`] — a serde-serializable description of
//!   a one-level series-or-parallel combination.
//!
//! ## Model
//!
//! Every resistor is treated as an ideal, real-valued, lumped DC
//! resistance (no reactance, no temperature drift, no tolerance, no
//! wire resistance). The governing relations are Ohm's law and
//! Kirchhoff's laws applied to these idealised elements:
//!
//! - Series: `R_eq = sum(Ri)`.
//! - Parallel: `1 / R_eq = sum(1/Ri)`.
//! - Voltage divider: `Vout = Vin * R2 / (R1 + R2)`.
//! - Current divider: `I1 = I_in * R2 / (R1 + R2)`.
//! - Wheatstone balance: `R1 / R2 = R3 / R4`.
//!
//! Resistances must be finite and strictly positive; source terms
//! (`Vin`, `I_in`, `Vex`) must be finite but may be zero or
//! negative. Out-of-domain inputs return an [`error::ResistorError`].
//!
//! ## Honest scope
//!
//! Research/educational grade. These are the standard closed-form /
//! numerical textbook models for ideal linear resistor networks.
//! This crate is NOT a clinical/medical tool and NOT a production
//! electrical-engineering tool: it ignores component tolerances,
//! power dissipation and thermal limits, parasitic inductance and
//! capacitance, AC/frequency effects, source impedance, and safety
//! margins. Do not use it to size real circuits or make
//! safety-relevant decisions.
//!
//! ## Example
//!
//! ```
//! use valenx_resistor_network::combination::{parallel, series};
//! use valenx_resistor_network::divider::voltage_divider;
//!
//! // 100 ohm || 100 ohm = 50 ohm, then in series with 150 ohm = 200 ohm.
//! let pair = parallel(&[100.0, 100.0]).unwrap();
//! let total = series(&[pair, 150.0]).unwrap();
//! assert!((total - 200.0).abs() < 1e-9);
//!
//! // A 12 V source across an equal-leg divider yields 6 V.
//! let vout = voltage_divider(12.0, 1000.0, 1000.0).unwrap();
//! assert!((vout - 6.0).abs() < 1e-9);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod bridge;
pub mod combination;
pub mod divider;
pub mod error;
pub mod network;

pub use bridge::{balancing_r4, detector_voltage, is_balanced, DEFAULT_BALANCE_TOL};
pub use combination::{parallel, parallel_pair, series};
pub use divider::{current_divider_i1, current_divider_i2, voltage_divider};
pub use error::{ErrorCategory, ResistorError};
pub use network::Combination;
