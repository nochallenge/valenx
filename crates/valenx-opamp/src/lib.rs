//! # valenx-opamp
//!
//! Closed-form operational-amplifier circuit models.
//!
//! ## What
//!
//! Small, dependency-light building blocks for the four canonical
//! op-amp topologies plus the single-pole gain-bandwidth relations:
//!
//! 1. [`Inverting`] — gain `G = -Rf / Rin`.
//! 2. [`NonInverting`] — gain `G = 1 + Rf / Rin` (always `>= 1`).
//! 3. [`VoltageFollower`] — the `Rf = 0` unity-gain buffer, `G = 1`.
//! 4. [`SummingAmplifier`] — `Vout = -Rf · Σ(Vᵢ / Rᵢ)`.
//! 5. [`Gbw`] — gain-bandwidth product, closed-loop bandwidth
//!    (`GBW / |gain|`), its dual the maximum gain for a required
//!    bandwidth (`GBW / bandwidth`), and unity-gain bandwidth (`= GBW`).
//!
//! Every numeric constructor validates its inputs and returns
//! [`OpAmpError`] rather than panicking.
//!
//! ## Model
//!
//! The amplifier topologies assume an **ideal op-amp**: infinite
//! open-loop gain, infinite input impedance, zero output impedance.
//! Negative feedback then forces `V+ = V-` and no current enters the
//! inputs, collapsing each circuit to a resistor ratio. The bandwidth
//! relations additionally assume the open-loop response is dominated by
//! a **single pole**, which makes the gain-bandwidth product constant.
//! All resistances are in ohms, voltages in volts, frequencies in
//! hertz.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are textbook closed-form and
//! numerical models of idealised circuits; they deliberately ignore
//! finite open-loop gain, input bias current and offset voltage,
//! multiple poles, slew-rate limiting, noise, output saturation and
//! temperature drift. This crate is NOT a clinical, medical, or
//! production electronic-engineering design tool, and its outputs must
//! not be used to qualify real hardware.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod bandwidth;
pub mod error;
pub mod ideal;

pub use bandwidth::Gbw;
pub use error::{OpAmpError, Result};
pub use ideal::{Inverting, NonInverting, SummingAmplifier, SummingInput, VoltageFollower};
