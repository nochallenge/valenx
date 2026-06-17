//! # valenx-mosfet
//!
//! Large-signal IV model for an n-channel enhancement-mode MOSFET, using
//! the classic square-law (Shockley level-1, long-channel) drain-current
//! equations with cutoff / triode / saturation region detection and the
//! small-signal transconductance `gm`.
//!
//! ## What
//!
//! A tiny, dependency-light device model: build a [`Mosfet`] from its
//! lumped transconductance parameter `k` (A/V²) and threshold voltage
//! `vth` (V), then evaluate, for any `(vgs, vds)` bias:
//!
//! - the operating [`Region`] ([`Mosfet::region`]),
//! - the drain current `Id` ([`Mosfet::drain_current`]),
//! - the transconductance `gm = dId/dVgs` ([`Mosfet::gm`]), and
//! - all three at once as an [`OperatingPoint`]
//!   ([`Mosfet::operating_point`]).
//!
//! Going the other way, the **bias-design inverse**
//! [`Mosfet::vgs_for_saturation_current`] (and its
//! [`Mosfet::overdrive_for_saturation_current`] core) returns the gate
//! voltage that carries a target saturation drain current,
//! `vov = sqrt(2·Id/k)`.
//!
//! ## Model
//!
//! With gate overdrive `vov = vgs − vth`:
//!
//! - **Cutoff** (`vov ≤ 0`): `Id = 0`.
//! - **Triode / linear** (`vds < vov`):
//!   `Id = k · (vov · vds − ½ · vds²)`.
//! - **Saturation** (`vds ≥ vov`): `Id = ½ · k · vov²`.
//! - **Transconductance** (saturation): `gm = k · vov` (and `0` in
//!   cutoff, where `Id ≡ 0`).
//! - **Saturation bias inverse**: `vov = sqrt(2 · Id / k)`,
//!   `vgs = vth + vov` — the gate bias for a target saturation `Id`.
//!
//! The triode and saturation expressions are continuous at the
//! pinch-off boundary `vds = vov`. Here `k = μ_n · C_ox · W / L` is the
//! process-and-geometry transconductance parameter. Voltages are in
//! volts, currents in amperes, `k` in A/V², and `gm` in siemens.
//!
//! ## Honest scope
//!
//! Research / educational grade. This crate implements the **idealized
//! closed-form long-channel square-law model only**. It deliberately
//! omits every second-order effect a real SPICE level-2/3 or BSIM model
//! includes: channel-length modulation (the `(1 + λ·vds)` factor), the
//! body / substrate-bias effect on `vth`, velocity saturation, mobility
//! degradation, subthreshold (weak-inversion) conduction, drain-induced
//! barrier lowering, gate leakage, temperature dependence, and all
//! short-channel corrections. It is **NOT** a clinical/medical or
//! production semiconductor-engineering tool and must not be used for
//! tape-out, circuit sign-off, or any safety-critical decision. Use it
//! to learn the square-law, sanity-check hand calculations, or seed a
//! richer model — not to replace one.
//!
//! ## Example
//!
//! ```
//! use valenx_mosfet::{Mosfet, Region};
//!
//! // A textbook NMOS: k = 0.5 mA/V^2, Vth = 1.0 V.
//! let m = Mosfet::nmos();
//!
//! // Below threshold: cut off, zero current.
//! assert_eq!(m.region(0.5, 2.0).unwrap(), Region::Cutoff);
//! assert!(m.drain_current(0.5, 2.0).unwrap().abs() < 1e-15);
//!
//! // Vgs = 3 V (overdrive 2 V), Vds = 5 V >= overdrive => saturation.
//! let op = m.operating_point(3.0, 5.0).unwrap();
//! assert_eq!(op.region, Region::Saturation);
//! // Id = 0.5 * k * vov^2 = 0.5 * 0.5e-3 * 4 = 1.0 mA.
//! assert!((op.id - 1.0e-3).abs() < 1e-12);
//! // gm = k * vov = 0.5e-3 * 2 = 1.0 mS.
//! assert!((op.gm - 1.0e-3).abs() < 1e-12);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod device;
pub mod error;

pub use device::{Mosfet, OperatingPoint, Region};
pub use error::{ErrorCategory, MosfetError, Result};
