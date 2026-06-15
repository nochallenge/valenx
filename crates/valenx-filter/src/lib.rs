//! # valenx-filter
//!
//! Closed-form **analog filter** design: the textbook RC and series-RLC
//! formulas that every circuits course starts from, packaged as
//! validated, queryable Rust types.
//!
//! ## What
//!
//! Two component-level filter models:
//!
//! - [`rc::RcFilter`] — a first-order **RC** section (a resistor and a
//!   capacitor). Gives the cutoff frequency, and the linear / decibel
//!   magnitude and phase of the low-pass (output across `C`) or
//!   high-pass (output across `R`) response at any frequency.
//! - [`rlc::RlcCircuit`] — a series **RLC** resonant circuit. Gives the
//!   resonant frequency, the quality factor `Q`, and the `-3 dB`
//!   bandwidth.
//!
//! Both return a [`response::Response`] (magnitude + phase) where a
//! per-frequency answer is wanted, and both reject non-physical inputs
//! up front via [`error::FilterError`].
//!
//! ```
//! use valenx_filter::{RcFilter, RlcCircuit};
//!
//! // A 1 kOhm + 1 uF low-pass: corner at ~159 Hz, -3 dB there.
//! let lp = RcFilter::low_pass(1_000.0, 1e-6).unwrap();
//! assert!((lp.cutoff_hz() - 159.154_9).abs() < 1e-3);
//! let mag = lp.magnitude(lp.cutoff_hz()).unwrap();
//! assert!((mag - 1.0 / 2.0_f64.sqrt()).abs() < 1e-9); // 0.7071
//!
//! // A series RLC: resonant frequency, Q, and bandwidth.
//! let tank = RlcCircuit::new(10.0, 10e-3, 1e-6).unwrap();
//! assert!((tank.quality_factor() - 10.0).abs() < 1e-9);
//! assert!((tank.bandwidth_hz() - tank.resonant_hz() / tank.quality_factor()).abs() < 1e-9);
//! ```
//!
//! ## Model
//!
//! The formulas are the standard linear-circuit results for ideal,
//! lumped, lossless reactive components:
//!
//! ```text
//! RC cutoff:        fc = 1 / (2 * pi * R * C)
//! RC low-pass:      |H(f)| = 1 / sqrt(1 + (f/fc)^2)
//! RC high-pass:     |H(f)| = (f/fc) / sqrt(1 + (f/fc)^2)
//! RLC resonance:    f0 = 1 / (2 * pi * sqrt(L * C))
//! RLC Q (series):   Q  = (1 / R) * sqrt(L / C)
//! RLC bandwidth:    BW = f0 / Q
//! ```
//!
//! At `f = fc` an RC section is at exactly `1/sqrt(2)` of its passband
//! gain (`-3 dB`, the half-power point). A higher RLC `Q` narrows the
//! bandwidth `BW = f0 / Q`, sharpening the resonance.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are **textbook closed-form**
//! models of ideal components, not a circuit simulator and **NOT a
//! clinical, medical, or production engineering tool**. Specifically:
//!
//! - Components are **ideal** — resistance is frequency-independent,
//!   inductors and capacitors are lossless and parasitic-free, and
//!   there is no source / load impedance, tolerance, or temperature
//!   drift. No SPICE-style netlist solve.
//! - Only **first-order RC** and **series-RLC** topologies are modelled.
//!   There is no cascading, no active (op-amp) filters, and no
//!   higher-order Butterworth / Chebyshev / elliptic synthesis.
//! - The RLC `-3 dB` edge frequencies use the symmetric narrow-band
//!   approximation `f0 +/- BW/2`, which is exact only in the
//!   high-`Q` limit. The resonant frequency, `Q`, and bandwidth
//!   themselves are exact.
//!
//! Within those bounds every number is the genuine analytic result and
//! is checked against the closed-form formula in the test suite.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod rc;
pub mod response;
pub mod rlc;

pub use error::{FilterError, Result};
pub use rc::{RcFilter, RcKind};
pub use response::Response;
pub use rlc::RlcCircuit;
