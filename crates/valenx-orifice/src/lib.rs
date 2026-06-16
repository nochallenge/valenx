//! # valenx-orifice
//!
//! Differential-pressure flow-meter sizing from the canonical
//! incompressible Bernoulli relation.
//!
//! ## What
//!
//! Given the geometry of a constriction in a pipe (a throat of diameter
//! `d` inside a pipe of diameter `D`), a fluid density `rho`, and a
//! discharge coefficient `Cd`, this crate computes the volumetric flow
//! rate `Q` from a measured pressure drop `dP` (the *forward* direction
//! a real meter performs) and, conversely, the `dP` needed to drive a
//! target `Q` (the *inverse* used to size the transmitter). It covers
//! orifice plates, flow nozzles, and Venturi tubes, each with a nominal
//! discharge coefficient, converts to mass flow, and reports the
//! permanent (unrecovered) pressure loss of a square-edged orifice plate
//! — the pumping penalty the meter imposes downstream.
//!
//! ## Model
//!
//! Bernoulli plus continuity for an incompressible stream, with all the
//! real-fluid losses lumped into the empirical discharge coefficient
//! `Cd`, give the working equation
//!
//! ```text
//! Q = Cd * A * sqrt( 2 * dP / (rho * (1 - beta^4)) )
//! ```
//!
//! with
//!
//! ```text
//! A    = pi * d^2 / 4          throat area              [m^2]
//! beta = d / D                 diameter ratio            [-]
//! E    = 1 / sqrt(1 - beta^4)  velocity-of-approach      [-]
//! ```
//!
//! The same relation inverted gives the exact pressure drop
//!
//! ```text
//! dP = (rho * (1 - beta^4) / 2) * ( Q / (Cd * A) )^2,
//! ```
//!
//! and mass flow follows as `mdot = rho * Q`. The forward and inverse
//! solutions round-trip to within floating-point precision.
//!
//! The measured differential `dP` is not all lost: downstream of the
//! plate the jet re-expands and recovers part of it. The fraction that
//! is *permanently* lost — the pumping penalty — follows the ISO 5167-2
//! square-edged-orifice relation
//!
//! ```text
//! dOmega / dP = ( sqrt(1 - beta^4 (1 - Cd^2)) - Cd beta^2 )
//!             / ( sqrt(1 - beta^4 (1 - Cd^2)) + Cd beta^2 ),
//! ```
//!
//! which lies in `(0, 1)`, falls from one toward zero as `beta` rises,
//! and collapses to the clean `(1 - beta^2) / (1 + beta^2)` when `Cd = 1`.
//!
//! Key behaviours the tests pin down: `Q` scales as `sqrt(dP)`; `Q` is
//! directly proportional to both throat area `A` and discharge
//! coefficient `Cd`; the diameter-ratio factor `1 / sqrt(1 - beta^4)`
//! always exceeds one; and a Venturi (`Cd ~ 0.98`) passes more flow than
//! a flow nozzle (`Cd ~ 0.97`), which passes more than an orifice plate
//! (`Cd ~ 0.61`) at the same geometry and pressure drop.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are textbook closed-form /
//! numerical models; this is NOT a clinical, medical, or production
//! engineering tool. `Cd` is treated as a caller-supplied constant and
//! the fluid is assumed incompressible (no gas-expansion `epsilon`
//! factor). The crate does NOT implement the Reynolds-number-dependent
//! Reader-Harris / Gallagher discharge-coefficient correlation, thermal
//! expansion, pressure-tapping geometry, or the installation rules and
//! uncertainty budgets of ISO 5167 / ASME MFC-3M, and it is no
//! substitute for an accredited flow calibration. The permanent-loss
//! relation is the square-edged orifice-plate formula; flow nozzles and
//! Venturis recover pressure in their diffusers and lose substantially
//! less, so it overestimates their loss. Use it to understand how a
//! differential-pressure meter scales, not to meter a real process.
//!
//! ## Example
//!
//! ```
//! use valenx_orifice::{Meter, MeterGeometry, MeterKind};
//!
//! // 50 mm orifice bore in a 100 mm pipe (beta = 0.5).
//! let geom = MeterGeometry::new(0.05, 0.10).unwrap();
//! let meter = Meter::with_typical_cd(geom, MeterKind::OrificePlate).unwrap();
//!
//! // Water, 50 kPa drop across the plate.
//! let q = meter.flow_rate(1000.0, 50_000.0).unwrap();
//! assert!(q > 0.0);
//!
//! // Inverting recovers the pressure drop exactly.
//! let dp = meter.pressure_drop(1000.0, q).unwrap();
//! assert!((dp - 50_000.0).abs() < 1e-6);
//!
//! // Only part of that 50 kPa is lost for good; the rest is recovered.
//! let lost = meter.permanent_pressure_loss(1000.0, q).unwrap();
//! assert!(lost < dp); // permanent loss is a fraction of the differential
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod meter;

pub use error::{ErrorCategory, OrificeError};
pub use meter::{Meter, MeterGeometry, MeterKind};
