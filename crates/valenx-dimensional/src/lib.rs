//! # valenx-dimensional
//!
//! Dimensionless groups for fluid mechanics and heat transfer.
//!
//! ## What
//!
//! A small, dependency-light library of the canonical dimensionless
//! numbers that recur throughout fluid mechanics, convective heat
//! transfer, and similitude analysis:
//!
//! - [`reynolds`] — Reynolds number `Re = rho v L / mu` (and a
//!   kinematic-viscosity form), plus a pipe-flow regime classifier.
//! - [`nusselt`] — Nusselt number `Nu = h L / k`.
//! - [`prandtl`] — Prandtl number `Pr = cp mu / k`.
//! - [`mach`] — Mach number `Ma = v / c` with a flow-speed regime
//!   classifier.
//! - [`froude`] — Froude number `Fr = v / sqrt(g L)` with an
//!   open-channel-flow regime classifier.
//! - [`biot`] — Biot number `Bi = h L / k` and the lumped-capacitance
//!   validity test.
//! - [`peclet`] — Peclet number, defined both directly
//!   (`Pe = rho cp v L / k`) and as the product `Pe = Re * Pr`.
//!
//! Each quantity is returned as a thin `f64` newtype (for example
//! [`reynolds::Reynolds`]) so the type system records which group a
//! value represents; every wrapper derives `serde` `Serialize` /
//! `Deserialize` and exposes its raw value via `value()`.
//!
//! ## Model
//!
//! These are the textbook closed-form definitions, nothing more. All
//! formulas assume a single **consistent system of units** — the doc
//! comments use SI (metre, second, kilogram, kelvin, pascal-second,
//! watt) but any coherent unit system works because every group is
//! dimensionless by construction. The library does **not** carry or
//! check physical units; it is the caller's responsibility to supply
//! consistent inputs.
//!
//! The regime thresholds (pipe-flow transition near `Re = 2300`, the
//! sonic boundary at `Ma = 1`, the critical open-channel boundary at
//! `Fr = 1`, the lumped-capacitance rule of thumb `Bi < 0.1`) are the
//! standard engineering rules of thumb. They are approximate, geometry-
//! and correlation-dependent in practice, and are provided only as
//! coarse classifiers — see each module's documentation.
//!
//! ## Honest scope
//!
//! Research / educational grade. This crate implements textbook
//! closed-form / numerical models with explicit, consistent-unit
//! conventions for learning, prototyping, and sanity-checking. It is
//! **not** a clinical, medical, or production engineering tool. Real
//! convective-transfer or compressible-flow design needs validated CFD
//! / thermal solvers, geometry-specific correlations, and certified
//! property data; the rule-of-thumb regime boundaries here are not a
//! substitute for engineering judgement or a qualified analysis.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod biot;
pub mod error;
pub mod froude;
pub mod mach;
pub mod nusselt;
pub mod peclet;
pub mod prandtl;
pub mod reynolds;

pub use biot::{Biot, LumpedCapacitance};
pub use error::DimensionlessError;
pub use froude::{ChannelRegime, Froude};
pub use mach::{Mach, SpeedRegime};
pub use nusselt::Nusselt;
pub use peclet::Peclet;
pub use prandtl::Prandtl;
pub use reynolds::{PipeRegime, Reynolds};
