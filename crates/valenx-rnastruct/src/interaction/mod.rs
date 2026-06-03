//! RNA-RNA interaction.
//!
//! Algorithms for how *two* RNA molecules behave together:
//!
//! - [`cofold`] — RNA-RNA cofolding: the joint minimum-free-energy
//!   structure of a two-strand complex (`RNAcofold`-class).
//! - [`accessibility`] — single-strand accessibility profiles: the
//!   per-base unpaired probability and window opening free energies.
//! - [`interaction`] — v1 seed-window IntaRNA-class interaction-site
//!   prediction combining a hybridisation duplex with the
//!   accessibility cost.
//! - [`intarna`] — full IntaRNA-class accessibility-aware interaction:
//!   seed + extension DP with internal loops / bulges, exposing
//!   intermolecular-pair lists and per-strand window-opening costs.

pub mod accessibility;
pub mod cofold;
pub mod intarna;
#[allow(clippy::module_inception)]
pub mod interaction;
