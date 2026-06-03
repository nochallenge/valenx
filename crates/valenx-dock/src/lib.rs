//! # valenx-dock
//!
//! Native Rust molecular-docking engine. Implements AutoDock Vina's
//! scoring function and search strategy in pure Rust so the Valenx
//! mesh editor can dock ligands without spawning the external `vina`
//! binary.
//!
//! The public entry point is [`runner::dock`]. Receptor + ligand are
//! parsed from PDBQT via [`valenx_bio::format::pdbqt`]; the search
//! returns a ranked list of [`pose::Pose`]s scored in Vina kcal/mol.
//!
//! Tunable defaults match Vina 1.2 ([`DEFAULT_EXHAUSTIVENESS`],
//! [`DEFAULT_NUM_MODES`], [`DEFAULT_ENERGY_RANGE`]).
//!
//! Reference: Trott & Olson, J. Comput. Chem. 31 (2010) 455-461.
//!
//! # Example
//!
//! ```no_run
//! use std::path::Path;
//! use nalgebra::Vector3;
//! use valenx_dock::{dock, DockConfig};
//!
//! let receptor = std::fs::read_to_string("receptor.pdbqt").unwrap();
//! let ligand   = std::fs::read_to_string("ligand.pdbqt").unwrap();
//! let mut cfg  = DockConfig::default();
//! cfg.center = Vector3::new(15.0, 54.0, 17.0);
//! cfg.size   = Vector3::new(22.5, 22.5, 22.5);
//! let poses = dock(&receptor, &ligand, &cfg, Path::new("out.pdbqt"), None).unwrap();
//! println!("best pose score = {:.3} kcal/mol", poses[0].1);
//! ```
//!
//! # Known Divergences from Upstream Vina
//!
//! * **No intra-ligand energy term.** Vina's full scoring function
//!   sums an inter-molecular term (implemented here) and an
//!   intra-ligand term over pairs separated by 4+ bonds (NOT
//!   implemented here). For sufficiently flexible ligands BFGS may
//!   drive torsions into internal clashes that real Vina would
//!   penalize. Top-1 RMSD parity on rigid-or-near-rigid ligands is
//!   unaffected; flexible-ligand parity is the open follow-on.
//! * **Fixed receptor only.** Flexible side chains (Vina supports
//!   them via the `flex` PDBQT input) are not modeled. The receptor
//!   PDBQT is loaded as rigid atoms.
//! * **No macrocycle sampling.** Vina 1.2 ships a macrocycle handler;
//!   valenx-dock treats all bonds as either rigid or fully torsional.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod atom_type;
pub mod cluster;
pub mod config;
pub mod dry_run;
pub mod error;
pub mod eval;
pub mod events;
pub mod grid;
pub mod ligand;
pub mod pose;
pub mod receptor;
pub mod runner;
pub mod schema;
pub mod score;
pub mod search;

pub use config::DockConfig;
pub use dry_run::{dock_dry_run, DryRunPlan};
pub use error::DockError;
pub use pose::Pose;
pub use runner::dock;

/// Vina 1.2 default Monte-Carlo restart count.
pub const DEFAULT_EXHAUSTIVENESS: u32 = 8;
/// Vina 1.2 default cap on poses written.
pub const DEFAULT_NUM_MODES: u32 = 9;
/// Vina 1.2 default kcal/mol cutoff above the best pose.
pub const DEFAULT_ENERGY_RANGE: f64 = 3.0;
/// Pair-interaction cutoff in Å (Vina internal: 8.0 Å).
pub const PAIR_CUTOFF: f64 = 8.0;
