//! Pose search algorithms.
//!
//! A docking search explores the space of ligand poses — translation,
//! orientation and per-bond torsion angles — looking for the global
//! energy minimum. This crate ships three classical search strategies
//! plus two drivers that wire them to a receptor:
//!
//! - [`ga`] — an AutoDock4-class **Lamarckian genetic algorithm**: a
//!   population of poses evolved by selection, crossover and mutation,
//!   with a local-search step whose result is written *back* into the
//!   genome (the "Lamarckian" part).
//! - [`mc`] — a **Monte-Carlo / simulated-annealing** pose search: a
//!   single walker that perturbs, scores and accepts by the Metropolis
//!   criterion under a cooling temperature schedule.
//! - [`ils`] — a Vina-class **iterated local search**: mutate, run a
//!   BFGS local optimisation, accept by Metropolis. Reuses
//!   [`valenx_dock`]'s BFGS minimiser.
//! - [`driver`] — [`driver::rigid_dock`] and
//!   [`driver::flexible_dock`], the two entry points
//!   that build affinity maps and run a chosen search.
//!
//! All three search algorithms score poses through a [`PoseObjective`]
//! — a closure-free wrapper over a [`valenx_dock::ligand::Ligand`] and
//! a precomputed [`crate::score::AffinityMapSet`] — so the search code
//! never re-evaluates receptor pairs in its inner loop.

pub mod driver;
pub mod flex_pose;
pub mod ga;
pub mod ils;
pub mod mc;
pub mod objective;
pub mod solis_wets;

pub use driver::{
    flexible_dock, induced_fit_dock_driver, rigid_dock, FlexibleSpec, SearchAlgorithm,
};
pub use flex_pose::{
    chi_from_axis_atoms, induced_fit_dock, induced_fit_solis_wets, ChiRotation, FlexPose,
    FlexPoseObjective, InducedFitResult,
};
pub use ga::{GaParams, LamarckianGa};
pub use ils::IlsParams;
pub use mc::{McParams, McSchedule};
pub use objective::PoseObjective;
pub use solis_wets::{solis_wets, SolisWetsParams};
