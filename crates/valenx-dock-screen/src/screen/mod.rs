//! Batch virtual screening and pose post-processing.
//!
//! Once a single ligand can be docked, virtual screening is the next
//! step: dock a whole library and rank it. This module covers the
//! screening pipeline and the pose analysis a screen needs:
//!
//! - [`batch`] — dock a library of ligands against a receptor and
//!   rank them by best score (parallelised across ligands).
//! - [`cluster`] — RMSD-based pose clustering and ranking, so a
//!   ligand's many poses collapse to a handful of distinct binding
//!   modes.
//! - [`ensemble`] — ensemble docking: dock against several receptor
//!   conformations and take the best (or a Boltzmann-weighted)
//!   per-ligand score.
//! - [`consensus`] — consensus scoring: combine several scoring
//!   functions by rank aggregation so no single function's bias
//!   dominates.
//! - [`rescore`] — an MM-GBSA-class generalized-Born rescoring of a
//!   docked pose, a more physical (and more expensive) energy than the
//!   docking score.

pub mod batch;
pub mod cluster;
pub mod consensus;
pub mod ensemble;
pub mod rescore;

pub use batch::{screen_library, LibraryEntry, ScreenResult};
pub use cluster::{cluster_poses, PoseCluster};
pub use consensus::{consensus_rank, ConsensusMethod};
pub use ensemble::{ensemble_dock, EnsembleMethod};
pub use rescore::{mmgbsa_rescore, MmGbsaTerms};
