//! Forward-in-time population simulation.
//!
//! This module group is the SLiM / fwdpy11 / simuPOP side of
//! `valenx-popgen`: an explicit generation-by-generation diploid
//! Wright-Fisher process.
//!
//! - [`wright_fisher`] — the [`WrightFisher`] simulator and its
//!   [`SimulationConfig`].
//! - [`selection`] — additive, multiplicative and epistatic
//!   [`SelectionModel`]s.
//! - [`mutation`] — infinite-sites and finite-sites [`MutationModel`]s.
//! - [`recombination`] — crossover and gene-conversion
//!   [`RecombinationModel`].
//! - [`migration`] — island and stepping-stone [`MigrationModel`]s for
//!   structured populations.
//! - [`demography`] — [`DemographicSchedule`] size histories and the
//!   closed-form [`Drift`] quantities.
//! - [`tree_recording`] — Wright-Fisher with `pyslim`/`tskit`-style
//!   genealogy recording.

pub mod demography;
pub mod migration;
pub mod mutation;
pub mod recombination;
pub mod selection;
pub mod tree_recording;
pub mod wright_fisher;

pub use demography::{DemographicSchedule, Drift};
pub use migration::MigrationModel;
pub use mutation::{MutationEvent, MutationModel};
pub use recombination::RecombinationModel;
pub use selection::{SelectionCoefficients, SelectionModel};
pub use tree_recording::{record_wright_fisher, RecordingConfig};
pub use wright_fisher::{SimulationConfig, WrightFisher};
