//! Stochastic simulation — features 14 through 17.
//!
//! The discrete-stochastic half of the simulation layer. Molecule
//! counts are integers and the chemical master equation is sampled
//! exactly or approximately.
//!
//! - [`rng`] — a small deterministic `splitmix64` PRNG so every
//!   simulation is bit-for-bit reproducible from its seed.
//! - [`ssa`] — the [`StochasticModel`] plus three algorithms: the
//!   Gillespie direct method, explicit tau-leaping with a
//!   negative-population guard, and the Gibson-Bruck next-reaction
//!   method (features 14, 15, 16).
//! - [`ensemble`] — many independent runs aggregated into per-species,
//!   per-time mean / variance / percentile statistics (feature 17).

pub mod ensemble;
pub mod rng;
pub mod ssa;

pub use ensemble::{run_ensemble, EnsembleStats, SsaMethod};
pub use rng::Rng;
pub use ssa::{StochasticModel, StochasticTrace};
