//! Phylogenetic simulation.
//!
//! Simulation generates trees and sequences under explicit models — for
//! testing inference methods, power analysis, and teaching.
//!
//! - [`coalescent`] — Kingman's coalescent: a genealogy of `n` sampled
//!   lineages running backward in time, under a constant or
//!   piecewise-constant population size.
//! - [`birthdeath`] — a forward birth-death branching process: lineages
//!   speciate at rate `λ` and go extinct at rate `μ`.
//! - [`seqgen`] — Seq-Gen-class sequence evolution: simulate an
//!   alignment by evolving a root sequence down a tree under a
//!   [substitution model](crate::likelihood::model).
//! - [`clock`] — root-to-tip molecular-clock dating: linear regression
//!   of each tip's root-to-tip distance against its sampling time
//!   yields a substitution rate and a root-age estimate.
//!
//! All randomness flows through the deterministic [`crate::rng::Rng`],
//! so a given seed reproduces a run exactly.

pub mod birthdeath;
pub mod clock;
pub mod coalescent;
pub mod seqgen;

pub use birthdeath::{simulate_birth_death, BirthDeathParams};
pub use clock::{root_to_tip_regression, ClockEstimate};
pub use coalescent::{simulate_coalescent, PopulationSize};
pub use seqgen::{simulate_sequences, SimulatedAlignment};
