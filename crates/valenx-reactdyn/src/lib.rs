//! # valenx-reactdyn
//!
//! Native **reaction-dynamics** core for Valenx — an engine-agnostic
//! molecular-dynamics shell plus an **ab-initio MD (AIMD)** backend.
//!
//! AIMD advances the nuclei with a velocity-Verlet integrator while the
//! forces come from **central finite differences of `valenx-qchem`'s
//! single-point energy** — qchem ships no analytic nuclear gradient
//! (`GeometryOptRequest::run` is a documented stub), so the gradient is
//! computed numerically. Everything runs in **atomic units** — qchem's
//! native bohr / hartree — so forces are hartree/bohr with no conversion;
//! see [`units`] for the two boundary conversions (timestep, mass).
//!
//! Phase 1 ships the shell + the AIMD backend; QM/MM and ReaxFF backends
//! plug in behind the same `ReactionEngine` trait later.

#![forbid(unsafe_code)]

pub mod engine;
pub mod error;
pub mod forces;
pub mod integrator;
pub mod kinetics;
pub mod mm;
pub mod qmmm;
pub mod reactive;
pub mod units;

pub use engine::{AimdEngine, Controls, Frame, ReactionEngine, System, Thermostat, Trajectory};
pub use error::{ReactDynError, Result};
pub use forces::Method;
pub use kinetics::{arrhenius_half_life_1st_order, arrhenius_rate, equilibrium_constant};
pub use mm::{classical_forces, Particle};
pub use qmmm::{Embedding, MmAtom, QmMmEngine, QmMmSystem};
pub use reactive::{morse_param, reactive_energy_forces, MorseParam, ReactiveEngine, ReactiveSystem};
