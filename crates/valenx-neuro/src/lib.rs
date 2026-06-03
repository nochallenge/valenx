//! # valenx-neuro
//!
//! Native **neural-interface / BCI** simulation for Valenx: the physics of
//! an implanted stimulating electrode acting on nearby neurons.
//!
//! Five coupled modules, each validated against a closed-form or textbook
//! result (see RFC 0011):
//!
//! - **Extracellular field** — `−∇·(σ∇φ)=I` in tissue, reusing
//!   `valenx-fem`'s steady-conduction solver (the operator is identical to
//!   `−∇·(k∇T)=q`).
//! - **Hodgkin–Huxley cable** — the membrane dynamics of axons.
//! - **Activating function** — the Rattay coupling from field to membrane.
//! - **Bioheat** — Pennes tissue heating from stimulation.
//! - **Electrode impedance** — access resistance + double-layer CPE.
//!
//! Scope is **research / education-grade** — quasi-static fields, idealized
//! geometry, and standard membrane models. It is not clinical or certified
//! software.

#![forbid(unsafe_code)]

pub mod activating;
pub mod bioheat;
pub mod cable;
pub mod engine;
pub mod error;
pub mod field;
pub mod impedance;
pub mod membrane;
pub mod myelinated;
pub mod units;

pub use activating::activating_along_x;
pub use bioheat::{analytic_point_heat_k, solve_point_heat, BioheatField};
pub use cable::{count_spikes, HhCable, HhCompartment, StimPulse};
pub use engine::{recruitment_curve, stimulate, Axon, Recruitment, Scene};
pub use error::{NeuroError, Result};
pub use field::{analytic_point_source_mv, ExtracellularField, TissueGrid};
pub use impedance::{Cpe, ElectrodeImpedance};
pub use membrane::{HhMembrane, ImplicitCable, Membrane};
pub use myelinated::MyelinatedFiber;
