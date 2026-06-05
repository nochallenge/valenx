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
pub mod aniso_field;
pub mod bioheat;
pub mod cable;
pub mod cable_theory;
pub mod chord_conductance;
pub mod conduction;
pub mod current_distance;
pub mod engine;
pub mod error;
pub mod field;
pub mod ghk;
pub mod impedance;
pub mod ionic;
pub mod membrane;
pub mod myelinated;
pub mod nernst;
pub mod recording;
pub mod safety;
pub mod steering;
pub mod strength_duration;
pub mod temperature;
pub mod units;

pub use activating::activating_along_x;
pub use aniso_field::{analytic_aniso_point_source_mv, AnisoTissue, Conductivity, SolvedField};
pub use bioheat::{analytic_point_heat_k, solve_point_heat, BioheatField};
pub use cable::{count_spikes, HhCable, HhCompartment, StimPulse};
pub use cable_theory::{
    charging_fraction, electrotonic_length, length_constant_cm, open_end_input_resistance,
    rall_equivalent_diameter, sealed_end_input_resistance, semi_infinite_input_resistance,
    steady_state_attenuation, time_constant_s, time_to_charge_fraction,
};
pub use chord_conductance::{chord_conductance_potential_mv, ConductanceChannel};
pub use conduction::{
    myelinated_conduction_velocity, unmyelinated_conduction_velocity, HURSH_FACTOR_M_PER_S_PER_UM,
};
pub use current_distance::{activation_radius, fit_constant, threshold_current};
pub use engine::{recruitment_curve, stimulate, Axon, Recruitment, Scene};
pub use error::{NeuroError, Result};
pub use field::{analytic_point_source_mv, ExtracellularField, TissueGrid};
pub use ghk::{ghk_potential_mv, GhkIon};
pub use impedance::{Cpe, ElectrodeImpedance};
pub use ionic::{driving_force_mv, ionic_current};
pub use membrane::{HhMembrane, ImplicitCable, Membrane};
pub use myelinated::MyelinatedFiber;
pub use nernst::{nernst_potential_mv, thermal_voltage_mv, BODY_TEMPERATURE_K};
pub use recording::{ExtracellularRecorder, Recording};
pub use safety::{charge_density, is_safe, max_safe_charge_per_phase, shannon_k};
pub use steering::ContactArray;
pub use strength_duration::{chronaxie, rheobase, threshold_amplitude};
pub use temperature::{q10_scale, TYPICAL_GATING_Q10};
