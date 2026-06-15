//! # valenx-thermoreg
//!
//! Human thermoregulation as a single-node heat balance:
//! `M - W = R + C + E + S`.
//!
//! ## What
//!
//! A small, self-contained model of whole-body human heat exchange. It
//! treats the body as one lumped thermal mass (a [`body::Body`]) sitting
//! in a thermal [`environment::Environment`], and closes the conceptual
//! physiology heat balance
//!
//! ```text
//!   M - W = R + C + E + S
//! ```
//!
//! — metabolic heat production `M`, external work `W`, and the
//! radiative / convective / evaporative losses `R` / `C` / `E`, with the
//! storage `S` as the residual that drives the core temperature. From
//! the stored heat it integrates the core temperature forward in time by
//! the calorimeter relation `dT = Q / (m c)`.
//!
//! Typical use:
//!
//! ```
//! use valenx_thermoreg::{Body, Environment, Metabolism, Sweat, heat_balance, step_core_temp};
//!
//! // A resting adult in a warm, still room, sweating lightly.
//! let body = Body::standard_adult(70.0, 170.0, 37.0, 33.0).unwrap();
//! let env = Environment::still_indoor(28.0).unwrap();
//! let met = Metabolism::resting();
//! let sweat = Sweat::from_rate(2.0e-5).unwrap(); // kg/s evaporated
//!
//! let balance = heat_balance(&body, &env, &met, &sweat);
//! // The books always close: M - W - (R + C + E) - S == 0.
//! assert!(balance.closure_residual_w().abs() < 1e-9);
//!
//! // Advance the core temperature by one minute.
//! let later = step_core_temp(&body, &env, &met, &sweat, 60.0).unwrap();
//! let _drift = later.core_temp_c - body.core_temp_c;
//! ```
//!
//! ## Model
//!
//! Every term is a textbook closed-form relation:
//!
//! - **Convection** — Newton's law of cooling,
//!   `C = h * A * (T_skin - T_air)` ([`Environment::convective_power`]).
//! - **Radiation** — Stefan-Boltzmann net long-wave exchange,
//!   `R = ε * σ * A * (T_skin^4 - T_radiant^4)` in kelvin
//!   ([`Environment::radiative_power`]).
//! - **Evaporation** — latent cooling,
//!   `E = m_dot_sweat * L` ([`Sweat::evaporative_power`]).
//! - **Storage** — the closing residual,
//!   `S = (M - W) - (R + C + E)` ([`heat_balance`]).
//! - **Core dynamics** — lumped capacitance,
//!   `dT = S * dt / (m * c)` ([`core_temp_change`] / [`step_core_temp`]).
//!
//! Body surface area can be estimated with the DuBois formula
//! ([`Body::dubois_surface_area`]); the default body and latent-heat
//! constants are the standard physiology values
//! ([`body::BODY_SPECIFIC_HEAT`], [`body::LATENT_HEAT_SWEAT`]).
//!
//! ## Honest scope
//!
//! Research / educational grade. These are the well-established,
//! textbook closed-form and lumped-parameter heat-transfer relations you
//! would find in a physiology or heat-transfer course — each term is the
//! genuine article and the balance closes exactly — but the model is
//! deliberately a one-node reduction:
//!
//! - **Single thermal node.** The body is one well-mixed mass at one
//!   temperature; there is no core-to-skin gradient, no blood-flow
//!   redistribution, and no per-segment resolution (it is not a
//!   Stolwijk / Fiala multi-node model).
//! - **No active control loop.** Sweat rate, skin temperature and
//!   metabolism are *inputs*; the crate does not model the
//!   hypothalamic set-point, shivering, vasomotor response, or sweat
//!   onset — you supply the control action, it computes the energetics.
//! - **Constant material properties** and skin temperature across a step;
//!   no humidity-limited evaporative ceiling, no clothing (`clo`)
//!   resistance, no respiratory heat loss.
//! - **Explicit-Euler integration** for the core temperature — fine for
//!   small steps, not a stiff-ODE solver.
//!
//! It is **NOT a clinical, medical, or production engineering tool**. Do
//! not use it for patient care, occupational heat-stress certification,
//! thermal-comfort compliance, or any safety-critical decision.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod balance;
pub mod body;
pub mod environment;
pub mod error;

pub use balance::{core_temp_change, heat_balance, step_core_temp, HeatBalance, Metabolism};
pub use body::{Body, Sweat, BODY_SPECIFIC_HEAT, LATENT_HEAT_SWEAT};
pub use environment::{Environment, KELVIN_OFFSET, SKIN_EMISSIVITY, STEFAN_BOLTZMANN};
pub use error::{ErrorCategory, ThermoregError};
