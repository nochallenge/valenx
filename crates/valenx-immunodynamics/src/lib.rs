//! # valenx-immunodynamics — within-host viral dynamics
//!
//! A small, self-contained simulator for the **target-cell-limited**
//! model of an acute within-host viral infection: how a population of
//! free virions rises, peaks, and falls as it exhausts the supply of
//! cells it can infect.
//!
//! ## What
//!
//! Describe the infection with four non-negative rate constants
//! ([`Parameters`]) and an initial [`State`] of three populations, then
//! [`simulate`] integrates the flight of the infection forward in time
//! and hands back a [`Trajectory`] of time / state samples. From the
//! trajectory you can read the engineering answer — when the viral load
//! peaks ([`Trajectory::peak_viral_load`]), how high, and what the
//! infection settles to — and from the parameters the basic
//! reproductive number [`Parameters::r0`].
//!
//! ```
//! use valenx_immunodynamics::{Parameters, State, simulate};
//!
//! // An acute-infection parameter set with R0 > 1.
//! let params = Parameters::new(2e-5, 1.0, 50.0, 5.0).expect("valid rates");
//! let start = State::new(1e5, 1.0, 0.0).expect("valid state");
//!
//! // R0 = beta * T0 * p / (delta * c).
//! let r0 = params.r0(start.target).expect("R0 defined");
//! assert!(r0 > 1.0);
//!
//! // Integrate to t = 30 with a fixed RK4 step, sampling every 10th step.
//! let traj = simulate(&params, &start, 30.0, 1e-3, 10).expect("valid run");
//! let (_, peak_time, peak_v) = traj.peak_viral_load().unwrap();
//! assert!(peak_v > 1.0);          // a real outbreak occurred
//! assert!(peak_time > 0.0);       // it peaked after the start
//! assert!(traj.all_non_negative(0.0)); // populations stay physical
//! ```
//!
//! ## Model
//!
//! Three compartments — uninfected **target cells** `T`, productively
//! **infected cells** `I`, and free **virions** `V` — evolve under the
//! standard three-equation system (Nowak & May 2000; Perelson 2002):
//!
//! ```text
//! dT/dt = -beta * T * V
//! dI/dt =  beta * T * V - delta * I
//! dV/dt =  p * I        - c * V
//! ```
//!
//! with infection rate `beta`, infected-cell death rate `delta`, virion
//! production rate `p`, and virion clearance rate `c`. Linearising about
//! the infection-free state `(T0, 0, 0)` gives the basic reproductive
//! number
//!
//! ```text
//! R0 = beta * T0 * p / (delta * c).
//! ```
//!
//! There is no target-cell replenishment, so the infection is *limited
//! by its target cells*: the finite pool `T0` depletes, the effective
//! reproductive number falls below one, and the viral load peaks and
//! then declines. Integration is a fixed-step explicit **RK4**
//! ([`integrate`]); every accepted step is clamped onto the non-negative
//! orthant so the discrete trajectory stays inside the model's physical
//! domain (populations are amounts, never negative).
//!
//! ## Honest scope
//!
//! Research / educational grade. Every piece here is the genuine
//! textbook article — the three-equation target-cell-limited model, the
//! closed-form `R0`, and a classic fourth-order Runge-Kutta integrator
//! that reproduces exact exponential decay in the degenerate limits and
//! the characteristic peak-then-decline of an acute infection. It is a
//! deliberate **v1** of a well-established, well-understood model, **not
//! a clinical, diagnostic, or production immunology tool**:
//!
//! - The model is the *basic* target-cell-limited form: no target-cell
//!   birth / death source term, no latently infected compartment, no
//!   explicit immune-effector (CTL / antibody) dynamics, no drug /
//!   pharmacokinetic terms, and no eclipse phase.
//! - It is deterministic and well-mixed (ODEs, not a stochastic or
//!   spatial agent model), so it describes mean behaviour, not the
//!   extinction / take-off randomness of a small founding population.
//! - The integrator is fixed-step RK4 — accurate and energy-faithful on
//!   these non-stiff kinetics, but with no adaptive error control or
//!   stiff solver.
//!
//! None of those omissions makes the result meaningless: the peak viral
//! load, its timing, the `R0` threshold, and the depletion-driven
//! turnover are all the real, named quantities of within-host viral
//! dynamics. Each omission is a documented, standard extension on the
//! way toward a fuller immunodynamics suite.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod integrate;
pub mod model;

pub use error::{ErrorCategory, ImmunoError};
pub use integrate::{rk4_step, simulate, Trajectory};
pub use model::{Parameters, State};
