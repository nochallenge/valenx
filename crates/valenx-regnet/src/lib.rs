//! # valenx-regnet — gene-regulatory-network dynamics
//!
//! A small, self-contained ODE model of transcriptional regulation, with
//! a fixed-step RK4 integrator and two canonical synthetic-biology
//! presets. State is a plain `Vec<f64>`; there is no linear-algebra
//! dependency.
//!
//! ## What
//!
//! - [`hill_activate`] / [`hill_repress`] — the sigmoidal Hill kinetics
//!   that turn a regulator concentration into a fractional transcription
//!   rate (the checked variants validate `k > 0`, `n > 0`).
//! - [`GeneRegulatoryNetwork`] — `N` genes, each with a production rate, a
//!   degradation rate, a basal leak, and a list of [`Regulator`] inputs;
//!   plus [`rk4_step`](GeneRegulatoryNetwork::rk4_step) and
//!   [`simulate`](GeneRegulatoryNetwork::simulate).
//! - [`Trajectory`] — the time series a run returns, with helpers to pull
//!   a single gene's [`series`](Trajectory::series) and to count its
//!   [`local maxima`](Trajectory::count_local_maxima) (the oscillation
//!   diagnostic).
//! - [`toggle_switch`] and [`repressilator`] presets (plus their
//!   `_default` parameterisations).
//!
//! ## Model
//!
//! Each gene `i` evolves by
//!
//! ```text
//! dx_i/dt = production_i * regulation_i(x) - degradation_i * x_i
//! ```
//!
//! where the dimensionless `regulation_i(x) ∈ [0, 1]` is the product of
//! the gene's Hill factors (activators via `x^n/(k^n+x^n)`, repressors via
//! `k^n/(k^n+x^n)`), lifted by a basal leak floor and capped at `1`. The
//! system is advanced with the **classical fourth-order Runge-Kutta**
//! method at a fixed step `dt`.
//!
//! The two presets reproduce the textbook qualitative behaviours: the
//! mutual-repression **toggle switch** is *bistable* (it settles into one
//! of two stable states depending on its initial condition), and the
//! cyclic-repression **repressilator** is a *sustained oscillator* (its
//! negative-feedback loop has no stable fixed point and settles into a
//! limit cycle).
//!
//! ## Honest scope
//!
//! **Research / educational grade.** These are textbook equations with
//! dimensionless example rates; the crate is **not a clinical, medical, or
//! production tool**. The presets illustrate the *qualitative* dynamics
//! (bistability, oscillation), not quantitative predictions for any real
//! gene, protein, or organism. There is no parameter fitting, no
//! stochasticity (it is a deterministic mean-field ODE — no intrinsic
//! noise, no single-molecule effects), no delays, and no spatial
//! structure. The crate writes no files: callers persist the returned
//! [`Trajectory`] themselves.
//!
//! ## Example
//!
//! ```rust
//! use valenx_regnet::{repressilator_default, GeneRegulatoryNetwork};
//!
//! // Three genes in a repression cycle -> a sustained oscillation.
//! let net: GeneRegulatoryNetwork = repressilator_default();
//!
//! // Integrate from an asymmetric start with RK4, dt = 0.01, 8000 steps.
//! let traj = net.simulate(&[1.0, 0.0, 0.0], 0.01, 8000).unwrap();
//!
//! // Gene 0's time series shows several peaks -> it oscillates.
//! let peaks = traj.count_local_maxima(0).unwrap();
//! assert!(peaks >= 3, "expected a sustained oscillation, saw {peaks} peaks");
//! ```

#![forbid(unsafe_code)]
// `missing_docs` is enforced via `[lints.rust]` in Cargo.toml (matching
// the sibling crates' convention); the crate-level attr is kept too so
// the intent is visible at the top of the file.
#![warn(missing_docs)]

pub mod error;
pub mod hill;
pub mod network;
pub mod presets;

pub use error::{RegnetError, Result};
pub use hill::{hill_activate, hill_activate_checked, hill_repress, hill_repress_checked};
pub use network::{Gene, GeneRegulatoryNetwork, Regulator, RegulatorKind, Trajectory};
pub use presets::{repressilator, repressilator_default, toggle_switch, toggle_switch_default};
