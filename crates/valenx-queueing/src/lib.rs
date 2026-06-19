//! # valenx-queueing — the M/M/1 queue
//!
//! A pure, dependency-light calculator for the classic **M/M/1**
//! queueing model: closed-form steady-state metrics for a single-server
//! queue with Poisson arrivals and exponential service.
//!
//! ## What
//!
//! Given an arrival rate `lambda` and a service rate `mu` (in the same
//! per-unit-time units), [`Mm1::new`] validates the pair and the
//! resulting [`Mm1`] exposes every standard steady-state quantity:
//!
//! - [`rho`](Mm1::rho) — traffic intensity / utilization
//!   `rho = lambda / mu`;
//! - [`p0`](Mm1::p0) — idle probability `P(0) = 1 - rho`;
//! - [`l`](Mm1::l) — mean number in system `L`;
//! - [`lq`](Mm1::lq) — mean number waiting `Lq`;
//! - [`w`](Mm1::w) — mean time in system `W`;
//! - [`wq`](Mm1::wq) — mean wait in queue `Wq`;
//! - [`state_probability`](Mm1::state_probability) — the stationary
//!   law `P(n) = (1 - rho) rho^n`, plus its [`cdf`](Mm1::cdf);
//! - [`metrics`](Mm1::metrics) — all of the above bundled into a
//!   serializable [`Metrics`] record;
//! - [`little_residual`](Mm1::little_residual) — the `L - lambda W`
//!   Little's-law self-consistency residual.
//!
//! Plus the two capacity-sizing inverses of `W = 1 / (mu - lambda)`:
//! [`service_rate_for_mean_response_time`](Mm1::service_rate_for_mean_response_time)
//! (`mu = lambda + 1 / W`), the service rate a single server needs to hit
//! a target mean response time at a given arrival rate, and its
//! complement
//! [`arrival_rate_for_mean_response_time`](Mm1::arrival_rate_for_mean_response_time)
//! (`lambda = mu - 1 / W`), the arrival rate a server of a given rate can
//! accept while holding that response time.
//!
//! ```
//! use valenx_queueing::Mm1;
//!
//! // lambda = 2, mu = 3  =>  rho = 2/3, L = 2, W = 1.
//! let q = Mm1::new(2.0, 3.0).unwrap();
//! assert!((q.rho() - 2.0 / 3.0).abs() < 1e-12);
//! assert!((q.l() - 2.0).abs() < 1e-12);
//! assert!((q.w() - 1.0).abs() < 1e-12);
//! // Little's law L = lambda W holds exactly.
//! assert!(q.little_residual().abs() < 1e-12);
//! ```
//!
//! ## Model
//!
//! M/M/1 (Kendall notation `A/S/c` with **M**arkovian arrivals,
//! **M**arkovian service, **1** server) is the continuous-time
//! birth-death chain on the customer count `n = 0, 1, 2, ...` with a
//! constant birth (arrival) rate `lambda` and a constant death
//! (service-completion) rate `mu`. Solving the global-balance equations
//! gives the geometric stationary distribution `P(n) = (1 - rho) rho^n`,
//! which exists iff the **stability condition** `rho = lambda / mu < 1`
//! holds. The mean-value formulas then follow:
//!
//! - `L  = rho / (1 - rho)`
//! - `Lq = rho^2 / (1 - rho) = L - rho`
//! - `W  = 1 / (mu - lambda)`
//! - `Wq = rho / (mu - lambda) = W - 1 / mu`
//!
//! and **Little's law** `L = lambda W` (likewise `Lq = lambda Wq`) ties
//! the count and time metrics together. As `rho -> 1` from below, all
//! four means diverge to `+inf`. The crate validates inputs up front, so
//! a constructed [`Mm1`] only ever yields finite results.
//!
//! ## Honest scope
//!
//! Research/educational grade. This is the **textbook closed-form
//! steady-state** M/M/1 model — exact rational formulas in `lambda` and
//! `mu` for the long-run *average* behaviour of an idealized queue
//! (single server, Poisson arrivals, exponential service, infinite
//! waiting room, work-conserving FCFS). It is **not** a discrete-event
//! simulator and says nothing about transient, finite-horizon, or
//! heavy-traffic-corrected behaviour, non-exponential (G/G/1) workloads,
//! multiple servers, finite buffers, priorities, or balking/reneging.
//! It is **not** a clinical/medical or production engineering /
//! capacity-planning tool — do not size real systems from it without a
//! model that matches their actual arrival and service distributions.
//!
//! ## Errors
//!
//! Every fallible function returns
//! [`Result<_, QueueingError>`](error::QueueingError); the error type
//! carries stable [`code`](error::QueueingError::code) and
//! [`category`](error::QueueingError::category) accessors for telemetry.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod mm1;

pub use error::{ErrorCategory, QueueingError, Result};
pub use mm1::{Metrics, Mm1};

// (intentionally no other re-exports — the crate's whole surface is the
// `Mm1` queue, its `Metrics` snapshot, and the error taxonomy.)
