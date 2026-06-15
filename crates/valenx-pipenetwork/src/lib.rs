//! # valenx-pipenetwork
//!
//! Steady-state hydraulic analysis of *looped* pipe networks by the
//! **Hardy-Cross** method.
//!
//! ## What
//!
//! Given a set of pipes — each with a resistance coefficient `k` and an
//! assumed flow `q` — grouped into one or more closed [`Loop`]s, this crate
//! iteratively corrects the flows until every loop is hydraulically
//! balanced. The result is the steady-state flow distribution that
//! simultaneously satisfies:
//!
//! - **energy / loop balance** — the signed head loss summed around each
//!   loop is zero, and
//! - **continuity / node balance** — flow into each junction equals flow
//!   out (preserved exactly from the initial guess; see below).
//!
//! The public surface is small:
//!
//! - [`Pipe`] — an edge with loss law `h = k q |q|`, plus the resistance
//!   helpers [`darcy_weisbach_k`] and [`hazen_williams_k`].
//! - [`Loop`] / [`LoopMember`] / [`Orientation`] — an oriented closed
//!   circuit and its Hardy-Cross [`Loop::correction`].
//! - [`Network`] — pipes + loops (+ optional [`Endpoints`] topology) and
//!   the driver [`Network::solve`], configured by [`SolveConfig`] and
//!   reporting a [`SolveReport`].
//! - [`NetworkError`] — a validated error taxonomy.
//!
//! ## Model
//!
//! Each pipe obeys the quadratic loss law
//!
//! ```text
//! h = k * q * |q|        (signed; magnitude k q^2)
//! ```
//!
//! For a loop, with member sign `s_i = +1` when the pipe is traversed in its
//! own flow direction and `s_i = -1` otherwise, the **Hardy-Cross loop
//! correction** is
//!
//! ```text
//! dQ = - sum_i ( s_i k_i q_i |q_i| )  /  sum_i ( 2 k_i |q_i| )
//! ```
//!
//! and is applied to every member as `q_i <- q_i + s_i dQ`. The numerator is
//! the net head loss around the loop; the denominator is the sum of the loss
//! derivatives `dh/dq = 2 k |q|`. This is one Newton step on the loop's head
//! balance, and because the correction only circulates flow around a closed
//! loop it never disturbs node continuity — so a continuity-satisfying
//! initial guess stays continuity-satisfying throughout, while the loop
//! sums converge to zero. [`Network::solve`] sweeps all loops Gauss-Seidel
//! style until the largest `|dQ|` falls below [`SolveConfig::tolerance`].
//!
//! The resistance coefficient `k` may be supplied directly or derived from
//! geometry: [`darcy_weisbach_k`] gives `k = 8 f L / (pi^2 g D^5)` for the
//! Darcy-Weisbach law, and [`hazen_williams_k`] linearises the empirical
//! Hazen-Williams law to the quadratic form at a chosen reference flow.
//!
//! ## Honest scope
//!
//! Research / educational grade. This crate implements **textbook
//! closed-form / numerical models** — the classical Hardy-Cross loop-balance
//! iteration with Darcy-Weisbach / (linearised) Hazen-Williams resistance.
//! It is **NOT a clinical / medical / production engineering tool**. It does
//! not model transients / water hammer, compressibility, pumps and
//! variable-head sources, two-phase or non-Newtonian flow, temperature
//! dependence, cavitation, or the full `q^1.852` Hazen-Williams exponent
//! (that law is fitted to a quadratic here). Results are only as good as the
//! caller-supplied resistance coefficients and the modelling assumptions
//! above; do not use them for real-world design, safety, or regulatory
//! decisions.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod loop_solver;
pub mod network;
pub mod pipe;

pub use error::{ErrorCategory, NetworkError};
pub use loop_solver::{Loop, LoopMember, Orientation};
pub use network::{Endpoints, Network, SolveConfig, SolveReport};
pub use pipe::{darcy_weisbach_k, hazen_williams_k, Pipe, GRAVITY_M_S2};
