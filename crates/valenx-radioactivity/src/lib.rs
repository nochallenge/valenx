//! # valenx-radioactivity
//!
//! Closed-form radioactive-decay kinematics: the textbook exponential
//! decay law for a single nuclide and the two-member Bateman parent ->
//! daughter chain.
//!
//! ## What
//!
//! Two small, fully analytic models:
//!
//! - [`Nuclide`] — a single radioactive species defined by its decay
//!   constant `lambda`. Gives the half-life, the mean (average) life, the
//!   number-remaining-vs-time law `N(t) = N0 * exp(-lambda * t)`, the
//!   activity `A = lambda * N` and its own decay `A(t) = A0 *
//!   exp(-lambda * t)`, the time to reach an arbitrary remaining
//!   fraction, and the number of decays over a time interval (the
//!   time-integral of activity, i.e. the cumulated-disintegration
//!   count). Constructible from any one of `lambda`, the half-life, or
//!   the mean life.
//! - [`DecayChain`] — a parent nuclide feeding a single radioactive
//!   daughter. Gives both populations and activities over time from the
//!   two-member Bateman solution, the daughter's grow-in-then-decay
//!   transient, the time of the daughter activity maximum, and the
//!   transient-/secular-equilibrium activity ratio.
//!
//! ## Model
//!
//! Radioactive decay is a first-order (exponential) process: every
//! surviving nucleus has a fixed probability per unit time `lambda` of
//! decaying, independent of age and of every other nucleus.
//!
//! - Single nuclide: `N(t) = N0 * exp(-lambda * t)`; half-life
//!   `t_half = ln(2) / lambda`; mean life `tau = 1 / lambda =
//!   t_half / ln(2)`; activity `A = lambda * N`, decaying as
//!   `A(t) = A0 * exp(-lambda * t)`. The number of decays in an interval
//!   is `N(t1) - N(t2)` — the time-integral of the activity, which over
//!   all time totals the initial population `N0`.
//! - Parent -> daughter (Bateman, with initial daughter `Nd0`):
//!
//!   ```text
//!   Np(t) = Np0 * exp(-lambda_p * t)
//!   Nd(t) = Np0 * lambda_p / (lambda_d - lambda_p)
//!               * (exp(-lambda_p * t) - exp(-lambda_d * t))
//!           + Nd0 * exp(-lambda_d * t)
//!   ```
//!
//!   The daughter activity peaks at
//!   `t_max = ln(lambda_d / lambda_p) / (lambda_d - lambda_p)`. For a
//!   longer-lived parent (`lambda_p < lambda_d`) the activities settle
//!   into transient equilibrium with ratio
//!   `A_d / A_p -> lambda_d / (lambda_d - lambda_p)`, tending to one
//!   (secular equilibrium) as `lambda_p / lambda_d -> 0`.
//!
//! All quantities are unit-agnostic: every time-valued output shares the
//! time unit in which `lambda` was supplied, and activities come out in
//! the reciprocal of that unit (becquerel when `lambda` is per-second).
//!
//! ## Honest scope
//!
//! Research/educational grade. These are textbook closed-form / numerical
//! models of idealised first-order kinetics: a single decay mode per
//! nuclide, no branching ratios, no nuclear data tables, no measurement
//! uncertainty, and only the first two members of a chain. It is NOT a
//! clinical/medical tool (no dosimetry, radiotherapy, or nuclear-medicine
//! decisions) and NOT a production engineering tool (no reactor,
//! shielding, criticality, or waste analysis). The output is the value of
//! a formula, not a measurement or a safety determination.
//!
//! ## Example
//!
//! ```
//! use valenx_radioactivity::Nuclide;
//!
//! // Iodine-131: half-life ~ 8.02 days.
//! let i131 = Nuclide::from_half_life(8.02).unwrap();
//! // After one half-life exactly half the sample remains.
//! let remaining = i131.remaining_fraction(i131.half_life()).unwrap();
//! assert!((remaining - 0.5).abs() < 1e-12);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod chain;
pub mod decay;
pub mod error;

pub use chain::DecayChain;
pub use decay::{Nuclide, LN_2};
pub use error::{ErrorCategory, RadioactivityError, Result};
