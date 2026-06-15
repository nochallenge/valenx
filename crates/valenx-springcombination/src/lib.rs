//! # valenx-springcombination
//!
//! Closed-form combinator for ideal linear (Hookean) springs.
//!
//! ## What
//!
//! Given one or more springs, each described by a single positive rate
//! `k` (newtons per metre), this crate computes:
//!
//! - the equivalent rate of springs in **parallel** ([`parallel_rate`]),
//! - the equivalent rate of springs in **series** ([`series_rate`]),
//! - either of those through one dispatcher ([`combine`] /
//!   [`Combination`]),
//! - the restoring **force** `F = k x` ([`Spring::force`],
//!   [`force_from_rate`]), and
//! - the stored elastic **energy** `U = 0.5 k x^2` ([`Spring::energy`],
//!   [`energy_from_rate`]).
//!
//! Every fallible call returns [`Result<_, SpringError>`]; the error type
//! exposes a stable [`code`](SpringError::code) and a coarse
//! [`category`](SpringError::category).
//!
//! ## Model
//!
//! A spring is reduced to its rate alone. The governing relations are the
//! textbook ones:
//!
//! ```text
//! F    = k * x                    (Hooke's law)
//! U    = 0.5 * k * x^2            (elastic potential energy)
//! parallel:  k_eq = sum(k_i)               (rates / stiffnesses add)
//! series:    1 / k_eq = sum(1 / k_i)       (compliances add)
//! ```
//!
//! Parallel arises when every spring shares the same displacement and the
//! forces add; series arises when every spring carries the same force and
//! the displacements add. The canonical checks fall straight out: two
//! equal springs of rate `k` give `2k` in parallel and `k/2` in series,
//! parallel is never softer than its stiffest member, and series is never
//! stiffer than its softest member.
//!
//! ```
//! use valenx_springcombination::{combine, Combination, Spring};
//!
//! let pair = [Spring::new(100.0).unwrap(), Spring::new(100.0).unwrap()];
//! let parallel = combine(Combination::Parallel, &pair).unwrap();
//! let series = combine(Combination::Series, &pair).unwrap();
//! assert!((parallel - 200.0).abs() < 1e-12); // 2k
//! assert!((series - 50.0).abs() < 1e-12); //   k/2
//! ```
//!
//! ## Honest scope
//!
//! Research / educational grade. These are the ideal-spring, single-
//! linear-regime formulas from a first-year mechanics course: textbook
//! closed-form / numerical models, **not** a clinical, medical, or
//! production engineering tool. They model a perfect Hookean element with
//! no preload, no end conditions, no buckling, no fatigue, no hysteresis,
//! no large-deflection non-linearity, no material creep and no
//! manufacturing tolerance, and they describe nothing about the geometry
//! or material of a physical coil. Do not size load-bearing or
//! safety-critical hardware from these numbers.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod springs;

pub use error::{ErrorCategory, SpringError};
pub use springs::{
    combine, energy_from_rate, force_from_rate, parallel_rate, series_rate, Combination, Spring,
};
