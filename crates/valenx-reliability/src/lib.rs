//! # valenx-reliability
//!
//! Closed-form reliability-engineering math: component lifetime models
//! and the algebra of building system reliability from component
//! reliabilities.
//!
//! ## What
//!
//! Three small, exact building blocks of classical reliability theory:
//!
//! - [`Exponential`] — the constant-hazard ("useful life") lifetime
//!   model. Reliability `R(t) = exp(-lambda t)`, mean time between
//!   failures `MTBF = 1 / lambda`, a flat hazard, and the memoryless
//!   property.
//! - [`Weibull`] — the two-parameter Weibull lifetime model
//!   `R(t) = exp(-(t / eta)^beta)`. A single shape parameter `beta`
//!   spans infant mortality (`beta < 1`), constant hazard (`beta = 1`,
//!   which reduces exactly to the exponential), and wear-out
//!   (`beta > 1`); the mean life `eta * Gamma(1 + 1 / beta)` is
//!   evaluated through a Lanczos [`gamma`] approximation. Its
//!   [`conditional_reliability`](Weibull::conditional_reliability) gives
//!   the mission reliability of an already-aged part — non-memoryless for
//!   `beta != 1`.
//! - [`system`] — reliability block diagrams: a [`system::series`]
//!   chain `prod(R_i)`, a [`system::parallel`] (active-redundant) block
//!   `1 - prod(1 - R_i)`, and the general [`system::k_out_of_n`]
//!   binomial structure.
//!
//! ```
//! use valenx_reliability::{Exponential, Weibull, system};
//!
//! // A part with MTBF = 1000 h: what fraction survives 500 h?
//! let part = Exponential::from_mtbf(1000.0).unwrap();
//! let r500 = part.reliability(500.0).unwrap();
//! assert!((r500 - (-0.5_f64).exp()).abs() < 1e-12);
//!
//! // Wear-out part (Weibull beta = 2) and its characteristic life.
//! let bearing = Weibull::new(2.0, 1500.0).unwrap();
//! assert!((bearing.reliability(1500.0).unwrap() - (-1.0_f64).exp()).abs() < 1e-12);
//!
//! // Two such parts in parallel (active redundancy) at 500 h.
//! let redundant = system::parallel(&[r500, r500]).unwrap();
//! assert!(redundant > r500);
//! ```
//!
//! ## Model
//!
//! Every quantity is a textbook closed form, evaluated directly — there
//! is no fitting, simulation, or sampling here:
//!
//! - Exponential: `R = exp(-lambda t)`, `F = 1 - R`,
//!   `f = lambda exp(-lambda t)`, `h = lambda`, `MTBF = 1 / lambda`,
//!   quantile `t = -ln(target) / lambda`.
//! - Weibull: `R = exp(-(t / eta)^beta)`,
//!   `h = (beta / eta) (t / eta)^(beta - 1)`, `f = h R`,
//!   `MTTF = eta Gamma(1 + 1 / beta)`, quantile
//!   `t = eta (-ln(target))^(1 / beta)`, conditional mission reliability
//!   `R(m | age) = R(age + m) / R(age)`.
//! - Systems (independent components): series `prod(R_i)`, parallel
//!   `1 - prod(1 - R_i)`, `k`-of-`n`
//!   `sum_{j=k}^{n} C(n, j) r^j (1 - r)^(n - j)`.
//!
//! All inputs are validated: rates, shapes and scales must be strictly
//! positive, times non-negative, and reliabilities finite probabilities
//! in `[0, 1]`. Invalid inputs return a [`ReliabilityError`] rather than
//! producing a meaningless number.
//!
//! ## Honest scope
//!
//! Research/educational grade. This crate implements textbook
//! closed-form and numerical reliability models for **studying and
//! teaching reliability theory** and for back-of-the-envelope analysis.
//! It is deliberately limited and is **not** a clinical, medical, or
//! production safety / reliability-certification engineering tool:
//!
//! - All components are assumed to fail **independently**; there is no
//!   common-cause failure, no load sharing, and no time-correlated
//!   degradation.
//! - The system models are **static** (reliabilities at a single
//!   instant). There is no repair, no maintenance / renewal process, no
//!   Markov availability model, and no minimal-cut / fault-tree
//!   enumeration of arbitrary topologies beyond series, parallel, and
//!   `k`-of-`n`.
//! - Parameters are taken as **exact and known**; there is no
//!   statistical estimation, censoring, confidence interval, or
//!   Bayesian updating from field data.
//! - It is **single-failure-mode** per component and covers only the
//!   exponential and two-parameter Weibull lifetime laws (no lognormal,
//!   gamma, mixture, or competing-risk models).
//!
//! Do not use it to certify, qualify, or make safety decisions about any
//! real system. Validate any number that matters against an established
//! reliability tool and the underlying assumptions.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod exponential;
pub mod system;
pub mod weibull;

pub use error::{ErrorCategory, ReliabilityError};
pub use exponential::Exponential;
pub use weibull::{gamma, Weibull};
