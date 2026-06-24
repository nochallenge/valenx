//! # valenx-uq — uncertainty quantification, surrogates, and sensitivity
//!
//! An **in-house, dependency-light** toolkit for *uncertainty quantification*
//! (UQ): the study of how uncertainty in a model's **inputs** propagates to its
//! **outputs**, and of which inputs matter most. It is the cross-cutting
//! numerical enabler underneath simulation work — it wraps **any** model that
//! maps a vector of real inputs to a vector of real outputs.
//!
//! ## What UQ is (and what this crate provides)
//!
//! Given a deterministic model `y = f(x)` and a description of how the inputs
//! `x` are uncertain (each input is a random variable with some
//! [`Distribution`]), UQ answers questions such as:
//!
//! * **Forward propagation** — what is the distribution of `y`? Its mean,
//!   variance, percentiles, confidence intervals? Estimated by sampling the
//!   inputs and pushing each sample through the model
//!   (`sampling` + `statistics`).
//! * **Sensitivity analysis** — which inputs drive the variance of `y`?
//!   Answered globally by **Sobol indices** (first-order `S_i` and total
//!   `S_Ti`) and **Morris elementary effects** (`mu_star`, `sigma`)
//!   (`sensitivity`).
//! * **Surrogate / metamodelling** — when `f` is expensive, fit a cheap
//!   approximation `f̂ ≈ f` from a modest sample and evaluate *that* instead
//!   (`surrogate`), or build a spectral **polynomial-chaos expansion** whose
//!   coefficients give the output mean and variance analytically (`pce`).
//! * **Quasi-random sampling** — a deterministic, low-discrepancy **Sobol'
//!   sequence** that fills `[0,1)^d` more evenly than pseudo-random points, for
//!   faster-converging quasi-Monte-Carlo integration (`sobol`).
//! * **Verification & validation** — order-of-accuracy fits, Roache's Grid
//!   Convergence Index, a Method-of-Manufactured-Solutions helper, and error
//!   norms for checking that a solver converges at its formal rate (`vnv`).
//!
//! The model itself is abstracted by the [`Model`] trait, so the same UQ
//! machinery applies to a one-line closure ([`FnModel`]) in a test or to any
//! valenx solver wrapped behind the trait.
//!
//! ## Determinism — no `rand`, seeded SplitMix64
//!
//! Reproducibility is a first-class requirement: an analysis that gives a
//! different answer on every run is hard to validate and impossible to
//! regression-test. This crate therefore takes **no `rand` dependency**.
//! All randomness comes from a tiny in-crate [`SplitMix64`] PRNG (the same
//! deterministic, seeded generator used in `valenx-photogrammetry`), with
//! standard-normal draws produced by the Box–Muller transform. Given the same
//! seed, every sample, every Sobol estimate, and every fitted surrogate is
//! bit-for-bit identical across runs and machines. The PRNG is **not** used
//! for any security purpose.
//!
//! ## Honesty / scope caveats
//!
//! These are textbook, research/educational-grade methods. Specifically:
//!
//! * **Monte-Carlo convergence is `O(1/√n)`.** Halving the standard error of a
//!   mean estimate needs roughly four times the samples. Latin-hypercube
//!   sampling reduces the *constant* (better space-filling, lower variance for
//!   smooth integrands) but not the asymptotic rate.
//! * **The Saltelli Sobol estimator has finite-sample variance.** First-order
//!   and total indices are *estimates*; at small `n` they can fall slightly
//!   outside `[0, 1]` or sum to slightly more/less than 1. They converge to
//!   the analytic indices as `n` grows. This crate uses the standard
//!   Saltelli A/B/AB cross-sampling design without any further bias
//!   correction.
//! * **The surrogate is a *global low-order polynomial* response surface**
//!   (degree ≤ 2) fitted by least squares. It captures smooth global trends,
//!   not sharp local features, and it does **not** interpolate the training
//!   points exactly (unless the data are exactly polynomial of that degree).
//!   A Gaussian-process / kriging surrogate — which interpolates and yields a
//!   predictive variance — is a documented future extension, **not** provided
//!   here.
//! * Every estimate is only as good as the input [`Distribution`]s supplied;
//!   garbage in, garbage out.
//!
//! ## Example
//!
//! ```
//! use valenx_uq::{FnModel, Model, Distribution, SplitMix64};
//! use valenx_uq::sampling::monte_carlo;
//! use valenx_uq::statistics;
//!
//! // Model: y = x0 + x1.
//! let model = FnModel::new(2, 1, |x| vec![x[0] + x[1]]);
//!
//! // Both inputs standard-normal-ish.
//! let dists = [
//!     Distribution::normal(0.0, 1.0).unwrap(),
//!     Distribution::normal(5.0, 2.0).unwrap(),
//! ];
//!
//! // Propagate 10_000 Monte-Carlo samples through the model (fixed seed).
//! let mut rng = SplitMix64::new(0x5EED_C0DE);
//! let inputs = monte_carlo(10_000, &dists, &mut rng);
//! let outputs: Vec<f64> = inputs.iter().map(|x| model.evaluate(x)[0]).collect();
//!
//! // Mean of (x0 + x1) ≈ 0 + 5 = 5.
//! let mean = statistics::mean(&outputs).unwrap();
//! assert!((mean - 5.0).abs() < 0.1);
//! ```

#![forbid(unsafe_code)]

pub mod distribution;
pub mod model;
pub mod pce;
pub mod reliability;
pub mod sampling;
pub mod sensitivity;
pub mod sobol;
pub mod statistics;
pub mod surrogate;
pub mod vnv;

mod error;
mod rng;

pub use distribution::Distribution;
pub use error::UqError;
pub use model::{FnModel, Model};
pub use pce::{Pce, PolyBasis};
pub use reliability::{
    form, pf_monte_carlo, sorm_breitung, FormConfig, FormResult, McResult, SormResult,
};
pub use rng::SplitMix64;
pub use sensitivity::{MorrisResult, SobolIndices};
pub use sobol::sobol_sequence;
pub use surrogate::PolynomialSurrogate;
pub use vnv::{l2_norm, linf_norm, rms_norm, ConvergenceStudy, Mms};
