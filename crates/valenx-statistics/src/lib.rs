//! # valenx-statistics
//!
//! Descriptive and inferential statistics over `f64` slices — the
//! small-sample, closed-form toolkit behind Valenx's plots and analyses.
//!
//! ## What
//!
//! - Location and dispersion ([`descriptive`]): [`descriptive::mean`],
//!   [`descriptive::sample_variance`] / [`descriptive::population_variance`]
//!   and their standard deviations, the [`descriptive::median`],
//!   [`descriptive::quartiles`] and the general [`descriptive::quantile`].
//! - Standardisation: the [`descriptive::z_score`] of one observation and
//!   [`descriptive::standardize`] of a whole sample.
//! - The standard normal distribution ([`normal`]): the density
//!   [`normal::pdf`], the cumulative [`normal::cdf`], the survival
//!   [`normal::sf`] and the [`normal::erf`] underneath.
//! - Inference ([`inference`]): the one-sample [`inference::t_statistic`],
//!   the inverse-normal quantile [`inference::probit`], and a known-sigma
//!   [`inference::normal_ci`] confidence interval for the mean.
//! - A validated error taxonomy ([`error`]) with the shared input checks in
//!   [`validate`].
//!
//! ## Model
//!
//! Every estimator is a textbook formula evaluated directly:
//!
//! - mean `x̄ = (1/n) Σ xᵢ`; the integers `1..=n` give the closed form
//!   `x̄ = (n + 1) / 2`;
//! - sample variance `s² = (1/(n−1)) Σ (xᵢ − x̄)²` (Bessel-corrected),
//!   population variance the same numerator over `n`; standard deviations
//!   are their square roots;
//! - quantiles by the **linear-interpolation type-7** rule (the NumPy / R
//!   default), so the median is the central order statistic(s);
//! - the standard normal via `φ(z) = e^{−z²/2} / √(2π)` and
//!   `Φ(z) = ½[1 + erf(z/√2)]` with the Abramowitz & Stegun 7.1.26 `erf`;
//! - the z-score `(x − μ) / σ`, the one-sample t-statistic
//!   `(x̄ − μ₀)/(s/√n)`, and the known-sigma interval `x̄ ± z·σ/√n` whose
//!   margin grows with `σ` and shrinks like `1/√n` with the sample size.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are textbook closed-form and
//! numerical models — exact for the descriptive estimators, and accurate to
//! the underlying special-function approximations elsewhere (the A&S `erf`
//! to `≈ 1.5 × 10⁻⁷`, the Acklam inverse-CDF to `≈ 1.15 × 10⁻⁹`). The
//! confidence interval is the *known-sigma z-interval*; there is no
//! Student-t quantile, so for unknown `σ` and small `n` it runs slightly
//! narrow. This crate is **NOT a clinical, medical, or production
//! engineering statistics tool** — for regulated or high-stakes work use a
//! validated, certified statistics package. Pure algorithms, no platform
//! dependencies.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod descriptive;
pub mod error;
pub mod inference;
pub mod normal;
pub mod validate;

pub use descriptive::{
    mean, median, population_std, population_variance, quantile, quartiles, sample_std,
    sample_variance, standardize, z_score, Quartiles,
};
pub use error::{Result, StatsError};
pub use inference::{normal_ci, probit, t_statistic, ConfidenceInterval};
pub use normal::{cdf, erf, erfc, pdf, sf};
