//! Approximate Bayesian Computation (ABC).
//!
//! ABC is the standard way to fit a population-genetic model whose
//! likelihood is intractable but which is cheap to *simulate*. The
//! rejection algorithm (Pritchard et al. 1999) is:
//!
//! 1. Draw a parameter vector `theta` from the prior.
//! 2. Simulate a dataset under `theta`.
//! 3. Reduce the dataset to a vector of summary statistics.
//! 4. **Accept** `theta` if its statistics are within tolerance of the
//!    *observed* statistics, **reject** otherwise.
//!
//! The accepted `theta`s approximate the posterior distribution.
//!
//! This module supplies the framework as generic building blocks:
//!
//! - [`Prior`] — a uniform prior over a parameter vector.
//! - [`AbcConfig`] — the number of simulations, tolerance and distance
//!   metric.
//! - [`abc_reject`] — runs the rejection sampler given user-supplied
//!   *simulate* and *summarise* closures, returning an
//!   [`AbcPosterior`].
//!
//! The caller plugs in any simulator (`valenx-popgen`'s Wright-Fisher
//! or coalescent are natural choices) by passing a closure that maps a
//! parameter vector to a summary-statistic vector.

use crate::error::{PopgenError, Result};
use crate::rng::Rng;

/// A uniform prior over a fixed-length parameter vector: each component
/// `i` is drawn uniformly from `[lower[i], upper[i])`.
#[derive(Clone, Debug, PartialEq)]
pub struct Prior {
    lower: Vec<f64>,
    upper: Vec<f64>,
}

impl Prior {
    /// Builds a uniform prior from per-parameter `[lower, upper]`
    /// bounds.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] on an empty bound list;
    /// [`PopgenError::Dimension`] if the two bound vectors differ in
    /// length; [`PopgenError::Invalid`] if any `lower >= upper`.
    pub fn uniform(lower: Vec<f64>, upper: Vec<f64>) -> Result<Self> {
        if lower.is_empty() {
            return Err(PopgenError::invalid("prior", "no parameters"));
        }
        if lower.len() != upper.len() {
            return Err(PopgenError::dimension(
                lower.len(),
                upper.len(),
                "prior bound vectors",
            ));
        }
        for (lo, hi) in lower.iter().zip(&upper) {
            if lo >= hi {
                return Err(PopgenError::invalid(
                    "prior",
                    "every lower bound must be below its upper bound",
                ));
            }
        }
        Ok(Prior { lower, upper })
    }

    /// Number of parameters.
    pub fn dimension(&self) -> usize {
        self.lower.len()
    }

    /// Draws one parameter vector from the prior.
    pub fn sample(&self, rng: &mut Rng) -> Vec<f64> {
        self.lower
            .iter()
            .zip(&self.upper)
            .map(|(&lo, &hi)| rng.uniform_range(lo, hi))
            .collect()
    }
}

/// Distance metric for comparing summary-statistic vectors.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Distance {
    /// Euclidean (`L2`) distance.
    Euclidean,
    /// Manhattan (`L1`) distance.
    Manhattan,
}

impl Distance {
    /// Computes the distance between two equal-length statistic
    /// vectors.
    fn between(&self, a: &[f64], b: &[f64]) -> f64 {
        match self {
            Distance::Euclidean => a
                .iter()
                .zip(b)
                .map(|(x, y)| (x - y).powi(2))
                .sum::<f64>()
                .sqrt(),
            Distance::Manhattan => {
                a.iter().zip(b).map(|(x, y)| (x - y).abs()).sum()
            }
        }
    }
}

/// Configuration for an ABC rejection run.
#[derive(Copy, Clone, Debug)]
pub struct AbcConfig {
    /// Number of prior draws / simulations to attempt.
    pub n_simulations: usize,
    /// Acceptance tolerance: a draw is accepted if its summary
    /// distance to the observed statistics is below this.
    pub tolerance: f64,
    /// Distance metric on the summary-statistic space.
    pub distance: Distance,
    /// RNG seed.
    pub seed: u64,
}

/// The posterior sample produced by ABC: the accepted parameter
/// vectors and their distances to the observed data.
#[derive(Clone, Debug, PartialEq)]
pub struct AbcPosterior {
    /// Accepted parameter vectors.
    accepted: Vec<Vec<f64>>,
    /// Distance of each accepted draw to the observed statistics.
    distances: Vec<f64>,
    /// Total number of simulations run.
    n_simulations: usize,
}

impl AbcPosterior {
    /// The accepted parameter vectors (the posterior sample).
    pub fn accepted(&self) -> &[Vec<f64>] {
        &self.accepted
    }

    /// Distances of the accepted draws to the observed data.
    pub fn distances(&self) -> &[f64] {
        &self.distances
    }

    /// Number of accepted draws.
    pub fn acceptance_count(&self) -> usize {
        self.accepted.len()
    }

    /// Acceptance rate: accepted / simulated.
    pub fn acceptance_rate(&self) -> f64 {
        if self.n_simulations == 0 {
            0.0
        } else {
            self.accepted.len() as f64 / self.n_simulations as f64
        }
    }

    /// Posterior mean of each parameter, or `None` if nothing was
    /// accepted.
    pub fn posterior_mean(&self) -> Option<Vec<f64>> {
        if self.accepted.is_empty() {
            return None;
        }
        let dim = self.accepted[0].len();
        let mut mean = vec![0.0; dim];
        for theta in &self.accepted {
            for (m, &v) in mean.iter_mut().zip(theta) {
                *m += v;
            }
        }
        for m in &mut mean {
            *m /= self.accepted.len() as f64;
        }
        Some(mean)
    }

    /// Posterior variance of each parameter, or `None` if fewer than
    /// two draws were accepted.
    pub fn posterior_variance(&self) -> Option<Vec<f64>> {
        if self.accepted.len() < 2 {
            return None;
        }
        let mean = self.posterior_mean()?;
        let dim = mean.len();
        let mut var = vec![0.0; dim];
        for theta in &self.accepted {
            for (k, v) in var.iter_mut().enumerate() {
                *v += (theta[k] - mean[k]).powi(2);
            }
        }
        for v in &mut var {
            *v /= self.accepted.len() as f64;
        }
        Some(var)
    }
}

/// Runs the ABC rejection algorithm.
///
/// `prior` supplies parameter draws; `simulate_and_summarise` maps a
/// parameter vector to a summary-statistic vector (this is where the
/// caller's simulator lives); `observed` is the summary of the real
/// data. A draw is accepted when its summary distance to `observed` is
/// below [`AbcConfig::tolerance`].
///
/// # Errors
/// [`PopgenError::Invalid`] on a zero simulation count or negative
/// tolerance; [`PopgenError::Dimension`] if a simulated summary vector
/// has a different length from `observed`.
pub fn abc_reject<F>(
    prior: &Prior,
    mut simulate_and_summarise: F,
    observed: &[f64],
    config: AbcConfig,
) -> Result<AbcPosterior>
where
    F: FnMut(&[f64]) -> Vec<f64>,
{
    if config.n_simulations == 0 {
        return Err(PopgenError::invalid(
            "n_simulations",
            "must be positive",
        ));
    }
    if config.tolerance < 0.0 {
        return Err(PopgenError::invalid("tolerance", "must be non-negative"));
    }
    if observed.is_empty() {
        return Err(PopgenError::invalid(
            "observed",
            "observed statistics vector is empty",
        ));
    }
    let mut rng = Rng::new(config.seed);
    let mut accepted = Vec::new();
    let mut distances = Vec::new();
    for _ in 0..config.n_simulations {
        let theta = prior.sample(&mut rng);
        let summary = simulate_and_summarise(&theta);
        if summary.len() != observed.len() {
            return Err(PopgenError::dimension(
                observed.len(),
                summary.len(),
                "simulated summary statistics",
            ));
        }
        let d = config.distance.between(&summary, observed);
        if d <= config.tolerance {
            accepted.push(theta);
            distances.push(d);
        }
    }
    Ok(AbcPosterior {
        accepted,
        distances,
        n_simulations: config.n_simulations,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn prior_samples_within_bounds() {
        let prior = Prior::uniform(vec![0.0, 10.0], vec![1.0, 20.0]).unwrap();
        let mut rng = Rng::new(1);
        for _ in 0..1000 {
            let theta = prior.sample(&mut rng);
            assert!((0.0..1.0).contains(&theta[0]));
            assert!((10.0..20.0).contains(&theta[1]));
        }
    }

    #[test]
    fn abc_recovers_a_known_parameter() {
        // Toy model: the "summary" is just the parameter itself, so
        // ABC should concentrate accepted draws near the observed
        // value.
        let prior = Prior::uniform(vec![0.0], vec![100.0]).unwrap();
        let observed = [42.0];
        let post = abc_reject(
            &prior,
            |theta| vec![theta[0]],
            &observed,
            AbcConfig {
                n_simulations: 20_000,
                tolerance: 2.0,
                distance: Distance::Euclidean,
                seed: 7,
            },
        )
        .unwrap();
        assert!(post.acceptance_count() > 0);
        let mean = post.posterior_mean().unwrap();
        // Accepted draws all lie within tolerance of 42.
        assert!((mean[0] - 42.0).abs() < 2.0, "posterior mean = {}", mean[0]);
    }

    #[test]
    fn tighter_tolerance_accepts_fewer() {
        let prior = Prior::uniform(vec![0.0], vec![100.0]).unwrap();
        let observed = [50.0];
        let run = |tol: f64| {
            abc_reject(
                &prior,
                |theta| vec![theta[0]],
                &observed,
                AbcConfig {
                    n_simulations: 10_000,
                    tolerance: tol,
                    distance: Distance::Euclidean,
                    seed: 3,
                },
            )
            .unwrap()
        };
        let loose = run(10.0);
        let tight = run(1.0);
        assert!(tight.acceptance_count() < loose.acceptance_count());
        assert!(tight.acceptance_rate() < loose.acceptance_rate());
    }

    #[test]
    fn posterior_mean_and_variance_are_consistent() {
        let prior = Prior::uniform(vec![0.0, 0.0], vec![10.0, 10.0]).unwrap();
        let observed = [5.0, 5.0];
        let post = abc_reject(
            &prior,
            |theta| vec![theta[0], theta[1]],
            &observed,
            AbcConfig {
                n_simulations: 5_000,
                tolerance: 3.0,
                distance: Distance::Manhattan,
                seed: 1,
            },
        )
        .unwrap();
        if post.acceptance_count() >= 2 {
            let mean = post.posterior_mean().unwrap();
            let var = post.posterior_variance().unwrap();
            assert_eq!(mean.len(), 2);
            assert_eq!(var.len(), 2);
            assert!(var.iter().all(|&v| v >= 0.0));
        }
    }

    #[test]
    fn abc_is_deterministic() {
        let prior = Prior::uniform(vec![0.0], vec![100.0]).unwrap();
        let cfg = AbcConfig {
            n_simulations: 2_000,
            tolerance: 5.0,
            distance: Distance::Euclidean,
            seed: 9,
        };
        let a = abc_reject(&prior, |t| vec![t[0]], &[30.0], cfg).unwrap();
        let b = abc_reject(&prior, |t| vec![t[0]], &[30.0], cfg).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn rejects_bad_input() {
        let prior = Prior::uniform(vec![0.0], vec![1.0]).unwrap();
        let cfg = AbcConfig {
            n_simulations: 0,
            tolerance: 1.0,
            distance: Distance::Euclidean,
            seed: 1,
        };
        assert!(abc_reject(&prior, |t| vec![t[0]], &[0.5], cfg).is_err());
        // A summary of the wrong length.
        let cfg2 = AbcConfig {
            n_simulations: 10,
            tolerance: 1.0,
            distance: Distance::Euclidean,
            seed: 1,
        };
        assert!(
            abc_reject(&prior, |_| vec![0.0, 0.0], &[0.5], cfg2).is_err()
        );
        // Bad prior bounds.
        assert!(Prior::uniform(vec![5.0], vec![1.0]).is_err());
    }
}
