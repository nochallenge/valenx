//! Priors for the Bayesian MCMC framework.
//!
//! In Bayesian phylogenetics the posterior is `P(tree, θ | data) ∝
//! P(data | tree, θ) · P(tree, θ)`. The first factor is the likelihood
//! (Felsenstein pruning, in [`crate::likelihood`]). This module encodes
//! the second — the priors.
//!
//! Three pieces, all evaluable in log space:
//!
//! - **Topology prior** — by default uniform over labelled rooted-tree
//!   topologies, the BEAST 2 default for a strict-clock analysis.
//! - **Branch-length prior** — independent **Exponential(λ)** per
//!   branch, the MrBayes default. The total log prior is the sum of
//!   independent log densities.
//! - **Model-parameter prior** — per substitution-model variant:
//!   * `Jc69` / `F81` — no parameters, log prior is 0.
//!   * `K80` / `Hky85` — `κ ~ Exp(λ_κ)` (positive, mean 1/λ).
//!   * `Gtr` — six exchangeabilities under a symmetric `Dirichlet(α)`
//!     after normalisation to sum to 1; equilibrium frequencies under a
//!     symmetric `Dirichlet(α_π)`. (The Dirichlet density is on the
//!     simplex of normalised rates; the Hastings ratio for the
//!     unnormalised Dirichlet move is handled in [`super::proposal`].)
//!
//! The values here are log densities, not log probabilities of discrete
//! events — for continuous parameters the topology / branch-length
//! prior ratios in the MH acceptance work because the *ratio* of two
//! prior densities is well-defined regardless of an absolute
//! normalising constant.

use crate::error::{PhyloError, Result};
use crate::likelihood::gamma::ln_gamma;
use crate::likelihood::model::SubstModel;
use crate::tree::Tree;

/// Joint prior over `(tree, branch lengths, substitution-model params)`.
#[derive(Debug, Clone)]
pub struct Prior {
    /// Rate `λ` of the per-branch Exponential prior on branch lengths.
    /// Mean branch length under the prior is `1/λ`. Must be positive.
    pub branch_length_rate: f64,
    /// Rate of the Exponential prior on the substitution-model `κ`
    /// (transition / transversion ratio) for K80 / HKY85. Must be
    /// positive.
    pub kappa_rate: f64,
    /// Concentration `α` of the symmetric Dirichlet prior on the six
    /// GTR exchangeabilities (after normalisation to sum to 1).
    pub gtr_rate_concentration: f64,
    /// Concentration `α_π` of the symmetric Dirichlet prior on the
    /// four equilibrium frequencies of F81 / HKY85 / GTR.
    pub frequency_concentration: f64,
    /// Rate `λ_α` of the Exponential prior on the gamma shape `α`
    /// (rate-heterogeneity parameter). Used when the chain optimises α
    /// jointly with the other parameters.
    pub gamma_alpha_rate: f64,
}

impl Default for Prior {
    /// Reasonable defaults: a diffuse-but-proper prior set close to the
    /// MrBayes defaults (`λ = 10`, `α = 1`).
    fn default() -> Self {
        Prior {
            branch_length_rate: 10.0,
            kappa_rate: 1.0,
            gtr_rate_concentration: 1.0,
            frequency_concentration: 1.0,
            gamma_alpha_rate: 1.0,
        }
    }
}

impl Prior {
    /// Validates the prior hyperparameters are positive.
    ///
    /// # Errors
    /// [`PhyloError::Invalid`] if any rate / concentration is `<= 0`.
    pub fn validate(&self) -> Result<()> {
        let positive = |name: &'static str, v: f64| -> Result<()> {
            if v.is_finite() && v > 0.0 {
                Ok(())
            } else {
                Err(PhyloError::invalid(name, "must be a positive finite real"))
            }
        };
        positive("branch_length_rate", self.branch_length_rate)?;
        positive("kappa_rate", self.kappa_rate)?;
        positive("gtr_rate_concentration", self.gtr_rate_concentration)?;
        positive("frequency_concentration", self.frequency_concentration)?;
        positive("gamma_alpha_rate", self.gamma_alpha_rate)?;
        Ok(())
    }

    /// Log prior density of the tree topology.
    ///
    /// Uniform over labelled binary topologies: `−log(N)` where `N` is
    /// the number of labelled trees on `n` leaves. The *ratio* of two
    /// topology priors is `1`, so the absolute value drops out of the
    /// MH acceptance — this returns `0.0`, the simplest correct value.
    pub fn log_topology(&self, _tree: &Tree) -> f64 {
        0.0
    }

    /// Log prior density of all branch lengths in `tree` under the
    /// per-branch Exponential prior.
    ///
    /// For each non-root node, `log p(t) = log λ − λ·t`.
    pub fn log_branch_lengths(&self, tree: &Tree) -> f64 {
        let lambda = self.branch_length_rate;
        let mut acc = 0.0;
        for id in 0..tree.node_count() {
            let node = tree.node(id);
            if node.parent.is_none() {
                continue;
            }
            let t = node.branch_length.unwrap_or(0.0).max(0.0);
            // log f(t; λ) = log λ - λ·t for t ≥ 0.
            acc += lambda.ln() - lambda * t;
        }
        acc
    }

    /// Log prior density of the substitution-model parameters.
    ///
    /// Variant-specific:
    /// - `Jc69` — `0.0` (no free parameters).
    /// - `K80 { kappa }` — `Exp(kappa_rate)`.
    /// - `F81 { freqs }` — symmetric `Dirichlet(α_π)` on frequencies.
    /// - `Hky85 { kappa, freqs }` — sum of the two above.
    /// - `Gtr { rates, freqs }` — Dirichlet on normalised rates +
    ///   Dirichlet on frequencies.
    pub fn log_model(&self, model: &SubstModel) -> f64 {
        match model {
            SubstModel::Jc69 => 0.0,
            SubstModel::K80 { kappa } => {
                if *kappa <= 0.0 {
                    f64::NEG_INFINITY
                } else {
                    self.kappa_rate.ln() - self.kappa_rate * kappa
                }
            }
            SubstModel::F81 { freqs } => {
                ln_dirichlet_pdf(freqs, self.frequency_concentration)
            }
            SubstModel::Hky85 { kappa, freqs } => {
                let k_lp = if *kappa <= 0.0 {
                    f64::NEG_INFINITY
                } else {
                    self.kappa_rate.ln() - self.kappa_rate * kappa
                };
                k_lp + ln_dirichlet_pdf(freqs, self.frequency_concentration)
            }
            SubstModel::Gtr { rates, freqs } => {
                let total: f64 = rates.iter().sum();
                if total <= 0.0 || rates.iter().any(|r| !r.is_finite() || *r <= 0.0) {
                    return f64::NEG_INFINITY;
                }
                let normalised: Vec<f64> = rates.iter().map(|r| r / total).collect();
                ln_dirichlet_pdf(&normalised, self.gtr_rate_concentration)
                    + ln_dirichlet_pdf(freqs, self.frequency_concentration)
            }
        }
    }

    /// Log prior density of a gamma rate-heterogeneity shape `α` under
    /// the Exponential prior — `log λ_α − λ_α·α`.
    pub fn log_gamma_alpha(&self, alpha: f64) -> f64 {
        if !alpha.is_finite() || alpha <= 0.0 {
            return f64::NEG_INFINITY;
        }
        self.gamma_alpha_rate.ln() - self.gamma_alpha_rate * alpha
    }
}

/// Log density of a Dirichlet(α, α, …, α) distribution at the simplex
/// point `xs` (must sum to 1, each component positive).
///
/// `log B(α₁, …, α_k) = Σ log Γ(α_i) − log Γ(Σ α_i)`, and
/// `log f(x; α) = − log B + Σ (α_i − 1) log x_i`. For a symmetric
/// Dirichlet `α_i = α` this simplifies; we evaluate the general form for
/// clarity. Returns `f64::NEG_INFINITY` on a degenerate point.
fn ln_dirichlet_pdf(xs: &[f64], alpha: f64) -> f64 {
    if alpha <= 0.0 {
        return f64::NEG_INFINITY;
    }
    let n = xs.len() as f64;
    if xs.iter().any(|&x| !x.is_finite() || x <= 0.0) {
        return f64::NEG_INFINITY;
    }
    let sum: f64 = xs.iter().sum();
    if (sum - 1.0).abs() > 1e-6 {
        return f64::NEG_INFINITY;
    }
    let log_b = n * ln_gamma(alpha) - ln_gamma(n * alpha);
    let mut acc = -log_b;
    for &x in xs {
        acc += (alpha - 1.0) * x.ln();
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::read_newick;

    #[test]
    fn default_prior_validates() {
        assert!(Prior::default().validate().is_ok());
    }

    #[test]
    fn bad_prior_rejected() {
        let p = Prior {
            branch_length_rate: 0.0,
            ..Prior::default()
        };
        assert!(p.validate().is_err());
        let p = Prior {
            kappa_rate: f64::NAN,
            ..Prior::default()
        };
        assert!(p.validate().is_err());
    }

    #[test]
    fn topology_log_prior_is_uniform() {
        let t1 = read_newick("((A,B),(C,D));").unwrap();
        let t2 = read_newick("((A,C),(B,D));").unwrap();
        let p = Prior::default();
        assert_eq!(p.log_topology(&t1), p.log_topology(&t2));
    }

    #[test]
    fn branch_length_prior_decreases_with_total_length() {
        let short = read_newick("((A:0.01,B:0.01):0.01,(C:0.01,D:0.01):0.01);").unwrap();
        let long = read_newick("((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0);").unwrap();
        let p = Prior::default();
        // Exp(10) puts most mass near 0; the long tree should have a
        // much lower prior density.
        assert!(p.log_branch_lengths(&short) > p.log_branch_lengths(&long));
    }

    #[test]
    fn model_log_prior_is_finite_for_jc69() {
        let p = Prior::default();
        assert_eq!(p.log_model(&SubstModel::Jc69), 0.0);
    }

    #[test]
    fn model_log_prior_finite_for_hky85() {
        let p = Prior::default();
        let lp = p.log_model(&SubstModel::Hky85 {
            kappa: 2.0,
            freqs: [0.25, 0.25, 0.25, 0.25],
        });
        assert!(lp.is_finite(), "lp = {lp}");
    }

    #[test]
    fn negative_kappa_has_zero_density() {
        let p = Prior::default();
        let lp = p.log_model(&SubstModel::K80 { kappa: -1.0 });
        assert_eq!(lp, f64::NEG_INFINITY);
    }

    #[test]
    fn gtr_log_prior_is_finite() {
        let p = Prior::default();
        let lp = p.log_model(&SubstModel::Gtr {
            rates: [1.0, 2.0, 0.5, 1.0, 2.0, 0.5],
            freqs: [0.25, 0.25, 0.25, 0.25],
        });
        assert!(lp.is_finite(), "lp = {lp}");
    }

    #[test]
    fn dirichlet_pdf_uniform_is_log_factorial() {
        // Dirichlet(1, 1, 1, 1) at uniform — log f = -log(3!) = -log 6
        // because the uniform Dirichlet density on the 3-simplex is
        // (k - 1)! everywhere.
        let xs = [0.25; 4];
        let lp = ln_dirichlet_pdf(&xs, 1.0);
        let want = 6.0_f64.ln();
        assert!((lp - want).abs() < 1e-9, "lp = {lp}, want = {want}");
    }

    #[test]
    fn dirichlet_pdf_off_simplex_is_neg_inf() {
        let xs = [0.1, 0.1, 0.1, 0.1]; // sums to 0.4, not 1
        assert_eq!(ln_dirichlet_pdf(&xs, 1.0), f64::NEG_INFINITY);
    }

    #[test]
    fn gamma_alpha_log_prior() {
        let p = Prior::default();
        let lp = p.log_gamma_alpha(0.5);
        assert!(lp.is_finite());
        assert_eq!(p.log_gamma_alpha(-1.0), f64::NEG_INFINITY);
    }
}
