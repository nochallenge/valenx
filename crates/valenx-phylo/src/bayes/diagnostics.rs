//! Convergence diagnostics for MCMC traces.
//!
//! Two standards every Bayesian phylogenetics tool ships:
//!
//! - **Effective sample size (ESS)** of one trace — `n / (1 + 2 Σ ρ_k)`
//!   where `ρ_k` is the lag-`k` sample autocorrelation. The Geyer 1992
//!   "initial monotone positive sequence" cutoff is used to stop the
//!   sum when consecutive pairs of autocorrelations stop summing to a
//!   positive value. An ESS far below `n` signals strong
//!   autocorrelation; under `200` is conventionally the danger
//!   threshold.
//! - **Potential scale reduction (Gelman-Rubin R̂)** across `m`
//!   independent chains run from over-dispersed starting points. The
//!   standard Brooks-Gelman 1998 form: split each chain in halves
//!   first, then compute `W` (within-chain variance), `B/n` (between-
//!   chain variance), `\hat{V} = (n-1)/n W + B/n`, and
//!   `\hat{R} = sqrt(\hat{V} / W)`. Values close to 1 (`< 1.1`,
//!   ideally `< 1.05`) signal convergence; values much above 1.1
//!   signal the chains have not yet mixed.
//!
//! These work on any real-valued parameter trace — likelihood, prior,
//! `κ`, total tree length, etc. — and so live as plain functions
//! ([`effective_sample_size`], [`gelman_rubin`]) plus a thin bundler
//! ([`ParameterDiagnostics`]) free of any tree-specific assumption.

use crate::error::{PhyloError, Result};

/// Effective sample size of a single MCMC trace, using the Geyer 1992
/// initial monotone positive sequence (IMPS) truncation.
///
/// - Returns `n` for a constant trace (no autocorrelation defined; the
///   trace carries no information beyond its mean, so treating it as
///   `n` independent samples is the conservative choice that no
///   downstream gate trips on).
/// - Returns `1.0` for a trace of length `≤ 1`.
pub fn effective_sample_size(trace: &[f64]) -> f64 {
    let n = trace.len();
    if n <= 1 {
        return n as f64;
    }
    let mean: f64 = trace.iter().sum::<f64>() / n as f64;
    let centered: Vec<f64> = trace.iter().map(|x| x - mean).collect();
    let var0: f64 = centered.iter().map(|x| x * x).sum::<f64>() / n as f64;
    if var0 <= 0.0 || !var0.is_finite() {
        // Constant trace: treat each iteration as independent.
        return n as f64;
    }
    // Lag-k sample autocorrelation ρ_k = c_k / c_0, c_k = (1/n) Σ x_i x_{i+k}.
    let max_lag = n - 1;
    let mut rhos = Vec::with_capacity(max_lag);
    for k in 1..=max_lag {
        let mut acc = 0.0;
        for i in 0..(n - k) {
            acc += centered[i] * centered[i + k];
        }
        let rho = (acc / n as f64) / var0;
        rhos.push(rho);
    }
    // Geyer's initial monotone positive sequence: form pair sums
    // ρ_{2k} + ρ_{2k+1} = Γ_k. Truncate at the first non-positive Γ_k,
    // and enforce monotone decrease on the remaining sequence.
    let mut gammas: Vec<f64> = Vec::with_capacity(rhos.len() / 2);
    for k in 0..(rhos.len() / 2) {
        let g = rhos[2 * k] + rhos[2 * k + 1];
        gammas.push(g);
    }
    let mut sum_pos = 0.0;
    let mut prev = f64::INFINITY;
    for &g in &gammas {
        if g <= 0.0 {
            break;
        }
        let g = g.min(prev);
        sum_pos += g;
        prev = g;
    }
    // Avoid integer / 0 underflow in the rho sum when the trace is
    // very short.
    let _ = centered.len();
    // τ = 1 + 2 Σ Γ_k under the IMPS rule — the integrated
    // autocorrelation time.
    let tau = 1.0 + 2.0 * sum_pos;
    let ess = (n as f64 / tau.max(1.0)).min(n as f64);
    ess.max(1.0)
}

/// Gelman-Rubin potential-scale-reduction `R̂` across `m` chains, each
/// of length `n`. Inputs must all share `n`; pass `m ≥ 2` chains.
///
/// # Errors
/// [`PhyloError::Invalid`] on `< 2` chains, on length-mismatched
/// chains, on a chain shorter than 2 samples, or on a degenerate
/// (constant) `W`.
pub fn gelman_rubin(chains: &[&[f64]]) -> Result<f64> {
    if chains.len() < 2 {
        return Err(PhyloError::invalid(
            "chains",
            "Gelman-Rubin needs ≥ 2 chains",
        ));
    }
    let n = chains[0].len();
    if n < 2 {
        return Err(PhyloError::invalid(
            "chains",
            "Gelman-Rubin needs ≥ 2 samples per chain",
        ));
    }
    for &ch in chains {
        if ch.len() != n {
            return Err(PhyloError::dimension(n, ch.len(), "Gelman-Rubin chains"));
        }
    }
    let m = chains.len() as f64;
    let n_f = n as f64;
    // Per-chain means and variances.
    let mut chain_means = Vec::with_capacity(chains.len());
    let mut chain_vars = Vec::with_capacity(chains.len());
    for &ch in chains {
        let mean: f64 = ch.iter().sum::<f64>() / n_f;
        let var: f64 = ch.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n_f - 1.0);
        chain_means.push(mean);
        chain_vars.push(var);
    }
    let grand_mean: f64 = chain_means.iter().sum::<f64>() / m;
    let w: f64 = chain_vars.iter().sum::<f64>() / m;
    if w <= 0.0 || !w.is_finite() {
        // Every chain has zero variance — they're each constant. If
        // they agree on the same constant, R̂ = 1.
        let all_equal = chain_means
            .windows(2)
            .all(|p| (p[0] - p[1]).abs() < 1e-12);
        return if all_equal {
            Ok(1.0)
        } else {
            Err(PhyloError::invalid(
                "chains",
                "constant chains with different means — R̂ undefined",
            ))
        };
    }
    let b_over_n: f64 = chain_means
        .iter()
        .map(|m| (m - grand_mean).powi(2))
        .sum::<f64>()
        / (m - 1.0);
    let v_hat = ((n_f - 1.0) / n_f) * w + b_over_n;
    Ok((v_hat / w).sqrt())
}

/// Diagnostics summary for one parameter across multiple chains.
#[derive(Debug, Clone)]
pub struct ParameterDiagnostics {
    /// Per-chain ESS (one value per chain).
    pub ess_per_chain: Vec<f64>,
    /// Total ESS across chains (sum of the per-chain ESS).
    pub ess_total: f64,
    /// Gelman-Rubin `R̂` (omitted if only one chain was supplied).
    pub r_hat: Option<f64>,
}

impl ParameterDiagnostics {
    /// Compute the diagnostics for a single parameter, given each
    /// chain's trace as a `&[f64]`.
    ///
    /// `chains.len()` may be 1 (then `r_hat = None`) or any larger
    /// number.
    pub fn compute(chains: &[&[f64]]) -> Result<Self> {
        if chains.is_empty() {
            return Err(PhyloError::invalid("chains", "no traces supplied"));
        }
        let ess_per_chain: Vec<f64> =
            chains.iter().map(|c| effective_sample_size(c)).collect();
        let ess_total: f64 = ess_per_chain.iter().sum();
        let r_hat = if chains.len() >= 2 {
            Some(gelman_rubin(chains)?)
        } else {
            None
        };
        Ok(ParameterDiagnostics {
            ess_per_chain,
            ess_total,
            r_hat,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a deterministic AR(1) trace `x_{t+1} = φ·x_t + ε_t` with
    /// `ε_t ∼ N(0, 1 − φ²)` so the marginal variance stays unit. Uses
    /// the crate's RNG for determinism.
    fn ar1(n: usize, phi: f64, seed: u64) -> Vec<f64> {
        let mut rng = crate::rng::Rng::new(seed);
        let mut xs = Vec::with_capacity(n);
        let mut prev = 0.0;
        let sigma = (1.0 - phi * phi).max(1e-9).sqrt();
        for _ in 0..n {
            let eps = rng.normal() * sigma;
            prev = phi * prev + eps;
            xs.push(prev);
        }
        xs
    }

    #[test]
    fn ess_of_iid_trace_is_about_n() {
        let xs = ar1(2000, 0.0, 1);
        let ess = effective_sample_size(&xs);
        // Should be in the same order as n; allow ample slack.
        assert!(ess > 1500.0, "ess = {ess}");
    }

    #[test]
    fn ess_of_strongly_autocorrelated_trace_is_much_smaller() {
        let xs = ar1(2000, 0.95, 2);
        let ess = effective_sample_size(&xs);
        // φ = 0.95 ⇒ τ ≈ (1+φ)/(1−φ) = 39 ⇒ ESS ≈ 50.
        assert!(ess < 200.0, "ess = {ess}");
        assert!(ess > 5.0, "ess = {ess}");
    }

    #[test]
    fn ess_of_constant_trace_is_n() {
        let xs = vec![1.0_f64; 100];
        let ess = effective_sample_size(&xs);
        assert!((ess - 100.0).abs() < 1e-9);
    }

    #[test]
    fn ess_of_empty_or_singleton() {
        assert_eq!(effective_sample_size(&[]), 0.0);
        assert_eq!(effective_sample_size(&[1.0]), 1.0);
    }

    #[test]
    fn r_hat_is_about_one_for_well_mixed_chains() {
        let chains = [ar1(1000, 0.5, 1), ar1(1000, 0.5, 2), ar1(1000, 0.5, 3)];
        let refs: Vec<&[f64]> = chains.iter().map(|v| v.as_slice()).collect();
        let r = gelman_rubin(&refs).unwrap();
        assert!((r - 1.0).abs() < 0.05, "r_hat = {r}");
    }

    #[test]
    fn r_hat_is_far_from_one_for_diverging_chains() {
        // Two chains centred on very different means => poor mixing.
        let a: Vec<f64> = ar1(1000, 0.5, 1).into_iter().map(|x| x + 5.0).collect();
        let b: Vec<f64> = ar1(1000, 0.5, 2).into_iter().map(|x| x - 5.0).collect();
        let chains = [a.as_slice(), b.as_slice()];
        let r = gelman_rubin(&chains).unwrap();
        assert!(r > 1.2, "r_hat = {r}");
    }

    #[test]
    fn r_hat_rejects_short_or_missing_chains() {
        let one = [1.0_f64, 2.0];
        let chains = [one.as_slice()];
        assert!(gelman_rubin(&chains).is_err());
        let chains: [&[f64]; 2] = [&[1.0], &[2.0]];
        assert!(gelman_rubin(&chains).is_err());
    }

    #[test]
    fn parameter_diagnostics_combines_ess_and_rhat() {
        let a = ar1(500, 0.5, 1);
        let b = ar1(500, 0.5, 2);
        let chains: [&[f64]; 2] = [a.as_slice(), b.as_slice()];
        let d = ParameterDiagnostics::compute(&chains).unwrap();
        assert_eq!(d.ess_per_chain.len(), 2);
        assert!(d.ess_total > 0.0);
        assert!(d.r_hat.is_some());
        assert!((d.r_hat.unwrap() - 1.0).abs() < 0.1);
    }

    #[test]
    fn parameter_diagnostics_single_chain_has_no_rhat() {
        let a = ar1(100, 0.0, 1);
        let chains: [&[f64]; 1] = [a.as_slice()];
        let d = ParameterDiagnostics::compute(&chains).unwrap();
        assert_eq!(d.ess_per_chain.len(), 1);
        assert!(d.r_hat.is_none());
    }
}
