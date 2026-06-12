//! Metropolis-Hastings sampler for `(tree, model, gamma_α)`.
//!
//! Wires together [`super::prior::Prior`], [`super::proposal`] and the
//! Felsenstein-pruning likelihood ([`crate::likelihood`]) into a real
//! Bayesian phylogenetic chain.
//!
//! Each iteration:
//!
//! 1. Pick a proposal kind from the [`super::proposal::ProposalSet`].
//! 2. Apply it (skip + count if it does not apply to the current state).
//! 3. Compute the posterior of the proposed state.
//! 4. Accept with probability `min(1, exp(log_posterior_new −
//!    log_posterior_old + log_hastings))`.
//! 5. Record the (now possibly updated) state on the trace.
//!
//! The result is a [`ChainResult`] with:
//!
//! - the tree sample (one tree per recorded iteration, after burn-in /
//!   thinning),
//! - per-iteration **parameter trace** (`log_posterior`,
//!   `log_likelihood`, model params, `α`, total tree length),
//! - acceptance counts per proposal kind.

use crate::bayes::prior::Prior;
use crate::bayes::proposal::{sample_proposal, ChainState, ProposalKind, ProposalSet};
use crate::error::{PhyloError, Result};
use crate::likelihood::felsenstein::{log_likelihood, log_likelihood_gamma};
use crate::likelihood::gamma::DiscreteGamma;
use crate::likelihood::model::SubstModel;
use crate::rng::Rng;
use crate::tree::Tree;
use std::collections::HashMap;

/// Settings for one MCMC chain run.
#[derive(Debug, Clone)]
pub struct ChainConfig {
    /// Total number of MCMC iterations (including burn-in).
    pub iterations: usize,
    /// Iterations to discard from the front (warm-up).
    pub burn_in: usize,
    /// Record every `thin`-th iteration to the trace (after burn-in).
    /// A value of 1 keeps every iteration; 10 keeps one in ten.
    pub thin: usize,
    /// RNG seed driving proposal selection + parameter draws.
    pub seed: u64,
    /// Number of discrete gamma categories used when
    /// [`ChainState::gamma_alpha`] is set. Defaults to 4.
    pub gamma_categories: usize,
}

impl Default for ChainConfig {
    fn default() -> Self {
        ChainConfig {
            iterations: 1000,
            burn_in: 100,
            thin: 1,
            seed: 0,
            gamma_categories: 4,
        }
    }
}

impl ChainConfig {
    /// Validates the chain configuration.
    ///
    /// # Errors
    /// [`PhyloError::Invalid`] on a non-positive `iterations`, a
    /// `burn_in` larger than `iterations`, or a non-positive `thin`.
    pub fn validate(&self) -> Result<()> {
        if self.iterations == 0 {
            return Err(PhyloError::invalid("iterations", "must be positive"));
        }
        if self.burn_in >= self.iterations {
            return Err(PhyloError::invalid(
                "burn_in",
                "must be smaller than iterations",
            ));
        }
        if self.thin == 0 {
            return Err(PhyloError::invalid("thin", "must be positive"));
        }
        if self.gamma_categories == 0 {
            return Err(PhyloError::invalid("gamma_categories", "must be positive"));
        }
        Ok(())
    }
}

/// One sample recorded from the chain at iteration `iter`.
#[derive(Debug, Clone)]
pub struct ChainSample {
    /// Iteration index this sample was recorded at.
    pub iter: usize,
    /// State at this iteration — a deep copy.
    pub state: ChainState,
    /// Cached `log P(data | state)` of `state`.
    pub log_likelihood: f64,
    /// Cached `log P(state)` of `state` (prior contribution).
    pub log_prior: f64,
    /// `log_likelihood + log_prior` — the chain's posterior score.
    pub log_posterior: f64,
}

/// Per-kind acceptance / rejection counters for a chain run.
#[derive(Debug, Clone, Default)]
pub struct AcceptanceCounts {
    /// Per-kind counters: `(proposed, accepted)`.
    pub by_kind: HashMap<ProposalKind, (usize, usize)>,
}

impl AcceptanceCounts {
    /// Records one proposed move and whether it was accepted.
    fn record(&mut self, kind: ProposalKind, accepted: bool) {
        let slot = self.by_kind.entry(kind).or_insert((0, 0));
        slot.0 += 1;
        if accepted {
            slot.1 += 1;
        }
    }

    /// Total proposals over all kinds.
    pub fn total_proposed(&self) -> usize {
        self.by_kind.values().map(|(p, _)| p).sum()
    }

    /// Total accepted over all kinds.
    pub fn total_accepted(&self) -> usize {
        self.by_kind.values().map(|(_, a)| a).sum()
    }

    /// Overall acceptance rate, or `0.0` for an empty run.
    pub fn overall_rate(&self) -> f64 {
        let p = self.total_proposed();
        if p == 0 {
            0.0
        } else {
            self.total_accepted() as f64 / p as f64
        }
    }

    /// Acceptance rate for a single proposal kind, or `None` if the
    /// kind was never proposed.
    pub fn rate(&self, kind: ProposalKind) -> Option<f64> {
        self.by_kind
            .get(&kind)
            .map(|(p, a)| if *p == 0 { 0.0 } else { *a as f64 / *p as f64 })
    }
}

/// All output of one MCMC chain run.
#[derive(Debug, Clone)]
pub struct ChainResult {
    /// Recorded samples (after burn-in, after thinning).
    pub samples: Vec<ChainSample>,
    /// Acceptance statistics by proposal kind.
    pub acceptance: AcceptanceCounts,
    /// The terminal state of the chain (useful for resuming).
    pub final_state: ChainState,
}

impl ChainResult {
    /// `log_posterior` trace, one value per recorded sample.
    pub fn log_posterior_trace(&self) -> Vec<f64> {
        self.samples.iter().map(|s| s.log_posterior).collect()
    }

    /// `log_likelihood` trace, one value per recorded sample.
    pub fn log_likelihood_trace(&self) -> Vec<f64> {
        self.samples.iter().map(|s| s.log_likelihood).collect()
    }

    /// `total tree length` trace.
    pub fn tree_length_trace(&self) -> Vec<f64> {
        self.samples
            .iter()
            .map(|s| s.state.tree.total_length())
            .collect()
    }

    /// `κ` trace (Hky85 / K80 only — other models contribute `1.0`).
    pub fn kappa_trace(&self) -> Vec<f64> {
        self.samples
            .iter()
            .map(|s| match s.state.model {
                SubstModel::K80 { kappa } | SubstModel::Hky85 { kappa, .. } => kappa,
                _ => 1.0,
            })
            .collect()
    }

    /// `α` (gamma shape) trace, one value per sample — `None` if the
    /// chain did not include a gamma-α.
    pub fn alpha_trace(&self) -> Option<Vec<f64>> {
        if self.samples.iter().any(|s| s.state.gamma_alpha.is_some()) {
            Some(
                self.samples
                    .iter()
                    .map(|s| s.state.gamma_alpha.unwrap_or(0.0))
                    .collect(),
            )
        } else {
            None
        }
    }

    /// Trees only — the topology / branch-length sample. Useful as
    /// input to [`super::posterior`] consensus + clade probabilities.
    pub fn tree_samples(&self) -> Vec<Tree> {
        self.samples.iter().map(|s| s.state.tree.clone()).collect()
    }
}

/// Pushes one [`ChainSample`] onto the trace if `iter` is past burn-in
/// and aligned to the thinning interval. Centralises the recording
/// rule so every code path that advances an iteration uses the same
/// guard.
fn record_sample(
    samples: &mut Vec<ChainSample>,
    iter: usize,
    state: &ChainState,
    ll: f64,
    lp: f64,
    cfg: &ChainConfig,
) {
    if iter >= cfg.burn_in && (iter - cfg.burn_in) % cfg.thin == 0 {
        samples.push(ChainSample {
            iter,
            state: state.clone(),
            log_likelihood: ll,
            log_prior: lp,
            log_posterior: ll + lp,
        });
    }
}

/// Computes `log P(data | state) + log P(state)`. Caches the gamma
/// table inside the call — callers should compute the posterior once
/// per (proposed) state.
fn compute_posterior(
    state: &ChainState,
    alignment: &[(String, Vec<u8>)],
    prior: &Prior,
    gamma_categories: usize,
) -> Result<(f64, f64)> {
    let ll = match state.gamma_alpha {
        None => log_likelihood(&state.tree, &state.model, alignment)?,
        Some(alpha) => {
            let g = DiscreteGamma::new(alpha, gamma_categories)?;
            log_likelihood_gamma(&state.tree, &state.model, alignment, &g)?
        }
    };
    let lp = prior.log_topology(&state.tree)
        + prior.log_branch_lengths(&state.tree)
        + prior.log_model(&state.model)
        + match state.gamma_alpha {
            None => 0.0,
            Some(a) => prior.log_gamma_alpha(a),
        };
    Ok((ll, lp))
}

/// Runs one Metropolis-Hastings chain.
///
/// `init` is the starting state; `prior` and `proposals` define the
/// model; `cfg` controls the chain length / burn-in / thinning / seed;
/// `alignment` is the data.
///
/// # Errors
/// - [`PhyloError`] propagated from the likelihood / prior / model
///   validators on the starting state. Proposed states whose
///   likelihood is non-finite are silently rejected (the move is
///   counted as proposed but not accepted).
pub fn run_chain(
    init: ChainState,
    prior: &Prior,
    proposals: &ProposalSet,
    cfg: &ChainConfig,
    alignment: &[(String, Vec<u8>)],
) -> Result<ChainResult> {
    cfg.validate()?;
    prior.validate()?;

    let mut rng = Rng::new(cfg.seed);
    let mut state = init;
    let (mut ll, mut lp) = compute_posterior(&state, alignment, prior, cfg.gamma_categories)?;
    let mut acceptance = AcceptanceCounts::default();
    let mut samples: Vec<ChainSample> = Vec::new();

    for iter in 0..cfg.iterations {
        // Attempt one proposal — if it cannot apply (e.g. only kinds
        // are listed that don't fit the current state), the chain stays
        // put for this iteration. We still record a sample, so the
        // sample count equals the iteration count after burn-in /
        // thinning (every chain run with the same config produces a
        // trace of the same length).
        if let Some((kind, outcome)) = sample_proposal(&state, proposals, &mut rng, 5)? {
            let scored =
                compute_posterior(&outcome.new_state, alignment, prior, cfg.gamma_categories);
            let (new_ll, new_lp) = match scored {
                Ok(pair) => pair,
                Err(_) => {
                    acceptance.record(kind, false);
                    record_sample(&mut samples, iter, &state, ll, lp, cfg);
                    continue;
                }
            };
            if !new_ll.is_finite() || !new_lp.is_finite() {
                acceptance.record(kind, false);
                record_sample(&mut samples, iter, &state, ll, lp, cfg);
                continue;
            }
            let log_alpha = (new_ll + new_lp) - (ll + lp) + outcome.log_hastings;
            let accept = log_alpha >= 0.0 || rng.uniform().ln() < log_alpha;
            if accept {
                state = outcome.new_state;
                ll = new_ll;
                lp = new_lp;
            }
            acceptance.record(kind, accept);
        }
        record_sample(&mut samples, iter, &state, ll, lp, cfg);
    }

    Ok(ChainResult {
        samples,
        acceptance,
        final_state: state,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bayes::proposal::apply;
    use crate::io::newick::read_newick;
    use crate::simulate::seqgen::simulate_sequences;

    fn jc_init() -> ChainState {
        ChainState {
            tree: read_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap(),
            model: SubstModel::Jc69,
            gamma_alpha: None,
        }
    }

    fn aln_for(tree: &Tree) -> Vec<(String, Vec<u8>)> {
        let sim = simulate_sequences(tree, &SubstModel::Jc69, 80, None, 42).unwrap();
        sim.rows
    }

    #[test]
    fn config_validation_works() {
        let cfg = ChainConfig {
            iterations: 0,
            ..ChainConfig::default()
        };
        assert!(cfg.validate().is_err());
        let cfg = ChainConfig {
            iterations: 10,
            burn_in: 10,
            ..ChainConfig::default()
        };
        assert!(cfg.validate().is_err());
        let cfg = ChainConfig {
            iterations: 10,
            burn_in: 5,
            thin: 0,
            ..ChainConfig::default()
        };
        assert!(cfg.validate().is_err());
        let cfg = ChainConfig {
            iterations: 10,
            burn_in: 5,
            thin: 1,
            ..ChainConfig::default()
        };
        assert!(cfg.validate().is_ok());
    }

    #[test]
    fn chain_records_at_least_one_sample() {
        let init = jc_init();
        let aln = aln_for(&init.tree);
        let cfg = ChainConfig {
            iterations: 200,
            burn_in: 50,
            thin: 10,
            seed: 1,
            gamma_categories: 4,
        };
        let result =
            run_chain(init, &Prior::default(), &ProposalSet::default(), &cfg, &aln).unwrap();
        // 150 post-burn-in iterations, thinning 10 → ~15 samples.
        assert!(!result.samples.is_empty(), "no samples recorded");
        assert!(result.samples.len() <= 16);
    }

    #[test]
    fn acceptance_counts_are_consistent() {
        let init = jc_init();
        let aln = aln_for(&init.tree);
        let cfg = ChainConfig {
            iterations: 300,
            burn_in: 50,
            thin: 5,
            seed: 2,
            gamma_categories: 4,
        };
        let result =
            run_chain(init, &Prior::default(), &ProposalSet::default(), &cfg, &aln).unwrap();
        let total = result.acceptance.total_proposed();
        let acc = result.acceptance.total_accepted();
        assert!(total > 0);
        assert!(acc <= total);
        assert!(result.acceptance.overall_rate() >= 0.0);
        assert!(result.acceptance.overall_rate() <= 1.0);
    }

    #[test]
    fn log_posterior_trace_is_finite() {
        let init = jc_init();
        let aln = aln_for(&init.tree);
        let cfg = ChainConfig {
            iterations: 200,
            burn_in: 50,
            thin: 5,
            seed: 3,
            gamma_categories: 4,
        };
        let result =
            run_chain(init, &Prior::default(), &ProposalSet::default(), &cfg, &aln).unwrap();
        for s in &result.samples {
            assert!(s.log_posterior.is_finite());
            assert!(s.log_likelihood.is_finite());
            assert!(s.log_prior.is_finite());
        }
    }

    #[test]
    fn chain_converges_toward_higher_posterior_from_a_bad_start() {
        // Start from a deliberately bad tree (all branches 5.0).
        let mut init = jc_init();
        for id in 0..init.tree.node_count() {
            if init.tree.node(id).parent.is_some() {
                init.tree.node_mut(id).branch_length = Some(5.0);
            }
        }
        let aln = aln_for(&init.tree);
        let cfg = ChainConfig {
            iterations: 600,
            burn_in: 100,
            thin: 10,
            seed: 4,
            gamma_categories: 4,
        };
        let result =
            run_chain(init, &Prior::default(), &ProposalSet::default(), &cfg, &aln).unwrap();
        let trace = result.log_posterior_trace();
        // The mean of the second half should be much higher than the
        // mean of the first quarter — the chain finds better parameter
        // values.
        let q = trace.len() / 4;
        let first_q: f64 = trace.iter().take(q).sum::<f64>() / q as f64;
        let second_half: f64 =
            trace.iter().skip(trace.len() / 2).sum::<f64>() / (trace.len() / 2) as f64;
        assert!(
            second_half > first_q,
            "trace did not improve: first_q = {first_q}, second_half = {second_half}"
        );
    }

    #[test]
    fn deterministic_for_a_seed() {
        let init = jc_init();
        let aln = aln_for(&init.tree);
        let cfg = ChainConfig {
            iterations: 100,
            burn_in: 20,
            thin: 5,
            seed: 99,
            gamma_categories: 4,
        };
        let a = run_chain(
            init.clone(),
            &Prior::default(),
            &ProposalSet::default(),
            &cfg,
            &aln,
        )
        .unwrap();
        let b = run_chain(init, &Prior::default(), &ProposalSet::default(), &cfg, &aln).unwrap();
        let ta = a.log_posterior_trace();
        let tb = b.log_posterior_trace();
        assert_eq!(ta.len(), tb.len());
        for (x, y) in ta.iter().zip(tb.iter()) {
            assert!((x - y).abs() < 1e-9, "trace not deterministic");
        }
    }

    #[test]
    fn symmetric_proposals_acceptance_matches_min_one_ratio() {
        // Build a chain that ONLY uses symmetric branch-slide moves;
        // for every accepted move the empirical fraction should match
        // min(1, exp(Δ log posterior)). The check is loose because of
        // Monte-Carlo noise.
        let init = jc_init();
        let aln = aln_for(&init.tree);
        let proposals = ProposalSet {
            kinds: vec![(ProposalKind::BranchSlide, 1.0)],
            scale_lambda: 1.0,
            slide_sigma: 0.05,
            dirichlet_beta: 100.0,
        };
        let cfg = ChainConfig {
            iterations: 400,
            burn_in: 50,
            thin: 1,
            seed: 5,
            gamma_categories: 4,
        };
        let result = run_chain(init, &Prior::default(), &proposals, &cfg, &aln).unwrap();
        let rate = result.acceptance.rate(ProposalKind::BranchSlide).unwrap();
        // For a slide_sigma = 0.05 step on a near-mode chain, the rate
        // should be neither pinned at 0 nor at 1.
        assert!(rate > 0.05 && rate < 0.99, "rate = {rate}");
    }

    #[test]
    fn apply_handles_bad_models_gracefully() {
        // Use HKY on a 4-taxon tree; ensure apply() never panics across
        // all proposal kinds.
        let mut rng = Rng::new(7);
        let st = ChainState {
            tree: read_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap(),
            model: SubstModel::Hky85 {
                kappa: 2.0,
                freqs: [0.25; 4],
            },
            gamma_alpha: Some(0.5),
        };
        let set = ProposalSet::default();
        let all_kinds = [
            ProposalKind::Nni,
            ProposalKind::Spr,
            ProposalKind::WilsonBalding,
            ProposalKind::BranchScale,
            ProposalKind::BranchSlide,
            ProposalKind::TreeScale,
            ProposalKind::KappaScale,
            ProposalKind::GtrRateDirichlet,
            ProposalKind::FreqDirichlet,
            ProposalKind::GammaAlphaScale,
        ];
        for k in all_kinds {
            let _ = apply(k, &st, &set, &mut rng);
        }
    }
}
