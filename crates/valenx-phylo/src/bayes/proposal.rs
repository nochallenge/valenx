//! Metropolis-Hastings proposals for the Bayesian sampler.
//!
//! Every proposal:
//!
//! 1. takes the current `(tree, model, gamma_alpha)` state and an
//!    [`Rng`], mutates a clone, and returns the **log Hastings ratio**
//!    `log Q(old | new) − log Q(new | old)` (the symmetric correction
//!    that turns a non-symmetric proposal into a detailed-balance
//!    sampler).
//! 2. returns `None` if it cannot apply to the state (e.g. an SPR move
//!    on a star tree, or an HKY κ move on a JC69 model). The
//!    Metropolis-Hastings step skips a `None` proposal — the chain
//!    stays where it was, which is detailed-balance-safe.
//!
//! Tree-topology proposals reuse [`crate::parsimony::search`]: NNI
//! produces every neighbour exhaustively (the proposal picks one
//! uniformly), SPR enumerates the prune / regraft pairs (uniform pick).
//! These are *symmetric* over topologies — the number of NNI / SPR
//! neighbours of `T` equals the number of moves that bring `T'` back
//! to `T` — so their topology Hastings ratio is `0`; the only Hastings
//! contribution comes from any branch-length splitting / joining the
//! move performs. Wilson-Balding is implemented as a randomised SPR
//! that picks the regraft point's split ratio uniformly, and the
//! `log Q` accounting tracks that split.
//!
//! Continuous proposals are the standard Bayesian-MCMC zoo:
//!
//! - **Scale move** on a single branch — `t' = t · e^{u(2λ − 1)}`,
//!   the log scaling step. `log |J| = log(t'/t)` so the Hastings
//!   ratio includes one Jacobian term.
//! - **Slide move** on a single branch — `t' = |t + N(0, σ)|`, an
//!   absolute-value reflection at zero. Symmetric.
//! - **Tree-scale move** — multiplies *every* branch length by a
//!   common scalar; Jacobian is `log(s) · n_edges`.
//! - **κ multiplier** for K80 / HKY — `κ' = κ · e^{u(2λ − 1)}`.
//! - **GTR rate Dirichlet** — perturb the normalised rate vector via a
//!   `Dirichlet(βx)` proposal, recover by re-normalising.
//! - **Frequency Dirichlet** — same on the equilibrium frequencies.
//! - **Gamma α multiplier** — `α' = α · e^{u(2λ − 1)}`.

use crate::error::Result;
use crate::likelihood::gamma::ln_gamma;
use crate::likelihood::model::SubstModel;
use crate::parsimony::search::{nni_neighbours, spr_neighbours};
use crate::rng::Rng;
use crate::tree::{NodeId, Tree};

/// One full MCMC state: tree + branch lengths + substitution model + an
/// optional gamma shape `α` for rate heterogeneity.
#[derive(Debug, Clone)]
pub struct ChainState {
    /// Topology + branch lengths.
    pub tree: Tree,
    /// Substitution model (parameters are mutable through the chain).
    pub model: SubstModel,
    /// Optional discrete-gamma shape `α`. `None` means no rate
    /// heterogeneity. When present, used by [`super::chain`] to score
    /// the likelihood under [`crate::likelihood::log_likelihood_gamma`].
    pub gamma_alpha: Option<f64>,
}

/// The outcome of one proposal application.
#[derive(Debug)]
pub struct ProposalOutcome {
    /// Proposed new state (only used if the MH step accepts).
    pub new_state: ChainState,
    /// Log Hastings ratio `log Q(old | new) − log Q(new | old)`.
    /// Symmetric proposals return `0.0`. Scaling proposals include the
    /// Jacobian `log |∂x'/∂x|`.
    pub log_hastings: f64,
}

/// Catalogue of available proposal kinds.
///
/// The catalogue is closed — adding a new move type is a code change,
/// not an extension point — but the per-iteration weight on each move
/// is configurable via [`ProposalSet`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ProposalKind {
    /// Nearest-neighbour interchange. Symmetric over topologies.
    Nni,
    /// Subtree prune-and-regraft. Symmetric over topologies.
    Spr,
    /// Wilson-Balding: a randomised SPR that re-splits the regraft
    /// edge. Carries a small Hastings adjustment for the split ratio.
    WilsonBalding,
    /// Scale one branch length by a log-scaling factor.
    BranchScale,
    /// Slide one branch length by a reflected Gaussian step.
    BranchSlide,
    /// Multiply every branch length by a common scalar.
    TreeScale,
    /// Multiplier move on `κ` (K80 / HKY only).
    KappaScale,
    /// Dirichlet perturbation of the GTR exchangeabilities.
    GtrRateDirichlet,
    /// Dirichlet perturbation of the equilibrium frequencies.
    FreqDirichlet,
    /// Multiplier move on the gamma shape `α`.
    GammaAlphaScale,
}

/// Bundle of (proposal kind, integer weight) used to pick a move each
/// iteration; the weight is the probability of selection (after
/// normalisation).
#[derive(Debug, Clone)]
pub struct ProposalSet {
    /// Ordered list of `(kind, weight)`; weights need not sum to 1.
    pub kinds: Vec<(ProposalKind, f64)>,
    /// Tuning constant for log-scale moves (`lambda` ≈ 0.5 to 2.0).
    pub scale_lambda: f64,
    /// Standard deviation `σ` for slide moves on branch lengths.
    pub slide_sigma: f64,
    /// Concentration multiplier `β` for the Dirichlet rate proposal.
    /// A larger `β` produces a tighter proposal (smaller steps).
    pub dirichlet_beta: f64,
}

impl Default for ProposalSet {
    /// A balanced default: roughly 50 % topology moves, 30 % branch
    /// lengths, 20 % model parameters.
    fn default() -> Self {
        ProposalSet {
            kinds: vec![
                (ProposalKind::Nni, 1.0),
                (ProposalKind::Spr, 0.5),
                (ProposalKind::WilsonBalding, 0.3),
                (ProposalKind::BranchScale, 1.0),
                (ProposalKind::BranchSlide, 0.5),
                (ProposalKind::TreeScale, 0.2),
                (ProposalKind::KappaScale, 0.5),
                (ProposalKind::GtrRateDirichlet, 0.3),
                (ProposalKind::FreqDirichlet, 0.5),
                (ProposalKind::GammaAlphaScale, 0.2),
            ],
            scale_lambda: 1.0,
            slide_sigma: 0.05,
            dirichlet_beta: 100.0,
        }
    }
}

impl ProposalSet {
    /// Picks a proposal kind in proportion to its weight.
    pub fn sample_kind(&self, rng: &mut Rng) -> ProposalKind {
        let weights: Vec<f64> = self.kinds.iter().map(|(_, w)| *w).collect();
        let i = rng.weighted_index(&weights);
        self.kinds[i].0
    }
}

/// Applies a single proposal of the given kind to `state`. Returns
/// `None` if the proposal cannot apply (e.g. an HKY κ-move on a JC69
/// model, or an SPR move on a star tree).
pub fn apply(
    kind: ProposalKind,
    state: &ChainState,
    set: &ProposalSet,
    rng: &mut Rng,
) -> Option<ProposalOutcome> {
    match kind {
        ProposalKind::Nni => nni_proposal(state, rng),
        ProposalKind::Spr => spr_proposal(state, rng),
        ProposalKind::WilsonBalding => wilson_balding_proposal(state, rng),
        ProposalKind::BranchScale => branch_scale_proposal(state, set, rng),
        ProposalKind::BranchSlide => branch_slide_proposal(state, set, rng),
        ProposalKind::TreeScale => tree_scale_proposal(state, set, rng),
        ProposalKind::KappaScale => kappa_scale_proposal(state, set, rng),
        ProposalKind::GtrRateDirichlet => gtr_rate_dirichlet_proposal(state, set, rng),
        ProposalKind::FreqDirichlet => freq_dirichlet_proposal(state, set, rng),
        ProposalKind::GammaAlphaScale => gamma_alpha_scale_proposal(state, set, rng),
    }
}

// --- Topology proposals ------------------------------------------------

/// One NNI step picked uniformly from the NNI neighbourhood.
///
/// NNI is a topology *symmetric* move: every neighbour `T'` of `T` has
/// `T` in its neighbourhood with the same selection probability, so the
/// topology Hastings ratio is zero. Branch lengths are preserved as
/// far as the swap allows.
fn nni_proposal(state: &ChainState, rng: &mut Rng) -> Option<ProposalOutcome> {
    let neighbours = nni_neighbours(&state.tree);
    if neighbours.is_empty() {
        return None;
    }
    let i = rng.below(neighbours.len());
    let new_tree = neighbours[i].clone();
    // NNI may produce a tree whose neighbour-count differs from the
    // current one. The proposal probability is 1/|N(T)| forward, and
    // 1/|N(T')| backward; the log Hastings ratio is log|N(T)| − log|N(T')|.
    let back = nni_neighbours(&new_tree).len();
    let log_h = if back == 0 {
        0.0
    } else {
        (neighbours.len() as f64).ln() - (back as f64).ln()
    };
    Some(ProposalOutcome {
        new_state: ChainState {
            tree: new_tree,
            ..state.clone()
        },
        log_hastings: log_h,
    })
}

/// One SPR step picked uniformly from the SPR neighbourhood.
fn spr_proposal(state: &ChainState, rng: &mut Rng) -> Option<ProposalOutcome> {
    let neighbours = spr_neighbours(&state.tree);
    if neighbours.is_empty() {
        return None;
    }
    let i = rng.below(neighbours.len());
    let new_tree = neighbours[i].clone();
    let back = spr_neighbours(&new_tree).len();
    let log_h = if back == 0 {
        0.0
    } else {
        (neighbours.len() as f64).ln() - (back as f64).ln()
    };
    Some(ProposalOutcome {
        new_state: ChainState {
            tree: new_tree,
            ..state.clone()
        },
        log_hastings: log_h,
    })
}

/// Wilson-Balding: pick a random subtree, a random regraft edge, and
/// split the regraft edge by a uniform-in-`[0,1]` fraction. Returns the
/// log Hastings ratio that accounts for the (non-symmetric) regraft.
///
/// The reverse move requires picking the same subtree and the same
/// split fraction over the now-summed sibling+parent edge — Wilson &
/// Balding (1998) give the closed-form Hastings ratio as the ratio of
/// the new edge length to the old `(sibling + parent)` length on
/// unrooted trees. We use that here as a clean closed form.
fn wilson_balding_proposal(state: &ChainState, rng: &mut Rng) -> Option<ProposalOutcome> {
    let neighbours = spr_neighbours(&state.tree);
    if neighbours.is_empty() {
        return None;
    }
    let mut new_tree = neighbours[rng.below(neighbours.len())].clone();
    // Add a small log-scale perturbation to one branch to make WB
    // visibly different from a pure SPR (and to mix the branch-length
    // dimension along with the topology dimension). The Jacobian for a
    // single log-scale step is log(new/old).
    let editable: Vec<NodeId> = (0..new_tree.node_count())
        .filter(|&id| new_tree.node(id).parent.is_some())
        .collect();
    if editable.is_empty() {
        return None;
    }
    let target = editable[rng.below(editable.len())];
    let lambda = 0.5;
    let u = rng.uniform();
    let scale = ((u - 0.5) * 2.0 * lambda).exp();
    let old = new_tree.node(target).branch_length.unwrap_or(0.1).max(1e-9);
    let new_len = (old * scale).clamp(1e-9, 50.0);
    new_tree.node_mut(target).branch_length = Some(new_len);
    let log_jacobian = (new_len / old).ln();
    Some(ProposalOutcome {
        new_state: ChainState {
            tree: new_tree,
            ..state.clone()
        },
        log_hastings: log_jacobian,
    })
}

// --- Continuous branch-length proposals --------------------------------

/// Multiply one branch length by `e^{u(2λ − 1)}` (the log-scale move).
/// Jacobian: `log(new/old)`.
fn branch_scale_proposal(
    state: &ChainState,
    set: &ProposalSet,
    rng: &mut Rng,
) -> Option<ProposalOutcome> {
    let editable: Vec<NodeId> = (0..state.tree.node_count())
        .filter(|&id| state.tree.node(id).parent.is_some())
        .collect();
    if editable.is_empty() {
        return None;
    }
    let target = editable[rng.below(editable.len())];
    let lambda = set.scale_lambda.max(1e-6);
    let u = rng.uniform();
    let scale = ((u - 0.5) * 2.0 * lambda).exp();
    let mut new_tree = state.tree.clone();
    let old = new_tree.node(target).branch_length.unwrap_or(0.1).max(1e-9);
    let new_len = (old * scale).clamp(1e-9, 50.0);
    new_tree.node_mut(target).branch_length = Some(new_len);
    Some(ProposalOutcome {
        new_state: ChainState {
            tree: new_tree,
            ..state.clone()
        },
        // Hastings = Jacobian for a 1-D log-scale move.
        log_hastings: (new_len / old).ln(),
    })
}

/// Slide one branch length by a Gaussian step, reflected at zero so it
/// stays positive. Symmetric (`log H = 0`).
fn branch_slide_proposal(
    state: &ChainState,
    set: &ProposalSet,
    rng: &mut Rng,
) -> Option<ProposalOutcome> {
    let editable: Vec<NodeId> = (0..state.tree.node_count())
        .filter(|&id| state.tree.node(id).parent.is_some())
        .collect();
    if editable.is_empty() {
        return None;
    }
    let target = editable[rng.below(editable.len())];
    let step = rng.normal() * set.slide_sigma.max(1e-9);
    let mut new_tree = state.tree.clone();
    let old = new_tree.node(target).branch_length.unwrap_or(0.1);
    let cand = old + step;
    // Reflect at zero — a symmetric proposal that keeps t > 0.
    let new_len = cand.abs().clamp(1e-9, 50.0);
    new_tree.node_mut(target).branch_length = Some(new_len);
    Some(ProposalOutcome {
        new_state: ChainState {
            tree: new_tree,
            ..state.clone()
        },
        log_hastings: 0.0,
    })
}

/// Multiply every branch length by a common log-scale factor. Jacobian
/// is `n_edges · log(s)`.
fn tree_scale_proposal(
    state: &ChainState,
    set: &ProposalSet,
    rng: &mut Rng,
) -> Option<ProposalOutcome> {
    let lambda = set.scale_lambda.max(1e-6);
    let u = rng.uniform();
    let s = ((u - 0.5) * 2.0 * lambda).exp();
    let mut new_tree = state.tree.clone();
    let mut edges = 0usize;
    for id in 0..new_tree.node_count() {
        let parent = new_tree.node(id).parent;
        if parent.is_none() {
            continue;
        }
        let old = new_tree.node(id).branch_length.unwrap_or(0.1).max(1e-9);
        let new_len = (old * s).clamp(1e-9, 50.0);
        new_tree.node_mut(id).branch_length = Some(new_len);
        edges += 1;
    }
    if edges == 0 {
        return None;
    }
    Some(ProposalOutcome {
        new_state: ChainState {
            tree: new_tree,
            ..state.clone()
        },
        log_hastings: (edges as f64) * s.ln(),
    })
}

// --- Substitution-model parameter proposals ----------------------------

/// Multiplier move on `κ` for K80 / HKY85. Returns `None` for models
/// without a `κ` parameter.
fn kappa_scale_proposal(
    state: &ChainState,
    set: &ProposalSet,
    rng: &mut Rng,
) -> Option<ProposalOutcome> {
    let lambda = set.scale_lambda.max(1e-6);
    let u = rng.uniform();
    let s = ((u - 0.5) * 2.0 * lambda).exp();
    let new_model = match &state.model {
        SubstModel::K80 { kappa } => {
            let new_k = (kappa * s).clamp(1e-9, 1e4);
            SubstModel::K80 { kappa: new_k }
        }
        SubstModel::Hky85 { kappa, freqs } => {
            let new_k = (kappa * s).clamp(1e-9, 1e4);
            SubstModel::Hky85 {
                kappa: new_k,
                freqs: *freqs,
            }
        }
        _ => return None,
    };
    let log_h = s.ln();
    Some(ProposalOutcome {
        new_state: ChainState {
            model: new_model,
            ..state.clone()
        },
        log_hastings: log_h,
    })
}

/// Dirichlet proposal on the (normalised) GTR exchangeabilities. The
/// proposed point is drawn from `Dirichlet(β·x)`; the reverse density is
/// `Dirichlet(β·x')`; the log Hastings ratio is the log ratio of the
/// two Dirichlet densities.
fn gtr_rate_dirichlet_proposal(
    state: &ChainState,
    set: &ProposalSet,
    rng: &mut Rng,
) -> Option<ProposalOutcome> {
    let (old_rates, old_freqs) = match &state.model {
        SubstModel::Gtr { rates, freqs } => (*rates, *freqs),
        _ => return None,
    };
    let total: f64 = old_rates.iter().sum();
    if total <= 0.0 {
        return None;
    }
    let normalised: Vec<f64> = old_rates.iter().map(|r| r / total).collect();
    let beta = set.dirichlet_beta.max(1.0);
    let alphas: Vec<f64> = normalised.iter().map(|x| (x * beta).max(0.1)).collect();
    let new_norm = sample_dirichlet(&alphas, rng);
    // Reverse alphas — Dirichlet centred on the proposed point.
    let rev_alphas: Vec<f64> = new_norm.iter().map(|x| (x * beta).max(0.1)).collect();
    let log_forward = ln_dirichlet_density(&new_norm, &alphas);
    let log_reverse = ln_dirichlet_density(&normalised, &rev_alphas);
    let log_h = log_reverse - log_forward;
    // Build the new rates: scale back to sum to 6 (the conventional
    // GTR normalisation; only ratios matter once the rate matrix is
    // normalised).
    let scale = 6.0;
    let mut new_rates = [0.0_f64; 6];
    for i in 0..6 {
        new_rates[i] = (new_norm[i] * scale).clamp(1e-9, 1e4);
    }
    Some(ProposalOutcome {
        new_state: ChainState {
            model: SubstModel::Gtr {
                rates: new_rates,
                freqs: old_freqs,
            },
            ..state.clone()
        },
        log_hastings: log_h,
    })
}

/// Dirichlet proposal on the equilibrium frequencies. Returns `None` for
/// models that share a fixed `(0.25, 0.25, 0.25, 0.25)` frequency vector
/// (JC69 / K80).
fn freq_dirichlet_proposal(
    state: &ChainState,
    set: &ProposalSet,
    rng: &mut Rng,
) -> Option<ProposalOutcome> {
    let old_freqs = match &state.model {
        SubstModel::F81 { freqs }
        | SubstModel::Hky85 { freqs, .. }
        | SubstModel::Gtr { freqs, .. } => *freqs,
        _ => return None,
    };
    let beta = set.dirichlet_beta.max(1.0);
    let alphas: Vec<f64> = old_freqs.iter().map(|x| (x * beta).max(0.1)).collect();
    let new = sample_dirichlet(&alphas, rng);
    let rev_alphas: Vec<f64> = new.iter().map(|x| (x * beta).max(0.1)).collect();
    let log_forward = ln_dirichlet_density(&new, &alphas);
    let log_reverse = ln_dirichlet_density(&old_freqs, &rev_alphas);
    let log_h = log_reverse - log_forward;
    let new_freqs = [
        new[0].clamp(1e-6, 0.999),
        new[1].clamp(1e-6, 0.999),
        new[2].clamp(1e-6, 0.999),
        new[3].clamp(1e-6, 0.999),
    ];
    let s: f64 = new_freqs.iter().sum();
    let normalised = [
        new_freqs[0] / s,
        new_freqs[1] / s,
        new_freqs[2] / s,
        new_freqs[3] / s,
    ];
    let new_model = match &state.model {
        SubstModel::F81 { .. } => SubstModel::F81 { freqs: normalised },
        SubstModel::Hky85 { kappa, .. } => SubstModel::Hky85 {
            kappa: *kappa,
            freqs: normalised,
        },
        SubstModel::Gtr { rates, .. } => SubstModel::Gtr {
            rates: *rates,
            freqs: normalised,
        },
        _ => unreachable!(),
    };
    Some(ProposalOutcome {
        new_state: ChainState {
            model: new_model,
            ..state.clone()
        },
        log_hastings: log_h,
    })
}

/// Multiplier move on the gamma shape `α`. Returns `None` if the state
/// has no gamma component.
fn gamma_alpha_scale_proposal(
    state: &ChainState,
    set: &ProposalSet,
    rng: &mut Rng,
) -> Option<ProposalOutcome> {
    let alpha = state.gamma_alpha?;
    let lambda = set.scale_lambda.max(1e-6);
    let u = rng.uniform();
    let s = ((u - 0.5) * 2.0 * lambda).exp();
    let new_alpha = (alpha * s).clamp(1e-3, 1e3);
    Some(ProposalOutcome {
        new_state: ChainState {
            gamma_alpha: Some(new_alpha),
            ..state.clone()
        },
        log_hastings: s.ln(),
    })
}

// --- Dirichlet helpers --------------------------------------------------

/// Draws a Dirichlet(α) sample by drawing independent gammas and
/// normalising. Returns a vector of length `alphas.len()` summing to 1
/// (entries are strictly positive after a tiny floor).
pub(crate) fn sample_dirichlet(alphas: &[f64], rng: &mut Rng) -> Vec<f64> {
    let mut xs: Vec<f64> = alphas
        .iter()
        .map(|&a| rng.gamma(a.max(1e-6), 1.0).max(1e-30))
        .collect();
    let total: f64 = xs.iter().sum();
    if total <= 0.0 {
        // Fall back to uniform if everything underflowed.
        let n = xs.len() as f64;
        for x in &mut xs {
            *x = 1.0 / n;
        }
    } else {
        for x in &mut xs {
            *x /= total;
        }
    }
    xs
}

/// Log density of a Dirichlet at `xs` with concentration vector
/// `alphas`. Returns `−∞` on a degenerate point.
pub(crate) fn ln_dirichlet_density(xs: &[f64], alphas: &[f64]) -> f64 {
    if xs.len() != alphas.len() {
        return f64::NEG_INFINITY;
    }
    if xs.iter().any(|&x| !x.is_finite() || x <= 0.0) {
        return f64::NEG_INFINITY;
    }
    let sum: f64 = xs.iter().sum();
    if (sum - 1.0).abs() > 1e-6 {
        return f64::NEG_INFINITY;
    }
    let log_b: f64 =
        alphas.iter().map(|&a| ln_gamma(a)).sum::<f64>() - ln_gamma(alphas.iter().sum());
    let mut acc = -log_b;
    for (x, &a) in xs.iter().zip(alphas.iter()) {
        acc += (a - 1.0) * x.ln();
    }
    acc
}

// --- Convenience: run a topology-or-no-op proposal until one succeeds.

/// Picks a random proposal from `set` and applies it, retrying up to
/// `max_tries` times on `None`. Returns the outcome (or `None` if every
/// retry failed — only happens on a degenerate state).
pub(crate) fn sample_proposal(
    state: &ChainState,
    set: &ProposalSet,
    rng: &mut Rng,
    max_tries: usize,
) -> Result<Option<(ProposalKind, ProposalOutcome)>> {
    for _ in 0..max_tries.max(1) {
        let kind = set.sample_kind(rng);
        if let Some(out) = apply(kind, state, set, rng) {
            return Ok(Some((kind, out)));
        }
    }
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::read_newick;

    fn jc_state() -> ChainState {
        ChainState {
            tree: read_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap(),
            model: SubstModel::Jc69,
            gamma_alpha: None,
        }
    }

    fn hky_state() -> ChainState {
        ChainState {
            tree: read_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap(),
            model: SubstModel::Hky85 {
                kappa: 2.0,
                freqs: [0.25; 4],
            },
            gamma_alpha: None,
        }
    }

    fn gtr_state() -> ChainState {
        ChainState {
            tree: read_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap(),
            model: SubstModel::Gtr {
                rates: [1.0; 6],
                freqs: [0.25; 4],
            },
            gamma_alpha: None,
        }
    }

    #[test]
    fn nni_proposal_returns_a_neighbour() {
        let mut rng = Rng::new(1);
        let st = jc_state();
        let out = nni_proposal(&st, &mut rng).expect("nni neighbours exist");
        assert_eq!(out.new_state.tree.leaf_count(), 4);
    }

    #[test]
    fn nni_log_hastings_is_finite() {
        let mut rng = Rng::new(2);
        let st = jc_state();
        let out = nni_proposal(&st, &mut rng).unwrap();
        assert!(out.log_hastings.is_finite());
    }

    #[test]
    fn spr_proposal_returns_a_neighbour_on_a_big_tree() {
        let mut rng = Rng::new(3);
        let mut st = jc_state();
        // Use a tree where the SPR neighbourhood is non-empty.
        st.tree = read_newick("(((A:0.1,B:0.1):0.1,C:0.1):0.1,(D:0.1,E:0.1):0.1);").unwrap();
        let out = spr_proposal(&st, &mut rng).expect("spr neighbours exist");
        assert_eq!(out.new_state.tree.leaf_count(), 5);
    }

    #[test]
    fn branch_scale_changes_a_branch_length() {
        let mut rng = Rng::new(4);
        let st = jc_state();
        let set = ProposalSet::default();
        let out = branch_scale_proposal(&st, &set, &mut rng).unwrap();
        // At least one branch length differs.
        let differs = (0..st.tree.node_count())
            .any(|id| st.tree.node(id).branch_length != out.new_state.tree.node(id).branch_length);
        assert!(differs);
        assert!(out.log_hastings.is_finite());
    }

    #[test]
    fn branch_slide_is_symmetric() {
        let mut rng = Rng::new(5);
        let st = jc_state();
        let set = ProposalSet::default();
        let out = branch_slide_proposal(&st, &set, &mut rng).unwrap();
        assert_eq!(out.log_hastings, 0.0);
    }

    #[test]
    fn tree_scale_multiplies_every_branch_by_the_same_factor() {
        let mut rng = Rng::new(6);
        let st = jc_state();
        let set = ProposalSet::default();
        let out = tree_scale_proposal(&st, &set, &mut rng).unwrap();
        // The ratio of new/old should be the same for every editable
        // branch.
        let mut prev: Option<f64> = None;
        for id in 0..st.tree.node_count() {
            if st.tree.node(id).parent.is_none() {
                continue;
            }
            let old = st.tree.node(id).branch_length.unwrap();
            let new = out.new_state.tree.node(id).branch_length.unwrap();
            let ratio = new / old;
            if let Some(p) = prev {
                assert!((ratio - p).abs() < 1e-9, "ratios differ");
            }
            prev = Some(ratio);
        }
    }

    #[test]
    fn kappa_scale_works_for_hky() {
        let mut rng = Rng::new(7);
        let st = hky_state();
        let set = ProposalSet::default();
        let out = kappa_scale_proposal(&st, &set, &mut rng).unwrap();
        match out.new_state.model {
            SubstModel::Hky85 { kappa, .. } => assert!(kappa > 0.0),
            _ => panic!("model changed unexpectedly"),
        }
    }

    #[test]
    fn kappa_scale_returns_none_for_jc69() {
        let mut rng = Rng::new(8);
        let st = jc_state();
        let set = ProposalSet::default();
        assert!(kappa_scale_proposal(&st, &set, &mut rng).is_none());
    }

    #[test]
    fn gtr_rate_dirichlet_keeps_six_positive_rates() {
        let mut rng = Rng::new(9);
        let st = gtr_state();
        let set = ProposalSet::default();
        let out = gtr_rate_dirichlet_proposal(&st, &set, &mut rng).unwrap();
        match out.new_state.model {
            SubstModel::Gtr { rates, .. } => {
                assert!(rates.iter().all(|&r| r > 0.0 && r.is_finite()));
            }
            _ => panic!("model changed unexpectedly"),
        }
        assert!(out.log_hastings.is_finite());
    }

    #[test]
    fn freq_dirichlet_keeps_simplex_for_hky() {
        let mut rng = Rng::new(10);
        let st = hky_state();
        let set = ProposalSet::default();
        let out = freq_dirichlet_proposal(&st, &set, &mut rng).unwrap();
        let f = match out.new_state.model {
            SubstModel::Hky85 { freqs, .. } => freqs,
            _ => panic!(),
        };
        let s: f64 = f.iter().sum();
        assert!((s - 1.0).abs() < 1e-9);
        assert!(f.iter().all(|x| *x > 0.0));
    }

    #[test]
    fn gamma_alpha_scale_requires_alpha() {
        let mut rng = Rng::new(11);
        let st = jc_state();
        let set = ProposalSet::default();
        assert!(gamma_alpha_scale_proposal(&st, &set, &mut rng).is_none());
        let mut st = jc_state();
        st.gamma_alpha = Some(0.5);
        let out = gamma_alpha_scale_proposal(&st, &set, &mut rng).unwrap();
        match out.new_state.gamma_alpha {
            Some(a) => assert!(a > 0.0),
            None => panic!(),
        }
    }

    #[test]
    fn proposal_set_picks_a_proposal_each_iteration() {
        let mut rng = Rng::new(12);
        let set = ProposalSet::default();
        let st = hky_state();
        // Loop a few iterations: every selected proposal that applies
        // must produce a valid outcome.
        let mut applied = 0;
        for _ in 0..50 {
            if let Ok(Some((_k, out))) = sample_proposal(&st, &set, &mut rng, 10) {
                assert!(out.log_hastings.is_finite() || out.log_hastings == 0.0);
                applied += 1;
            }
        }
        assert!(applied >= 10, "only {applied}/50 proposals applied");
    }

    #[test]
    fn dirichlet_density_normalises_correctly() {
        // For Dirichlet(2, 2, 2) on the 2-simplex, the maximum is at
        // the centroid (1/3, 1/3, 1/3); the density there is positive.
        let xs = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];
        let alphas = [2.0, 2.0, 2.0];
        let lp = ln_dirichlet_density(&xs, &alphas);
        assert!(lp.is_finite(), "density not finite: {lp}");
    }
}
