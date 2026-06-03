//! Bayesian phylogenetic inference — a real Metropolis-Hastings MCMC
//! sampler over `(tree, branch lengths, substitution-model parameters)`.
//!
//! This is the [valenx-phylo] commercial-depth pass on Bayesian
//! inference: the framework that BEAST 2 / MrBayes / RevBayes ship as
//! the modern standard for tree inference with uncertainty.
//!
//! Modules:
//!
//! - [`prior`] — joint prior `P(tree, θ)`. Uniform topology +
//!   Exponential per-branch + per-model-parameter priors (Exponential
//!   on `κ`, symmetric Dirichlet on GTR rates / equilibrium
//!   frequencies, Exponential on gamma `α`).
//! - [`proposal`] — the Metropolis-Hastings proposal zoo: NNI / SPR /
//!   Wilson-Balding topology proposals, branch-length scale / slide /
//!   tree-scale, `κ` multiplier, GTR rate Dirichlet, frequency
//!   Dirichlet, gamma `α` multiplier — each returns its log Hastings
//!   ratio so the MH acceptance is correct.
//! - [`chain`] — [`run_chain`] runs one Metropolis-Hastings chain and
//!   produces a [`ChainResult`] (samples + per-iteration parameter
//!   trace + per-kind acceptance counts).
//! - [`diagnostics`] — Effective sample size (Geyer IMPS) and
//!   Gelman-Rubin `R̂` for cross-chain convergence checks.
//! - [`posterior`] — Posterior summaries: majority-rule consensus
//!   tree, MAP tree, per-clade posterior probabilities.
//!
//! Re-exports the most commonly used entry points so callers can write
//! `use valenx_phylo::bayes::{run_chain, ChainConfig, Prior};`.

pub mod chain;
pub mod diagnostics;
pub mod posterior;
pub mod prior;
pub mod proposal;

pub use chain::{run_chain, AcceptanceCounts, ChainConfig, ChainResult, ChainSample};
pub use diagnostics::{
    effective_sample_size, gelman_rubin, ParameterDiagnostics,
};
pub use posterior::{
    clade_posterior_table, clade_probability, summarise_posterior, CladePosterior,
    PosteriorSummary,
};
pub use prior::Prior;
pub use proposal::{ChainState, ProposalKind, ProposalSet};
