//! Maximum-likelihood phylogenetics.
//!
//! Likelihood methods score a tree by the probability of the observed
//! alignment given the tree, a [substitution model](model) and branch
//! lengths. They are the statistical backbone of PhyML, RAxML-NG and
//! IQ-TREE.
//!
//! - [`model`] — the standard nucleotide substitution models (JC69,
//!   K80, F81, HKY85, GTR): each builds a rate matrix `Q` and, by
//!   eigendecomposition, a transition-probability matrix `P(t) = e^{Qt}`.
//! - [`gamma`] — discrete-gamma rate heterogeneity (Yang 1994):
//!   site-to-site rate variation approximated by a handful of rate
//!   categories.
//! - [`felsenstein`] — Felsenstein's pruning algorithm (1981): the
//!   log-likelihood of a tree + alignment + model, computed in one
//!   postorder traversal.
//! - [`optimize`] — maximum-likelihood branch-length optimisation
//!   (golden-section per branch) and an NNI hill-climb on the
//!   likelihood score.
//! - [`ancestral`] — marginal ancestral-state reconstruction: the
//!   posterior state distribution at each internal node.

pub mod ancestral;
pub mod felsenstein;
pub mod gamma;
pub mod model;
pub mod optimize;

pub use ancestral::{ancestral_states, AncestralResult};
pub use felsenstein::{log_likelihood, log_likelihood_gamma, LikelihoodModel};
pub use gamma::DiscreteGamma;
pub use model::{SubstModel, TransitionMatrix};
pub use optimize::{
    optimize_branch_lengths, optimize_topology_ml, optimize_topology_ml_multistart,
    optimize_topology_ml_spr, MlSearchReport,
};
