//! # valenx-phylo — phylogenetics and molecular evolution
//!
//! A native-Rust replacement
//! for the inference core of BEAST 2, MrBayes, IQ-TREE, RAxML-NG,
//! PhyML, FastTree, RevBayes, Seq-Gen, FigTree and Dendroscope —
//! pure algorithms, no neural-network weights and no external
//! processes.
//!
//! It builds on [`valenx_bioseq`] (Block 6.1) and [`valenx_align`]
//! (Block 6.2): a [`valenx_align::Msa`] is the natural input for the
//! [`distance`] and [`likelihood`] modules.
//!
//! ## What it does
//!
//! - **Model** ([`tree`]) — the rooted / unrooted [`Tree`] arena with
//!   branch lengths, leaf and internal labels, preorder / postorder
//!   traversal and LCA / patristic-distance queries.
//! - **I/O & rendering** ([`io`], [`render`]) — Newick and NEXUS
//!   read + write, PhyloXML / NeXML writers, and ASCII-cladogram /
//!   SVG-phylogram rendering.
//! - **Distance methods** ([`distance`]) — p-distance, Jukes-Cantor,
//!   Kimura-2P and Tamura-Nei distance matrices, then UPGMA, WPGMA,
//!   neighbor-joining and BIONJ clustering.
//! - **Parsimony** ([`parsimony`]) — Fitch small-parsimony, Sankoff
//!   weighted parsimony and an NNI / SPR hill-climbing large-parsimony
//!   search.
//! - **Likelihood** ([`likelihood`]) — the JC69 / K80 / F81 / HKY85 /
//!   GTR substitution models, Felsenstein pruning log-likelihood,
//!   discrete-gamma rate heterogeneity, branch-length and topology
//!   optimisation, and marginal ancestral-state reconstruction.
//! - **Support & comparison** ([`compare`]) — bootstrap support,
//!   Robinson-Foulds and quartet distances, majority-rule / strict
//!   consensus, and tree rerooting / ladderising / pruning.
//! - **Simulation** ([`simulate`]) — coalescent and birth-death tree
//!   simulation, Seq-Gen-class sequence evolution along a tree, and
//!   root-to-tip molecular-clock dating.
//!
//! ## Errors
//!
//! Every fallible function returns
//! [`Result<_, PhyloError>`](error::PhyloError); the error type
//! carries stable [`code`](error::PhyloError::code) and
//! [`category`](error::PhyloError::category) accessors for telemetry.
//!
//! ## v1 scope
//!
//! This is a real working v1, not production parity with the reference
//! tools. Each module documents its own simplifications inline; the
//! current scope:
//!
//! - **Bayesian MCMC** ([`bayes`]) ships as a real
//!   Metropolis-Hastings sampler — NNI / SPR / Wilson-Balding
//!   topology proposals, branch-length scale / slide / tree-scale,
//!   `κ` multiplier, GTR rate + frequency Dirichlet, gamma `α`
//!   multiplier, with Geyer-IMPS effective sample size and the
//!   Gelman-Rubin `R̂` cross-chain diagnostic — but it is **not**
//!   BEAST 2 / MrBayes parity: no relaxed-clock or tip-dating models,
//!   no reversible-jump model selection, no operator-tuning auto-
//!   adaptation, no BEAUTi-style XML, no MC³ / coupled chains.
//! - **Maximum-likelihood topology search** is now an **NNI + SPR**
//!   hill-climb ([`likelihood::optimize_topology_ml_spr`]) with
//!   [`likelihood::optimize_topology_ml_multistart`] for restarting
//!   from multiple starting trees, closer to the IQ-TREE / RAxML-NG
//!   default but without their model selection, ultrafast bootstrap,
//!   or perturbation kernels.
//! - The random-number generator is a small deterministic PCG so
//!   simulations are reproducible but not cryptographic; and the
//!   substitution models are the standard nucleotide family only (no
//!   codon or amino-acid models).

#![forbid(unsafe_code)]
#![warn(missing_docs)]
// Surface future `&str` byte-offset slicing in clippy review — this
// crate parses untrusted text, where non-char-boundary slices panic.
// WARN (not deny): most existing slices are safe ASCII; this only flags
// NEW ones.
#![allow(clippy::string_slice, reason = "parsers slice ASCII fixed-format records at byte offsets from find() or constant ASCII prefixes, always valid char boundaries")]

pub mod bayes;
pub mod compare;
pub mod distance;
pub mod error;
pub mod io;
pub mod likelihood;
pub mod parsimony;
pub mod render;
pub mod rng;
pub mod simulate;
pub mod tree;

// --- Convenience re-exports of the most-used types --------------------

pub use error::{ErrorCategory, PhyloError, Result};
pub use tree::{Node, NodeId, Tree};

pub use bayes::{
    run_chain, summarise_posterior, ChainConfig, ChainResult, ChainState,
    ParameterDiagnostics, PosteriorSummary, Prior, ProposalSet,
};
pub use distance::{distance_matrix, DistMatrix, DistanceModel};
pub use io::{read_newick, read_nexus, write_newick, write_nexus};
pub use likelihood::{log_likelihood, SubstModel};
pub use parsimony::{fitch_parsimony, sankoff_parsimony};
pub use render::{render_ascii, render_svg};

#[cfg(test)]
mod tests {
    use super::*;

    /// End-to-end: parse a tree, round-trip it through Newick, and
    /// confirm the topology survives.
    #[test]
    fn newick_round_trip_end_to_end() {
        let src = "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);";
        let tree = read_newick(src).unwrap();
        assert_eq!(tree.leaf_count(), 4);
        let again = read_newick(&write_newick(&tree)).unwrap();
        assert_eq!(tree.leaf_labels(), again.leaf_labels());
    }
}
