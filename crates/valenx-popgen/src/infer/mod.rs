//! Inference and the genotype-matrix exchange type.
//!
//! - [`matrix`] — the [`GenotypeMatrix`], the `n_samples x n_sites`
//!   biallelic matrix that every statistic, exporter and inference
//!   routine consumes. A [`crate::model::Population`], a coalescent
//!   genealogy and a [`crate::coalescent::TreeSequence`] can all be
//!   projected into one.
//! - [`abc`] — the Approximate Bayesian Computation rejection
//!   framework: [`Prior`], [`AbcConfig`], [`abc_reject`] and the
//!   [`AbcPosterior`].

pub mod abc;
pub mod matrix;

pub use abc::{abc_reject, AbcConfig, AbcPosterior, Distance, Prior};
pub use matrix::GenotypeMatrix;
