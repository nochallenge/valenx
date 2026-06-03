//! Boltzmann-ensemble analysis.
//!
//! Where [`crate::fold`] finds the single best structure, this module
//! characterises the whole *ensemble* of structures a sequence can
//! adopt:
//!
//! - [`partition`] — the McCaskill partition function, the ensemble
//!   free energy, and the base-pair probability matrix.
//! - [`linear_partition`] — LinearPartition: a linear-time beam-search
//!   partition function + approximate base-pair probabilities for long
//!   sequences.
//! - [`centroid`] — the centroid structure (`p > 0.5` pairs) and the
//!   maximum-expected-accuracy (MEA) structure.
//! - [`suboptimal`] — Zuker-Stiegler enumeration of all structures
//!   within ΔE of the MFE.
//! - [`sampling`] — Boltzmann stochastic structure sampling
//!   (stochastic traceback) with the deterministic [`rng`].
//! - [`melting`] — a temperature-dependent melting curve.
//! - [`kinetics`] — Kinfold-class kinetic folding: stochastic
//!   Monte-Carlo Metropolis / Kawasaki walks through the structure
//!   ensemble, trajectory output, and ensemble first-passage statistics.

pub mod centroid;
pub mod kinetics;
pub mod linear_partition;
pub mod melting;
pub mod partition;
pub mod rng;
pub mod sampling;
pub mod suboptimal;
