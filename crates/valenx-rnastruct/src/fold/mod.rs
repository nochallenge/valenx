//! Secondary-structure folding.
//!
//! This module groups the folding algorithms — turning an RNA
//! *sequence* into a *structure*:
//!
//! - [`turner2004`] — the complete published Turner-2004 nearest-
//!   neighbor parameter tables (stacking, loop-length, terminal
//!   mismatches, dangles, small-loop special cases, the multiloop
//!   model, the coaxial-stacking term). Pure data.
//! - [`energy`] — the loop-energy layer on top of [`turner2004`]: it
//!   assembles a table lookup, the length term, the terminal penalty,
//!   the mismatch bonus and any special-case override into the single
//!   per-loop free energy every energy-based folder draws on.
//! - [`coaxial`] — the coaxial-stacking correction: the stabilising
//!   energy when two helices in a loop lie end to end. Makes
//!   multi-helix folds match ViennaRNA `RNAfold -d2`.
//! - [`nussinov`] — the Nussinov-Jacobson maximum-base-pairing DP,
//!   the simplest folder (no energy model).
//! - [`eval`] — free-energy evaluation of a given (sequence,
//!   structure) pair by loop decomposition.
//! - [`zuker`] — the Zuker-Stiegler minimum-free-energy folding DP
//!   under the Turner model — the workhorse folder.
//! - [`linear`] — LinearFold: linear-time beam-search MFE folding for
//!   long sequences (Huang et al. 2019).
//! - [`constraint`] — hard structural constraints and soft
//!   pseudo-energies (the mechanism behind SHAPE folding).
//! - [`shape`] — SHAPE-reactivity-directed folding (Deigan model).

pub mod coaxial;
pub mod constraint;
pub mod energy;
pub mod eval;
pub mod linear;
pub mod nussinov;
pub mod shape;
pub mod turner2004;
pub mod zuker;
