//! # valenx-select
//!
//! The **selection funnel**: turn many scored candidates into a ranked,
//! diverse, defensible shortlist the operator can cut to any `N`. The value
//! here is the *cut* — from thousands to a handful that "fail differently" —
//! not volume.
//!
//! ## What
//!
//! - [`consensus`] — rank candidates by *agreement* across orthogonal scoring
//!   methods (docking, MM-GBSA, interface confidence …) via Borda-style
//!   fractional-rank aggregation, and surface where methods **disagree** as a
//!   per-candidate low-confidence signal.
//! - [`diversity`] — pick a structurally/feature-diverse subset:
//!   [`diversity::farthest_point_select`] (MaxMin spread) and
//!   [`diversity::sphere_exclusion_select`] (score-ordered cluster-and-top), so
//!   the chosen `N` are not near-duplicates.
//! - [`funnel`] — [`funnel::select_shortlist`]: consensus → diversify → top-`N`,
//!   carrying each candidate's calibrated confidence and safety flags through to
//!   the [`funnel::Shortlist`]. `N` is just where the operator cuts.
//!
//! This crate is **standalone**: it takes generic per-method scores and feature
//! vectors, so it needs none of the upstream engines to build or test. Wire it
//! to `valenx-score` (ComparableScore), `valenx-calibrate` (confidence) and
//! `valenx-offtarget` (safety flags) at the call site.
//!
//! ## Model
//!
//! Consensus uses fractional ranks in `[0, 1]` (1 = best) per method; the
//! consensus is their mean and the disagreement their standard deviation.
//! Diversity selection is greedy: MaxMin repeatedly adds the point farthest
//! from those already chosen; sphere exclusion walks candidates best-score-first
//! and accepts one only if it lies beyond `radius` of every accepted point
//! (Butina-style cluster-and-top).
//!
//! ## Honest scope
//!
//! Research/educational grade. These are standard, transparent **selection
//! heuristics** — rank aggregation and diversity picking. They improve the odds
//! that a shortlist is robust and non-redundant, but they are not a guarantee:
//! a diverse, high-consensus candidate can still fail in the wet lab, and the
//! disagreement signal flags uncertainty, it does not resolve it. Nothing here
//! is a validated success predictor.
//!
//! ## Example
//!
//! ```
//! use valenx_select::consensus::consensus_borda;
//!
//! // Two methods that agree: candidate 0 is best in both.
//! let res = consensus_borda(&[vec![3.0, 1.0, 2.0], vec![9.0, 1.0, 5.0]]).unwrap();
//! assert!(res.consensus[0] > res.consensus[2]);  // 0 ranks above 2
//! assert!(res.disagreement[0].abs() < 1e-12);    // they agree -> no disagreement
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod consensus;
pub mod diversity;
pub mod error;
pub mod funnel;

pub use consensus::{consensus_borda, ConsensusResult};
pub use diversity::{farthest_point_select, sphere_exclusion_select};
pub use error::SelectError;
pub use funnel::{select_shortlist, Candidate, Shortlist, ShortlistEntry};
