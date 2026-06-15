//! # valenx-epitope-map
//!
//! Linear **B-cell epitope** propensity by the classic sliding-window
//! hydrophilicity method: hydrophilic, surface-exposed stretches of a protein
//! are the ones antibodies tend to recognise, so a windowed hydrophilicity
//! profile with its peaks flagged is a first-pass map of candidate linear
//! epitopes. Complements the T-cell (MHC) screen in `valenx-immuno`.
//!
//! ## What
//!
//! - [`scale::PropensityScale`] — a per-residue propensity scale;
//!   [`scale::hydrophilicity_kd`] is the built-in default.
//! - [`map::propensity_profile`] — the windowed mean propensity along the
//!   sequence.
//! - [`map::linear_epitope_regions`] — contiguous spans whose windowed
//!   propensity stays at or above a threshold (the predicted epitopes).
//!
//! ## Model
//!
//! This is the Parker / Hopp-Woods approach: average a residue propensity over
//! a sliding window (classically 7 residues) and call the peaks. **Scale note,
//! stated honestly:** the original method uses the Parker (1986) or Hopp-Woods
//! (1981) hydrophilicity scales. To avoid shipping numeric constants this crate
//! has not independently verified, the built-in [`scale::hydrophilicity_kd`]
//! scale is the **negated Kyte-Doolittle** hydropathy (a hydrophilicity scale
//! whose values this codebase does use and test elsewhere). Supply your own
//! [`scale::PropensityScale`] (e.g. the verified Parker values) for method
//! fidelity.
//!
//! ## Honest scope
//!
//! Research/educational. Sliding-window hydrophilicity is a transparent,
//! decades-old heuristic, not a validated epitope predictor (modern tools —
//! BepiPred, DiscoTope — use structure and machine learning). A flagged region
//! is a hypothesis for experimental epitope mapping, nothing more.
//!
//! ## Example
//!
//! ```
//! use valenx_epitope_map::map::linear_epitope_regions;
//! use valenx_epitope_map::scale::hydrophilicity_kd;
//!
//! // A hydrophilic charged stretch is flagged as a candidate epitope.
//! let regions = linear_epitope_regions("DEKRDEKR", &hydrophilicity_kd(), 3, 0.0).unwrap();
//! assert!(!regions.is_empty());
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod map;
pub mod scale;

pub use error::EpitopeError;
pub use map::{linear_epitope_regions, propensity_profile};
pub use scale::{hydrophilicity_kd, PropensityScale};
