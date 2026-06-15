//! # valenx-developability
//!
//! Sequence-level **developability** flags for a designed protein or biologic —
//! the biophysical "will this be a nightmare to express, purify, and formulate?"
//! triage that complements binding and safety screens.
//!
//! ## What
//!
//! - [`hydrophobicity::gravy`] — the grand average of hydropathy (mean
//!   Kyte-Doolittle value); [`hydrophobicity::hydropathy_profile`] — a windowed
//!   profile; [`hydrophobicity::aggregation_prone_regions`] — windows whose mean
//!   hydropathy exceeds a threshold (an illustrative aggregation proxy).
//! - [`charge::net_charge_at_ph`] — net charge from the ionizable residues and
//!   termini; [`charge::isoelectric_point`] — the pH of zero net charge.
//! - [`report::assess`] / [`report::DevelopabilityReport`] — GRAVY, pI, net
//!   charge at physiological pH, aggregation-prone-region count, and flags.
//!
//! ## Model
//!
//! Hydrophobicity uses the published **Kyte-Doolittle** hydropathy scale (1982).
//! Charge uses the Henderson-Hasselbalch equation over the ionizable side chains
//! (D, E, C, Y, H, K, R) and the N-/C-termini with a standard **EMBOSS** pKa
//! set; the isoelectric point is found by bisection on the (monotone) net-charge
//! curve.
//!
//! ## Honest scope
//!
//! Research/educational. These are standard, transparent textbook calculations,
//! but real developability assessment (aggregation, viscosity, chemical
//! liabilities, expression) uses experimentally-trained tools and wet-lab
//! measurement. The aggregation-prone-region call here is a hydrophobicity
//! heuristic, not a validated predictor (TANGO / AGGRESCAN / Camsol), and pI
//! depends on the pKa set chosen. Treat every output as a flag for follow-up,
//! not a verdict.
//!
//! ## Example
//!
//! ```
//! use valenx_developability::hydrophobicity::gravy;
//!
//! // GRAVY of "IV" is the mean of the Kyte-Doolittle values 4.5 and 4.2.
//! assert!((gravy("IV").unwrap() - 4.35).abs() < 1e-9);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod aa;
pub mod charge;
pub mod error;
pub mod hydrophobicity;
pub mod report;

pub use charge::{isoelectric_point, net_charge_at_ph};
pub use error::DevelopabilityError;
pub use hydrophobicity::{aggregation_prone_regions, gravy, hydropathy_profile};
pub use report::{assess, DevelopabilityReport};
