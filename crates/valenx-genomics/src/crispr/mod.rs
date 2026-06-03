//! CRISPR — guide design, off-target enumeration, edit-outcome analysis.
//!
//! The CRISPR-tooling block — a native replacement for CHOPCHOP,
//! Cas-OFFinder and CRISPResso2:
//!
//! - [`guide`] — PAM scanning for any nuclease ([`guide::PamSpec`]),
//!   both-strand protospacer discovery, and a Doench-Rule-Set-2-style
//!   on-target efficiency score.
//! - [`offtarget`] — a mismatch-tolerant both-strand genome scan
//!   ([`offtarget::enumerate_off_targets`]) and a CFD-style off-target
//!   activity score with a guide-specificity aggregate.
//! - [`editing`] — CRISPResso-class amplicon edit-outcome analysis:
//!   align edited reads to the reference, classify the indel spectrum
//!   and report editing efficiency.
//!
//! ## v1 scope
//!
//! The on-target and off-target scores are **transparent
//! feature-weighted heuristics** in the spirit of Doench Rule-Set-2
//! and the CFD matrix — the right shape and ranking, not the published
//! trained coefficients (a trained model the project's "no trained
//! weights" rule keeps out). The off-target scan is exhaustive within
//! a mismatch budget but not seed-index-accelerated; bulges are not
//! modelled. See each module's own note.

pub mod editing;
pub mod guide;
pub mod offtarget;
