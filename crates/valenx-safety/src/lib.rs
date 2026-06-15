//! # valenx-safety
//!
//! The **safety gate's consolidator**: take a candidate's independent screens —
//! off-target / cross-reactivity, immunogenicity, CRISPR off-target — and merge
//! them into **one per-candidate [`RiskReport`]**: severity-ranked [`RiskFlag`]s
//! with their source and evidence, plus an aggregate severity.
//!
//! ## What
//!
//! - [`Severity`] — `Info < Low < Moderate < High < Critical`.
//! - [`RiskFlag`] — one flag: source, severity, detail, optional numeric
//!   evidence. Build directly, or from a screen via [`flag::offtarget_flag`],
//!   [`flag::immunogenicity_flag`], [`flag::crispr_offtarget_flag`] (each returns
//!   `None` when the screen is below threshold).
//! - [`RiskReport`] — the consolidated record for one candidate
//!   ([`report::consolidate`]); [`RiskReport::aggregate_severity`] is the worst
//!   flag, [`RiskReport::render`] is the human report.
//!
//! ## The one rule it enforces by construction
//!
//! **It never says "safe".** There is no safe/approved/green-light field.
//! [`RiskReport::requires_human_review`] is always `true`. An empty flag list
//! means *no screen raised a flag* — which is **not** the same as safe, and the
//! report says so. A high-identity off-target (e.g. GDF-11 against an
//! anti-myostatin design) is surfaced as a flag for a human to weigh, never an
//! automatic pass or block.
//!
//! ## Honest scope
//!
//! Research/educational. The screens it consolidates are heuristics, and the
//! severity bands here are illustrative cutoffs, not validated thresholds. This
//! is a triage aid that organises evidence for a qualified reviewer — it is
//! **not** a validated safety assessment and makes no safety guarantee.
//!
//! ## Example
//!
//! ```
//! use valenx_safety::flag::offtarget_flag;
//! use valenx_safety::{RiskReport, Severity};
//!
//! // A ~90% off-target hit (GDF-11 vs an anti-myostatin design) is a serious flag.
//! let flag = offtarget_flag("GDF-11", 0.899, 0.80).unwrap().unwrap();
//! let report = RiskReport::new("design_A").unwrap().with_flag(flag);
//! assert!(report.aggregate_severity().unwrap() >= Severity::High);
//! assert!(report.requires_human_review()); // always
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod flag;
pub mod report;
pub mod severity;

pub use error::SafetyError;
pub use flag::RiskFlag;
pub use report::{consolidate, RiskReport};
pub use severity::Severity;
