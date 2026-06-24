//! # valenx-orchestrator
//!
//! The **single-command design funnel**. It chains the three standalone
//! pipeline crates into one run so an operator can go from many scored
//! candidates to one signed, fingerprintable dossier with a single call:
//!
//! 1. [`valenx_select`] — consensus-rank candidates across orthogonal scoring
//!    methods, then take a feature-diverse top-`N`.
//! 2. [`valenx_safety`] — consolidate each shortlisted candidate's off-target,
//!    immunogenicity and CRISPR screens into one risk report.
//! 3. [`valenx_dossier`] — assemble goal, ranked candidates, calibration
//!    provenance and software list into one content-hashed [`RunDossier`].
//!
//! [`run_funnel`] is the whole pipeline; [`FunnelOutcome`] carries the dossier
//! plus every intermediate product. [`run_funnel_seqs`] is the sequence-driven
//! entry point: give it candidate amino-acid sequences and a reference panel and
//! it *runs the screens itself* — off-target identity ([`valenx_offtarget`]),
//! T-cell epitope density ([`valenx_immuno`]), developability
//! ([`valenx_developability`]) and linear B-cell epitopes
//! ([`valenx_epitope_map`]) — deriving diversity features from amino-acid
//! composition, then funnels through the same selection → safety → dossier path.
//!
//! ## What it does *not* do
//!
//! - **Never fabricates.** The de-novo *generate*, structure-based *dock* and
//!   physics-based *score* stages need gated resources (trained weights, a
//!   docking engine, an experimental structure, a GPU or a license). When one is
//!   declared unavailable it is recorded as `BLOCKED: <dep>` and skipped — the
//!   core still runs on whatever real candidate scores were supplied. See
//!   [`GatedStage`] / [`BlockedStage`].
//! - **Never marks a candidate safe.** A screen with no input is recorded as
//!   `Info: not run`, not as a pass. Every run
//!   [`requires_human_signoff`](FunnelOutcome::requires_human_signoff).
//!
//! ## Digital-engineering spine
//!
//! Alongside the biologic funnel, the [`digeng`] module is a standalone
//! **systems-engineering / MBSE** layer — the digital thread that ties design
//! parameters → requirements → trade studies into one auditable record:
//! [`digeng::Requirement`]s with signed-margin [`digeng::Verdict`]s,
//! [`digeng::DesignPoint`] compliance, a full-factorial
//! [`digeng::TradeStudy`] driver with Pareto-front extraction, and
//! requirements-↔-metrics [`digeng::coverage`]. It is pure systems-engineering
//! book-keeping (zero dual-use concern) and shares no state with the funnel;
//! sampling-based / UQ design-of-experiments is a documented future hook, not a
//! dependency.
//!
//! ## Honest scope
//!
//! Research/educational grade. This crate is *plumbing* — it composes
//! transparent selection and safety heuristics and records provenance. It adds
//! no new predictive power and is not a validated design or safety pipeline.
//!
//! ## Example
//!
//! ```
//! use valenx_orchestrator::{run_funnel, FunnelCandidate, FunnelConfig, GatedStage};
//! use valenx_dossier::CalibrationStatus;
//!
//! let candidates = vec![
//!     FunnelCandidate {
//!         id: "design_A".into(),
//!         method_scores: vec![3.0, 9.0], // best in both methods
//!         features: vec![0.0, 0.0],
//!         calibrated_confidence: None,
//!         offtarget: None,
//!         immunogenicity: None,
//!         crispr_offtarget_sites: None,
//!         developability_flags: Vec::new(),
//!         bcell_epitope_regions: None,
//!     },
//!     FunnelCandidate {
//!         id: "design_B".into(),
//!         method_scores: vec![1.0, 1.0],
//!         features: vec![1.0, 1.0],
//!         calibrated_confidence: None,
//!         offtarget: None,
//!         immunogenicity: None,
//!         crispr_offtarget_sites: None,
//!         developability_flags: Vec::new(),
//!         bcell_epitope_regions: None,
//!     },
//! ];
//!
//! let config = FunnelConfig {
//!     top_n: 2,
//!     diversity_radius: 0.5,
//!     offtarget_threshold: 0.8,
//!     immunogenicity_threshold: 0.1,
//!     crispr_threshold: 1,
//!     calibration: CalibrationStatus::blocked("no held-out ground truth"),
//!     blocked_stages: vec![GatedStage::Generate, GatedStage::Dock, GatedStage::Score],
//! };
//!
//! let out = run_funnel("inhibit target X", &candidates, &config).unwrap();
//! assert_eq!(out.dossier.ranked()[0].id, "design_A"); // best consensus
//! assert!(out.requires_human_signoff());              // always
//! assert_eq!(out.blocked.len(), 3);                   // gated stages skipped, not faked
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod digeng;
pub mod error;
pub mod funnel;
pub mod seq;

pub use error::OrchestratorError;
pub use funnel::{
    run_funnel, BlockedStage, FunnelCandidate, FunnelConfig, FunnelOutcome, GatedStage,
    OfftargetEvidence,
};
pub use seq::{run_funnel_seqs, ScreenConfig, SeqCandidate};

// Digital-engineering / MBSE spine (see [`digeng`]).
pub use digeng::{
    coverage, pareto_front, Comparator, ComplianceReport, CoverageReport, DesignPoint, Direction,
    Objective, Parameter, ParameterSweep, Requirement, TradeResult, TradeStudy, TradeStudyOutcome,
    Verdict,
};

// Re-export the dossier type the public API returns inside [`FunnelOutcome`].
pub use valenx_dossier::RunDossier;
