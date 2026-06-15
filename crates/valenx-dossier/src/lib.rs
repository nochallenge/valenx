//! # valenx-dossier
//!
//! Assemble a biologic-design run into **one signed, reproducible dossier** — the
//! artifact an operator reads to cut a shortlist, and an auditor reads to verify
//! a run. It records the goal, the ranked candidates with their broken-out
//! scores, their calibrated confidence *or* an honest "blocked" note when no
//! calibration ground truth exists, and their safety flags — then content-hashes
//! the whole thing through [`valenx_repro`] so the run is tamper-evident and
//! reproducible.
//!
//! ## What
//!
//! - [`ScoredCandidate`] — one candidate: id, comparable score, broken-out score
//!   components, optional calibrated confidence, and safety flags.
//! - [`CalibrationStatus`] — either [`CalibrationStatus::Calibrated`] (method +
//!   calibration set recorded) or [`CalibrationStatus::Blocked`] (e.g. *"SKEMPI
//!   dataset"*), so a missing benchmark is stated, never papered over with a fake
//!   confidence.
//! - [`RunDossier`] — the run: goal, candidates, calibration status, software
//!   manifest. [`RunDossier::ranked`] sorts by score, [`RunDossier::fingerprint`]
//!   is the reproducible SHA-256 hash, [`RunDossier::render`] is the human report.
//!
//! ## Two rules it enforces by construction
//!
//! 1. **Never "safe".** There is no "safe"/"approved" field.
//!    [`RunDossier::requires_human_signoff`] is always `true`; the rendered report
//!    states that an empty flag list is *not* a safety assertion and that
//!    promotion to wet-lab testing is a human decision.
//! 2. **Never fake calibration.** With no labelled benchmark the dossier carries
//!    [`CalibrationStatus::Blocked`] and reports the missing dependency, rather
//!    than emitting an uncalibrated number dressed as a probability.
//!
//! ## Honest scope
//!
//! Research/educational. This crate is the connective record-keeping layer — it
//! does not itself score, calibrate, or screen; it faithfully assembles and
//! hashes whatever the upstream stages produced. A reproducible fingerprint
//! proves a run is *unchanged*, not that its predictions are *correct*.
//!
//! ## Example
//!
//! ```
//! use valenx_dossier::{CalibrationStatus, RunDossier, ScoredCandidate};
//!
//! let d = RunDossier::new("Inhibit myostatin (GDF-8)", CalibrationStatus::blocked("SKEMPI dataset"))
//!     .unwrap()
//!     .with_candidate(ScoredCandidate::new("design_A", 0.82).unwrap())
//!     .with_candidate(ScoredCandidate::new("design_B", 0.91).unwrap());
//! // ranked best-first; nothing is ever marked "safe"
//! assert_eq!(d.ranked()[0].id, "design_B");
//! assert!(d.requires_human_signoff());
//! assert_eq!(d.fingerprint().unwrap().len(), 64);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod calibration;
pub mod candidate;
pub mod dossier;
pub mod error;

pub use calibration::CalibrationStatus;
pub use candidate::ScoredCandidate;
pub use dossier::RunDossier;
pub use error::DossierError;
