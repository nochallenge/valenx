//! # valenx-repro
//!
//! **Reproducibility-bundle packaging** — turns a computational study into a
//! single, verifiable package so that someone else (a skeptical reviewer, a
//! collaborator, your future self) can check exactly what was run and confirm
//! the results were not tampered with.
//!
//! ## What goes in a bundle
//!
//! A [`ReproBundle`] records:
//!
//! - **inputs** and **outputs** as [`Artifact`]s — each carrying a SHA-256 of
//!   its bytes, so the data is content-addressed,
//! - the **parameters** ([`Parameter`]) used,
//! - the ordered **workflow provenance** ([`ProvenanceStep`]) — which tool, at
//!   which version, with which arguments, in what order,
//! - the **software manifest** ([`SoftwareRef`]).
//!
//! From that it derives a single deterministic **[`fingerprint`]** (a SHA-256
//! over the canonicalised contents): any change to any artifact, parameter,
//! software entry or step flips it, so the fingerprint is a tamper-evidence
//! root you can publish alongside a finding. It also generates a templated
//! [`methods_scaffold`] and [`abstract_scaffold`] — plain fill-in-the-blanks
//! text built from the recorded facts, **no language model involved**.
//!
//! ```
//! use valenx_repro::{Artifact, ArtifactRole, ProvenanceStep, ReproBundle, SoftwareRef};
//!
//! let bundle = ReproBundle::new("MSTN translate", "Translate the MSTN CDS")
//!     .unwrap()
//!     .with_software(SoftwareRef::new("valenx-bioseq", "0.1.0"))
//!     .with_artifact(Artifact::from_bytes("mstn.fasta", ArtifactRole::Input, b"ATGCAA"))
//!     .with_step(ProvenanceStep::new(1, "translate", "0.1.0", "frame +1"))
//!     .unwrap();
//! assert_eq!(bundle.fingerprint().len(), 64); // a SHA-256, in hex
//! ```
//!
//! [`fingerprint`]: ReproBundle::fingerprint
//!
//! ## Honest scope
//!
//! This crate is **export-only**: it writes no files and performs **no network
//! action** — no auto-submission to a preprint server, no e-mailing, nothing
//! sent anywhere. It builds the package and computes the checksums; persisting
//! and sharing it is the caller's (and the human's) job. And a reproducible
//! bundle makes an **in-silico hypothesis** auditable — it is **not** a
//! substitute for wet-lab validation or peer review.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod bundle;
pub mod error;
pub mod scaffold;

pub use bundle::{Artifact, ArtifactRole, Parameter, ProvenanceStep, ReproBundle, SoftwareRef};
pub use error::ReproError;
pub use scaffold::{abstract_scaffold, methods_scaffold};
