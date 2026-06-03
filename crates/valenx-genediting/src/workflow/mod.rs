//! Workflow — Group F (features 27–30).
//!
//! The orchestration layer that ties Groups A–E into end-to-end
//! design workflows:
//!
//! - [`advisor`] — feature 27: the edit-strategy advisor — given a
//!   desired genomic change, compares nuclease vs base vs prime editing
//!   and recommends one with a rationale.
//! - [`variant`] — feature 28: the variant-correction planner — given
//!   a pathogenic variant, designs the editing strategy that reverts
//!   it to wild type.
//! - [`batch`] — feature 29: the batch editing-design driver over many
//!   targets, with a result table and captured per-target failures.
//! - [`driver`] — feature 30: the top-level
//!   [`driver::run_editing_design`] driver, the [`driver::EditingReport`]
//!   bundle (and, for mRNA, [`crate::mrna::design::MrnaDesignReport`]),
//!   and a typed [`driver::GeneditingRequest`] /
//!   [`driver::GeneditingResponse`] surface an external LLM can drive
//!   over an MCP tool.
//!
//! ## v1 scope
//!
//! The advisor and the planner use transparent decision rubrics, not
//! trained recommenders. The driver chooses one design per request (a
//! dispatcher, not an exhaustive search). The LLM surface captures
//! errors into a serialisable response rather than propagating them.

pub mod advisor;
pub mod batch;
pub mod driver;
pub mod variant;
