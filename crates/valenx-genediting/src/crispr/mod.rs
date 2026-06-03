//! CRISPR nuclease editing — Group A (features 1–5).
//!
//! The nuclease-editing block: designing double-strand-break editing
//! experiments with Cas9- and Cas12-class nucleases. It is the
//! workflow layer on top of [`valenx_genomics`]' CRISPR primitives —
//! the PAM scanner, the on-target score and the off-target enumerator
//! are reused directly, never duplicated.
//!
//! - [`nuclease`] — feature 1: the nuclease database (SpCas9,
//!   SpCas9-NG, SaCas9, Cas12a / Cpf1, Cas12f, Cas13, xCas9) with PAM,
//!   guide length, cut-site offset and end chemistry.
//! - [`guide_design`] — feature 2: the guide-RNA design workflow for a
//!   target region, with on-target ranking and optional off-target
//!   scoring.
//! - [`knockout`] — feature 3: NHEJ knockout strategy design —
//!   early-exon / functional-domain targeting and frameshift-optimal
//!   guide placement.
//! - [`donor`] — feature 4: HDR knock-in donor-template design —
//!   homology arms plus silent PAM- / seed-blocking mutations to stop
//!   re-cutting.
//! - [`donor_opt`] — commercial-depth HDR donor optimisation:
//!   multiple stacked silent mutations across the seed and PAM, with
//!   CAI-aware codon-swap selection and cryptic-splice-site avoidance.
//! - [`multiplex`] — feature 5: multiplex editing and gRNA-array
//!   design.
//! - [`offtarget_fm`] — commercial-depth genome-wide off-target search
//!   over a `valenx-align` SA-IS FM-index. Reuses the production
//!   index + the standard BWA-style seed-and-extend mismatch-tolerant
//!   search; the per-window scan in [`guide_design`] stays the default
//!   path because it ranges over a *user-supplied* genome window.
//!
//! ## v1 scope
//!
//! Every efficiency / quality score in this block is a transparent
//! feature-weighted heuristic, documented as such — there is no
//! trained model anywhere in it (the project's "no trained-weights"
//! rule). Donor and knockout design operate on a linear reference /
//! gene-model window the caller supplies; they do not fetch a genome
//! or resolve isoforms.

pub mod donor;
pub mod donor_opt;
pub mod guide_design;
pub mod knockout;
pub mod multiplex;
pub mod nuclease;
pub mod offtarget_fm;
