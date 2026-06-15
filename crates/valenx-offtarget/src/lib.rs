//! # valenx-offtarget
//!
//! Off-target / cross-reactivity screening for a designed protein candidate —
//! the second arm of a biologic-design safety gate, alongside immunogenicity.
//! It asks: "does this candidate look dangerously similar to something it should
//! not bind — a human self-protein, an essential off-target?".
//!
//! ## What
//!
//! Given a candidate sequence and a *reference set* (e.g. a panel of
//! self-proteins or essential targets), it computes two complementary
//! similarity measures against every reference, flags the references that
//! exceed a threshold, and reduces them to a single risk summary:
//!
//! - [`similarity::best_ungapped_identity`] — the best fractional residue
//!   identity over any ungapped overlap of the two sequences (a positional
//!   look-alike score).
//! - [`similarity::kmer_jaccard`] — Jaccard overlap of the two sequences'
//!   k-mer (k-peptide) sets (a composition / shared-motif score).
//! - [`screen::screen`] / [`screen::OffTargetReport`] — run both over the
//!   reference set and report the flagged hits and the maximum identity.
//!
//! ## Model
//!
//! Sequence similarity is the first-line, fully-transparent screen for
//! cross-reactivity: a candidate that is a near-identical local match to a
//! self-protein is an off-target liability. `best_ungapped_identity` slides the
//! shorter sequence across the longer and takes the maximum fraction of
//! matching positions; `kmer_jaccard` compares the sets of length-`k`
//! sub-peptides. Both are in `[0, 1]`; 1.0 means identical content.
//!
//! ## Honest scope
//!
//! Research/educational grade. Sequence-similarity screening is genuinely how
//! off-target risk is first triaged, but it is a **heuristic**: it does not
//! model structure, binding, expression, or tissue context, and ungapped
//! identity misses gapped homologies that a full Smith–Waterman / BLAST search
//! would find. This is **not** a validated off-target or cross-reactivity
//! safety predictor; a flag (or its absence) must be confirmed with proper
//! homology search and experiment. Do not use for clinical, regulatory, or
//! safety-of-use decisions.
//!
//! ## Example
//!
//! ```
//! use valenx_offtarget::screen::screen;
//!
//! let candidate = "MKTAYIAKQR";
//! let refs = [
//!     ("self_protein_1", "MKTAYIAKQR"), // identical -> flagged
//!     ("unrelated",      "WWWWWWWWWW"),
//! ];
//! let report = screen(candidate, &refs, 3, 0.8).unwrap();
//! assert_eq!(report.flagged.len(), 1);
//! assert_eq!(report.flagged[0].reference_id, "self_protein_1");
//! assert!((report.max_identity - 1.0).abs() < 1e-12);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod aa;
pub mod error;
pub mod screen;
pub mod similarity;

pub use error::OffTargetError;
pub use screen::{screen, OffTargetHit, OffTargetReport};
pub use similarity::{best_ungapped_identity, kmer_jaccard};
