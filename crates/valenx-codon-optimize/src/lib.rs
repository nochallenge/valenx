//! # valenx-codon-optimize
//!
//! Turn a designed protein into a **synthesis-ready DNA coding sequence** and
//! score it — the "make it orderable from a gene synthesis vendor" step of a
//! biologic-design pipeline.
//!
//! ## What
//!
//! - [`code::translate_codon`] / [`code::synonymous_codons`] — the standard
//!   genetic code (NCBI translation table 1).
//! - [`usage::CodonUsage`] — relative-adaptiveness weights for a host, built
//!   from codon frequencies or supplied directly; [`usage::CodonUsage::optimal_codon`]
//!   gives the most-adapted synonymous codon.
//! - [`optimize::reverse_translate_optimal`] — protein → DNA picking the
//!   highest-adaptiveness codon per residue.
//! - [`optimize::cai`] — the Codon Adaptation Index (Sharp & Li, 1987); the
//!   geometric mean of relative adaptiveness over the coding sequence.
//! - [`optimize::gc_content`] — GC fraction.
//!
//! ## Model
//!
//! The genetic code is the canonical NCBI table 1 (verifiable, not a parameter).
//! The Codon Adaptation Index is the classic Sharp & Li geometric-mean measure
//! of how well a sequence's codons match a reference host's preferences;
//! "optimization" greedily picks the maximal-adaptiveness synonymous codon per
//! residue. By convention CAI excludes Met and Trp (single-codon) and stop
//! codons.
//!
//! ## Honest scope
//!
//! Research/educational. The genetic code is canonical. The **host codon-usage
//! weights are caller-supplied** — the included [`usage::illustrative_weights`]
//! is an *illustrative* set, **not** a specific organism's verified table; for
//! real work supply a measured table (e.g. from the Kazusa codon-usage database
//! for your expression host) and record its source. Greedy max-CAI optimization
//! also ignores real-world constraints (restriction sites, mRNA secondary
//! structure, repeats, GC windows), so the output is a starting point for a
//! synthesis vendor's checks, not a finished construct.
//!
//! ## Example
//!
//! ```
//! use valenx_codon_optimize::code::{synonymous_codons, translate_codon};
//!
//! assert_eq!(translate_codon("ATG"), Some('M'));
//! assert_eq!(translate_codon("TGG"), Some('W'));
//! assert_eq!(translate_codon("TAA"), Some('*')); // stop
//! assert_eq!(synonymous_codons('M'), vec!["ATG"]); // Met is single-codon
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod code;
pub mod error;
pub mod optimize;
pub mod usage;

pub use error::CodonError;
pub use optimize::{cai, gc_content, reverse_translate_optimal};
pub use usage::{illustrative_weights, CodonUsage};
