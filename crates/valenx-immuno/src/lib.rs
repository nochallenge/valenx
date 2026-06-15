//! # valenx-immuno
//!
//! T-cell epitope / MHC-binding **immunogenicity screening** by the classic
//! position-specific scoring matrix (PSSM) method — the safety-screen arm of a
//! biologic-design pipeline that asks "could this designed protein provoke an
//! unwanted immune response?".
//!
//! ## What
//!
//! Given an amino-acid sequence (a candidate binder, enzyme, or any protein),
//! this crate slides a fixed-length window across it, scores every peptide
//! against a [`matrix::Pssm`] for a chosen MHC allele, and reports the
//! high-scoring windows as predicted T-cell epitopes. From those it derives a
//! single [`scan::epitope_density`] summary — a coarse immunogenicity flag.
//!
//! Core pieces:
//!
//! - [`matrix::Pssm`] — a per-position, per-residue weight matrix; `score`
//!   sums the weights of a peptide's residues.
//! - [`scan::scan`] / [`scan::scan_threshold`] / [`scan::top_n`] — slide the
//!   matrix over a protein and rank the windows.
//! - [`scan::epitope_density`] — fraction of windows above a threshold.
//! - [`library::illustrative_hla_a0201`] — one clearly-labelled *illustrative*
//!   class-I anchor motif so the engine is usable out of the box.
//!
//! ## Model
//!
//! The matrix method (Parker-, SYFPEITHI-style) treats MHC binding as additive
//! over peptide positions: `score(p) = sum_i W[i][p_i]`, where `W[i][a]` is the
//! contribution of residue `a` at position `i`. Class-I binding is dominated by
//! a few *anchor* positions (for many alleles, P2 and the C-terminus P9), so a
//! matrix that rewards the favoured anchor residues ranks plausible binders
//! above non-binders. This is a transparent, fully-auditable baseline — every
//! number that produces a score is inspectable.
//!
//! ## Honest scope
//!
//! Research/educational grade. The matrix method is a real, citable technique,
//! but it is a **baseline**: it ignores peptide processing, TAP transport,
//! proteasomal cleavage, MHC-allele coverage beyond the supplied matrix, and
//! the non-additive effects that modern experimentally-trained predictors
//! (IEDB, NetMHCpan) capture. The matrix shipped in [`library`] is an
//! **illustrative** encoding of documented anchor preferences, *not* a
//! quantitative experimental matrix. Nothing here is a validated clinical
//! immunogenicity assessment: any flag must be confirmed experimentally and
//! interpreted by a qualified immunologist. Do not use it to make clinical,
//! regulatory, or safety-of-use decisions.
//!
//! ## Example
//!
//! ```
//! use valenx_immuno::library::illustrative_hla_a0201;
//! use valenx_immuno::scan::{scan, top_n};
//!
//! let pssm = illustrative_hla_a0201();
//! // A short stretch of protein to screen (>= 9 residues).
//! let hits = scan(&pssm, "GILGFVFTLKYAS").unwrap();
//! assert_eq!(hits.len(), "GILGFVFTLKYAS".len() - 9 + 1);
//! let best = top_n(hits, 1);
//! assert_eq!(best[0].peptide.len(), 9);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod aa;
pub mod error;
pub mod library;
pub mod matrix;
pub mod scan;

pub use error::ImmunoError;
pub use matrix::Pssm;
pub use scan::{epitope_density, scan, scan_threshold, top_n, EpitopeHit};
