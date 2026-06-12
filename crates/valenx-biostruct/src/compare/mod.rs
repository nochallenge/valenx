//! Structure comparison: pairwise alignment, TM-score and a
//! FoldSeek-class structural-search descriptor.
//!
//! - [`align`] — sequence-anchored iterative-superposition structure
//!   alignment (the v1 path), plus the TM-score formula.
//! - [`tmalign`] — TM-align-class **sequence-independent**
//!   iterative-DP structure aligner and a CE-style aligned-fragment-pair
//!   variant. The production default for sequence-divergent /
//!   structurally-similar pairs.
//! - [`foldseek`] — a 3Di-like structural alphabet for fast
//!   fold-similarity screening.

pub mod align;
pub mod foldseek;
pub mod tmalign;

pub use align::{align_chains, tm_d0, tm_score, StructureAlignment};
pub use foldseek::{
    best_ungapped_similarity, descriptor_identity, structural_descriptor, StructuralDescriptor,
};
pub use tmalign::{align_chains as align_chains_tm, align_chains_ce, TmAlignment, TmSeedKind};
