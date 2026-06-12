//! Utilities — edit distance, dot-plots, SAM conversion, read mapping.
//!
//! The miscellaneous tools that round out the crate:
//!
//! - [`editdist`] — [`editdist::levenshtein`] edit distance and
//!   Myers' [`bit-parallel`](editdist::myers_bit_parallel) algorithm.
//! - [`dotplot`] — [`dotplot::dot_plot`] windowed-identity
//!   self-alignment matrices.
//! - [`sam`] — SAM [`crate::pairwise::Cigar`] ⇄
//!   [`crate::pairwise::Alignment`] conversion and the
//!   [`sam::SamRecord`] type.
//! - [`mapper`] — [`mapper::ReadMapper`], a seed-and-extend short-read
//!   aligner that emits SAM.

pub mod dotplot;
pub mod editdist;
pub mod mapper;
pub mod sam;

pub use dotplot::{dot_plot, self_dot_plot, DiagonalRun, DotPlot};
pub use editdist::{levenshtein, levenshtein_bounded, myers_bit_parallel};
pub use mapper::{
    mapping_quality, InsertSizeModel, MapperParams, MappingResult, PairedMappingResult, ReadMapper,
    Strand,
};
pub use sam::{alignment_from_cigar, alignment_to_cigar, cigar_to_rows, SamRecord};
