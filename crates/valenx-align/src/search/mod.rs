//! Homology search and indexing.
//!
//! The search family covers the indexing structures and heuristics
//! that replace BLAST+, BWA, Bowtie2 and minimap2:
//!
//! - [`kmer`] — a [`kmer::KmerIndex`] inverted index from k-mers to
//!   positions, over one or many sequences.
//! - [`seed`] — [`seed::SeedSearch`], the BLAST-class seed-and-extend
//!   heuristic with diagonal binning and X-drop extension, reporting
//!   scored [`seed::Hsp`] segment pairs.
//! - [`stats`] — [`stats::KarlinAltschul`] E-value and bit-score
//!   statistics.
//! - [`fmindex`] — [`fmindex::FmIndex`], a Burrows-Wheeler / FM-index
//!   for exact substring search (the BWA-class core).
//! - [`minimizer`] — `(k, w)`-[`minimizer::Minimizer`] sketches (the
//!   minimap2-class seed sample).
//! - [`chain`] — colinear DP [`chain::chain_anchors`] chaining of seed
//!   anchors.

pub mod chain;
pub mod fmindex;
pub mod kmer;
pub mod minimizer;
pub mod seed;
pub mod stats;

pub use chain::{chain_anchors, Anchor, Chain, ChainParams};
pub use fmindex::{FmIndex, Smem};
pub use kmer::{KmerHit, KmerIndex};
pub use minimizer::{minimizer_sketch, Minimizer};
pub use seed::{Hsp, SeedParams, SeedSearch};
pub use stats::KarlinAltschul;
