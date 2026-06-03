//! Read processing — FASTQ QC, trimming, filtering, dedup and depth.
//!
//! Everything that happens to sequencing reads between the instrument
//! and the variant caller:
//!
//! - [`qcstats`] — a FastQC-class [`qcstats::FastqcReport`]: per-base
//!   quality quartiles, per-read quality distribution, base
//!   composition, GC and length histograms.
//! - [`trim`] — 3′ adapter trimming with a mismatch-tolerant overlap
//!   scan, plus leading/trailing and sliding-window quality trimming.
//! - [`filter`] — whole-read removal on length, mean quality, `N`
//!   fraction and a linguistic-complexity score.
//! - [`dedup`] — coordinate-based PCR-duplicate marking (sets the SAM
//!   `DUPLICATE` flag) and sequence-identical de-duplication.
//! - [`coverage`] — per-base [`coverage::DepthProfile`] across a
//!   reference, with mean / median / breadth statistics.

pub mod coverage;
pub mod dedup;
pub mod filter;
pub mod qcstats;
pub mod trim;
