//! Assembly — contig statistics and two genome assemblers.
//!
//! - [`stats`] — QUAST-class assembly statistics
//!   ([`stats::assembly_stats`]): N50 / L50 / N75 / N90, total size,
//!   GC, plus a canonical k-mer spectrum for genome-size estimation.
//! - [`debruijn`] — a De Bruijn graph assembler
//!   ([`debruijn::assemble`]) for short reads: build the graph, clip
//!   tips, pop simple bubbles, emit unitig contigs.
//! - [`olc`] — an overlap-layout-consensus mini-assembler
//!   ([`olc::assemble_olc`]) for long reads: all-pairs suffix-prefix
//!   overlap, greedy layout, majority consensus.
//!
//! ## v1 scope
//!
//! Both assemblers are correct, deterministic graph-algorithm v1s —
//! not SPAdes / hifiasm at genome scale. The De Bruijn assembler
//! builds its graph in memory and handles simple bubbles only; the OLC
//! assembler's overlap step is an O(n²) all-pairs scan. See each
//! module's own note.

pub mod debruijn;
pub mod olc;
pub mod stats;
