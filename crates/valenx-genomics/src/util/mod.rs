//! Cross-cutting utilities for `valenx-genomics`.
//!
//! - [`rng`] — a small deterministic `SplitMix64` pseudo-random
//!   generator, the reproducible randomness behind the read simulators
//!   and the subsampling helpers (no `rand`-crate dependency).
//! - [`subsample`] — random and seeded FASTQ / FASTA subsampling and
//!   downsampling.

pub mod rng;
pub mod subsample;
