//! Read simulation — Illumina, long-read and paired-end generators.
//!
//! Synthetic reads with known truth are the backbone of benchmarking a
//! mapper, a variant caller or an assembler. This module replaces the
//! daily-driver read simulators:
//!
//! - [`illumina`] — an ART-class Illumina short-read simulator with a
//!   position-specific substitution-error model and a per-cycle
//!   quality profile.
//! - [`longread`] — a pbsim / Badread-class long-read simulator with
//!   an indel-heavy PacBio / Nanopore error model and a broad
//!   read-length distribution.
//! - [`paired`] — paired-end generation with a configurable
//!   [`paired::InsertSizeModel`] in Illumina FR orientation.
//!
//! Every routine is seeded through [`crate::util::rng`] — the same
//! `(seed, input)` always produces the same reads.

pub mod illumina;
pub mod longread;
pub mod paired;
