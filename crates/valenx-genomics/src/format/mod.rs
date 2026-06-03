//! Genomic file formats — parsers, writers and the pileup engine.
//!
//! Every NGS workflow is glued together by a handful of text formats.
//! This module implements the daily-driver set, each as a typed
//! data-model plus a round-trippable parser and writer:
//!
//! - [`sam`] — SAM alignment records ([`sam::SamFile`]) with full
//!   CIGAR, flag-bit and optional-tag handling.
//! - [`vcf`] — VCF variant records ([`vcf::VcfFile`]) with typed
//!   `##INFO` / `##FORMAT` metadata and per-sample genotype columns.
//! - [`bed`] — BED intervals ([`bed::BedFile`]) across the full BED12
//!   column set, with merge / overlap utilities.
//! - [`gff`] — GFF3 and GTF feature annotation ([`gff::GffFile`]) with
//!   dialect auto-detection.
//! - [`pileup`] — turns aligned [`sam::SamRecord`]s into per-position
//!   [`pileup::PileupColumn`]s and renders the samtools text-pileup
//!   format.
//!
//! ## v1 scope
//!
//! All five formats are handled in their **text** encodings. The
//! binary siblings — BAM and CRAM (alignment) and BCF (variants) — are
//! out of v1 scope: they need a BGZF / DEFLATE codec, and the Round-6
//! budget deliberately keeps a compression dependency off the
//! workspace tree. The text parsers cover the formats every analysis
//! pipeline reads and writes day to day.

pub mod bed;
pub mod gff;
pub mod pileup;
pub mod sam;
pub mod vcf;
