//! # valenx-bio
//!
//! Canonical types + format I/O for Valenx's biology + biotech
//! adapters. Three canonical types live here:
//!
//! - [`sequence::Sequence`] — DNA / RNA / protein sequences with
//!   IUPAC alphabet validation.
//! - [`structure::Structure`] — atomic / residue / chain hierarchy
//!   for proteins, nucleic acids, small molecules.
//! - [`trajectory::Trajectory`] — per-frame atomic coordinates from
//!   MD output (DCD / XTC / TRR).
//!
//! Format readers + writers live under [`mod@format`]: FASTA, PDB,
//! mmCIF, DCD. Each adapter pulls just the formats it produces or
//! consumes, keeping per-adapter dep weight small.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod alignment;
pub mod alphabet;
pub mod fastq;
pub mod format;
pub mod sequence;
pub mod structure;
pub mod trajectory;
pub mod vcf;

pub use alignment::{Alignment, AlignmentError};
pub use fastq::FastqRecord;
pub use sequence::{Alphabet, Sequence};
pub use structure::{Atom, Chain, Residue, Structure};
pub use trajectory::Trajectory;
pub use vcf::{Vcf, VcfError, VcfRecord};
