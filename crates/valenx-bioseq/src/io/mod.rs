//! Sequence file I/O: FASTA, FASTQ, GenBank and EMBL readers/writers,
//! the Phred/Solexa quality codec, and the shared INSDC
//! location-string parser.

pub mod embl;
pub mod fasta;
pub mod fastq;
pub mod genbank;
pub mod locstr;
pub mod quality;
