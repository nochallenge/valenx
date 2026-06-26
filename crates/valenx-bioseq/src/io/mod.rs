//! Sequence file I/O: FASTA, FASTQ, GenBank and EMBL readers/writers,
//! the Phred/Solexa quality codec, the shared INSDC location-string
//! parser, and native VCF/SAM parsing (via the pure-Rust `noodles`
//! crates, no C htslib).

pub mod embl;
pub mod fasta;
pub mod fastq;
pub mod genbank;
pub mod locstr;
pub mod noodles_formats;
pub mod quality;
