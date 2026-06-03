//! # valenx-genomics — NGS / variant tooling
//!
//! A native-Rust replacement
//! for the daily-driver portion of samtools, bcftools, VCFtools, the
//! alignment-based core of GATK and FreeBayes, the read simulators ART
//! / wgsim / DWGSIM, and the CRISPR tools CHOPCHOP / Cas-OFFinder /
//! CRISPResso2 — pure deterministic algorithms, no neural-network
//! weights and no external processes.
//!
//! It builds on [`valenx_bioseq`] (Block 6.1) for the
//! [`Seq`](valenx_bioseq::Seq) type and the FASTQ quality codec, and
//! on [`valenx_align`] (Block 6.2) for CIGAR conversion and the
//! seed-and-extend read mapper used by the CRISPR edit-outcome
//! analyser.
//!
//! ## What it does
//!
//! - **Formats** ([`mod@format`]) — SAM, VCF, BED, GFF3 / GTF parsers
//!   and writers, and a pileup engine that turns aligned reads into
//!   per-position evidence.
//! - **Read processing** ([`reads`]) — a FastQC-class QC report,
//!   adapter / quality trimming, read filtering, duplicate marking and
//!   per-base depth computation.
//! - **Variant calling** ([`variant`]) — a v1 pileup-based SNV and
//!   short-indel caller with a Bayesian diploid genotype-likelihood
//!   model, **plus** a GATK HaplotypeCaller-class local-haplotype-
//!   reassembly caller (active-region detection → De-Bruijn local
//!   assembly → quality-aware PairHMM → diploid marginalisation —
//!   see [`variant::haplotype`]); variant filtering / annotation,
//!   allele-frequency statistics and VCF normalisation.
//! - **Read simulation** ([`simulate`]) — an Illumina-style short-read
//!   simulator with a position-specific error model, a long-read
//!   (PacBio / Nanopore-style) simulator and paired-end generation.
//! - **CRISPR** ([`crispr`]) — guide-RNA design with PAM scanning and
//!   a Doench-style on-target score, off-target enumeration with a
//!   CFD-style score, and CRISPResso-class edit-outcome analysis.
//! - **Assembly** ([`assembly`]) — k-mer assembly statistics, a De
//!   Bruijn graph assembler and an overlap-layout-consensus
//!   mini-assembler.
//! - **Utilities** ([`util`]) — a deterministic PRNG and seeded
//!   FASTA / FASTQ subsampling.
//! - **Pipeline** ([`pipeline`]) — batch helpers and a bundled
//!   [`pipeline::GenomicsReport`].
//!
//! ## Errors
//!
//! Every fallible function returns
//! [`Result<_, GenomicsError>`](error::GenomicsError). The error type
//! carries stable [`code`](error::GenomicsError::code) and
//! [`category`](error::GenomicsError::category) accessors for telemetry.
//!
//! ## v1 scope
//!
//! This is a real working v1, not production parity with the 30-year
//! reference tools. Each module documents its own simplifications
//! inline; the notable ones are: the file formats are handled in their
//! **text** encodings only (no binary BAM / CRAM / BCF — that needs a
//! BGZF codec the Round-6 budget keeps off the dependency tree); the
//! v1 pileup variant caller is a real per-site Bayesian model and the
//! new haplotype caller in [`variant::haplotype`] closes the
//! GATK-class local-reassembly gap (single-sample, biallelic per
//! locus — multi-sample joint calling stays a documented gap); the
//! CRISPR on-target and off-target scores are the *published feature
//! weights* of Doench Rule-Set-2 and the CFD matrix, not
//! re-trained models; and the assemblers are correct graph-algorithm
//! v1s, not SPAdes / hifiasm at scale.

#![forbid(unsafe_code)]

pub mod assembly;
pub mod crispr;
pub mod error;
pub mod format;
pub mod pipeline;
pub mod reads;
pub mod simulate;
pub mod util;
pub mod variant;

// --- Convenience re-exports of the most-used types --------------------

pub use error::{ErrorCategory, GenomicsError, Result};
pub use format::bed::{BedFile, BedRecord};
pub use format::gff::{GffFile, GffRecord};
pub use format::pileup::{build_pileup, PileupColumn, Reference};
pub use format::sam::{Cigar, SamFile, SamRecord};
pub use format::vcf::{VcfFile, VcfRecord};
pub use pipeline::GenomicsReport;
pub use variant::call::{call_variants, Variant};
pub use variant::haplotype::{
    call_haplotype_variants, ActiveRegion, ActiveRegionParams, HaplotypeCallParams,
    LocalAssemblyParams, PairHmmParams, VariantCallMethod,
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crate_error_reexported() {
        let e = GenomicsError::not_yet("x");
        assert_eq!(e.category_enum(), ErrorCategory::Capability);
    }
}
