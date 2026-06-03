//! Variant calling — pileup caller, haplotype caller, genotyping,
//! filtering, statistics.
//!
//! The variant-discovery half of an NGS pipeline:
//!
//! - [`genotype`] — a Bayesian diploid genotype-likelihood model
//!   ([`genotype::genotype_site`]) over per-base error probabilities.
//! - [`call`] — the v1 pileup-based SNV and short-indel caller
//!   ([`call::call_variants`]) that tallies alleles, applies depth /
//!   fraction / quality gates and genotypes each passing site.
//! - [`haplotype`] — a **GATK HaplotypeCaller-class** local-
//!   reassembly variant caller
//!   ([`haplotype::call_haplotype_variants`]) — active-region
//!   detection, local De-Bruijn-graph haplotype assembly,
//!   quality-aware PairHMM read-vs-haplotype likelihoods, diploid
//!   genotype marginalisation. The new high-stakes default; the
//!   pileup caller stays available via
//!   [`haplotype::VariantCallMethod`].
//! - [`filter`] — quality-gate filtering and annotation, including a
//!   strand-bias score, plus conversion to [`crate::format::vcf`]
//!   records.
//! - [`stats`] — allele-frequency, transition/transversion,
//!   genotype-class and Hardy-Weinberg statistics over a VCF.
//! - [`normalize`] — VCF normalisation: decompose multiallelics, trim
//!   redundant bases and left-align indels.
//!
//! ## v1 scope
//!
//! The pileup caller is a real per-site model — pileup tally,
//! threshold, Bayesian genotype likelihood — and stays the lightest
//! option; the haplotype caller closes the GATK-class gap (local
//! De-Bruijn assembly + PairHMM + marginalisation). Single-sample
//! biallelic calling is supported in both; multi-sample joint calling
//! and proper multi-allelic representation stay documented gaps. See
//! the inline notes in each module.

pub mod call;
pub mod filter;
pub mod genotype;
pub mod haplotype;
pub mod normalize;
pub mod stats;
