//! Tests, frequencies and quantitative genetics on diploid
//! [`crate::model::Population`]s.
//!
//! Where [`crate::stats`] operates on a haplotype genotype matrix,
//! `popstats` works on diploid *individuals* — the genotype counts,
//! equilibrium tests and quantitative-genetics quantities that need
//! genotype (not just allele) information.
//!
//! - [`hardy_weinberg`] — the chi-square and exact Hardy-Weinberg
//!   equilibrium tests.
//! - [`frequencies`] — allele / genotype frequencies, observed and
//!   expected heterozygosity, the inbreeding coefficient `F`.
//! - [`quant_genetics`] — breeding values, additive genetic variance,
//!   heritability and the breeder's equation.
//! - [`microsatellite`] — the stepwise mutation model and
//!   microsatellite summary statistics.

pub mod frequencies;
pub mod hardy_weinberg;
pub mod microsatellite;
pub mod quant_genetics;

pub use frequencies::{
    genotype_counts, locus_stats, mean_heterozygosity, mean_inbreeding, LocusStats,
};
pub use hardy_weinberg::{
    hwe_chi_square, hwe_exact, GenotypeCounts, HweResult,
};
pub use microsatellite::{
    allele_size_variance, m_ratio, mean_allele_size, microsat_heterozygosity,
    StepwiseMutationModel,
};
pub use quant_genetics::{
    narrow_sense_heritability, response_to_selection, AdditiveTrait,
};
