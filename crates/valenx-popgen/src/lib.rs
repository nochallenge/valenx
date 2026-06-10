//! # valenx-popgen — population genetics and evolution simulators
//!
//! A native-Rust replacement for
//! SLiM, msprime / tskit, fwdpy11, simuPOP, Nemo, discoal, `ms` and
//! stdpopsim — pure algorithms, no neural-network weights and no
//! external processes.
//!
//! It builds on [`valenx_bioseq`] (Block 6.1) for the nucleotide
//! alphabet and on [`valenx_phylo`] (Block 6.3) for the genealogy data
//! model: a coalescent genealogy *is* a [`valenx_phylo::Tree`], so the
//! whole of `valenx-phylo`'s rendering and comparison machinery applies
//! to a simulated tree for free.
//!
//! ## What it does
//!
//! - **Model** ([`model`], [`rng`]) — the [`Population`] /
//!   [`Individual`] / [`Genome`] data model and a deterministic
//!   seedable PCG random-number generator.
//! - **Forward simulation** ([`forward`]) — a discrete-generation
//!   diploid Wright-Fisher simulator with additive / multiplicative /
//!   epistatic selection, infinite- and finite-sites mutation,
//!   crossover and gene-conversion recombination, island and
//!   stepping-stone migration, piecewise demography, and
//!   `pyslim`/`tskit`-style tree recording.
//! - **Coalescent** ([`mod@coalescent`]) — the backward Kingman
//!   coalescent, the structured (multi-deme) coalescent, the
//!   coalescent with recombination (an ancestral recombination
//!   graph), the succinct `tskit`-class [`coalescent::TreeSequence`]
//!   and mutation overlay onto a genealogy.
//! - **Summary statistics** ([`stats`]) — the site-frequency spectrum,
//!   nucleotide diversity, Watterson's theta, Tajima's D, Fu & Li's D,
//!   Fay & Wu's H, Hudson and Weir-Cockerham Fst, linkage
//!   disequilibrium and EHH / iHS selection scans.
//! - **Tests & quant-gen** ([`popstats`]) — the Hardy-Weinberg
//!   equilibrium tests, allele / genotype frequencies, heterozygosity
//!   and inbreeding, quantitative-genetics breeding values and
//!   heritability, and the stepwise microsatellite mutation model.
//! - **Inference & catalog** ([`infer`], [`catalog`], [`io`]) — the
//!   Approximate Bayesian Computation rejection framework, a
//!   `stdpopsim`-class species catalog, and VCF / `ms` / Newick
//!   export of simulations.
//!
//! ## Errors
//!
//! Every fallible function returns
//! [`Result<_, PopgenError>`](error::PopgenError); the error type
//! carries stable [`code`](error::PopgenError::code) and
//! [`category`](error::PopgenError::category) accessors for telemetry.
//!
//! ## Scope
//!
//! Each module documents its own simplifications inline. The
//! Wright-Fisher simulator uses discrete non-overlapping generations
//! and viability-only selection (no age structure, no separate
//! sexes). The forward-time tree-recorder is `pyslim`-class — it
//! records every crossover per meiosis as a separate edge and runs a
//! drop-unreachable + unary-chain-squash simplify pass on the table.
//! The ARG runs Hudson's (1983) algorithm with a segment-list
//! ancestral-material representation and a piecewise-constant
//! recombination-rate map ([`coalescent::RecombinationMap`]) for hot-
//! and cold-spot modelling. The summary-statistic module exposes both
//! site-mode (mutation-based) and branch-mode (`tskit`-class)
//! per-window estimators via [`stats::tree_stats`]. Fst's
//! Weir-Cockerham estimator falls back to expected heterozygosity
//! when individual phase is unavailable; iHS is reported
//! un-genome-standardised; the `stdpopsim`-class catalog uses
//! representative round parameters rather than the precise curated
//! published estimates; and the random-number generator is a small
//! deterministic PCG so simulations are reproducible but not
//! cryptographic.

#![forbid(unsafe_code)]

pub mod catalog;
pub mod coalescent;
pub mod error;
pub mod forward;
pub mod infer;
pub mod io;
pub mod model;
pub mod popstats;
pub mod rng;
pub mod stats;

// --- Convenience re-exports of the most-used types --------------------

pub use error::{ErrorCategory, PopgenError, Result};
pub use model::{Genome, Individual, Ploidy, Population, Site};
pub use rng::Rng;

pub use catalog::{Catalog, DemographicModel, SpeciesModel};
pub use coalescent::{
    coalescent, overlay_mutations, overlay_on_tree, simulate_arg,
    structured_coalescent, ArgParams, PopHistory, RecombinationMap, TreeSequence,
};
pub use forward::{
    record_wright_fisher, DemographicSchedule, Drift, MigrationModel,
    MutationModel, RecombinationModel, SelectionModel, SimulationConfig,
    WrightFisher,
};
pub use infer::{abc_reject, AbcConfig, AbcPosterior, GenotypeMatrix, Prior};
pub use io::{read_ms, write_ms, write_newick_genealogy, write_vcf};
pub use stats::{
    branch_diversity, branch_divergence, equal_windows, expected_heterozygosity, fst_hudson,
    fst_weir_cockerham, genotype_concordance, minor_allele_frequency, nucleotide_diversity,
    site_frequency_spectrum,
    tajimas_d, wattersons_theta, windowed_branch_diversity,
    windowed_segregating_sites, windowed_site_divergence,
    windowed_site_diversity, WindowedStats,
};

#[cfg(test)]
mod tests {
    use super::*;

    /// End-to-end: simulate a coalescent genealogy, overlay mutations,
    /// and recover diversity statistics from the resulting matrix.
    #[test]
    fn coalescent_to_statistics_end_to_end() {
        let labels: Vec<String> = (0..20).map(|i| format!("s{i}")).collect();
        let tree = coalescent(&labels, &PopHistory::Constant(1000.0), 42).unwrap();
        assert_eq!(tree.leaf_count(), 20);

        let gm = overlay_on_tree(&tree, 2e-3, 7).unwrap();
        assert!(gm.n_sites() > 0, "no segregating sites generated");

        // Diversity statistics are all finite and sensible.
        let pi = nucleotide_diversity(&gm).unwrap();
        let theta = wattersons_theta(&gm).unwrap();
        let d = tajimas_d(&gm).unwrap();
        assert!(pi >= 0.0);
        assert!(theta >= 0.0);
        assert!(d.is_finite());

        // The SFS sums to the segregating-site count.
        let sfs = site_frequency_spectrum(&gm).unwrap();
        assert_eq!(sfs.segregating_sites(), gm.segregating_sites());
    }

    /// End-to-end: a forward Wright-Fisher run yields a population
    /// whose genotype matrix exports to a valid VCF.
    #[test]
    fn forward_to_vcf_end_to_end() {
        let cfg =
            SimulationConfig::neutral(20, 30, 3e-3, 1e-4, 1000.0, 1).unwrap();
        let pop = WrightFisher::new(cfg).unwrap().run().unwrap();
        let gm = pop.genotype_matrix();
        // 20 diploid individuals -> 40 haplotype rows -> valid VCF.
        let vcf = write_vcf(&gm, "sim").unwrap();
        assert!(vcf.contains("##fileformat=VCFv4.2"));
    }
}
