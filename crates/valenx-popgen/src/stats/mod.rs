//! Summary statistics over a [`crate::infer::GenotypeMatrix`] or a
//! [`crate::coalescent::TreeSequence`].
//!
//! These are the population-genetic summaries `VCFtools`,
//! `scikit-allel`, `libsequence` and `tskit` compute. Every function in
//! the site-mode submodules takes a genotype matrix (or the [`Sfs`]
//! derived from one), so it applies equally to a forward simulation, a
//! coalescent genealogy or imported VCF data. The
//! [`mod@tree_stats`] submodule adds the **tskit statistics framework**:
//! site-mode, branch-mode and per-window variants computed directly on
//! a tree sequence.
//!
//! - [`sfs`] — the folded / unfolded site-frequency spectrum.
//! - [`diversity`] — nucleotide diversity `pi`, Watterson's `theta`,
//!   Tajima's D, Fu & Li's D, Fay & Wu's H.
//! - [`fst`] — Hudson and Weir-Cockerham Fst (population
//!   differentiation).
//! - [`ld`] — linkage disequilibrium `D`, `D'`, `r^2` and the pairwise
//!   LD matrix.
//! - [`selection_scan`] — EHH, integrated EHH and iHS haplotype
//!   selection scans.
//! - [`tree_stats`] — windowed site/branch π, divergence, segregating
//!   sites over a [`crate::coalescent::TreeSequence`].

pub mod concordance;
pub mod diversity;
pub mod fst;
pub mod ld;
pub mod selection_scan;
pub mod sfs;
pub mod tree_stats;

pub use concordance::genotype_concordance;
pub use diversity::{
    expected_heterozygosity, fay_wu_h, fu_li_d, minor_allele_frequency, nucleotide_diversity,
    pairwise_differences, tajimas_d, wattersons_theta,
};
pub use fst::{fst_hudson, fst_weir_cockerham};
pub use ld::{ld_d, ld_d_prime, ld_matrix, ld_pair, ld_r_squared, LdStats};
pub use selection_scan::{ehh, ihs, integrated_ehh};
pub use sfs::{folded_spectrum, site_frequency_spectrum, Sfs};
pub use tree_stats::{
    branch_diversity, branch_divergence, equal_windows, windowed_branch_diversity,
    windowed_segregating_sites, windowed_site_divergence, windowed_site_diversity,
    WindowedStats,
};
