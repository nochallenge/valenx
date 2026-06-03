//! Allele / genotype frequencies, heterozygosity and inbreeding.
//!
//! Given a diploid [`crate::model::Population`] this module computes the
//! everyday descriptive statistics of a locus:
//!
//! - **Allele frequencies** — the derived-allele frequency `p`.
//! - **Genotype frequencies** — the observed `AA : Aa : aa`
//!   proportions.
//! - **Observed heterozygosity `Ho`** — the fraction of individuals
//!   that are heterozygous.
//! - **Expected heterozygosity `He`** (gene diversity) — `2p(1-p)`, the
//!   heterozygosity expected under Hardy-Weinberg.
//! - **Inbreeding coefficient `F`** — `(He - Ho) / He`, the relative
//!   heterozygote deficit; positive under inbreeding, negative under an
//!   excess of heterozygotes.
//!
//! Multi-locus averages ([`mean_heterozygosity`]) summarise a whole
//! population.

use crate::error::{PopgenError, Result};
use crate::model::Population;
use crate::popstats::hardy_weinberg::GenotypeCounts;

/// Descriptive statistics of one biallelic locus.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct LocusStats {
    /// Derived-allele frequency `p`.
    pub derived_frequency: f64,
    /// Observed heterozygosity `Ho`.
    pub observed_heterozygosity: f64,
    /// Expected heterozygosity `He = 2p(1-p)`.
    pub expected_heterozygosity: f64,
    /// Inbreeding coefficient `F = (He - Ho) / He`.
    pub inbreeding_coefficient: f64,
}

/// Tallies the diploid genotype counts at `site` in a population.
///
/// # Errors
/// [`PopgenError::Invalid`] if `site` is out of range or the population
/// is not diploid.
pub fn genotype_counts(pop: &Population, site: usize) -> Result<GenotypeCounts> {
    if pop.ploidy() != crate::model::Ploidy::Diploid {
        return Err(PopgenError::invalid(
            "population",
            "genotype counts require a diploid population",
        ));
    }
    if site >= pop.site_count() {
        return Err(PopgenError::invalid("site", "index out of range"));
    }
    let mut counts = GenotypeCounts {
        aa: 0,
        ab: 0,
        bb: 0,
    };
    for ind in pop.individuals() {
        match ind.genotype(site) {
            0 => counts.aa += 1,
            1 => counts.ab += 1,
            _ => counts.bb += 1,
        }
    }
    Ok(counts)
}

/// Computes the [`LocusStats`] for `site` in a population.
///
/// # Errors
/// [`PopgenError::Invalid`] if `site` is out of range, the population
/// is not diploid, or it is empty.
pub fn locus_stats(pop: &Population, site: usize) -> Result<LocusStats> {
    let counts = genotype_counts(pop, site)?;
    let n = counts.total();
    if n == 0 {
        return Err(PopgenError::invalid("population", "no individuals"));
    }
    let p = counts.derived_frequency();
    let ho = counts.ab as f64 / n as f64;
    let he = 2.0 * p * (1.0 - p);
    let f = if he < 1e-12 {
        0.0 // monomorphic locus: F is undefined, report 0
    } else {
        (he - ho) / he
    };
    Ok(LocusStats {
        derived_frequency: p,
        observed_heterozygosity: ho,
        expected_heterozygosity: he,
        inbreeding_coefficient: f,
    })
}

/// Mean observed and expected heterozygosity across every site of a
/// population: `(mean_Ho, mean_He)`.
///
/// A population with no sites returns `(0, 0)`.
///
/// # Errors
/// [`PopgenError::Invalid`] if the population is not diploid or empty.
pub fn mean_heterozygosity(pop: &Population) -> Result<(f64, f64)> {
    if pop.ploidy() != crate::model::Ploidy::Diploid {
        return Err(PopgenError::invalid(
            "population",
            "heterozygosity requires a diploid population",
        ));
    }
    if pop.size() == 0 {
        return Err(PopgenError::invalid("population", "no individuals"));
    }
    let s = pop.site_count();
    if s == 0 {
        return Ok((0.0, 0.0));
    }
    let mut sum_ho = 0.0;
    let mut sum_he = 0.0;
    for site in 0..s {
        let st = locus_stats(pop, site)?;
        sum_ho += st.observed_heterozygosity;
        sum_he += st.expected_heterozygosity;
    }
    Ok((sum_ho / s as f64, sum_he / s as f64))
}

/// Population-wide mean inbreeding coefficient `F`, averaged over the
/// polymorphic sites only (monomorphic sites have undefined `F`).
///
/// Returns `0.0` if there are no polymorphic sites.
///
/// # Errors
/// [`PopgenError::Invalid`] if the population is not diploid or empty.
pub fn mean_inbreeding(pop: &Population) -> Result<f64> {
    if pop.ploidy() != crate::model::Ploidy::Diploid {
        return Err(PopgenError::invalid(
            "population",
            "inbreeding requires a diploid population",
        ));
    }
    if pop.size() == 0 {
        return Err(PopgenError::invalid("population", "no individuals"));
    }
    let mut sum = 0.0;
    let mut polymorphic = 0usize;
    for site in 0..pop.site_count() {
        let st = locus_stats(pop, site)?;
        if st.expected_heterozygosity > 1e-12 {
            sum += st.inbreeding_coefficient;
            polymorphic += 1;
        }
    }
    Ok(if polymorphic == 0 {
        0.0
    } else {
        sum / polymorphic as f64
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Ploidy, Population, Site};

    /// Builds a diploid population with one site and the given diploid
    /// genotypes (`0`/`1`/`2` derived copies per individual).
    fn pop_with_genotypes(genotypes: &[u8]) -> Population {
        let mut pop =
            Population::founder(genotypes.len(), Ploidy::Diploid, 100.0).unwrap();
        let s = pop.add_site(Site::at(50.0));
        for (i, &g) in genotypes.iter().enumerate() {
            let ind = &mut pop.individuals_mut()[i];
            if g >= 1 {
                ind.genomes[0].set_derived(s);
            }
            if g >= 2 {
                ind.genomes[1].set_derived(s);
            }
        }
        pop
    }

    #[test]
    fn genotype_counts_tally_correctly() {
        // 2 AA, 3 Aa, 1 aa.
        let pop = pop_with_genotypes(&[0, 0, 1, 1, 1, 2]);
        let c = genotype_counts(&pop, 0).unwrap();
        assert_eq!(c.aa, 2);
        assert_eq!(c.ab, 3);
        assert_eq!(c.bb, 1);
    }

    #[test]
    fn hwe_population_has_f_near_zero() {
        // p = 0.5, HWE genotypes 1 AA : 2 Aa : 1 aa.
        let pop = pop_with_genotypes(&[0, 1, 1, 2]);
        let st = locus_stats(&pop, 0).unwrap();
        assert!((st.derived_frequency - 0.5).abs() < 1e-12);
        // Ho = 2/4 = 0.5; He = 2*0.5*0.5 = 0.5; F = 0.
        assert!((st.observed_heterozygosity - 0.5).abs() < 1e-12);
        assert!((st.expected_heterozygosity - 0.5).abs() < 1e-12);
        assert!(st.inbreeding_coefficient.abs() < 1e-12);
    }

    #[test]
    fn heterozygote_deficit_gives_positive_f() {
        // p = 0.5 but only homozygotes -> Ho = 0 -> F = 1.
        let pop = pop_with_genotypes(&[0, 0, 2, 2]);
        let st = locus_stats(&pop, 0).unwrap();
        assert!((st.observed_heterozygosity).abs() < 1e-12);
        assert!((st.inbreeding_coefficient - 1.0).abs() < 1e-12);
    }

    #[test]
    fn heterozygote_excess_gives_negative_f() {
        // p = 0.5, all heterozygous -> Ho = 1 > He = 0.5 -> F < 0.
        let pop = pop_with_genotypes(&[1, 1, 1, 1]);
        let st = locus_stats(&pop, 0).unwrap();
        assert!(st.inbreeding_coefficient < 0.0);
    }

    #[test]
    fn mean_heterozygosity_averages_sites() {
        let pop = pop_with_genotypes(&[0, 1, 1, 2]);
        let (ho, he) = mean_heterozygosity(&pop).unwrap();
        assert!((ho - 0.5).abs() < 1e-12);
        assert!((he - 0.5).abs() < 1e-12);
    }

    #[test]
    fn monomorphic_site_has_zero_diversity() {
        let pop = pop_with_genotypes(&[0, 0, 0, 0]);
        let st = locus_stats(&pop, 0).unwrap();
        assert_eq!(st.expected_heterozygosity, 0.0);
        assert_eq!(st.inbreeding_coefficient, 0.0);
        // Monomorphic sites are skipped by mean_inbreeding.
        assert_eq!(mean_inbreeding(&pop).unwrap(), 0.0);
    }

    #[test]
    fn rejects_out_of_range_site() {
        let pop = pop_with_genotypes(&[0, 1]);
        assert!(locus_stats(&pop, 9).is_err());
    }
}
