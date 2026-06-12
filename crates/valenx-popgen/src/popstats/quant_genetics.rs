//! Quantitative genetics: breeding values, genetic variance,
//! heritability.
//!
//! Most traits of interest are *quantitative* — controlled by many
//! loci of small effect plus environment. This module supplies the
//! core quantitative-genetics machinery on top of `valenx-popgen`'s
//! diploid genotypes.
//!
//! - **Breeding value** ([`AdditiveTrait::breeding_value`]) — the sum
//!   of an individual's additive allelic effects, the heritable part
//!   of its genotype. With per-allele effect `a` at a locus, the three
//!   genotypes contribute `0`, `a`, `2a`.
//! - **Additive genetic variance `Va`**
//!   ([`AdditiveTrait::additive_variance`]) — the variance of breeding
//!   values across the population: `Va = sum_l 2 p_l q_l a_l^2`.
//! - **Phenotype** ([`AdditiveTrait::phenotype`]) — breeding value
//!   plus a Gaussian environmental deviation.
//! - **Heritability `h^2`** ([`narrow_sense_heritability`]) — the
//!   fraction of phenotypic variance that is additive-genetic,
//!   `Va / (Va + Ve)`.
//! - **Breeder's equation** ([`response_to_selection`]) — predicted
//!   response `R = h^2 * S` to a selection differential `S`.

use crate::error::{PopgenError, Result};
use crate::model::Population;
use crate::rng::Rng;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// An additive quantitative trait: per-locus allelic effects plus an
/// environmental-variance parameter.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct AdditiveTrait {
    /// site index -> additive effect `a` of one derived-allele copy.
    effects: BTreeMap<usize, f64>,
    /// Environmental variance `Ve` added to the breeding value to give
    /// a phenotype.
    environmental_variance: f64,
}

impl AdditiveTrait {
    /// A trait with no loci and the given environmental variance.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] if `environmental_variance` is
    /// negative.
    pub fn new(environmental_variance: f64) -> Result<Self> {
        if environmental_variance < 0.0 {
            return Err(PopgenError::invalid(
                "environmental_variance",
                "must be non-negative",
            ));
        }
        Ok(AdditiveTrait {
            effects: BTreeMap::new(),
            environmental_variance,
        })
    }

    /// Sets the additive effect of one derived-allele copy at `site`.
    pub fn set_effect(&mut self, site: usize, a: f64) -> &mut Self {
        self.effects.insert(site, a);
        self
    }

    /// The environmental variance `Ve`.
    pub fn environmental_variance(&self) -> f64 {
        self.environmental_variance
    }

    /// The breeding value of one individual: the sum over loci of
    /// `effect * derived-copy-count`.
    pub fn breeding_value(&self, individual: &crate::model::Individual) -> f64 {
        self.effects
            .iter()
            .map(|(&site, &a)| a * individual.genotype(site) as f64)
            .sum()
    }

    /// Mean breeding value across a population.
    pub fn mean_breeding_value(&self, pop: &Population) -> f64 {
        if pop.size() == 0 {
            return 0.0;
        }
        let sum: f64 = pop
            .individuals()
            .iter()
            .map(|ind| self.breeding_value(ind))
            .sum();
        sum / pop.size() as f64
    }

    /// The additive genetic variance `Va` of the trait in a
    /// population, computed directly as the variance of breeding
    /// values.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] if the population is empty.
    pub fn additive_variance(&self, pop: &Population) -> Result<f64> {
        if pop.size() == 0 {
            return Err(PopgenError::invalid("population", "no individuals"));
        }
        let mean = self.mean_breeding_value(pop);
        let var: f64 = pop
            .individuals()
            .iter()
            .map(|ind| (self.breeding_value(ind) - mean).powi(2))
            .sum::<f64>()
            / pop.size() as f64;
        Ok(var)
    }

    /// The theoretical additive variance under Hardy-Weinberg:
    /// `Va = sum_l 2 p_l q_l a_l^2`.
    ///
    /// This is the textbook decomposition; it equals
    /// [`additive_variance`](Self::additive_variance) when the
    /// population is in HWE and loci are unlinked.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] if the population is empty.
    pub fn theoretical_additive_variance(&self, pop: &Population) -> Result<f64> {
        if pop.size() == 0 {
            return Err(PopgenError::invalid("population", "no individuals"));
        }
        let mut va = 0.0;
        for (&site, &a) in &self.effects {
            if site < pop.site_count() {
                let p = pop.allele_frequency(site)?;
                va += 2.0 * p * (1.0 - p) * a * a;
            }
        }
        Ok(va)
    }

    /// A phenotype: the individual's breeding value plus a Gaussian
    /// environmental deviation with variance `Ve`.
    pub fn phenotype(&self, individual: &crate::model::Individual, rng: &mut Rng) -> f64 {
        let env = rng.normal() * self.environmental_variance.sqrt();
        self.breeding_value(individual) + env
    }
}

/// Narrow-sense heritability `h^2 = Va / (Va + Ve)`.
///
/// `0` means the trait variance is entirely environmental, `1` entirely
/// additive-genetic. Returns `0` when total variance is zero.
///
/// # Errors
/// [`PopgenError::Invalid`] on a negative variance.
pub fn narrow_sense_heritability(va: f64, ve: f64) -> Result<f64> {
    if va < 0.0 || ve < 0.0 {
        return Err(PopgenError::invalid("variance", "must be non-negative"));
    }
    let total = va + ve;
    Ok(if total < 1e-12 { 0.0 } else { va / total })
}

/// The breeder's equation: predicted response to selection
/// `R = h^2 * S` for a selection differential `S`.
///
/// # Errors
/// [`PopgenError::Invalid`] if `h2` is outside `[0, 1]`.
pub fn response_to_selection(h2: f64, selection_differential: f64) -> Result<f64> {
    if !(0.0..=1.0).contains(&h2) {
        return Err(PopgenError::invalid(
            "h2",
            "heritability must lie in [0, 1]",
        ));
    }
    Ok(h2 * selection_differential)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Ploidy, Population, Site};

    /// Diploid population with one site and the given diploid genotype
    /// dosages (`0`/`1`/`2`).
    fn pop(genotypes: &[u8]) -> Population {
        let mut p = Population::founder(genotypes.len(), Ploidy::Diploid, 100.0).unwrap();
        let s = p.add_site(Site::at(50.0));
        for (i, &g) in genotypes.iter().enumerate() {
            let ind = &mut p.individuals_mut()[i];
            if g >= 1 {
                ind.genomes[0].set_derived(s);
            }
            if g >= 2 {
                ind.genomes[1].set_derived(s);
            }
        }
        p
    }

    #[test]
    fn breeding_value_sums_allelic_effects() {
        let population = pop(&[0, 1, 2]);
        let mut trait_def = AdditiveTrait::new(0.0).unwrap();
        trait_def.set_effect(0, 2.5);
        // genotype 0 -> 0; 1 -> 2.5; 2 -> 5.0.
        assert_eq!(trait_def.breeding_value(&population.individuals()[0]), 0.0);
        assert!((trait_def.breeding_value(&population.individuals()[1]) - 2.5).abs() < 1e-12);
        assert!((trait_def.breeding_value(&population.individuals()[2]) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn additive_variance_matches_the_theoretical_formula_under_hwe() {
        // HWE genotypes for p = 0.5: 1 AA : 2 Aa : 1 aa.
        let population = pop(&[0, 1, 1, 2]);
        let mut trait_def = AdditiveTrait::new(0.0).unwrap();
        trait_def.set_effect(0, 1.0);
        let empirical = trait_def.additive_variance(&population).unwrap();
        let theoretical = trait_def
            .theoretical_additive_variance(&population)
            .unwrap();
        // 2pq a^2 = 2*0.5*0.5*1 = 0.5.
        assert!((theoretical - 0.5).abs() < 1e-12);
        assert!(
            (empirical - theoretical).abs() < 1e-9,
            "empirical {empirical} vs theoretical {theoretical}"
        );
    }

    #[test]
    fn monomorphic_locus_has_zero_additive_variance() {
        let population = pop(&[0, 0, 0, 0]);
        let mut trait_def = AdditiveTrait::new(0.0).unwrap();
        trait_def.set_effect(0, 3.0);
        assert_eq!(
            trait_def
                .theoretical_additive_variance(&population)
                .unwrap(),
            0.0
        );
    }

    #[test]
    fn heritability_is_a_proportion() {
        assert!((narrow_sense_heritability(1.0, 1.0).unwrap() - 0.5).abs() < 1e-12);
        assert!((narrow_sense_heritability(3.0, 1.0).unwrap() - 0.75).abs() < 1e-12);
        assert_eq!(narrow_sense_heritability(0.0, 0.0).unwrap(), 0.0);
        assert!(narrow_sense_heritability(-1.0, 1.0).is_err());
    }

    #[test]
    fn breeders_equation_predicts_response() {
        // h^2 = 0.4, S = 10 -> R = 4.
        assert!((response_to_selection(0.4, 10.0).unwrap() - 4.0).abs() < 1e-12);
        assert!(response_to_selection(1.5, 10.0).is_err());
    }

    #[test]
    fn phenotype_equals_breeding_value_with_zero_environment() {
        let population = pop(&[2]);
        let mut trait_def = AdditiveTrait::new(0.0).unwrap();
        trait_def.set_effect(0, 1.0);
        let mut rng = Rng::new(1);
        let ph = trait_def.phenotype(&population.individuals()[0], &mut rng);
        assert!((ph - 2.0).abs() < 1e-12);
    }

    #[test]
    fn rejects_negative_environmental_variance() {
        assert!(AdditiveTrait::new(-1.0).is_err());
    }
}
