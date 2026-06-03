//! The stepwise mutation model — microsatellite evolution.
//!
//! Microsatellites (short tandem repeats) do not follow the
//! infinite-sites model: a mutation typically adds or removes a single
//! repeat unit, so alleles are *integers* (repeat counts) and mutation
//! takes a step of +-1 (Ohta & Kimura 1973, the **stepwise mutation
//! model**, SMM).
//!
//! [`StepwiseMutationModel`] simulates a population of repeat-count
//! alleles forward in time under drift + stepwise mutation, and the
//! module provides the microsatellite-specific summary statistics:
//!
//! - **Allele-size variance** — the variance of repeat counts, the
//!   quantity SMM-based distances (`(delta-mu)^2`) are built on.
//! - **Expected heterozygosity** under SMM equilibrium.
//! - **Garza-Williamson's M ratio** ([`m_ratio`]) — the number of
//!   distinct alleles over the allele-size range; a bottleneck drops M
//!   because rare size classes are lost faster than the range shrinks.

use crate::error::{PopgenError, Result};
use crate::rng::Rng;

/// A stepwise (single-step) microsatellite mutation model.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct StepwiseMutationModel {
    /// Per-allele per-generation mutation rate.
    pub mutation_rate: f64,
    /// Smallest allowed repeat count (alleles cannot mutate below it).
    pub min_repeats: i32,
    /// Largest allowed repeat count.
    pub max_repeats: i32,
}

impl StepwiseMutationModel {
    /// A model with the given rate and `[min, max]` repeat bounds.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] on a rate outside `[0, 1]` or
    /// `min >= max`.
    pub fn new(mutation_rate: f64, min_repeats: i32, max_repeats: i32) -> Result<Self> {
        if !(0.0..=1.0).contains(&mutation_rate) {
            return Err(PopgenError::invalid(
                "mutation_rate",
                "must lie in [0, 1]",
            ));
        }
        if min_repeats >= max_repeats {
            return Err(PopgenError::invalid(
                "repeat bounds",
                "min_repeats must be below max_repeats",
            ));
        }
        Ok(StepwiseMutationModel {
            mutation_rate,
            min_repeats,
            max_repeats,
        })
    }

    /// Applies one mutation step to a single allele: with probability
    /// `mutation_rate` it steps +-1 (each direction equally likely),
    /// clamped to the `[min, max]` bounds.
    pub fn mutate_allele(&self, allele: i32, rng: &mut Rng) -> i32 {
        if !rng.bernoulli(self.mutation_rate) {
            return allele;
        }
        let step = if rng.bernoulli(0.5) { 1 } else { -1 };
        (allele + step).clamp(self.min_repeats, self.max_repeats)
    }

    /// Simulates a Wright-Fisher population of `n` microsatellite
    /// alleles for `generations`, starting all at `start_repeats`.
    ///
    /// Each generation: resample `n` alleles with replacement from the
    /// previous generation (drift), then apply stepwise mutation. The
    /// returned vector is the final generation's repeat counts.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] if `n == 0` or `start_repeats` is
    /// outside the bounds.
    pub fn simulate(
        &self,
        n: usize,
        generations: usize,
        start_repeats: i32,
        seed: u64,
    ) -> Result<Vec<i32>> {
        if n == 0 {
            return Err(PopgenError::invalid("n", "population must be non-empty"));
        }
        if start_repeats < self.min_repeats || start_repeats > self.max_repeats {
            return Err(PopgenError::invalid(
                "start_repeats",
                "outside the [min, max] bounds",
            ));
        }
        let mut rng = Rng::new(seed);
        let mut alleles = vec![start_repeats; n];
        for _ in 0..generations {
            let mut next = Vec::with_capacity(n);
            for _ in 0..n {
                let parent = alleles[rng.below(n)];
                next.push(self.mutate_allele(parent, &mut rng));
            }
            alleles = next;
        }
        Ok(alleles)
    }
}

/// Mean repeat count of an allele sample.
pub fn mean_allele_size(alleles: &[i32]) -> f64 {
    if alleles.is_empty() {
        return 0.0;
    }
    alleles.iter().map(|&a| a as f64).sum::<f64>() / alleles.len() as f64
}

/// Variance of repeat counts — the basis of SMM genetic distances.
pub fn allele_size_variance(alleles: &[i32]) -> f64 {
    if alleles.is_empty() {
        return 0.0;
    }
    let mean = mean_allele_size(alleles);
    alleles
        .iter()
        .map(|&a| (a as f64 - mean).powi(2))
        .sum::<f64>()
        / alleles.len() as f64
}

/// Expected heterozygosity of an allele sample: `1 - sum p_i^2` over
/// the distinct repeat-count classes (Nei's gene diversity).
pub fn microsat_heterozygosity(alleles: &[i32]) -> f64 {
    use std::collections::HashMap;
    if alleles.is_empty() {
        return 0.0;
    }
    let n = alleles.len() as f64;
    let mut freq: HashMap<i32, usize> = HashMap::new();
    for &a in alleles {
        *freq.entry(a).or_insert(0) += 1;
    }
    let homozygosity: f64 = freq
        .values()
        .map(|&c| {
            let p = c as f64 / n;
            p * p
        })
        .sum();
    1.0 - homozygosity
}

/// The Garza-Williamson M ratio: the number of distinct alleles
/// divided by the overall allele-size range (`max - min + 1`).
///
/// A recent bottleneck drives M below ~0.68 because intermediate
/// size classes are lost while the size range persists.
///
/// Returns `0.0` for an empty sample.
pub fn m_ratio(alleles: &[i32]) -> f64 {
    use std::collections::HashSet;
    if alleles.is_empty() {
        return 0.0;
    }
    let distinct: HashSet<i32> = alleles.iter().copied().collect();
    let min = *alleles.iter().min().unwrap();
    let max = *alleles.iter().max().unwrap();
    let range = (max - min + 1) as f64;
    distinct.len() as f64 / range
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mutate_allele_steps_by_one() {
        let model = StepwiseMutationModel::new(1.0, 0, 100).unwrap();
        let mut rng = Rng::new(1);
        // Rate 1.0 -> always mutates; step is +-1.
        for _ in 0..200 {
            let a = model.mutate_allele(50, &mut rng);
            assert!(a == 49 || a == 51, "stepped further than 1: {a}");
        }
    }

    #[test]
    fn mutation_respects_bounds() {
        let model = StepwiseMutationModel::new(1.0, 5, 10).unwrap();
        let mut rng = Rng::new(2);
        // At the lower bound a -1 step is clamped.
        for _ in 0..100 {
            let a = model.mutate_allele(5, &mut rng);
            assert!((5..=10).contains(&a));
        }
        // At the upper bound a +1 step is clamped.
        for _ in 0..100 {
            let a = model.mutate_allele(10, &mut rng);
            assert!((5..=10).contains(&a));
        }
    }

    #[test]
    fn simulation_spreads_allele_sizes() {
        // Starting all identical, drift + stepwise mutation builds
        // allele-size variance over time.
        let model = StepwiseMutationModel::new(0.05, 0, 200).unwrap();
        let early = model.simulate(100, 5, 100, 42).unwrap();
        let late = model.simulate(100, 200, 100, 42).unwrap();
        assert!(
            allele_size_variance(&late) > allele_size_variance(&early),
            "variance did not grow"
        );
    }

    #[test]
    fn simulation_is_deterministic() {
        let model = StepwiseMutationModel::new(0.05, 0, 200).unwrap();
        let a = model.simulate(50, 30, 100, 7).unwrap();
        let b = model.simulate(50, 30, 100, 7).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn heterozygosity_of_a_monomorphic_sample_is_zero() {
        assert_eq!(microsat_heterozygosity(&[10, 10, 10, 10]), 0.0);
        // A maximally diverse sample approaches 1.
        let diverse: Vec<i32> = (0..100).collect();
        assert!(microsat_heterozygosity(&diverse) > 0.98);
    }

    #[test]
    fn m_ratio_is_one_for_a_full_contiguous_range() {
        // Alleles 5..=10 with every class present -> M = 6/6 = 1.
        let alleles = vec![5, 6, 7, 8, 9, 10];
        assert!((m_ratio(&alleles) - 1.0).abs() < 1e-12);
        // Gaps in the range drop M below 1.
        let gapped = vec![5, 5, 10, 10];
        assert!(m_ratio(&gapped) < 1.0);
    }

    #[test]
    fn allele_size_variance_of_constant_sample_is_zero() {
        assert_eq!(allele_size_variance(&[7, 7, 7]), 0.0);
        assert!((mean_allele_size(&[4, 6, 8]) - 6.0).abs() < 1e-12);
    }

    #[test]
    fn model_rejects_bad_params() {
        assert!(StepwiseMutationModel::new(1.5, 0, 10).is_err());
        assert!(StepwiseMutationModel::new(0.1, 10, 5).is_err());
        let model = StepwiseMutationModel::new(0.1, 0, 10).unwrap();
        assert!(model.simulate(0, 5, 5, 1).is_err());
        assert!(model.simulate(10, 5, 99, 1).is_err()); // start out of bounds
    }
}
