//! Mutation models for the forward simulator.
//!
//! A [`MutationModel`] decides, each generation, where new mutations
//! fall on the offspring genomes.
//!
//! - **Infinite-sites** ([`MutationModel::InfiniteSites`]) — every
//!   mutation creates a *brand-new* segregating site at a previously
//!   unoccupied real-valued position. No site is ever hit twice, so a
//!   derived allele is unambiguously "the mutation". This is the model
//!   behind `ms` and the coalescent. The per-genome number of new
//!   mutations is `Poisson(mu * L)`, where `mu` is the per-base rate
//!   and `L` the sequence length.
//! - **Finite-sites** ([`MutationModel::FiniteSites`]) — the genome has
//!   a fixed number of sites each mutating at a per-site rate; a
//!   mutation *flips* the allele, so recurrent and back mutation are
//!   possible. This is the model behind finite-locus simulators and
//!   the one that makes the infinite-sites assumption falsifiable.
//!
//! The model owns no RNG; [`MutationModel::draw_mutations`] takes a
//! `&mut Rng` so the caller controls the seed.

use crate::error::{PopgenError, Result};
use crate::rng::Rng;
use serde::{Deserialize, Serialize};

/// A mutation event to apply to one genome.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum MutationEvent {
    /// Create a new infinite-sites variant at this genomic position.
    /// The simulator allocates the site index and records the position.
    NewSite {
        /// Real-valued genomic coordinate of the new variant.
        position: f64,
    },
    /// Flip the allele at an existing finite-sites locus.
    Flip {
        /// Site index to flip (ancestral <-> derived).
        site: usize,
    },
}

/// How and where new mutations arise.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum MutationModel {
    /// Infinite-sites: every mutation is a new variant at a fresh
    /// position. `rate` is the per-base-pair, per-generation mutation
    /// rate.
    InfiniteSites {
        /// Per-base-pair per-generation mutation rate (`mu`).
        rate: f64,
    },
    /// Finite-sites: a fixed `n_sites`-locus genome, each locus
    /// mutating at `rate` per generation; mutation flips the allele.
    FiniteSites {
        /// Number of fixed loci.
        n_sites: usize,
        /// Per-site per-generation mutation rate.
        rate: f64,
    },
}

impl MutationModel {
    /// Validates the model parameters.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] on a negative rate, a rate above 1
    /// (per-site rates are probabilities), or zero finite-sites loci.
    pub fn validate(&self) -> Result<()> {
        match self {
            MutationModel::InfiniteSites { rate } => {
                if *rate < 0.0 {
                    return Err(PopgenError::invalid(
                        "mutation_rate",
                        "must be non-negative",
                    ));
                }
            }
            MutationModel::FiniteSites { n_sites, rate } => {
                if *n_sites == 0 {
                    return Err(PopgenError::invalid(
                        "n_sites",
                        "finite-sites model needs at least one locus",
                    ));
                }
                if !(0.0..=1.0).contains(rate) {
                    return Err(PopgenError::invalid(
                        "mutation_rate",
                        "per-site rate must lie in [0, 1]",
                    ));
                }
            }
        }
        Ok(())
    }

    /// Draws the mutation events that fall on a single genome this
    /// generation.
    ///
    /// For [`MutationModel::InfiniteSites`] the count is
    /// `Poisson(rate * sequence_length)` and each event gets a uniform
    /// position on `[0, sequence_length)`. For
    /// [`MutationModel::FiniteSites`] each locus independently flips
    /// with probability `rate` (`sequence_length` is ignored).
    pub fn draw_mutations(&self, rng: &mut Rng, sequence_length: f64) -> Vec<MutationEvent> {
        match self {
            MutationModel::InfiniteSites { rate } => {
                let expected = rate * sequence_length.max(0.0);
                let n = rng.poisson(expected);
                (0..n)
                    .map(|_| MutationEvent::NewSite {
                        position: rng.uniform() * sequence_length.max(f64::MIN_POSITIVE),
                    })
                    .collect()
            }
            MutationModel::FiniteSites { n_sites, rate } => (0..*n_sites)
                .filter(|_| rng.bernoulli(*rate))
                .map(|site| MutationEvent::Flip { site })
                .collect(),
        }
    }

    /// The fixed locus count for a finite-sites model, or `None` for
    /// infinite-sites.
    pub fn fixed_site_count(&self) -> Option<usize> {
        match self {
            MutationModel::InfiniteSites { .. } => None,
            MutationModel::FiniteSites { n_sites, .. } => Some(*n_sites),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn infinite_sites_count_tracks_rate_and_length() {
        let model = MutationModel::InfiniteSites { rate: 1e-3 };
        let mut rng = Rng::new(1);
        let trials = 20_000;
        let total: usize = (0..trials)
            .map(|_| model.draw_mutations(&mut rng, 1000.0).len())
            .sum();
        // Expected per-genome mutations = mu * L = 1.0.
        let mean = total as f64 / trials as f64;
        assert!((mean - 1.0).abs() < 0.05, "mean = {mean}");
    }

    #[test]
    fn infinite_sites_positions_in_range() {
        let model = MutationModel::InfiniteSites { rate: 0.1 };
        let mut rng = Rng::new(2);
        for _ in 0..200 {
            for ev in model.draw_mutations(&mut rng, 500.0) {
                match ev {
                    MutationEvent::NewSite { position } => {
                        assert!((0.0..500.0).contains(&position));
                    }
                    _ => panic!("infinite-sites produced a Flip"),
                }
            }
        }
    }

    #[test]
    fn finite_sites_flip_rate() {
        let model = MutationModel::FiniteSites {
            n_sites: 100,
            rate: 0.05,
        };
        let mut rng = Rng::new(3);
        let trials = 5_000;
        let total: usize = (0..trials)
            .map(|_| model.draw_mutations(&mut rng, 0.0).len())
            .sum();
        // Expected flips per genome = n_sites * rate = 5.
        let mean = total as f64 / trials as f64;
        assert!((mean - 5.0).abs() < 0.1, "mean = {mean}");
    }

    #[test]
    fn validate_catches_bad_params() {
        assert!(MutationModel::InfiniteSites { rate: -1.0 }
            .validate()
            .is_err());
        assert!(MutationModel::FiniteSites {
            n_sites: 0,
            rate: 0.1
        }
        .validate()
        .is_err());
        assert!(MutationModel::FiniteSites {
            n_sites: 10,
            rate: 1.5
        }
        .validate()
        .is_err());
        assert!(MutationModel::InfiniteSites { rate: 1e-8 }
            .validate()
            .is_ok());
    }

    #[test]
    fn fixed_site_count_reflects_model() {
        assert_eq!(
            MutationModel::InfiniteSites { rate: 1e-8 }.fixed_site_count(),
            None
        );
        assert_eq!(
            MutationModel::FiniteSites {
                n_sites: 42,
                rate: 0.01
            }
            .fixed_site_count(),
            Some(42)
        );
    }
}
