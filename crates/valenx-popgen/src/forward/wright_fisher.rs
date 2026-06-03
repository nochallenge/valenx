//! The Wright-Fisher forward simulator.
//!
//! [`WrightFisher`] runs the canonical discrete-generation, diploid
//! Wright-Fisher process forward in time. Each generation:
//!
//! 1. The next generation's census size is read from the
//!    [`DemographicSchedule`].
//! 2. For each offspring, two parents are drawn from the *current*
//!    generation with probability proportional to their relative
//!    fitness ([`SelectionModel`]) — viability selection.
//! 3. Each parent contributes a gamete formed by recombination
//!    ([`RecombinationModel`]) of its two homologs.
//! 4. New mutations ([`MutationModel`]) are layered onto the offspring
//!    genomes.
//! 5. With a structured population the offspring's deme is resolved
//!    through the [`MigrationModel`].
//!
//! After `generations` steps the final [`Population`] is returned. The
//! whole run is deterministic in the seed — re-running with the same
//! [`SimulationConfig`] reproduces the population bit-for-bit.
//!
//! This is the engine SLiM, fwdpy11 and simuPOP provide; the v1
//! simplifications are documented on [`SimulationConfig`].

use crate::error::{PopgenError, Result};
use crate::forward::demography::DemographicSchedule;
use crate::forward::migration::MigrationModel;
use crate::forward::mutation::{MutationEvent, MutationModel};
use crate::forward::recombination::RecombinationModel;
use crate::forward::selection::SelectionModel;
use crate::model::{Genome, Individual, Population, Site};
use crate::rng::Rng;

/// Everything needed to run a Wright-Fisher simulation.
///
/// ## v1 simplifications
///
/// - Generations are **discrete and non-overlapping** — there is no
///   age structure or overlapping cohorts.
/// - Selection is **viability** selection only (it weights parent
///   choice); there is no separate fecundity or sexual-selection term.
/// - Mating is random within (post-migration) demes; there is no
///   assortative mating and no explicit separate sexes.
/// - The genealogy is *not* recorded here — use
///   [`crate::forward::tree_recording`] for that.
#[derive(Clone, Debug)]
pub struct SimulationConfig {
    /// Forward-in-time population-size schedule.
    pub demography: DemographicSchedule,
    /// Number of generations to simulate.
    pub generations: usize,
    /// Fitness model (use [`SelectionModel::neutral`] for drift only).
    pub selection: SelectionModel,
    /// Mutation model.
    pub mutation: MutationModel,
    /// Recombination model.
    pub recombination: RecombinationModel,
    /// Optional structured-population migration. `None` is panmixia.
    pub migration: Option<MigrationModel>,
    /// Length of the simulated genomic segment.
    pub sequence_length: f64,
    /// RNG seed — fixes the entire run.
    pub seed: u64,
}

impl SimulationConfig {
    /// A minimal neutral panmictic configuration: constant size,
    /// infinite-sites mutation, crossover recombination.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] if `n == 0` or a rate is out of range.
    pub fn neutral(
        n: usize,
        generations: usize,
        mutation_rate: f64,
        recombination_rate: f64,
        sequence_length: f64,
        seed: u64,
    ) -> Result<Self> {
        Ok(SimulationConfig {
            demography: DemographicSchedule::constant(n)?,
            generations,
            selection: SelectionModel::neutral(),
            mutation: MutationModel::InfiniteSites {
                rate: mutation_rate,
            },
            recombination: RecombinationModel::crossover_only(recombination_rate),
            migration: None,
            sequence_length,
            seed,
        })
    }

    /// Validates every component of the configuration.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] / [`PopgenError::Dimension`] from any
    /// component's own validator.
    pub fn validate(&self) -> Result<()> {
        if self.sequence_length <= 0.0 {
            return Err(PopgenError::invalid(
                "sequence_length",
                "must be positive",
            ));
        }
        self.selection.validate()?;
        self.mutation.validate()?;
        self.recombination.validate()?;
        Ok(())
    }
}

/// The Wright-Fisher forward simulator.
#[derive(Debug)]
pub struct WrightFisher {
    config: SimulationConfig,
}

impl WrightFisher {
    /// Builds a simulator from a configuration.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] if the configuration fails validation.
    pub fn new(config: SimulationConfig) -> Result<Self> {
        config.validate()?;
        Ok(WrightFisher { config })
    }

    /// Runs the simulation and returns the final [`Population`].
    ///
    /// The founding generation is all-ancestral; mutations accumulate
    /// from generation 1 onward.
    ///
    /// # Errors
    /// [`PopgenError`] propagated from population construction.
    pub fn run(&self) -> Result<Population> {
        let mut rng = Rng::new(self.config.seed);
        let n0 = self.config.demography.initial_size();
        let mut pop = Population::founder(
            n0,
            crate::model::Ploidy::Diploid,
            self.config.sequence_length,
        )?;

        // Mutable site-position map; infinite-sites mutation appends.
        let mut positions: Vec<f64> = Vec::new();
        // Finite-sites models pre-allocate their fixed loci.
        if let Some(k) = self.config.mutation.fixed_site_count() {
            for i in 0..k {
                let pos = (i as f64 + 0.5) / k as f64 * self.config.sequence_length;
                positions.push(pos);
                pop.add_site(Site::nucleotide(pos, 'A', 'T'));
            }
        }

        for gen in 1..=self.config.generations {
            let next_n = self.config.demography.size_at(gen);
            pop = self.step(&pop, next_n, &mut positions, &mut rng)?;
        }
        Ok(pop)
    }

    /// Advances the population by exactly one generation.
    fn step(
        &self,
        current: &Population,
        next_n: usize,
        positions: &mut Vec<f64>,
        rng: &mut Rng,
    ) -> Result<Population> {
        // Pre-compute parental fitnesses for proportional selection.
        let fitness: Vec<f64> = current
            .individuals()
            .iter()
            .map(|ind| {
                self.config
                    .selection
                    .fitness(&ind.genomes[0], &ind.genomes[1])
            })
            .collect();
        let any_fit = fitness.iter().any(|&w| w > 0.0);

        let n_demes = self
            .config
            .migration
            .as_ref()
            .map(|m| m.deme_count())
            .unwrap_or(1);

        let mut offspring = Vec::with_capacity(next_n);
        for child_id in 0..next_n {
            // Resolve the offspring's deme (panmixia => deme 0).
            let deme = if n_demes > 1 {
                child_id % n_demes
            } else {
                0
            };
            // Backward migration: the parents live in the source deme.
            let parent_deme = match &self.config.migration {
                Some(m) => m.sample_source(deme, rng),
                None => 0,
            };

            let p1 = self.pick_parent(current, &fitness, any_fit, parent_deme, rng);
            let p2 = self.pick_parent(current, &fitness, any_fit, parent_deme, rng);

            let g1 = self.gamete(&current.individuals()[p1], positions, rng);
            let g2 = self.gamete(&current.individuals()[p2], positions, rng);

            let mut child = Individual {
                id: child_id,
                genomes: vec![g1, g2],
                deme,
            };
            // Layer on new mutations.
            self.mutate(&mut child, positions, rng);
            offspring.push(child);
        }

        // Rebuild the site map: infinite-sites mutation may have grown
        // `positions`.
        let sites: Vec<Site> = positions
            .iter()
            .map(|&p| Site::nucleotide(p, 'A', 'T'))
            .collect();
        Population::new(
            sites,
            offspring,
            crate::model::Ploidy::Diploid,
            self.config.sequence_length,
        )
    }

    /// Picks one parent index, weighting by fitness within a deme.
    fn pick_parent(
        &self,
        current: &Population,
        fitness: &[f64],
        any_fit: bool,
        deme: usize,
        rng: &mut Rng,
    ) -> usize {
        // Candidate indices: members of the source deme.
        let candidates: Vec<usize> = current
            .individuals()
            .iter()
            .enumerate()
            .filter(|(_, ind)| ind.deme == deme || self.config.migration.is_none())
            .map(|(i, _)| i)
            .collect();
        if candidates.is_empty() {
            return rng.below(current.size());
        }
        if !any_fit || self.config.selection.is_neutral() {
            return candidates[rng.below(candidates.len())];
        }
        let weights: Vec<f64> = candidates.iter().map(|&i| fitness[i]).collect();
        candidates[rng.weighted_index(&weights)]
    }

    /// Forms one gamete from a diploid parent via recombination.
    fn gamete(&self, parent: &Individual, positions: &[f64], rng: &mut Rng) -> Genome {
        self.config.recombination.recombine(
            &parent.genomes[0],
            &parent.genomes[1],
            positions,
            self.config.sequence_length,
            rng,
        )
    }

    /// Applies new mutations to both of an offspring's genomes.
    fn mutate(&self, child: &mut Individual, positions: &mut Vec<f64>, rng: &mut Rng) {
        for genome_idx in 0..child.genomes.len() {
            let events = self
                .config
                .mutation
                .draw_mutations(rng, self.config.sequence_length);
            for ev in events {
                match ev {
                    MutationEvent::NewSite { position } => {
                        // Infinite-sites: allocate a brand-new site.
                        let site = positions.len();
                        positions.push(position);
                        child.genomes[genome_idx].set_derived(site);
                    }
                    MutationEvent::Flip { site } => {
                        // Finite-sites: flip an existing locus.
                        if site < positions.len() {
                            child.genomes[genome_idx].flip(site);
                        }
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forward::selection::SelectionCoefficients;

    #[test]
    fn neutral_run_produces_a_population() {
        let cfg =
            SimulationConfig::neutral(20, 30, 1e-3, 1e-4, 1000.0, 42).unwrap();
        let sim = WrightFisher::new(cfg).unwrap();
        let pop = sim.run().unwrap();
        assert_eq!(pop.size(), 20);
        assert_eq!(pop.ploidy(), crate::model::Ploidy::Diploid);
    }

    #[test]
    fn run_is_deterministic_in_the_seed() {
        let cfg = SimulationConfig::neutral(15, 20, 2e-3, 1e-4, 500.0, 7).unwrap();
        let a = WrightFisher::new(cfg.clone()).unwrap().run().unwrap();
        let b = WrightFisher::new(cfg).unwrap().run().unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn mutation_accumulates_segregating_sites() {
        // A reasonably high mutation rate over many generations must
        // leave segregating variation behind.
        let cfg =
            SimulationConfig::neutral(30, 40, 5e-3, 1e-4, 2000.0, 99).unwrap();
        let pop = WrightFisher::new(cfg).unwrap().run().unwrap();
        assert!(pop.site_count() > 0, "no mutations accumulated");
    }

    #[test]
    fn strong_positive_selection_raises_allele_frequency() {
        // Seed a beneficial allele model and confirm that, started
        // from a finite-sites locus, selection pushes it up relative
        // to neutral drift on average across seeds.
        let mut coeffs = SelectionCoefficients::neutral();
        coeffs.set(0, 0.5, 0.5); // strongly beneficial, codominant
        let make = |selection: SelectionModel, seed: u64| -> Population {
            let cfg = SimulationConfig {
                demography: DemographicSchedule::constant(60).unwrap(),
                generations: 30,
                selection,
                mutation: MutationModel::FiniteSites {
                    n_sites: 1,
                    rate: 0.02,
                },
                recombination: RecombinationModel::none(),
                migration: None,
                sequence_length: 100.0,
                seed,
            };
            WrightFisher::new(cfg).unwrap().run().unwrap()
        };
        // Average derived-allele frequency at the selected locus.
        let mean_freq = |selection: SelectionModel| -> f64 {
            let mut acc = 0.0;
            for seed in 0..12 {
                let pop = make(selection.clone(), seed);
                acc += pop.allele_frequency(0).unwrap();
            }
            acc / 12.0
        };
        let selected = mean_freq(SelectionModel::Additive(coeffs));
        let neutral = mean_freq(SelectionModel::neutral());
        assert!(
            selected > neutral,
            "selection ({selected}) did not beat drift ({neutral})"
        );
    }

    #[test]
    fn structured_population_assigns_demes() {
        let cfg = SimulationConfig {
            demography: DemographicSchedule::constant(40).unwrap(),
            generations: 10,
            selection: SelectionModel::neutral(),
            mutation: MutationModel::InfiniteSites { rate: 1e-3 },
            recombination: RecombinationModel::crossover_only(1e-4),
            migration: Some(MigrationModel::island(4, 0.1).unwrap()),
            sequence_length: 500.0,
            seed: 5,
        };
        let pop = WrightFisher::new(cfg).unwrap().run().unwrap();
        // Every deme index should appear.
        let demes: std::collections::HashSet<usize> =
            pop.individuals().iter().map(|i| i.deme).collect();
        assert_eq!(demes.len(), 4);
    }

    #[test]
    fn bottleneck_changes_the_population_size() {
        let cfg = SimulationConfig {
            demography: DemographicSchedule::bottleneck(50, 5, 3, 8).unwrap(),
            generations: 20,
            selection: SelectionModel::neutral(),
            mutation: MutationModel::InfiniteSites { rate: 1e-3 },
            recombination: RecombinationModel::none(),
            migration: None,
            sequence_length: 500.0,
            seed: 1,
        };
        let pop = WrightFisher::new(cfg).unwrap().run().unwrap();
        // After generation 8 the schedule is back to 50.
        assert_eq!(pop.size(), 50);
    }

    #[test]
    fn invalid_config_is_rejected() {
        let mut cfg =
            SimulationConfig::neutral(10, 5, 1e-3, 1e-4, 100.0, 1).unwrap();
        cfg.sequence_length = -1.0;
        assert!(WrightFisher::new(cfg).is_err());
    }
}
