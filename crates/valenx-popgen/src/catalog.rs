//! A `stdpopsim`-class species and demographic-model catalog.
//!
//! `stdpopsim` is a community catalog of *standard* population-genetic
//! models — for each species a reference effective size, mutation and
//! recombination rates, a genetic map, and named published demographic
//! histories. Researchers simulate from the catalog so results are
//! comparable across studies.
//!
//! This module ships a small built-in catalog: a handful of model
//! species, each with a [`SpeciesModel`] (genome-wide parameters) and
//! one or more named [`DemographicModel`]s. A catalog entry plugs
//! straight into [`fn@crate::coalescent::coalescent`] (via its
//! [`crate::coalescent::PopHistory`]) or the Wright-Fisher simulator.
//!
//! ## v1 scope
//!
//! The catalog is deliberately small and the parameters are *round
//! representative values*, not the precise published estimates the
//! real `stdpopsim` curates — it demonstrates the mechanism. The
//! genetic map is a single uniform recombination rate per species
//! rather than a per-chromosome map.

use crate::coalescent::kingman::PopHistory;
use crate::error::{PopgenError, Result};

/// Genome-wide parameters of a catalogued species.
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesModel {
    /// Catalog identifier, e.g. `"HomSap"`.
    pub id: &'static str,
    /// Common name.
    pub common_name: &'static str,
    /// Reference diploid effective population size.
    pub effective_size: f64,
    /// Per-base-pair per-generation mutation rate.
    pub mutation_rate: f64,
    /// Per-base-pair per-generation recombination rate (a single
    /// uniform genetic-map rate for this v1 catalog).
    pub recombination_rate: f64,
    /// Generation time in years.
    pub generation_time: f64,
}

/// A named demographic model for a species.
#[derive(Clone, Debug, PartialEq)]
pub struct DemographicModel {
    /// Model identifier, e.g. `"Constant"` or `"OutOfAfrica"`.
    pub id: &'static str,
    /// Short human-readable description.
    pub description: &'static str,
    /// The size history as `(duration_in_generations, size)` segments
    /// ordered from the present backward.
    pub epochs: Vec<(f64, f64)>,
}

impl DemographicModel {
    /// Converts the model to a [`PopHistory`] usable by the coalescent.
    ///
    /// A single-epoch model collapses to [`PopHistory::Constant`].
    pub fn to_pop_history(&self) -> PopHistory {
        if self.epochs.len() == 1 {
            PopHistory::Constant(self.epochs[0].1)
        } else {
            PopHistory::Piecewise(self.epochs.clone())
        }
    }
}

/// The built-in species + demographic-model catalog.
#[derive(Clone, Debug)]
pub struct Catalog {
    species: Vec<(SpeciesModel, Vec<DemographicModel>)>,
}

impl Catalog {
    /// Builds the standard built-in catalog.
    pub fn standard() -> Self {
        let species = vec![
            // --- Homo sapiens --------------------------------------
            (
                SpeciesModel {
                    id: "HomSap",
                    common_name: "human",
                    effective_size: 10_000.0,
                    mutation_rate: 1.29e-8,
                    recombination_rate: 1.0e-8,
                    generation_time: 29.0,
                },
                vec![
                    DemographicModel {
                        id: "Constant",
                        description: "constant-size panmictic population",
                        epochs: vec![(f64::INFINITY, 10_000.0)],
                    },
                    DemographicModel {
                        id: "OutOfAfricaExpansion",
                        description: "an ancestral population, a bottleneck at the \
                             out-of-Africa migration, then recent expansion",
                        epochs: vec![
                            // Recent: a large, recently expanded population.
                            (500.0, 50_000.0),
                            // The out-of-Africa bottleneck.
                            (1_500.0, 2_000.0),
                            // The larger ancestral African population.
                            (f64::INFINITY, 12_000.0),
                        ],
                    },
                ],
            ),
            // --- Drosophila melanogaster ---------------------------
            (
                SpeciesModel {
                    id: "DroMel",
                    common_name: "fruit fly",
                    effective_size: 1_700_000.0,
                    mutation_rate: 5.49e-9,
                    recombination_rate: 2.4e-8,
                    generation_time: 0.1,
                },
                vec![DemographicModel {
                    id: "Constant",
                    description: "constant-size panmictic population",
                    epochs: vec![(f64::INFINITY, 1_700_000.0)],
                }],
            ),
            // --- Arabidopsis thaliana ------------------------------
            (
                SpeciesModel {
                    id: "AraTha",
                    common_name: "thale cress",
                    effective_size: 250_000.0,
                    mutation_rate: 7.0e-9,
                    recombination_rate: 3.0e-8,
                    generation_time: 1.0,
                },
                vec![
                    DemographicModel {
                        id: "Constant",
                        description: "constant-size panmictic population",
                        epochs: vec![(f64::INFINITY, 250_000.0)],
                    },
                    DemographicModel {
                        id: "SouthMiddleAtlas",
                        description: "a recent contraction following a larger \
                             ancestral population",
                        epochs: vec![(10_000.0, 25_000.0), (f64::INFINITY, 250_000.0)],
                    },
                ],
            ),
            // --- Escherichia coli ----------------------------------
            (
                SpeciesModel {
                    id: "EscCol",
                    common_name: "E. coli",
                    effective_size: 1_800_000.0,
                    mutation_rate: 1.0e-9,
                    recombination_rate: 0.0,
                    generation_time: 0.00057,
                },
                vec![DemographicModel {
                    id: "Constant",
                    description: "constant-size clonal population",
                    epochs: vec![(f64::INFINITY, 1_800_000.0)],
                }],
            ),
        ];

        Catalog { species }
    }

    /// All catalogued species ids.
    pub fn species_ids(&self) -> Vec<&'static str> {
        self.species.iter().map(|(s, _)| s.id).collect()
    }

    /// Number of catalogued species.
    pub fn species_count(&self) -> usize {
        self.species.len()
    }

    /// Looks up a species' genome-wide model by catalog id.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] if no species has that id.
    pub fn species(&self, id: &str) -> Result<&SpeciesModel> {
        self.species
            .iter()
            .find(|(s, _)| s.id == id)
            .map(|(s, _)| s)
            .ok_or_else(|| PopgenError::invalid("species", format!("no catalog entry `{id}`")))
    }

    /// The named demographic models available for a species.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] if no species has that id.
    pub fn demographic_models(&self, species_id: &str) -> Result<&[DemographicModel]> {
        self.species
            .iter()
            .find(|(s, _)| s.id == species_id)
            .map(|(_, m)| m.as_slice())
            .ok_or_else(|| {
                PopgenError::invalid("species", format!("no catalog entry `{species_id}`"))
            })
    }

    /// Looks up a specific demographic model.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] if the species or model id is unknown.
    pub fn demographic_model(&self, species_id: &str, model_id: &str) -> Result<&DemographicModel> {
        self.demographic_models(species_id)?
            .iter()
            .find(|m| m.id == model_id)
            .ok_or_else(|| {
                PopgenError::invalid(
                    "model",
                    format!("species `{species_id}` has no model `{model_id}`"),
                )
            })
    }
}

impl Default for Catalog {
    fn default() -> Self {
        Catalog::standard()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_catalog_has_several_species() {
        let cat = Catalog::standard();
        assert!(cat.species_count() >= 4);
        let ids = cat.species_ids();
        assert!(ids.contains(&"HomSap"));
        assert!(ids.contains(&"DroMel"));
    }

    #[test]
    fn species_lookup_returns_parameters() {
        let cat = Catalog::standard();
        let human = cat.species("HomSap").unwrap();
        assert_eq!(human.common_name, "human");
        assert!(human.effective_size > 0.0);
        assert!(human.mutation_rate > 0.0);
        assert!(cat.species("NoSuchSpecies").is_err());
    }

    #[test]
    fn demographic_models_are_listed() {
        let cat = Catalog::standard();
        let models = cat.demographic_models("HomSap").unwrap();
        assert!(models.len() >= 2);
        // Every species has a Constant model.
        for id in cat.species_ids() {
            let ms = cat.demographic_models(id).unwrap();
            assert!(ms.iter().any(|m| m.id == "Constant"));
        }
    }

    #[test]
    fn constant_model_becomes_a_constant_pop_history() {
        let cat = Catalog::standard();
        let model = cat.demographic_model("HomSap", "Constant").unwrap();
        match model.to_pop_history() {
            PopHistory::Constant(n) => assert!((n - 10_000.0).abs() < 1e-6),
            _ => panic!("constant model did not collapse"),
        }
    }

    #[test]
    fn multi_epoch_model_becomes_piecewise() {
        let cat = Catalog::standard();
        let model = cat
            .demographic_model("HomSap", "OutOfAfricaExpansion")
            .unwrap();
        match model.to_pop_history() {
            PopHistory::Piecewise(segs) => assert!(segs.len() >= 2),
            _ => panic!("multi-epoch model did not become piecewise"),
        }
    }

    #[test]
    fn catalog_pop_history_drives_the_coalescent() {
        use crate::coalescent::coalescent;
        let cat = Catalog::standard();
        let model = cat.demographic_model("HomSap", "Constant").unwrap();
        let labels: Vec<String> = (0..6).map(|i| format!("s{i}")).collect();
        let tree = coalescent(&labels, &model.to_pop_history(), 42).unwrap();
        assert_eq!(tree.leaf_count(), 6);
    }

    #[test]
    fn unknown_model_is_an_error() {
        let cat = Catalog::standard();
        assert!(cat.demographic_model("HomSap", "Nonexistent").is_err());
        assert!(cat.demographic_model("Nonexistent", "Constant").is_err());
    }
}
