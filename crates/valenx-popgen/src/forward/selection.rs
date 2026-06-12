//! Fitness models for the forward simulator.
//!
//! A [`SelectionModel`] turns a diploid genome pair into a relative
//! fitness `w >= 0`. The Wright-Fisher simulator
//! ([`crate::forward::WrightFisher`]) draws parents in proportion to
//! `w`, so a value of `1.0` is neutral, `> 1` is favoured and `< 1` is
//! purged.
//!
//! Three classic models are provided:
//!
//! - **Additive** — each derived allele copy adds `s` (with a
//!   dominance coefficient `h` controlling the heterozygote effect):
//!   `w = 1 + (#derived copies scaled by h, s)`. This is the standard
//!   single-locus selection model `1, 1 + hs, 1 + 2s`. Effects sum
//!   across loci, so total fitness is `1 + sum_i effect_i`.
//! - **Multiplicative** — locus fitnesses *multiply* rather than add:
//!   `w = prod_i (1 + effect_i)`. The natural model when loci act
//!   independently on survival.
//! - **Epistatic** — a pairwise interaction term: on top of the
//!   additive main effects, every pair of selected loci that are both
//!   homozygous-derived contributes an extra `epsilon`. A positive
//!   `epsilon` is synergistic, a negative one is antagonistic.

use crate::error::{PopgenError, Result};
use crate::model::Genome;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Per-site selection coefficients keyed by site index.
///
/// `s` is the selection coefficient and `h` the dominance coefficient
/// (`h = 0` recessive, `0.5` additive/codominant, `1` dominant). A site
/// absent from the map is neutral.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct SelectionCoefficients {
    /// site index -> (s, h).
    effects: BTreeMap<usize, (f64, f64)>,
}

impl SelectionCoefficients {
    /// An empty (entirely neutral) coefficient set.
    pub fn neutral() -> Self {
        SelectionCoefficients::default()
    }

    /// Records a selection coefficient `s` and dominance `h` at `site`.
    pub fn set(&mut self, site: usize, s: f64, h: f64) -> &mut Self {
        self.effects.insert(site, (s, h));
        self
    }

    /// The `(s, h)` pair at `site`, or `(0, 0.5)` (neutral) if unset.
    pub fn get(&self, site: usize) -> (f64, f64) {
        self.effects.get(&site).copied().unwrap_or((0.0, 0.5))
    }

    /// `true` if no site carries a non-zero coefficient.
    pub fn is_neutral(&self) -> bool {
        self.effects.values().all(|&(s, _)| s == 0.0)
    }

    /// The selected site indices, ascending.
    pub fn sites(&self) -> Vec<usize> {
        self.effects.keys().copied().collect()
    }
}

/// How per-locus fitness effects combine into a genome-pair fitness.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum SelectionModel {
    /// Locus effects **sum**: `w = 1 + sum_i effect_i`, clamped at 0.
    Additive(SelectionCoefficients),
    /// Locus effects **multiply**: `w = prod_i (1 + effect_i)`.
    Multiplicative(SelectionCoefficients),
    /// Additive main effects plus a pairwise epistasis term `epsilon`
    /// added for every pair of loci both homozygous-derived.
    Epistatic {
        /// Per-site additive main effects.
        coefficients: SelectionCoefficients,
        /// Pairwise interaction coefficient.
        epsilon: f64,
    },
}

impl SelectionModel {
    /// A neutral model (every genome pair has fitness `1.0`).
    pub fn neutral() -> Self {
        SelectionModel::Additive(SelectionCoefficients::neutral())
    }

    /// `true` if this model assigns fitness `1.0` to every genome pair.
    pub fn is_neutral(&self) -> bool {
        match self {
            SelectionModel::Additive(c) | SelectionModel::Multiplicative(c) => c.is_neutral(),
            SelectionModel::Epistatic {
                coefficients,
                epsilon,
            } => coefficients.is_neutral() && *epsilon == 0.0,
        }
    }

    /// Per-locus fitness effect of a diploid genotype: `0` derived
    /// copies -> `0`, `1` -> `h*s`, `2` -> `s`.
    fn locus_effect(copies: u8, s: f64, h: f64) -> f64 {
        match copies {
            0 => 0.0,
            1 => h * s,
            _ => s,
        }
    }

    /// Computes the relative fitness of a diploid `(maternal,
    /// paternal)` genome pair. The result is clamped to be non-negative
    /// (a fitness below 0 is meaningless — it means lethal, i.e. 0).
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] if the two genomes are not a diploid
    /// pair is *not* checked here — the caller guarantees diploidy;
    /// instead this returns the computed fitness directly.
    pub fn fitness(&self, maternal: &Genome, paternal: &Genome) -> f64 {
        let copies = |site: usize| maternal.allele(site) + paternal.allele(site);
        let w = match self {
            SelectionModel::Additive(c) => {
                let mut acc = 1.0;
                for site in c.sites() {
                    let (s, h) = c.get(site);
                    acc += Self::locus_effect(copies(site), s, h);
                }
                acc
            }
            SelectionModel::Multiplicative(c) => {
                let mut acc = 1.0;
                for site in c.sites() {
                    let (s, h) = c.get(site);
                    acc *= 1.0 + Self::locus_effect(copies(site), s, h);
                }
                acc
            }
            SelectionModel::Epistatic {
                coefficients,
                epsilon,
            } => {
                let mut acc = 1.0;
                let sites = coefficients.sites();
                for &site in &sites {
                    let (s, h) = coefficients.get(site);
                    acc += Self::locus_effect(copies(site), s, h);
                }
                // Pairwise epistasis: both loci homozygous-derived.
                for i in 0..sites.len() {
                    for j in (i + 1)..sites.len() {
                        if copies(sites[i]) == 2 && copies(sites[j]) == 2 {
                            acc += epsilon;
                        }
                    }
                }
                acc
            }
        };
        w.max(0.0)
    }

    /// Validates that every coefficient keeps fitness in a sane range
    /// (`s >= -1`, `0 <= h <= 1`).
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] on an out-of-domain coefficient.
    pub fn validate(&self) -> Result<()> {
        let check = |c: &SelectionCoefficients| -> Result<()> {
            for site in c.sites() {
                let (s, h) = c.get(site);
                if s < -1.0 {
                    return Err(PopgenError::invalid(
                        "selection_coefficient",
                        "s must be >= -1 (s = -1 is lethal)",
                    ));
                }
                if !(0.0..=1.0).contains(&h) {
                    return Err(PopgenError::invalid(
                        "dominance_coefficient",
                        "h must lie in [0, 1]",
                    ));
                }
            }
            Ok(())
        };
        match self {
            SelectionModel::Additive(c) | SelectionModel::Multiplicative(c) => check(c),
            SelectionModel::Epistatic { coefficients, .. } => check(coefficients),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn neutral_model_is_flat() {
        let m = SelectionModel::neutral();
        assert!(m.is_neutral());
        let g = Genome::from_derived(vec![0, 1, 2]);
        assert!((m.fitness(&g, &g) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn additive_with_dominance() {
        let mut c = SelectionCoefficients::neutral();
        c.set(0, 0.1, 0.5); // s = 0.1, codominant
        let m = SelectionModel::Additive(c);
        let anc = Genome::ancestral();
        let het = Genome::from_derived(vec![0]);
        // 0 derived copies: w = 1.
        assert!((m.fitness(&anc, &anc) - 1.0).abs() < 1e-12);
        // 1 copy (heterozygote): w = 1 + h*s = 1.05.
        assert!((m.fitness(&het, &anc) - 1.05).abs() < 1e-12);
        // 2 copies (homozygote): w = 1 + s = 1.1.
        assert!((m.fitness(&het, &het) - 1.1).abs() < 1e-12);
    }

    #[test]
    fn recessive_deleterious_hides_in_heterozygote() {
        let mut c = SelectionCoefficients::neutral();
        c.set(0, -0.2, 0.0); // recessive lethal-ish, h = 0
        let m = SelectionModel::Additive(c);
        let anc = Genome::ancestral();
        let der = Genome::from_derived(vec![0]);
        // Heterozygote: h*s = 0, so fitness stays 1.
        assert!((m.fitness(&der, &anc) - 1.0).abs() < 1e-12);
        // Homozygote derived: w = 1 - 0.2 = 0.8.
        assert!((m.fitness(&der, &der) - 0.8).abs() < 1e-12);
    }

    #[test]
    fn multiplicative_differs_from_additive_at_two_loci() {
        let mut c = SelectionCoefficients::neutral();
        c.set(0, 0.1, 0.5);
        c.set(1, 0.1, 0.5);
        let hom = Genome::from_derived(vec![0, 1]);
        let add = SelectionModel::Additive(c.clone());
        let mul = SelectionModel::Multiplicative(c);
        // Additive: 1 + 0.1 + 0.1 = 1.2.
        assert!((add.fitness(&hom, &hom) - 1.2).abs() < 1e-12);
        // Multiplicative: 1.1 * 1.1 = 1.21.
        assert!((mul.fitness(&hom, &hom) - 1.21).abs() < 1e-12);
    }

    #[test]
    fn epistasis_adds_an_interaction_term() {
        let mut c = SelectionCoefficients::neutral();
        c.set(0, 0.05, 0.5);
        c.set(1, 0.05, 0.5);
        let m = SelectionModel::Epistatic {
            coefficients: c,
            epsilon: 0.3,
        };
        let hom = Genome::from_derived(vec![0, 1]);
        // Main effects 1 + 0.05 + 0.05 = 1.1, plus epsilon 0.3 = 1.4.
        assert!((m.fitness(&hom, &hom) - 1.4).abs() < 1e-12);
        // A heterozygote at one locus -> no homozygous pair -> no
        // epistasis term.
        let het = Genome::ancestral();
        assert!(m.fitness(&hom, &het) < 1.2);
    }

    #[test]
    fn fitness_clamps_at_zero() {
        let mut c = SelectionCoefficients::neutral();
        c.set(0, -1.0, 1.0); // lethal
        let m = SelectionModel::Additive(c);
        let der = Genome::from_derived(vec![0]);
        assert_eq!(m.fitness(&der, &der), 0.0);
    }

    #[test]
    fn validate_rejects_out_of_domain() {
        let mut c = SelectionCoefficients::neutral();
        c.set(0, -2.0, 0.5);
        assert!(SelectionModel::Additive(c).validate().is_err());
        let mut c2 = SelectionCoefficients::neutral();
        c2.set(0, 0.1, 1.5);
        assert!(SelectionModel::Additive(c2).validate().is_err());
    }
}
