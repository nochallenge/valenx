//! The population-genetics data model.
//!
//! The central abstraction is a [`Population`] of [`Individual`]s, each
//! carrying one or two [`Genome`]s (a *haplotype* — an ordered list of
//! the derived alleles it carries at a set of [`Site`]s).
//!
//! ## Allele representation
//!
//! `valenx-popgen` works with **biallelic, derived-allele-counted**
//! genomes throughout. A [`Site`] has a real-valued genomic
//! `position`; a [`Genome`] stores the *sorted set of site indices* at
//! which it carries the derived (`1`) allele. Every other site is
//! ancestral (`0`). This sparse encoding is what `msprime` / `tskit`
//! use and it keeps neutral-locus simulations cheap: a 0.001-frequency
//! variant costs one `usize` per carrier rather than a bit per genome.
//!
//! A [`Ploidy`] of [`Ploidy::Haploid`] gives every individual one
//! genome; [`Ploidy::Diploid`] gives two. The Wright-Fisher forward
//! simulator ([`crate::forward`]) is diploid; the coalescent
//! ([`mod@crate::coalescent`]) samples haploid lineages.

use crate::error::{PopgenError, Result};
use serde::{Deserialize, Serialize};

/// Number of genome copies per individual.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Ploidy {
    /// One genome copy per individual.
    Haploid,
    /// Two genome copies per individual.
    Diploid,
}

impl Ploidy {
    /// Number of genome copies (`1` or `2`).
    pub fn copies(self) -> usize {
        match self {
            Ploidy::Haploid => 1,
            Ploidy::Diploid => 2,
        }
    }
}

/// A single segregating site on the genome.
///
/// Sites are biallelic: an ancestral allele (`0`) and a derived allele
/// (`1`). `position` is a real-valued coordinate on `[0, sequence
/// length)` — recombination (crossover, gene conversion) acts on this
/// continuous coordinate.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Site {
    /// Genomic coordinate of the site (real-valued, increasing).
    pub position: f64,
    /// Optional ancestral-allele label (e.g. a nucleotide); `None`
    /// for an abstract biallelic locus.
    pub ancestral: Option<char>,
    /// Optional derived-allele label.
    pub derived: Option<char>,
}

impl Site {
    /// A bare biallelic site at `position` with no nucleotide labels.
    pub fn at(position: f64) -> Self {
        Site {
            position,
            ancestral: None,
            derived: None,
        }
    }

    /// A site with explicit ancestral / derived nucleotide labels.
    pub fn nucleotide(position: f64, ancestral: char, derived: char) -> Self {
        Site {
            position,
            ancestral: Some(ancestral),
            derived: Some(derived),
        }
    }
}

/// A single haplotype: the sorted set of site indices at which this
/// genome carries the **derived** allele.
///
/// Any site index *not* in [`derived`](Genome::derived) is ancestral.
/// The vector is kept sorted and deduplicated by every mutating method
/// so set operations (genotype lookup, recombination) stay O(n).
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Genome {
    /// Sorted, deduplicated site indices carrying the derived allele.
    derived: Vec<usize>,
}

impl Genome {
    /// An all-ancestral genome (no derived alleles).
    pub fn ancestral() -> Self {
        Genome::default()
    }

    /// Builds a genome from a list of derived-allele site indices.
    /// The list is sorted and deduplicated.
    pub fn from_derived(mut sites: Vec<usize>) -> Self {
        sites.sort_unstable();
        sites.dedup();
        Genome { derived: sites }
    }

    /// The sorted derived-allele site indices.
    pub fn derived(&self) -> &[usize] {
        &self.derived
    }

    /// Number of derived alleles this genome carries.
    pub fn derived_count(&self) -> usize {
        self.derived.len()
    }

    /// Allelic state at `site`: `1` if derived, `0` if ancestral.
    pub fn allele(&self, site: usize) -> u8 {
        u8::from(self.derived.binary_search(&site).is_ok())
    }

    /// Sets the derived allele at `site` (idempotent).
    pub fn set_derived(&mut self, site: usize) {
        if let Err(pos) = self.derived.binary_search(&site) {
            self.derived.insert(pos, site);
        }
    }

    /// Clears the derived allele at `site`, reverting it to ancestral.
    pub fn clear_derived(&mut self, site: usize) {
        if let Ok(pos) = self.derived.binary_search(&site) {
            self.derived.remove(pos);
        }
    }

    /// Flips the allele at `site` (ancestral <-> derived). Used by the
    /// finite-sites recurrent-mutation model.
    pub fn flip(&mut self, site: usize) {
        match self.derived.binary_search(&site) {
            Ok(pos) => {
                self.derived.remove(pos);
            }
            Err(pos) => self.derived.insert(pos, site),
        }
    }
}

/// One member of a [`Population`]: an `id` plus its genome copies.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Individual {
    /// Stable identifier within the population.
    pub id: usize,
    /// Genome copies: one if haploid, two if diploid.
    pub genomes: Vec<Genome>,
    /// Optional deme index for structured / migration models.
    pub deme: usize,
}

impl Individual {
    /// A new individual with all-ancestral genomes.
    pub fn ancestral(id: usize, ploidy: Ploidy, deme: usize) -> Self {
        Individual {
            id,
            genomes: (0..ploidy.copies()).map(|_| Genome::ancestral()).collect(),
            deme,
        }
    }

    /// Diploid genotype at `site`: the number of derived copies
    /// (`0`, `1` or `2`). For a haploid individual returns `0` or `1`.
    pub fn genotype(&self, site: usize) -> u8 {
        self.genomes.iter().map(|g| g.allele(site)).sum()
    }

    /// `true` if the individual is heterozygous at `site` (diploid,
    /// exactly one derived copy).
    pub fn is_heterozygous(&self, site: usize) -> bool {
        self.genotype(site) == 1 && self.genomes.len() == 2
    }
}

/// A panmictic (or structured) collection of [`Individual`]s sharing a
/// common list of [`Site`]s.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Population {
    /// The site map shared by every genome in the population.
    sites: Vec<Site>,
    /// The individuals, in arbitrary order.
    individuals: Vec<Individual>,
    /// Ploidy of every individual.
    ploidy: Ploidy,
    /// Length of the simulated genomic segment.
    sequence_length: f64,
}

impl Population {
    /// Builds a population of `n` all-ancestral individuals over a
    /// genome of length `sequence_length` with no segregating sites.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] if `n == 0` or `sequence_length <= 0`.
    pub fn founder(n: usize, ploidy: Ploidy, sequence_length: f64) -> Result<Self> {
        if n == 0 {
            return Err(PopgenError::invalid("n", "population must be non-empty"));
        }
        if sequence_length <= 0.0 {
            return Err(PopgenError::invalid(
                "sequence_length",
                "must be positive",
            ));
        }
        Ok(Population {
            sites: Vec::new(),
            individuals: (0..n)
                .map(|i| Individual::ancestral(i, ploidy, 0))
                .collect(),
            ploidy,
            sequence_length,
        })
    }

    /// Builds a population directly from explicit parts.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] on an empty individual list;
    /// [`PopgenError::Dimension`] if any individual's ploidy disagrees
    /// with `ploidy`.
    pub fn new(
        sites: Vec<Site>,
        individuals: Vec<Individual>,
        ploidy: Ploidy,
        sequence_length: f64,
    ) -> Result<Self> {
        if individuals.is_empty() {
            return Err(PopgenError::invalid(
                "individuals",
                "population must be non-empty",
            ));
        }
        if sequence_length <= 0.0 {
            return Err(PopgenError::invalid(
                "sequence_length",
                "must be positive",
            ));
        }
        for ind in &individuals {
            if ind.genomes.len() != ploidy.copies() {
                return Err(PopgenError::dimension(
                    ploidy.copies(),
                    ind.genomes.len(),
                    "individual genome copies",
                ));
            }
        }
        Ok(Population {
            sites,
            individuals,
            ploidy,
            sequence_length,
        })
    }

    /// Number of individuals in the population (the census size).
    pub fn size(&self) -> usize {
        self.individuals.len()
    }

    /// Ploidy of the population.
    pub fn ploidy(&self) -> Ploidy {
        self.ploidy
    }

    /// Length of the simulated genomic segment.
    pub fn sequence_length(&self) -> f64 {
        self.sequence_length
    }

    /// The shared site map.
    pub fn sites(&self) -> &[Site] {
        &self.sites
    }

    /// Number of segregating sites currently tracked.
    pub fn site_count(&self) -> usize {
        self.sites.len()
    }

    /// The individuals.
    pub fn individuals(&self) -> &[Individual] {
        &self.individuals
    }

    /// Mutable access to the individuals — used by the forward
    /// simulator to overwrite a generation in place.
    pub fn individuals_mut(&mut self) -> &mut Vec<Individual> {
        &mut self.individuals
    }

    /// Total number of genome copies (`size * ploidy`) — the
    /// chromosome sample size `n` for allele-frequency statistics.
    pub fn chromosome_count(&self) -> usize {
        self.individuals.iter().map(|i| i.genomes.len()).sum()
    }

    /// Appends a new site to the shared map and returns its index.
    /// Existing genomes are unaffected (they remain ancestral there).
    pub fn add_site(&mut self, site: Site) -> usize {
        self.sites.push(site);
        self.sites.len() - 1
    }

    /// Count of derived alleles at `site`, summed over every genome
    /// copy in the population.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] if `site` is out of range.
    pub fn derived_count(&self, site: usize) -> Result<usize> {
        if site >= self.sites.len() {
            return Err(PopgenError::invalid("site", "index out of range"));
        }
        Ok(self
            .individuals
            .iter()
            .flat_map(|i| i.genomes.iter())
            .map(|g| g.allele(site) as usize)
            .sum())
    }

    /// Derived-allele frequency at `site` in `[0, 1]`.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] if `site` is out of range.
    pub fn allele_frequency(&self, site: usize) -> Result<f64> {
        let derived = self.derived_count(site)?;
        Ok(derived as f64 / self.chromosome_count().max(1) as f64)
    }

    /// Builds a borrowed genotype matrix: one row per chromosome copy,
    /// one column per site, each entry `0` (ancestral) or `1`
    /// (derived). The natural input for [`crate::stats`].
    pub fn genotype_matrix(&self) -> crate::infer::GenotypeMatrix {
        let positions: Vec<f64> = self.sites.iter().map(|s| s.position).collect();
        let mut rows = Vec::with_capacity(self.chromosome_count());
        for ind in &self.individuals {
            for g in &ind.genomes {
                let row: Vec<u8> =
                    (0..self.sites.len()).map(|s| g.allele(s)).collect();
                rows.push(row);
            }
        }
        crate::infer::GenotypeMatrix::from_rows(rows, positions)
            .expect("population rows are uniform by construction")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn genome_set_and_allele() {
        let mut g = Genome::ancestral();
        assert_eq!(g.allele(3), 0);
        g.set_derived(3);
        g.set_derived(1);
        g.set_derived(3); // idempotent
        assert_eq!(g.derived(), &[1, 3]);
        assert_eq!(g.allele(3), 1);
        assert_eq!(g.allele(1), 1);
        assert_eq!(g.allele(2), 0);
        g.clear_derived(1);
        assert_eq!(g.allele(1), 0);
        g.flip(5);
        assert_eq!(g.allele(5), 1);
        g.flip(5);
        assert_eq!(g.allele(5), 0);
    }

    #[test]
    fn from_derived_sorts_and_dedups() {
        let g = Genome::from_derived(vec![5, 1, 5, 3, 1]);
        assert_eq!(g.derived(), &[1, 3, 5]);
        assert_eq!(g.derived_count(), 3);
    }

    #[test]
    fn individual_genotype_counts_copies() {
        let mut ind = Individual::ancestral(0, Ploidy::Diploid, 0);
        ind.genomes[0].set_derived(2);
        assert_eq!(ind.genotype(2), 1);
        assert!(ind.is_heterozygous(2));
        ind.genomes[1].set_derived(2);
        assert_eq!(ind.genotype(2), 2);
        assert!(!ind.is_heterozygous(2));
    }

    #[test]
    fn founder_population_is_all_ancestral() {
        let pop = Population::founder(10, Ploidy::Diploid, 1000.0).unwrap();
        assert_eq!(pop.size(), 10);
        assert_eq!(pop.chromosome_count(), 20);
        assert_eq!(pop.site_count(), 0);
    }

    #[test]
    fn founder_rejects_bad_args() {
        assert!(Population::founder(0, Ploidy::Haploid, 1.0).is_err());
        assert!(Population::founder(5, Ploidy::Haploid, 0.0).is_err());
    }

    #[test]
    fn allele_frequency_after_adding_a_variant() {
        let mut pop = Population::founder(4, Ploidy::Diploid, 100.0).unwrap();
        let s = pop.add_site(Site::at(50.0));
        // Make 3 of 8 chromosomes derived.
        pop.individuals_mut()[0].genomes[0].set_derived(s);
        pop.individuals_mut()[1].genomes[0].set_derived(s);
        pop.individuals_mut()[1].genomes[1].set_derived(s);
        assert_eq!(pop.derived_count(s).unwrap(), 3);
        assert!((pop.allele_frequency(s).unwrap() - 3.0 / 8.0).abs() < 1e-12);
        assert!(pop.derived_count(99).is_err());
    }

    #[test]
    fn genotype_matrix_shape() {
        let mut pop = Population::founder(3, Ploidy::Diploid, 100.0).unwrap();
        pop.add_site(Site::at(10.0));
        pop.add_site(Site::at(20.0));
        let gm = pop.genotype_matrix();
        assert_eq!(gm.n_samples(), 6);
        assert_eq!(gm.n_sites(), 2);
    }

    #[test]
    fn new_rejects_ploidy_mismatch() {
        let bad = Individual {
            id: 0,
            genomes: vec![Genome::ancestral()], // haploid genome...
            deme: 0,
        };
        // ...declared in a diploid population.
        let r = Population::new(vec![], vec![bad], Ploidy::Diploid, 10.0);
        assert!(r.is_err());
    }
}
