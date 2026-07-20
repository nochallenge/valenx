//! Expected genotype frequencies under Hardy-Weinberg equilibrium.
//!
//! Under random mating, no selection, no mutation, no migration, and an
//! infinite population, a di-allelic locus with allele frequencies `p`
//! and `q = 1 - p` settles in one generation to the genotype
//! frequencies given by the binomial expansion of `(p + q)^2`:
//!
//! ```text
//! freq(AA) = p^2
//! freq(Aa) = 2 * p * q
//! freq(aa) = q^2
//! ```
//!
//! These three frequencies sum to `(p + q)^2 = 1`. The canonical
//! ground-truth case is `p = q = 0.5`, giving `0.25 / 0.50 / 0.25`.

use crate::allele::AlleleFreq;
use crate::error::Result;

/// Expected Hardy-Weinberg genotype frequencies for a di-allelic locus.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct GenotypeFreq {
    /// Expected frequency of the `AA` homozygote, `p^2`.
    pub aa_hom: f64,
    /// Expected frequency of the `Aa` heterozygote, `2pq`.
    pub het: f64,
    /// Expected frequency of the `aa` homozygote, `q^2`.
    pub aa: f64,
}

impl GenotypeFreq {
    /// Compute the equilibrium genotype frequencies from a validated
    /// [`AlleleFreq`].
    ///
    /// Always succeeds because `AlleleFreq` already guarantees
    /// `p, q in [0, 1]`; the [`Result`] is kept for signature symmetry
    /// with the rest of the crate.
    ///
    /// # Errors
    ///
    /// Never fails for a frequency pair produced by this crate; the
    /// fallible return type mirrors the other model functions.
    pub fn from_alleles(freq: &AlleleFreq) -> Result<Self> {
        let p = freq.p;
        let q = freq.q;
        Ok(GenotypeFreq {
            aa_hom: p * p,
            het: 2.0 * p * q,
            aa: q * q,
        })
    }

    /// Compute equilibrium genotype frequencies directly from `p`
    /// (with `q = 1 - p`).
    ///
    /// # Errors
    ///
    /// Returns an error if `p` is non-finite or outside `[0, 1]`.
    pub fn from_p(p: f64) -> Result<Self> {
        let freq = AlleleFreq::from_p(p)?;
        GenotypeFreq::from_alleles(&freq)
    }

    /// Sum of the three genotype frequencies. Equals `1` (up to
    /// rounding) for any valid allele-frequency pair.
    pub fn total(&self) -> f64 {
        self.aa_hom + self.het + self.aa
    }

    /// Expected genotype *counts* for a sample of `n` diploid
    /// individuals: each frequency scaled by `n`. Useful as the
    /// "expected" vector of a chi-square goodness-of-fit test.
    pub fn expected_counts(&self, n: f64) -> [f64; 3] {
        [self.aa_hom * n, self.het * n, self.aa * n]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    #[test]
    fn canonical_half_half_case() {
        // p = q = 0.5 => 0.25 / 0.50 / 0.25 (the textbook ground truth).
        let g = GenotypeFreq::from_p(0.5).unwrap();
        assert!((g.aa_hom - 0.25).abs() < EPS, "AA = {}", g.aa_hom);
        assert!((g.het - 0.50).abs() < EPS, "Aa = {}", g.het);
        assert!((g.aa - 0.25).abs() < EPS, "aa = {}", g.aa);
    }

    #[test]
    fn frequencies_sum_to_one() {
        for &p in &[0.0, 0.1, 0.37, 0.5, 0.82, 1.0] {
            let g = GenotypeFreq::from_p(p).unwrap();
            assert!(
                (g.total() - 1.0).abs() < EPS,
                "sum at p={p} = {}",
                g.total()
            );
        }
    }

    #[test]
    fn asymmetric_worked_example() {
        // p = 0.7, q = 0.3 => 0.49 / 0.42 / 0.09 (hand-computed).
        let g = GenotypeFreq::from_p(0.7).unwrap();
        assert!((g.aa_hom - 0.49).abs() < EPS, "AA = {}", g.aa_hom);
        assert!((g.het - 0.42).abs() < EPS, "Aa = {}", g.het);
        assert!((g.aa - 0.09).abs() < EPS, "aa = {}", g.aa);
    }

    #[test]
    fn fixed_allele_gives_single_genotype() {
        // p = 1 => everyone AA.
        let g = GenotypeFreq::from_p(1.0).unwrap();
        assert!((g.aa_hom - 1.0).abs() < EPS);
        assert!(g.het.abs() < EPS);
        assert!(g.aa.abs() < EPS);
    }

    #[test]
    fn heterozygosity_is_maximal_at_half() {
        // 2pq is maximised at p = 0.5, where it equals 0.5.
        let g_half = GenotypeFreq::from_p(0.5).unwrap();
        let g_other = GenotypeFreq::from_p(0.3).unwrap();
        assert!(g_half.het > g_other.het);
        assert!((g_half.het - 0.5).abs() < EPS);
    }

    #[test]
    fn expected_counts_scale_with_sample_size() {
        // p = q = 0.5, n = 100 => 25 / 50 / 25.
        let g = GenotypeFreq::from_p(0.5).unwrap();
        let counts = g.expected_counts(100.0);
        assert!((counts[0] - 25.0).abs() < EPS);
        assert!((counts[1] - 50.0).abs() < EPS);
        assert!((counts[2] - 25.0).abs() < EPS);
    }

    #[test]
    fn from_alleles_matches_from_p() {
        let freq = AlleleFreq::from_p(0.42).unwrap();
        let a = GenotypeFreq::from_alleles(&freq).unwrap();
        let b = GenotypeFreq::from_p(0.42).unwrap();
        assert!((a.aa_hom - b.aa_hom).abs() < EPS);
        assert!((a.het - b.het).abs() < EPS);
        assert!((a.aa - b.aa).abs() < EPS);
    }

    #[test]
    fn from_p_rejects_bad_input() {
        assert!(GenotypeFreq::from_p(1.2).is_err());
        assert!(GenotypeFreq::from_p(f64::NAN).is_err());
    }
}
