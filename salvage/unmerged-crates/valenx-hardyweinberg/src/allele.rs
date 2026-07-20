//! Allele-frequency estimation for a single di-allelic locus.
//!
//! For a locus with two alleles `A` and `a`, a diploid individual is
//! one of three genotypes: homozygous `AA`, heterozygous `Aa`, or
//! homozygous `aa`. Given counts of each genotype in a sample of `N`
//! individuals (so `2N` allele copies), the maximum-likelihood allele
//! frequencies follow by gene counting:
//!
//! ```text
//! p = freq(A) = (2 * n_AA + n_Aa) / (2N)
//! q = freq(a) = (2 * n_aa + n_Aa) / (2N)
//! ```
//!
//! and by construction `p + q = 1` exactly (every allele copy is either
//! `A` or `a`).

use crate::error::{HweError, Result};

/// A di-allelic locus summarised by allele frequencies.
///
/// Invariant on a value returned by this crate: `p` and `q` are each in
/// `[0, 1]` and `p + q == 1` up to floating-point rounding.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct AlleleFreq {
    /// Frequency of the reference allele `A`.
    pub p: f64,
    /// Frequency of the alternate allele `a`.
    pub q: f64,
}

impl AlleleFreq {
    /// Build allele frequencies from genotype counts at a single
    /// di-allelic locus.
    ///
    /// `n_aa_hom` is the count of `AA` homozygotes, `n_het` the count of
    /// `Aa` heterozygotes, and `n_aa` the count of `aa` homozygotes.
    /// Counts are `f64` so weighted / expected counts are accepted, but
    /// each must be finite and non-negative and their sum must be
    /// strictly positive.
    ///
    /// # Errors
    ///
    /// Returns [`HweError::NotFinite`] if any count is `NaN` / infinite,
    /// and [`HweError::OutOfDomain`] if a count is negative or the total
    /// number of individuals is zero.
    pub fn from_genotype_counts(n_aa_hom: f64, n_het: f64, n_aa: f64) -> Result<Self> {
        let counts = [("n_AA", n_aa_hom), ("n_Aa", n_het), ("n_aa", n_aa)];
        for (label, c) in counts {
            if !c.is_finite() {
                return Err(HweError::not_finite(label));
            }
            if c < 0.0 {
                return Err(HweError::out_of_domain(
                    label,
                    format!("genotype count {c} is negative"),
                ));
            }
        }

        let total_individuals = n_aa_hom + n_het + n_aa;
        if total_individuals <= 0.0 {
            return Err(HweError::out_of_domain(
                "total",
                "sample has zero individuals; cannot estimate allele frequencies",
            ));
        }

        let total_alleles = 2.0 * total_individuals;
        let a_copies = 2.0 * n_aa_hom + n_het;
        let p = a_copies / total_alleles;
        let q = 1.0 - p;
        Ok(AlleleFreq { p, q })
    }

    /// Build allele frequencies directly from the reference-allele
    /// frequency `p`. The alternate frequency is `q = 1 - p`.
    ///
    /// # Errors
    ///
    /// Returns an error if `p` is non-finite or outside `[0, 1]`.
    pub fn from_p(p: f64) -> Result<Self> {
        let p = crate::error::unit_interval("p", p)?;
        Ok(AlleleFreq { p, q: 1.0 - p })
    }
}

/// Frequency of the major (more common) allele.
///
/// Convenience wrapper: returns `max(p, q)` for an already-validated
/// [`AlleleFreq`]. Always in `[0.5, 1]`.
pub fn major_allele_frequency(freq: &AlleleFreq) -> f64 {
    freq.p.max(freq.q)
}

/// Frequency of the minor (less common) allele — the classic "MAF".
///
/// Returns `min(p, q)`. Always in `[0, 0.5]`.
pub fn minor_allele_frequency(freq: &AlleleFreq) -> f64 {
    freq.p.min(freq.q)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    #[test]
    fn equal_homozygotes_give_half_half() {
        // 25 AA, 50 Aa, 25 aa: A copies = 2*25 + 50 = 100, total = 200
        // => p = 0.5, ground-truth allele symmetry.
        let f = AlleleFreq::from_genotype_counts(25.0, 50.0, 25.0).unwrap();
        assert!((f.p - 0.5).abs() < EPS, "p = {}", f.p);
        assert!((f.q - 0.5).abs() < EPS, "q = {}", f.q);
    }

    #[test]
    fn textbook_worked_example() {
        // Hartl & Clark style worked case: 35 AA, 20 Aa, 45 aa (N = 100).
        // A copies = 2*35 + 20 = 90 of 200 => p = 0.45, q = 0.55.
        let f = AlleleFreq::from_genotype_counts(35.0, 20.0, 45.0).unwrap();
        assert!((f.p - 0.45).abs() < EPS, "p = {}", f.p);
        assert!((f.q - 0.55).abs() < EPS, "q = {}", f.q);
    }

    #[test]
    fn p_and_q_always_sum_to_one() {
        let f = AlleleFreq::from_genotype_counts(13.0, 7.0, 31.0).unwrap();
        assert!((f.p + f.q - 1.0).abs() < EPS, "p + q = {}", f.p + f.q);
    }

    #[test]
    fn fixed_allele_limit() {
        // No a alleles anywhere => p = 1, q = 0.
        let f = AlleleFreq::from_genotype_counts(40.0, 0.0, 0.0).unwrap();
        assert!((f.p - 1.0).abs() < EPS);
        assert!(f.q.abs() < EPS);
    }

    #[test]
    fn from_p_round_trips() {
        let f = AlleleFreq::from_p(0.3).unwrap();
        assert!((f.p - 0.3).abs() < EPS);
        assert!((f.q - 0.7).abs() < EPS);
    }

    #[test]
    fn major_and_minor_split_correctly() {
        let f = AlleleFreq::from_genotype_counts(35.0, 20.0, 45.0).unwrap(); // p=0.45,q=0.55
        assert!((major_allele_frequency(&f) - 0.55).abs() < EPS);
        assert!((minor_allele_frequency(&f) - 0.45).abs() < EPS);
    }

    #[test]
    fn rejects_negative_count() {
        let err = AlleleFreq::from_genotype_counts(-1.0, 5.0, 5.0).unwrap_err();
        assert_eq!(err.code(), "hwe.out_of_domain");
    }

    #[test]
    fn rejects_nan_count() {
        let err = AlleleFreq::from_genotype_counts(f64::NAN, 5.0, 5.0).unwrap_err();
        assert_eq!(err.code(), "hwe.not_finite");
    }

    #[test]
    fn rejects_empty_sample() {
        let err = AlleleFreq::from_genotype_counts(0.0, 0.0, 0.0).unwrap_err();
        assert_eq!(err.code(), "hwe.out_of_domain");
    }

    #[test]
    fn from_p_rejects_out_of_range() {
        assert!(AlleleFreq::from_p(1.5).is_err());
        assert!(AlleleFreq::from_p(-0.1).is_err());
        assert!(AlleleFreq::from_p(f64::INFINITY).is_err());
    }
}
