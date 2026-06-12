//! Fst — population differentiation.
//!
//! Fst measures how much of a sample's genetic variation is *between*
//! sub-populations rather than within them. `Fst = 0` means panmixia,
//! `Fst = 1` complete isolation (fixed for different alleles).
//!
//! Two standard estimators are provided, both taking two
//! [`GenotypeMatrix`] sub-samples (one per population):
//!
//! - **Hudson's estimator** ([`fst_hudson`]) — `Fst = 1 - Hw / Hb`,
//!   where `Hw` is the mean within-population heterozygosity and `Hb`
//!   the between-population heterozygosity. It is nearly unbiased with
//!   respect to sample size and is the recommended ratio-of-averages
//!   estimator (Bhatia et al. 2013).
//! - **Weir & Cockerham's estimator** ([`fst_weir_cockerham`]) — the
//!   1984 analysis-of-variance estimator `theta`, which partitions
//!   allele-frequency variance into among-population (`a`) and
//!   within-population (`b`, `c`) components.
//!
//! Both return a *genome-wide* estimate as a ratio of summed numerators
//! to summed denominators — the statistically sound way to average Fst.

use crate::error::{PopgenError, Result};
use crate::infer::GenotypeMatrix;

/// Hudson's Fst estimator for two populations.
///
/// For each shared site, with derived frequencies `p1`, `p2` and sample
/// sizes `n1`, `n2`:
///
/// - numerator `= (p1 - p2)^2 - p1(1-p1)/(n1-1) - p2(1-p2)/(n2-1)`
/// - denominator `= p1(1-p2) + p2(1-p1)`
///
/// `Fst` is `sum(numerator) / sum(denominator)`.
///
/// # Errors
/// [`PopgenError::Invalid`] if either population has fewer than two
/// samples; [`PopgenError::Dimension`] if the two matrices differ in
/// site count.
pub fn fst_hudson(pop1: &GenotypeMatrix, pop2: &GenotypeMatrix) -> Result<f64> {
    let (n1, n2) = check_pair(pop1, pop2)?;
    let mut num = 0.0;
    let mut den = 0.0;
    for col in 0..pop1.n_sites() {
        let p1 = pop1.frequency(col)?;
        let p2 = pop2.frequency(col)?;
        let term = (p1 - p2).powi(2)
            - p1 * (1.0 - p1) / (n1 as f64 - 1.0)
            - p2 * (1.0 - p2) / (n2 as f64 - 1.0);
        num += term;
        den += p1 * (1.0 - p2) + p2 * (1.0 - p1);
    }
    Ok(if den.abs() < 1e-12 { 0.0 } else { num / den })
}

/// Weir & Cockerham's (1984) Fst estimator `theta` for two populations.
///
/// For each shared site this accumulates the among-population variance
/// component `a` and the within-population components `b` and `c`;
/// `theta = sum(a) / sum(a + b + c)`.
///
/// # Errors
/// [`PopgenError::Invalid`] if either population has fewer than two
/// samples; [`PopgenError::Dimension`] if the matrices differ in site
/// count.
pub fn fst_weir_cockerham(pop1: &GenotypeMatrix, pop2: &GenotypeMatrix) -> Result<f64> {
    let (n1u, n2u) = check_pair(pop1, pop2)?;
    let n1 = n1u as f64;
    let n2 = n2u as f64;
    let r = 2.0; // number of populations
    let n_bar = (n1 + n2) / r;
    // Squared coefficient of variation of sample sizes.
    let nc = (n1 + n2 - (n1 * n1 + n2 * n2) / (n1 + n2)) / (r - 1.0);

    let mut sum_a = 0.0;
    let mut sum_abc = 0.0;
    for col in 0..pop1.n_sites() {
        let p1 = pop1.frequency(col)?;
        let p2 = pop2.frequency(col)?;
        // Sample-size-weighted average allele frequency.
        let p_bar = (n1 * p1 + n2 * p2) / (n1 + n2);
        // Variance of allele frequency among populations (weighted).
        let s2 = (n1 * (p1 - p_bar).powi(2) + n2 * (p2 - p_bar).powi(2)) / ((r - 1.0) * n_bar);
        // Average observed heterozygosity. With only allele counts
        // (haplotype rows) the observed heterozygosity is taken as the
        // expected 2p(1-p) under random union — the standard fallback
        // when individual genotype phase is unavailable.
        let h_bar = (n1 * 2.0 * p1 * (1.0 - p1) + n2 * 2.0 * p2 * (1.0 - p2)) / (n1 + n2);
        // Weir & Cockerham 1984 variance components.
        let a = (n_bar / nc)
            * (s2
                - (1.0 / (n_bar - 1.0))
                    * (p_bar * (1.0 - p_bar) - ((r - 1.0) / r) * s2 - 0.25 * h_bar));
        let b = (n_bar / (n_bar - 1.0))
            * (p_bar * (1.0 - p_bar)
                - ((r - 1.0) / r) * s2
                - ((2.0 * n_bar - 1.0) / (4.0 * n_bar)) * h_bar);
        let c = 0.5 * h_bar;
        sum_a += a;
        sum_abc += a + b + c;
    }
    Ok(if sum_abc.abs() < 1e-12 {
        0.0
    } else {
        sum_a / sum_abc
    })
}

/// Validates a population pair and returns their sample sizes.
fn check_pair(pop1: &GenotypeMatrix, pop2: &GenotypeMatrix) -> Result<(usize, usize)> {
    let n1 = pop1.n_samples();
    let n2 = pop2.n_samples();
    if n1 < 2 || n2 < 2 {
        return Err(PopgenError::invalid(
            "population",
            "each population needs at least two samples for Fst",
        ));
    }
    if pop1.n_sites() != pop2.n_sites() {
        return Err(PopgenError::dimension(
            pop1.n_sites(),
            pop2.n_sites(),
            "Fst population site counts",
        ));
    }
    Ok((n1, n2))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matrix(rows: Vec<Vec<u8>>) -> GenotypeMatrix {
        let cols = rows[0].len();
        let pos: Vec<f64> = (0..cols).map(|i| i as f64).collect();
        GenotypeMatrix::from_rows(rows, pos).unwrap()
    }

    #[test]
    fn identical_populations_have_near_zero_fst() {
        // Two populations with the same allele frequencies. Hudson's
        // estimator carries an explicit sample-size bias correction
        // p(1-p)/(n-1), so on a *tiny* sample it can come out
        // appreciably negative for identical populations — that is
        // correct behaviour of the unbiased estimator, not a bug. A
        // meaningful "near zero" test therefore needs a sample large
        // enough for the bias term to be small: 20 samples per
        // population, each at allele frequency 0.5.
        let pop = |n: usize| {
            let mut rows: Vec<Vec<u8>> = Vec::new();
            for r in 0..n {
                // half the samples carry the derived allele at each site
                let bit = if r < n / 2 { 1 } else { 0 };
                rows.push(vec![bit, 1 - bit]);
            }
            matrix(rows)
        };
        let p1 = pop(20);
        let p2 = pop(20);
        let fst_h = fst_hudson(&p1, &p2).unwrap();
        let fst_wc = fst_weir_cockerham(&p1, &p2).unwrap();
        assert!(fst_h.abs() < 0.1, "Hudson Fst = {fst_h}");
        assert!(fst_wc.abs() < 0.1, "WC Fst = {fst_wc}");
    }

    #[test]
    fn fixed_for_different_alleles_gives_high_fst() {
        // Population 1 all derived, population 2 all ancestral.
        let p1 = matrix(vec![vec![1], vec![1], vec![1], vec![1]]);
        let p2 = matrix(vec![vec![0], vec![0], vec![0], vec![0]]);
        let fst_h = fst_hudson(&p1, &p2).unwrap();
        assert!(fst_h > 0.8, "expected high Fst, got {fst_h}");
        let fst_wc = fst_weir_cockerham(&p1, &p2).unwrap();
        assert!(fst_wc > 0.8, "expected high WC Fst, got {fst_wc}");
    }

    #[test]
    fn partial_differentiation_is_intermediate() {
        // Pop 1 mostly derived, pop 2 mostly ancestral.
        let p1 = matrix(vec![vec![1], vec![1], vec![1], vec![0]]);
        let p2 = matrix(vec![vec![0], vec![0], vec![0], vec![1]]);
        let fst = fst_hudson(&p1, &p2).unwrap();
        assert!(fst > 0.0 && fst < 1.0, "Fst = {fst}");
    }

    #[test]
    fn rejects_mismatched_site_counts() {
        let p1 = matrix(vec![vec![1, 0], vec![0, 1]]);
        let p2 = matrix(vec![vec![1], vec![0]]);
        assert!(fst_hudson(&p1, &p2).is_err());
        assert!(fst_weir_cockerham(&p1, &p2).is_err());
    }

    #[test]
    fn rejects_singleton_populations() {
        let p1 = matrix(vec![vec![1]]);
        let p2 = matrix(vec![vec![0]]);
        assert!(fst_hudson(&p1, &p2).is_err());
    }

    #[test]
    fn fst_estimators_match_hand_computed_values() {
        // GROUND TRUTH. Two 4-sample populations, one shared site: pop1 has
        // derived frequency p1 = 3/4, pop2 p2 = 1/4. The existing Fst tests
        // only check sign/range (≈0 / >0.8 / in-(0,1)); this pins BOTH
        // estimators to exact hand-computed values.
        //
        // Hudson (Bhatia 2013): num = (p1−p2)² − p1(1−p1)/(n1−1) − p2(1−p2)/(n2−1)
        //   = 0.25 − 0.1875/3 − 0.1875/3 = 0.25 − 0.0625 − 0.0625 = 0.125;
        //   den = p1(1−p2) + p2(1−p1) = 0.5625 + 0.0625 = 0.625;
        //   Fst = 0.125/0.625 = 0.2.
        // Weir & Cockerham (1984), r=2, n̄=4, nc=4, p̄=0.5, s²=0.125, h̄=0.375
        //   → a = 0.09375, b = 0.03125, c = 0.1875;
        //   θ = a/(a+b+c) = 0.09375/0.3125 = 0.3.
        let p1 = matrix(vec![vec![1], vec![1], vec![1], vec![0]]); // p = 3/4
        let p2 = matrix(vec![vec![1], vec![0], vec![0], vec![0]]); // p = 1/4
        let fst_h = fst_hudson(&p1, &p2).unwrap();
        assert!((fst_h - 0.2).abs() < 1e-12, "Hudson Fst {fst_h} != 0.2");
        let fst_wc = fst_weir_cockerham(&p1, &p2).unwrap();
        assert!(
            (fst_wc - 0.3).abs() < 1e-9,
            "Weir-Cockerham Fst {fst_wc} != 0.3"
        );
    }
}
