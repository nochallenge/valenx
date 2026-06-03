//! The Hardy-Weinberg equilibrium test.
//!
//! Under Hardy-Weinberg equilibrium (HWE) — random mating, no
//! selection, drift, mutation or migration at the locus — a biallelic
//! locus with derived-allele frequency `p` has genotype frequencies
//! `(1-p)^2 : 2p(1-p) : p^2` for the AA / Aa / aa classes.
//!
//! Given the *observed* diploid genotype counts `(n_AA, n_Aa, n_aa)`
//! this module tests whether they are consistent with HWE:
//!
//! - [`hwe_chi_square`] — Pearson's chi-square goodness-of-fit against
//!   the expected counts (1 degree of freedom), with the asymptotic
//!   p-value. Fast, standard, but unreliable for small samples or rare
//!   alleles.
//! - [`hwe_exact`] — the exact test (Wigginton, Cutler & Abecasis
//!   2005): the probability of the observed heterozygote count, plus
//!   every *less probable* configuration with the same allele count.
//!   Correct for any sample size — the recommended test.

use crate::error::{PopgenError, Result};

/// Observed diploid genotype counts at a biallelic locus.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct GenotypeCounts {
    /// Homozygous-ancestral (`AA`) individuals.
    pub aa: usize,
    /// Heterozygous (`Aa`) individuals.
    pub ab: usize,
    /// Homozygous-derived (`bb`) individuals.
    pub bb: usize,
}

impl GenotypeCounts {
    /// Total number of genotyped individuals.
    pub fn total(&self) -> usize {
        self.aa + self.ab + self.bb
    }

    /// Derived-allele frequency `p`.
    pub fn derived_frequency(&self) -> f64 {
        let alleles = 2 * self.total();
        if alleles == 0 {
            return 0.0;
        }
        (2 * self.bb + self.ab) as f64 / alleles as f64
    }
}

/// Result of a Hardy-Weinberg test.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct HweResult {
    /// The test statistic (chi-square value, or the observed-config
    /// probability for the exact test).
    pub statistic: f64,
    /// The p-value: probability of a result this extreme under HWE.
    pub p_value: f64,
}

impl HweResult {
    /// `true` if HWE is rejected at significance level `alpha`.
    pub fn rejects_at(&self, alpha: f64) -> bool {
        self.p_value < alpha
    }
}

/// Pearson chi-square goodness-of-fit test for HWE.
///
/// The expected counts are derived from the observed allele frequency;
/// the statistic has 1 degree of freedom. The p-value is the upper
/// tail of the chi-square(1) distribution.
///
/// # Errors
/// [`PopgenError::Invalid`] if no individuals are supplied.
pub fn hwe_chi_square(counts: GenotypeCounts) -> Result<HweResult> {
    let n = counts.total();
    if n == 0 {
        return Err(PopgenError::invalid("counts", "no genotyped individuals"));
    }
    let p = counts.derived_frequency();
    let q = 1.0 - p;
    let nf = n as f64;
    let exp_aa = q * q * nf;
    let exp_ab = 2.0 * p * q * nf;
    let exp_bb = p * p * nf;

    let term = |obs: usize, exp: f64| -> f64 {
        if exp < 1e-12 {
            0.0
        } else {
            let d = obs as f64 - exp;
            d * d / exp
        }
    };
    let chi2 =
        term(counts.aa, exp_aa) + term(counts.ab, exp_ab) + term(counts.bb, exp_bb);
    Ok(HweResult {
        statistic: chi2,
        p_value: chi_square_sf_1df(chi2),
    })
}

/// The exact Hardy-Weinberg test (Wigginton-Cutler-Abecasis 2005).
///
/// Conditioning on the observed allele counts, this enumerates every
/// possible heterozygote count, computes the probability of each
/// configuration, and sums the probabilities of all configurations no
/// more probable than the observed one. That sum is the exact p-value.
///
/// # Errors
/// [`PopgenError::Invalid`] if no individuals are supplied.
pub fn hwe_exact(counts: GenotypeCounts) -> Result<HweResult> {
    let n = counts.total();
    if n == 0 {
        return Err(PopgenError::invalid("counts", "no genotyped individuals"));
    }
    // Minor- and major-allele counts.
    let derived = 2 * counts.bb + counts.ab;
    let ancestral = 2 * counts.aa + counts.ab;
    let n_rare = derived.min(ancestral);
    let n_common = derived.max(ancestral);
    let obs_het = counts.ab;

    // Probability (up to a shared constant) of `het` heterozygotes
    // given the allele counts, via the recurrence in WCA 2005.
    // het and the rare-homozygote count have the same parity as
    // n_rare; iterate het over that grid.
    let mut probs: Vec<(usize, f64)> = Vec::new();
    // Start from het = n_rare mod 2 and step by 2.
    let mut het = n_rare % 2;
    // Use log-space weights to avoid overflow for large n.
    while het <= n_rare {
        let rare_hom = (n_rare - het) / 2;
        let common_hom = (n_common - het) / 2;
        // log multinomial weight: it suffices to compare relative
        // weights, so use ln of the standard HWE configuration count.
        let log_w = ln_factorial(n)
            - ln_factorial(rare_hom)
            - ln_factorial(het)
            - ln_factorial(common_hom)
            + het as f64 * std::f64::consts::LN_2;
        probs.push((het, log_w));
        het += 2;
    }
    // Convert log-weights to normalised probabilities.
    let max_log = probs
        .iter()
        .map(|&(_, w)| w)
        .fold(f64::NEG_INFINITY, f64::max);
    let weights: Vec<(usize, f64)> = probs
        .iter()
        .map(|&(h, w)| (h, (w - max_log).exp()))
        .collect();
    let total: f64 = weights.iter().map(|&(_, w)| w).sum();
    let obs_prob = weights
        .iter()
        .find(|&&(h, _)| h == obs_het)
        .map(|&(_, w)| w / total)
        .unwrap_or(0.0);
    // p-value: sum of probabilities of configurations no more likely
    // than the observed one.
    let p_value: f64 = weights
        .iter()
        .map(|&(_, w)| w / total)
        .filter(|&pr| pr <= obs_prob + 1e-12)
        .sum();
    Ok(HweResult {
        statistic: obs_prob,
        p_value: p_value.min(1.0),
    })
}

/// Survival function (upper tail) of the chi-square distribution with
/// one degree of freedom.
///
/// For 1 d.f., `chi2 = z^2`, so `P(X > chi2) = 2 * (1 - Phi(sqrt(chi2)))`
/// = `erfc(sqrt(chi2 / 2))`.
fn chi_square_sf_1df(chi2: f64) -> f64 {
    if chi2 <= 0.0 {
        return 1.0;
    }
    erfc((chi2 / 2.0).sqrt())
}

/// Complementary error function via Abramowitz & Stegun 7.1.26 — a
/// rational approximation accurate to ~1e-7.
fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.327_591_1 * z);
    let poly = t
        * (0.254_829_592
            + t * (-0.284_496_736
                + t * (1.421_413_741
                    + t * (-1.453_152_027 + t * 1.061_405_429))));
    let approx = poly * (-z * z).exp();
    if x >= 0.0 {
        approx
    } else {
        2.0 - approx
    }
}

/// Natural log of `n!` via `ln_gamma(n + 1)`.
fn ln_factorial(n: usize) -> f64 {
    ln_gamma(n as f64 + 1.0)
}

/// Lanczos approximation to `ln(Gamma(x))` for `x > 0`.
fn ln_gamma(x: f64) -> f64 {
    // Lanczos g = 7, n = 9 coefficients.
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    let g = 7.0;
    if x < 0.5 {
        // Reflection formula.
        std::f64::consts::PI.ln()
            - (std::f64::consts::PI * x).sin().abs().ln()
            - ln_gamma(1.0 - x)
    } else {
        let x = x - 1.0;
        let mut a = C[0];
        let t = x + g + 0.5;
        for (i, &c) in C.iter().enumerate().skip(1) {
            a += c / (x + i as f64);
        }
        0.5 * (2.0 * std::f64::consts::PI).ln()
            + (x + 0.5) * t.ln()
            - t
            + a.ln()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perfect_hwe_has_a_high_p_value() {
        // p = 0.5, n = 100: expected (25, 50, 25). Feed exactly that.
        let counts = GenotypeCounts {
            aa: 25,
            ab: 50,
            bb: 25,
        };
        let chi = hwe_chi_square(counts).unwrap();
        assert!(chi.statistic < 1e-6, "chi2 = {}", chi.statistic);
        assert!(chi.p_value > 0.99, "p = {}", chi.p_value);
        let exact = hwe_exact(counts).unwrap();
        assert!(exact.p_value > 0.5, "exact p = {}", exact.p_value);
    }

    #[test]
    fn strong_heterozygote_deficit_is_rejected() {
        // p = 0.5 but almost no heterozygotes -> strong HWE departure.
        let counts = GenotypeCounts {
            aa: 50,
            ab: 0,
            bb: 50,
        };
        let chi = hwe_chi_square(counts).unwrap();
        assert!(chi.p_value < 0.001, "p = {}", chi.p_value);
        assert!(chi.rejects_at(0.05));
        let exact = hwe_exact(counts).unwrap();
        assert!(exact.p_value < 0.001, "exact p = {}", exact.p_value);
    }

    #[test]
    fn heterozygote_excess_is_also_rejected() {
        // All heterozygotes -> impossible under HWE for p = 0.5.
        let counts = GenotypeCounts {
            aa: 0,
            ab: 100,
            bb: 0,
        };
        let chi = hwe_chi_square(counts).unwrap();
        assert!(chi.p_value < 0.01, "p = {}", chi.p_value);
    }

    #[test]
    fn derived_frequency_is_correct() {
        let counts = GenotypeCounts {
            aa: 10,
            ab: 4,
            bb: 6,
        };
        // derived alleles = 2*6 + 4 = 16, total alleles = 40 -> 0.4.
        assert!((counts.derived_frequency() - 0.4).abs() < 1e-12);
        assert_eq!(counts.total(), 20);
    }

    #[test]
    fn exact_p_value_is_a_probability() {
        for counts in [
            GenotypeCounts {
                aa: 5,
                ab: 10,
                bb: 5,
            },
            GenotypeCounts {
                aa: 1,
                ab: 2,
                bb: 7,
            },
            GenotypeCounts {
                aa: 30,
                ab: 1,
                bb: 0,
            },
        ] {
            let r = hwe_exact(counts).unwrap();
            assert!(
                (0.0..=1.0).contains(&r.p_value),
                "p out of range: {}",
                r.p_value
            );
        }
    }

    #[test]
    fn rejects_empty_input() {
        let empty = GenotypeCounts {
            aa: 0,
            ab: 0,
            bb: 0,
        };
        assert!(hwe_chi_square(empty).is_err());
        assert!(hwe_exact(empty).is_err());
    }

    #[test]
    fn ln_gamma_matches_known_factorials() {
        // ln(5!) = ln(120).
        assert!((ln_factorial(5) - 120.0_f64.ln()).abs() < 1e-6);
        // ln(0!) = 0.
        assert!(ln_factorial(0).abs() < 1e-6);
    }
}
