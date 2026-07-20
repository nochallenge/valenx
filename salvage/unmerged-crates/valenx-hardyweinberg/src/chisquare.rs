//! Chi-square goodness-of-fit test for Hardy-Weinberg equilibrium.
//!
//! Given the *observed* genotype counts `O = (n_AA, n_Aa, n_aa)` from a
//! sample, the test asks whether they are consistent with the
//! Hardy-Weinberg *expected* counts `E` computed from the allele
//! frequencies estimated from those same data. The Pearson statistic is
//!
//! ```text
//! chi2 = sum_i (O_i - E_i)^2 / E_i
//! ```
//!
//! summed over the three genotype classes. Because the allele frequency
//! `p` is estimated from the data, one degree of freedom is spent on
//! it, leaving
//!
//! ```text
//! df = (classes - 1) - (params estimated) = (3 - 1) - 1 = 1
//! ```
//!
//! for a single di-allelic locus. The upper-tail p-value is the
//! probability that a chi-square random variable with `df` degrees of
//! freedom exceeds the observed statistic; for `df = 1` it has the
//! closed form `p = erfc(sqrt(chi2 / 2))`.

use crate::allele::AlleleFreq;
use crate::error::{HweError, Result};
use crate::genotype::GenotypeFreq;

/// Result of a Hardy-Weinberg chi-square goodness-of-fit test.
#[derive(Clone, Debug, PartialEq)]
pub struct HweTest {
    /// Pearson chi-square statistic, `sum (O - E)^2 / E`.
    pub chi_square: f64,
    /// Degrees of freedom (`1` for a di-allelic locus with `p`
    /// estimated from the data).
    pub degrees_of_freedom: u32,
    /// Upper-tail p-value: `P(X > chi_square)` for `X ~ chi-square(df)`.
    pub p_value: f64,
    /// Observed genotype counts in the order `[AA, Aa, aa]`.
    pub observed: [f64; 3],
    /// Hardy-Weinberg expected counts in the order `[AA, Aa, aa]`.
    pub expected: [f64; 3],
}

impl HweTest {
    /// Whether the null hypothesis of Hardy-Weinberg equilibrium is
    /// *rejected* at significance level `alpha` (i.e. `p_value < alpha`).
    ///
    /// A `true` result means the observed genotype proportions deviate
    /// from equilibrium more than sampling alone would plausibly
    /// explain. Common choices for `alpha` are `0.05` or `0.01`.
    pub fn rejects_hwe(&self, alpha: f64) -> bool {
        self.p_value < alpha
    }
}

/// Run the Hardy-Weinberg chi-square goodness-of-fit test on observed
/// genotype counts for one di-allelic locus.
///
/// `obs_aa_hom`, `obs_het`, and `obs_aa` are the observed counts of the
/// `AA`, `Aa`, and `aa` genotypes. The allele frequency `p` is
/// estimated from these counts; the expected counts are then `N * p^2`,
/// `N * 2pq`, `N * q^2`.
///
/// # Errors
///
/// Returns [`HweError::NotFinite`] / [`HweError::OutOfDomain`] for
/// non-finite or negative counts or an empty sample (via
/// [`AlleleFreq::from_genotype_counts`]), and
/// [`HweError::NoDegreesOfFreedom`] when one allele is fixed so an
/// expected cell is zero and the statistic is undefined.
pub fn hwe_chi_square_test(obs_aa_hom: f64, obs_het: f64, obs_aa: f64) -> Result<HweTest> {
    let freq = AlleleFreq::from_genotype_counts(obs_aa_hom, obs_het, obs_aa)?;
    let n = obs_aa_hom + obs_het + obs_aa;
    let expected = GenotypeFreq::from_alleles(&freq)?.expected_counts(n);
    let observed = [obs_aa_hom, obs_het, obs_aa];

    // With a fixed allele (p = 0 or p = 1) two expected cells are zero;
    // the locus is monomorphic and the goodness-of-fit test has no
    // information to test, so df collapses. Reject explicitly rather
    // than dividing by zero.
    if freq.p <= 0.0 || freq.q <= 0.0 {
        return Err(HweError::no_degrees_of_freedom(
            "one allele is fixed (monomorphic locus); HWE chi-square is undefined",
        ));
    }

    let mut chi_square = 0.0;
    for i in 0..3 {
        let e = expected[i];
        let diff = observed[i] - e;
        chi_square += diff * diff / e;
    }

    let degrees_of_freedom = 1u32;
    let p_value = chi_square_upper_tail(chi_square, degrees_of_freedom)?;

    Ok(HweTest {
        chi_square,
        degrees_of_freedom,
        p_value,
        observed,
        expected,
    })
}

/// Upper-tail probability `P(X > x)` for `X ~ chi-square(df)`.
///
/// Equals the regularised upper incomplete gamma function
/// `Q(df/2, x/2)`. Exposed publicly so callers can convert their own
/// chi-square statistics (with arbitrary positive `df`) into p-values.
///
/// # Errors
///
/// Returns [`HweError::NotFinite`] if `x` is `NaN` / infinite,
/// [`HweError::OutOfDomain`] if `x` is negative, and
/// [`HweError::NoDegreesOfFreedom`] if `df == 0`.
pub fn chi_square_upper_tail(x: f64, df: u32) -> Result<f64> {
    let x = crate::error::finite("chi_square", x)?;
    if x < 0.0 {
        return Err(HweError::out_of_domain(
            "chi_square",
            format!("statistic {x} is negative"),
        ));
    }
    if df == 0 {
        return Err(HweError::no_degrees_of_freedom("df must be at least 1"));
    }
    if x == 0.0 {
        return Ok(1.0);
    }
    Ok(gammq(0.5 * f64::from(df), 0.5 * x))
}

/// Regularised upper incomplete gamma function `Q(a, x) = 1 - P(a, x)`.
///
/// Numerical Recipes recipe: series expansion for `x < a + 1`,
/// continued fraction otherwise. Valid for `a > 0`, `x >= 0`.
fn gammq(a: f64, x: f64) -> f64 {
    if x < a + 1.0 {
        1.0 - gamma_series(a, x)
    } else {
        gamma_continued_fraction(a, x)
    }
}

/// Lower regularised incomplete gamma `P(a, x)` via its series.
fn gamma_series(a: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    let ln_gamma_a = ln_gamma(a);
    let mut ap = a;
    let mut sum = 1.0 / a;
    let mut del = sum;
    for _ in 0..1000 {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if del.abs() < sum.abs() * 1e-15 {
            break;
        }
    }
    sum * (-x + a * x.ln() - ln_gamma_a).exp()
}

/// Upper regularised incomplete gamma `Q(a, x)` via its continued
/// fraction (Lentz's algorithm).
fn gamma_continued_fraction(a: f64, x: f64) -> f64 {
    let ln_gamma_a = ln_gamma(a);
    let tiny = 1e-300;
    let mut b = x + 1.0 - a;
    let mut c = 1.0 / tiny;
    let mut d = 1.0 / b;
    let mut h = d;
    for i in 1..1000 {
        let an = -(i as f64) * (i as f64 - a);
        b += 2.0;
        d = an * d + b;
        if d.abs() < tiny {
            d = tiny;
        }
        c = b + an / c;
        if c.abs() < tiny {
            c = tiny;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-15 {
            break;
        }
    }
    (-x + a * x.ln() - ln_gamma_a).exp() * h
}

/// Natural log of the gamma function, `ln Gamma(z)`, via the Lanczos
/// approximation (`g = 7`, 9 coefficients). Accurate to ~1e-15 for the
/// `z > 0` range used here.
fn ln_gamma(z: f64) -> f64 {
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
    // Reflection is unnecessary: all callers pass z = df/2 >= 0.5.
    let z = z - 1.0;
    let mut a = C[0];
    let t = z + g + 0.5;
    for (i, &coeff) in C.iter().enumerate().skip(1) {
        a += coeff / (z + i as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (z + 0.5) * t.ln() - t + a.ln()
}

/// Complementary error function `erfc(x)`, used only by tests as an
/// independent ground-truth for the `df = 1` p-value
/// (`P(chi2 > x) = erfc(sqrt(x / 2))`). Abramowitz & Stegun 7.1.26
/// rational approximation; max abs error ~1.5e-7.
#[cfg(test)]
fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let tau = t
        * (-z * z - 1.265_512_23
            + t * (1.000_023_68
                + t * (0.374_091_96
                    + t * (0.096_784_18
                        + t * (-0.186_288_06
                            + t * (0.278_868_07
                                + t * (-1.135_203_98
                                    + t * (1.488_515_87
                                        + t * (-0.822_152_23 + t * 0.170_872_77)))))))))
            .exp();
    if x >= 0.0 {
        tau
    } else {
        2.0 - tau
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    #[test]
    fn perfect_equilibrium_gives_zero_statistic() {
        // Observed exactly equals expected for p = q = 0.5:
        // 25 / 50 / 25. Then chi2 = 0 and p-value = 1.
        let t = hwe_chi_square_test(25.0, 50.0, 25.0).unwrap();
        assert!(t.chi_square.abs() < 1e-9, "chi2 = {}", t.chi_square);
        assert!((t.p_value - 1.0).abs() < EPS, "p = {}", t.p_value);
        assert_eq!(t.degrees_of_freedom, 1);
        assert!(!t.rejects_hwe(0.05));
    }

    #[test]
    fn expected_counts_are_hwe_proportions() {
        // 25/50/25 already at equilibrium -> expected == observed.
        let t = hwe_chi_square_test(25.0, 50.0, 25.0).unwrap();
        for i in 0..3 {
            assert!((t.expected[i] - t.observed[i]).abs() < 1e-9);
        }
    }

    #[test]
    fn hand_worked_chi_square_value() {
        // Classic textbook deviation case. Observed 30 AA, 30 Aa, 40 aa
        // (N = 100). A copies = 2*30 + 30 = 90 / 200 => p = 0.45, q = 0.55.
        // Expected: AA = 100*0.2025 = 20.25, Aa = 100*0.495 = 49.5,
        //           aa = 100*0.3025 = 30.25.
        // chi2 = (30-20.25)^2/20.25 + (30-49.5)^2/49.5 + (40-30.25)^2/30.25
        //      = 4.6944... + 7.6818... + 3.1425...
        //      = 15.5188 (hand-computed).
        let t = hwe_chi_square_test(30.0, 30.0, 40.0).unwrap();
        assert!(
            (t.expected[0] - 20.25).abs() < EPS,
            "E_AA = {}",
            t.expected[0]
        );
        assert!(
            (t.expected[1] - 49.5).abs() < EPS,
            "E_Aa = {}",
            t.expected[1]
        );
        assert!(
            (t.expected[2] - 30.25).abs() < EPS,
            "E_aa = {}",
            t.expected[2]
        );
        assert!(
            (t.chi_square - 15.518_793).abs() < 1e-4,
            "chi2 = {}",
            t.chi_square
        );
        // chi2 ~ 15.5 at df=1 is far in the tail -> reject HWE.
        assert!(t.rejects_hwe(0.05));
        assert!(t.rejects_hwe(0.001));
    }

    #[test]
    fn excess_heterozygosity_is_detected() {
        // 5 AA, 90 Aa, 5 aa (N = 100): p = q = 0.5, expected 25/50/25,
        // but huge heterozygote excess. chi2 = (5-25)^2/25 * 2 +
        // (90-50)^2/50 = 16 + 16 + 32 = 64 (hand-computed).
        let t = hwe_chi_square_test(5.0, 90.0, 5.0).unwrap();
        assert!(
            (t.chi_square - 64.0).abs() < 1e-6,
            "chi2 = {}",
            t.chi_square
        );
        assert!(t.rejects_hwe(0.001));
    }

    #[test]
    fn p_value_matches_erfc_closed_form_df1() {
        // For df = 1, P(chi2 > x) = erfc(sqrt(x/2)). Cross-check the
        // incomplete-gamma path against the independent erfc series.
        for &x in &[0.5_f64, 1.0, 2.0, 3.841_459, 6.635, 15.5] {
            let viagamma = chi_square_upper_tail(x, 1).unwrap();
            let viaerfc = erfc((x / 2.0).sqrt());
            assert!(
                (viagamma - viaerfc).abs() < 1e-6,
                "x={x}: gamma={viagamma}, erfc={viaerfc}"
            );
        }
    }

    #[test]
    fn critical_value_p_is_five_percent() {
        // The df=1 chi-square 0.05 critical value is 3.841459; its
        // upper-tail probability must be 0.05.
        let p = chi_square_upper_tail(3.841_459, 1).unwrap();
        assert!((p - 0.05).abs() < 1e-4, "p = {p}");
    }

    #[test]
    fn upper_tail_is_monotone_decreasing() {
        let p_small = chi_square_upper_tail(1.0, 1).unwrap();
        let p_large = chi_square_upper_tail(10.0, 1).unwrap();
        assert!(p_large < p_small);
        assert!((0.0..=1.0).contains(&p_small));
        assert!((0.0..=1.0).contains(&p_large));
    }

    #[test]
    fn zero_statistic_has_unit_p_value() {
        assert!((chi_square_upper_tail(0.0, 1).unwrap() - 1.0).abs() < 1e-12);
        assert!((chi_square_upper_tail(0.0, 3).unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn df_two_mean_check() {
        // A chi-square(2) variable has median ln(4) ~ 1.3862944, where
        // the upper-tail probability is exactly 0.5 (df=2 is exponential
        // with mean 2: P(X > x) = exp(-x/2), so x = 2 ln 2).
        let p = chi_square_upper_tail(2.0 * std::f64::consts::LN_2, 2).unwrap();
        assert!((p - 0.5).abs() < 1e-6, "p = {p}");
    }

    #[test]
    fn fixed_locus_is_rejected_as_no_df() {
        // All AA: p = 1, monomorphic -> undefined statistic.
        let err = hwe_chi_square_test(50.0, 0.0, 0.0).unwrap_err();
        assert_eq!(err.code(), "hwe.no_degrees_of_freedom");
    }

    #[test]
    fn upper_tail_rejects_bad_input() {
        assert!(chi_square_upper_tail(-1.0, 1).is_err());
        assert!(chi_square_upper_tail(f64::NAN, 1).is_err());
        assert!(chi_square_upper_tail(1.0, 0).is_err());
    }

    #[test]
    fn test_propagates_count_errors() {
        assert!(hwe_chi_square_test(-1.0, 5.0, 5.0).is_err());
        assert!(hwe_chi_square_test(0.0, 0.0, 0.0).is_err());
    }
}
