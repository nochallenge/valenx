//! Ground-truth analytic tests for `valenx-statistics`.
//!
//! Each test pins an estimator to a value derived by hand or from a known
//! closed form (textbook reference data, the integers `1..=n`, the standard
//! normal's documented properties). All float comparisons use an absolute
//! tolerance, never `==`.

use valenx_statistics::descriptive::{
    mean, median, population_std, population_variance, quantile, quartiles, sample_std,
    sample_variance, standardize, z_score,
};
use valenx_statistics::error::StatsError;
use valenx_statistics::inference::{normal_ci, probit, t_statistic};
use valenx_statistics::normal::{cdf, erf, pdf, sf};

/// Absolute tolerance for "exact" closed-form arithmetic.
const EXACT: f64 = 1e-12;
/// Absolute tolerance for the A&S `erf` / Acklam probit approximations.
const APPROX: f64 = 1e-6;

// ---------------------------------------------------------------------------
// Mean
// ---------------------------------------------------------------------------

#[test]
fn mean_of_one_to_n_is_n_plus_1_over_2() {
    // x̄(1..=n) = (n+1)/2 for several n.
    for n in [1usize, 2, 5, 10, 100, 1000] {
        let data: Vec<f64> = (1..=n).map(|i| i as f64).collect();
        let expected = (n as f64 + 1.0) / 2.0;
        let got = mean(&data).unwrap();
        assert!(
            (got - expected).abs() < EXACT,
            "n={n}: mean {got} != {expected}"
        );
    }
}

#[test]
fn mean_of_constant_is_the_constant() {
    let data = [7.5, 7.5, 7.5, 7.5];
    assert!((mean(&data).unwrap() - 7.5).abs() < EXACT);
}

#[test]
fn mean_rejects_empty_and_nonfinite() {
    assert!(matches!(
        mean(&[]),
        Err(StatsError::EmptySample { estimator: "mean" })
    ));
    assert!(matches!(
        mean(&[1.0, f64::NAN]),
        Err(StatsError::NonFinite { name: "sample" })
    ));
    assert!(matches!(
        mean(&[1.0, f64::INFINITY]),
        Err(StatsError::NonFinite { .. })
    ));
}

// ---------------------------------------------------------------------------
// Variance / standard deviation
// ---------------------------------------------------------------------------

/// Textbook worked example: {2,4,4,4,5,5,7,9} has mean 5,
/// population variance 4 (σ=2), sample variance 32/7.
const TEXTBOOK: [f64; 8] = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];

#[test]
fn population_variance_textbook() {
    assert!((mean(&TEXTBOOK).unwrap() - 5.0).abs() < EXACT);
    assert!((population_variance(&TEXTBOOK).unwrap() - 4.0).abs() < EXACT);
    assert!((population_std(&TEXTBOOK).unwrap() - 2.0).abs() < EXACT);
}

#[test]
fn sample_variance_uses_n_minus_1() {
    // Sum of squared deviations is 32; sample divides by n-1=7.
    let sample_v = sample_variance(&TEXTBOOK).unwrap();
    let pop_v = population_variance(&TEXTBOOK).unwrap();
    assert!((sample_v - 32.0 / 7.0).abs() < EXACT);
    // The Bessel relation: s² = σ² · n/(n-1).
    let n = TEXTBOOK.len() as f64;
    assert!((sample_v - pop_v * n / (n - 1.0)).abs() < EXACT);
    // Sample variance strictly exceeds population variance for n>1.
    assert!(sample_v > pop_v);
}

#[test]
fn std_is_sqrt_of_variance() {
    assert!(
        (sample_std(&TEXTBOOK).unwrap() - sample_variance(&TEXTBOOK).unwrap().sqrt()).abs() < EXACT
    );
    assert!(
        (population_std(&TEXTBOOK).unwrap() - population_variance(&TEXTBOOK).unwrap().sqrt()).abs()
            < EXACT
    );
}

#[test]
fn variance_of_constant_is_zero() {
    let data = [3.0, 3.0, 3.0, 3.0, 3.0];
    assert!(population_variance(&data).unwrap().abs() < EXACT);
    assert!(sample_variance(&data).unwrap().abs() < EXACT);
}

#[test]
fn sample_variance_needs_two_points() {
    match sample_variance(&[42.0]) {
        Err(StatsError::TooFewObservations {
            estimator: "sample_variance",
            needed: 2,
            got: 1,
        }) => {}
        other => panic!("expected TooFewObservations, got {other:?}"),
    }
    // Population variance, by contrast, is fine with a single point (zero).
    assert!(population_variance(&[42.0]).unwrap().abs() < EXACT);
}

// ---------------------------------------------------------------------------
// Median / quartiles / quantile
// ---------------------------------------------------------------------------

#[test]
fn median_odd_and_even() {
    // Odd count: middle order statistic (also checks it sorts).
    assert!((median(&[3.0, 1.0, 2.0]).unwrap() - 2.0).abs() < EXACT);
    // Even count: mean of the two central order statistics.
    assert!((median(&[1.0, 2.0, 3.0, 4.0]).unwrap() - 2.5).abs() < EXACT);
    // Single point: the point itself.
    assert!((median(&[9.0]).unwrap() - 9.0).abs() < EXACT);
}

#[test]
fn quartiles_type7() {
    // 1..=5: type-7 quartiles are 2, 3, 4 (NumPy/R default).
    let q = quartiles(&[1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
    assert!((q.q1 - 2.0).abs() < EXACT);
    assert!((q.q2 - 3.0).abs() < EXACT);
    assert!((q.q3 - 4.0).abs() < EXACT);
    assert!((q.iqr() - 2.0).abs() < EXACT);
    // Q2 equals the median.
    assert!((q.q2 - median(&[1.0, 2.0, 3.0, 4.0, 5.0]).unwrap()).abs() < EXACT);
}

#[test]
fn quantile_endpoints_and_interpolation() {
    let data = [1.0, 2.0, 3.0, 4.0];
    // q=0 -> min, q=1 -> max.
    assert!((quantile(&data, 0.0).unwrap() - 1.0).abs() < EXACT);
    assert!((quantile(&data, 1.0).unwrap() - 4.0).abs() < EXACT);
    // h = 0.5*(4-1) = 1.5 -> halfway between x[1]=2 and x[2]=3 -> 2.5.
    assert!((quantile(&data, 0.5).unwrap() - 2.5).abs() < EXACT);
    // h = 0.25*3 = 0.75 -> x[0] + 0.75*(x[1]-x[0]) = 1.75.
    assert!((quantile(&data, 0.25).unwrap() - 1.75).abs() < EXACT);
}

#[test]
fn quantile_is_monotone_nondecreasing_in_q() {
    let data = [5.0, 1.0, 3.0, 9.0, 7.0, 2.0];
    let mut prev = f64::NEG_INFINITY;
    for k in 0..=10 {
        let q = k as f64 / 10.0;
        let v = quantile(&data, q).unwrap();
        assert!(
            v >= prev - EXACT,
            "quantile not monotone at q={q}: {v} < {prev}"
        );
        prev = v;
    }
}

#[test]
fn quantile_rejects_out_of_range_q() {
    match quantile(&[1.0, 2.0], 1.5) {
        Err(StatsError::OutOfRange {
            name: "q", value, ..
        }) => {
            assert!((value - 1.5).abs() < EXACT);
        }
        other => panic!("expected OutOfRange, got {other:?}"),
    }
    assert!(matches!(
        quantile(&[1.0, 2.0], -0.1),
        Err(StatsError::OutOfRange { name: "q", .. })
    ));
    assert!(matches!(
        quantile(&[], 0.5),
        Err(StatsError::EmptySample { .. })
    ));
}

// ---------------------------------------------------------------------------
// Standard normal: pdf / cdf
// ---------------------------------------------------------------------------

#[test]
fn cdf_at_zero_is_one_half() {
    assert!((cdf(0.0) - 0.5).abs() < EXACT);
}

#[test]
fn pdf_peak_and_symmetry() {
    // Peak height 1/sqrt(2*pi) at the origin.
    let peak = 1.0 / (2.0 * std::f64::consts::PI).sqrt();
    assert!((pdf(0.0) - peak).abs() < EXACT);
    // Even function: pdf(-z) == pdf(z).
    for z in [0.3_f64, 1.0, 2.5, 4.0] {
        assert!((pdf(-z) - pdf(z)).abs() < EXACT, "pdf not symmetric at {z}");
    }
}

#[test]
fn cdf_reflection_and_monotone() {
    // Reflection identity Phi(-z) = 1 - Phi(z).
    for z in [0.2_f64, 0.7, 1.5, 2.0, 3.0] {
        assert!(
            (cdf(-z) - (1.0 - cdf(z))).abs() < APPROX,
            "reflection fails at {z}"
        );
    }
    // Monotone increasing.
    let mut prev = -1.0;
    for k in -50..=50 {
        let z = k as f64 / 10.0;
        let c = cdf(z);
        assert!(c >= prev - APPROX, "cdf not monotone at z={z}");
        prev = c;
    }
}

#[test]
fn cdf_known_quantiles_and_empirical_rule() {
    // Standard normal landmark values (7-figure references).
    assert!((cdf(1.0) - 0.841_344_7).abs() < APPROX);
    assert!((cdf(1.959_963_985) - 0.975).abs() < APPROX);
    assert!((cdf(2.575_829_304) - 0.995).abs() < APPROX);
    // 68-95-99.7 empirical rule.
    assert!((cdf(1.0) - cdf(-1.0) - 0.682_689_5).abs() < APPROX);
    assert!((cdf(2.0) - cdf(-2.0) - 0.954_499_7).abs() < APPROX);
    assert!((cdf(3.0) - cdf(-3.0) - 0.997_300_2).abs() < APPROX);
}

#[test]
fn sf_complements_cdf() {
    for z in [-1.5_f64, 0.0, 0.8, 2.3] {
        assert!((sf(z) + cdf(z) - 1.0).abs() < APPROX, "sf+cdf != 1 at {z}");
    }
}

#[test]
fn erf_odd_and_pinned() {
    assert!(erf(0.0).abs() < EXACT);
    assert!((erf(-1.2) + erf(1.2)).abs() < EXACT);
    // erf(1) ≈ 0.842700793.
    assert!((erf(1.0) - 0.842_700_793).abs() < APPROX);
}

// ---------------------------------------------------------------------------
// z-score
// ---------------------------------------------------------------------------

#[test]
fn z_score_formula() {
    // z = (x - mean)/std.
    assert!((z_score(130.0, 100.0, 15.0).unwrap() - 2.0).abs() < EXACT);
    assert!((z_score(85.0, 100.0, 15.0).unwrap() - (-1.0)).abs() < EXACT);
    // At the mean the score is zero.
    assert!(z_score(50.0, 50.0, 4.0).unwrap().abs() < EXACT);
}

#[test]
fn z_score_rejects_nonpositive_std() {
    assert!(matches!(
        z_score(1.0, 0.0, 0.0),
        Err(StatsError::NonPositiveScale { name: "std", .. })
    ));
    assert!(matches!(
        z_score(1.0, 0.0, -2.0),
        Err(StatsError::NonPositiveScale { .. })
    ));
}

#[test]
fn standardize_has_zero_mean_unit_sample_std() {
    let data = [10.0, 12.0, 23.0, 23.0, 16.0, 23.0, 21.0, 16.0];
    let z = standardize(&data).unwrap();
    // Standardised data: sample mean 0, sample std 1.
    assert!(mean(&z).unwrap().abs() < EXACT);
    assert!((sample_std(&z).unwrap() - 1.0).abs() < EXACT);
}

// ---------------------------------------------------------------------------
// t-statistic
// ---------------------------------------------------------------------------

#[test]
fn t_statistic_zero_when_mean_matches() {
    let data = [4.0, 5.0, 6.0];
    assert!(t_statistic(&data, 5.0).unwrap().abs() < EXACT);
}

#[test]
fn t_statistic_known_value() {
    // data {1,2,3,4,5}: mean=3, s=sqrt(2.5), n=5, SE=sqrt(2.5)/sqrt(5)=sqrt(0.5).
    // Test mu0=0: t = 3 / sqrt(0.5) = 3*sqrt(2) ≈ 4.242640687.
    let data = [1.0, 2.0, 3.0, 4.0, 5.0];
    let t = t_statistic(&data, 0.0).unwrap();
    assert!((t - 3.0 * 2.0_f64.sqrt()).abs() < EXACT, "t={t}");
    // Sign flips below the hypothesised mean.
    assert!(t_statistic(&data, 10.0).unwrap() < 0.0);
}

#[test]
fn t_statistic_rejects_small_or_constant() {
    assert!(matches!(
        t_statistic(&[1.0], 0.0),
        Err(StatsError::TooFewObservations { .. })
    ));
    assert!(matches!(
        t_statistic(&[5.0, 5.0, 5.0], 0.0),
        Err(StatsError::NonPositiveScale { .. })
    ));
}

// ---------------------------------------------------------------------------
// probit (inverse normal CDF)
// ---------------------------------------------------------------------------

#[test]
fn probit_landmarks() {
    assert!(probit(0.5).unwrap().abs() < 1e-9);
    assert!((probit(0.975).unwrap() - 1.959_963_985).abs() < APPROX);
    assert!((probit(0.995).unwrap() - 2.575_829_304).abs() < APPROX);
    // Odd about the median.
    assert!((probit(0.1).unwrap() + probit(0.9).unwrap()).abs() < APPROX);
}

#[test]
fn probit_inverts_cdf() {
    // cdf(probit(p)) ≈ p across the interior.
    for k in 1..20 {
        let p = k as f64 / 20.0;
        let z = probit(p).unwrap();
        assert!((cdf(z) - p).abs() < APPROX, "cdf(probit({p}))={}", cdf(z));
    }
}

#[test]
fn probit_rejects_endpoints() {
    assert!(matches!(
        probit(0.0),
        Err(StatsError::OutOfRange { name: "p", .. })
    ));
    assert!(matches!(
        probit(1.0),
        Err(StatsError::OutOfRange { name: "p", .. })
    ));
}

// ---------------------------------------------------------------------------
// Confidence interval
// ---------------------------------------------------------------------------

#[test]
fn ci_centered_on_mean_with_known_margin() {
    // n=4, sigma=2, 95%: margin = z_{0.975} * 2/sqrt(4) = 1.959964 * 1.
    let data = [10.0, 10.0, 10.0, 10.0];
    let ci = normal_ci(&data, 2.0, 0.95).unwrap();
    assert!((ci.point - 10.0).abs() < EXACT);
    let expected_margin = 1.959_963_985 * 2.0 / 4.0_f64.sqrt();
    assert!((ci.margin - expected_margin).abs() < APPROX);
    assert!((ci.lower - (10.0 - expected_margin)).abs() < APPROX);
    assert!((ci.upper - (10.0 + expected_margin)).abs() < APPROX);
    assert!((ci.width() - 2.0 * expected_margin).abs() < APPROX);
    assert!(ci.contains(10.0));
}

#[test]
fn ci_grows_with_sigma() {
    let data = [0.0, 0.0, 0.0, 0.0, 0.0];
    let narrow = normal_ci(&data, 1.0, 0.95).unwrap();
    let wide = normal_ci(&data, 5.0, 0.95).unwrap();
    assert!(
        wide.margin > narrow.margin,
        "larger sigma must widen the interval: {} vs {}",
        wide.margin,
        narrow.margin
    );
    // Linear in sigma: 5x the spread -> 5x the margin.
    assert!((wide.margin - 5.0 * narrow.margin).abs() < APPROX);
}

#[test]
fn ci_shrinks_with_n() {
    // Same sigma and level, more data -> tighter interval (~1/sqrt(n)).
    let small: Vec<f64> = vec![0.0; 4];
    let large: Vec<f64> = vec![0.0; 100];
    let m_small = normal_ci(&small, 3.0, 0.95).unwrap().margin;
    let m_large = normal_ci(&large, 3.0, 0.95).unwrap().margin;
    assert!(m_large < m_small, "more data must shrink the interval");
    // n goes 4 -> 100 (25x), so margin scales by sqrt(4/100) = 1/5.
    assert!((m_large - m_small / 5.0).abs() < APPROX);
}

#[test]
fn ci_widens_with_confidence() {
    let data = [1.0, 2.0, 3.0, 4.0];
    let c90 = normal_ci(&data, 2.0, 0.90).unwrap().margin;
    let c95 = normal_ci(&data, 2.0, 0.95).unwrap().margin;
    let c99 = normal_ci(&data, 2.0, 0.99).unwrap().margin;
    assert!(
        c90 < c95 && c95 < c99,
        "higher confidence must widen: {c90} {c95} {c99}"
    );
}

#[test]
fn ci_rejects_bad_params() {
    let data = [1.0, 2.0, 3.0];
    assert!(matches!(
        normal_ci(&data, 0.0, 0.95),
        Err(StatsError::NonPositiveScale { name: "sigma", .. })
    ));
    assert!(matches!(
        normal_ci(&data, 1.0, 1.0),
        Err(StatsError::OutOfRange {
            name: "confidence",
            ..
        })
    ));
    assert!(matches!(
        normal_ci(&data, 1.0, 0.0),
        Err(StatsError::OutOfRange {
            name: "confidence",
            ..
        })
    ));
    assert!(matches!(
        normal_ci(&[], 1.0, 0.95),
        Err(StatsError::EmptySample { .. })
    ));
}
