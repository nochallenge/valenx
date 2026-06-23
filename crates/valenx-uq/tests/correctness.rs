//! Correctness gate for `valenx-uq`.
//!
//! Each test pins a fixed PRNG seed so the assertions are deterministic and
//! reproducible across runs and machines. The five groups mirror the crate's
//! required-test contract:
//!
//! 1. Monte-Carlo recovers the mean/variance of a known additive model.
//! 2. Latin-hypercube is truly stratified, and beats MC on a smooth function.
//! 3. Sobol indices match the analytic Ishigami values.
//! 4. The polynomial surrogate fits a quadratic exactly and noisy data well.
//! 5. Validation / degenerate inputs are handled (no panic, no NaN).

use std::f64::consts::PI;
use valenx_uq::sampling::{latin_hypercube, monte_carlo, stratum_indices};
use valenx_uq::sensitivity::{morris, sobol_indices};
use valenx_uq::statistics;
use valenx_uq::surrogate::PolynomialSurrogate;
use valenx_uq::{Distribution, FnModel, Model, SplitMix64, UqError};

// --- 1. Monte-Carlo forward propagation -----------------------------------

#[test]
fn monte_carlo_recovers_mean_and_variance_of_sum() {
    // Model: y = x0 + x1.
    let model = FnModel::new(2, 1, |x| vec![x[0] + x[1]]);

    // x0 ~ N(1, 2), x1 ~ N(-3, 1).  For independent inputs:
    //   E[y]   = 1 + (-3)       = -2
    //   Var[y] = 2² + 1²        = 5
    let dists = [
        Distribution::normal(1.0, 2.0).unwrap(),
        Distribution::normal(-3.0, 1.0).unwrap(),
    ];

    let mut rng = SplitMix64::new(0x1234_5678);
    let inputs = monte_carlo(200_000, &dists, &mut rng);
    let outputs: Vec<f64> = inputs.iter().map(|x| model.evaluate(x)[0]).collect();

    let mean = statistics::mean(&outputs).unwrap();
    let var = statistics::variance(&outputs).unwrap();

    assert!(
        (mean - (-2.0)).abs() < 0.05,
        "MC mean {mean} should be ~ -2"
    );
    assert!((var - 5.0).abs() < 0.1, "MC variance {var} should be ~ 5");
}

// --- 2. Latin-hypercube stratification + variance -------------------------

#[test]
fn lhs_is_a_true_stratification() {
    // For n samples in d dims, each dimension's stratum indices must be a
    // permutation of 0..n (every band hit exactly once).
    let n = 64;
    let dists = [
        Distribution::uniform(-1.0, 1.0).unwrap(),
        Distribution::normal(0.0, 1.0).unwrap(),
        Distribution::triangular(0.0, 2.0, 5.0).unwrap(),
    ];

    let mut rng = SplitMix64::new(0xABCD_0001);
    let samples = latin_hypercube(n, &dists, &mut rng);
    assert_eq!(samples.len(), n);

    let table = stratum_indices(&samples, &dists);
    for (j, strata) in table.iter().enumerate() {
        let mut sorted = strata.clone();
        sorted.sort_unstable();
        let expected: Vec<usize> = (0..n).collect();
        assert_eq!(
            sorted, expected,
            "dimension {j}: stratum indices must be a permutation of 0..{n}"
        );
    }
}

#[test]
fn lhs_mean_is_at_least_as_accurate_as_mc() {
    // On a smooth function, LHS should estimate the mean at least as well as
    // plain MC at equal n (lower-variance space-filling design).
    //
    // f(x) = sin(x0) * exp(-0.5 * x1^2)  with x ~ U(-2,2)^2.
    let f = |x: &[f64]| (x[0]).sin() * (-0.5 * x[1] * x[1]).exp();
    let model = FnModel::new(2, 1, move |x| vec![f(x)]);
    let dists = [
        Distribution::uniform(-2.0, 2.0).unwrap(),
        Distribution::uniform(-2.0, 2.0).unwrap(),
    ];

    // Ground truth by a fine deterministic grid (no PRNG): average over the
    // box. f is odd in x0 over a symmetric box, so the true mean is 0.
    let truth = grid_mean(&f, -2.0, 2.0, -2.0, 2.0, 400);

    let n = 256;
    let mut mc_rng = SplitMix64::new(7);
    let mut lhs_rng = SplitMix64::new(7);

    let mc = monte_carlo(n, &dists, &mut mc_rng);
    let lhs = latin_hypercube(n, &dists, &mut lhs_rng);

    let mc_mean =
        statistics::mean(&mc.iter().map(|x| model.evaluate(x)[0]).collect::<Vec<_>>()).unwrap();
    let lhs_mean =
        statistics::mean(&lhs.iter().map(|x| model.evaluate(x)[0]).collect::<Vec<_>>()).unwrap();

    let mc_err = (mc_mean - truth).abs();
    let lhs_err = (lhs_mean - truth).abs();

    assert!(
        lhs_err <= mc_err + 1e-9,
        "LHS error {lhs_err} should be <= MC error {mc_err} (truth={truth})"
    );
}

// --- 3. Sobol on the Ishigami function ------------------------------------

#[test]
fn sobol_matches_ishigami_analytic_indices() {
    // Ishigami: f = sin(x0) + a sin²(x1) + b x2⁴ sin(x0), x_i ~ U(-π, π).
    const A: f64 = 7.0;
    const B: f64 = 0.1;
    let model = FnModel::new(3, 1, |x| {
        let f = x[0].sin() + A * x[1].sin().powi(2) + B * x[2].powi(4) * x[0].sin();
        vec![f]
    });
    let dists = [
        Distribution::uniform(-PI, PI).unwrap(),
        Distribution::uniform(-PI, PI).unwrap(),
        Distribution::uniform(-PI, PI).unwrap(),
    ];

    // Analytic Ishigami variances (standard reference values):
    //   D  = a²/8 + b π⁴/5 + b² π⁸/18 + 1/2
    //   D1 = b π⁴/5 + b² π⁸/50 + 1/2
    //   D2 = a²/8
    //   D3 = 0
    //   DT1 = D1 + b² π⁸ (1/18 - 1/50)   (interaction of x0,x2)
    //   DT3 = b² π⁸ (1/18 - 1/50)
    let pi4 = PI.powi(4);
    let pi8 = PI.powi(8);
    let d_tot = A * A / 8.0 + B * pi4 / 5.0 + B * B * pi8 / 18.0 + 0.5;
    let d1 = B * pi4 / 5.0 + B * B * pi8 / 50.0 + 0.5;
    let d2 = A * A / 8.0;
    let interaction = B * B * pi8 * (1.0 / 18.0 - 1.0 / 50.0);

    let s1_true = d1 / d_tot;
    let s2_true = d2 / d_tot;
    let s3_true = 0.0;
    let st1_true = (d1 + interaction) / d_tot;
    let st2_true = d2 / d_tot; // x1 has no interactions ⇒ total == first order
    let st3_true = interaction / d_tot;

    let mut rng = SplitMix64::new(0x5A1_7E11);
    let res = sobol_indices(&model, &dists, 200_000, 0, &mut rng).expect("sobol ok");

    // Saltelli estimator at finite N: ~0.03 absolute tolerance is comfortable.
    let tol = 0.03;
    assert!(
        (res.first_order[0] - s1_true).abs() < tol,
        "S1 {} vs {s1_true}",
        res.first_order[0]
    );
    assert!(
        (res.first_order[1] - s2_true).abs() < tol,
        "S2 {} vs {s2_true}",
        res.first_order[1]
    );
    assert!(
        (res.first_order[2] - s3_true).abs() < tol,
        "S3 {} vs {s3_true}",
        res.first_order[2]
    );
    assert!(
        (res.total[0] - st1_true).abs() < tol,
        "ST1 {} vs {st1_true}",
        res.total[0]
    );
    assert!(
        (res.total[1] - st2_true).abs() < tol,
        "ST2 {} vs {st2_true}",
        res.total[1]
    );
    assert!(
        (res.total[2] - st3_true).abs() < tol,
        "ST3 {} vs {st3_true}",
        res.total[2]
    );
}

#[test]
fn morris_ranks_inputs_by_elementary_effect() {
    // Morris is a screening method based on *elementary effects* (local
    // gradients), which is distinct from variance share. On a LINEAR model the
    // elementary effect of each input is exactly its coefficient, with zero
    // spread — a clean, analytically-known ground truth for the screening.
    //
    //   y = 10*x0 + 1*x1 + 0*x2   ⇒   mu_star ≈ [10, 1, 0], sigma ≈ [0, 0, 0].
    let model = FnModel::new(3, 1, |x| vec![10.0 * x[0] + 1.0 * x[1] + 0.0 * x[2]]);
    let dists = [
        Distribution::uniform(-1.0, 1.0).unwrap(),
        Distribution::uniform(-1.0, 1.0).unwrap(),
        Distribution::uniform(-1.0, 1.0).unwrap(),
    ];

    let mut rng = SplitMix64::new(0x_0033_5510);
    let res = morris(&model, &dists, 400, 8, 0, &mut rng).expect("morris ok");

    // Strict, correct ranking x0 > x1 > x2, with values near the coefficients.
    assert!(
        res.mu_star[0] > res.mu_star[1] && res.mu_star[1] > res.mu_star[2],
        "mu_star ranking {:?} should be strictly decreasing",
        res.mu_star
    );
    assert!(
        (res.mu_star[0] - 10.0).abs() < 1e-6,
        "x0 mu* {}",
        res.mu_star[0]
    );
    assert!(
        (res.mu_star[1] - 1.0).abs() < 1e-6,
        "x1 mu* {}",
        res.mu_star[1]
    );
    assert!(res.mu_star[2].abs() < 1e-9, "x2 mu* {}", res.mu_star[2]);
    // A linear model has constant elementary effects ⇒ ~zero spread.
    for (i, &s) in res.sigma.iter().enumerate() {
        assert!(s < 1e-6, "linear sigma[{i}] {s} should be ~0");
    }

    // And on a non-linear model, sigma must become non-trivial: the
    // elementary effects vary across the input space.
    let nonlinear = FnModel::new(2, 1, |x| vec![x[0] * x[0] + x[0] * x[1]]);
    let d2 = [
        Distribution::uniform(-2.0, 2.0).unwrap(),
        Distribution::uniform(-2.0, 2.0).unwrap(),
    ];
    let mut rng2 = SplitMix64::new(0x_00C0_FFEE);
    let nl = morris(&nonlinear, &d2, 400, 8, 0, &mut rng2).expect("morris ok");
    assert!(
        nl.sigma[0] > 0.1,
        "non-linear sigma[0] {} should be clearly > 0",
        nl.sigma[0]
    );
}

// --- 4. Polynomial surrogate ----------------------------------------------

#[test]
fn surrogate_fits_quadratic_exactly() {
    // True model: y = 1 + 2*x0 - 3*x1 + 0.5*x0² + x0*x1.
    let truth = |x: &[f64]| 1.0 + 2.0 * x[0] - 3.0 * x[1] + 0.5 * x[0] * x[0] + x[0] * x[1];

    // A deterministic grid of training points (no PRNG needed).
    let mut samples = Vec::new();
    let mut values = Vec::new();
    for i in 0..7 {
        for j in 0..7 {
            let x0 = -3.0 + i as f64;
            let x1 = -3.0 + j as f64;
            let s = vec![x0, x1];
            values.push(truth(&s));
            samples.push(s);
        }
    }

    let surrogate = PolynomialSurrogate::fit(&samples, &values, 2).expect("fit ok");
    assert!(
        surrogate.r_squared() > 1.0 - 1e-9,
        "R² {} should be ~1 on exactly-polynomial data",
        surrogate.r_squared()
    );

    // Predictions should match the truth at held-out points.
    for &(x0, x1) in &[(0.25_f64, -1.5_f64), (2.7, 2.1), (-2.2, 0.4)] {
        let pred = surrogate.predict(&[x0, x1]).unwrap();
        let exact = truth(&[x0, x1]);
        assert!(
            (pred - exact).abs() < 1e-6,
            "predict {pred} vs exact {exact} at ({x0},{x1})"
        );
    }
}

#[test]
fn surrogate_handles_noisy_data_with_high_but_imperfect_fit() {
    // Linear truth y = 3 + 2*x with small additive deterministic-PRNG noise.
    let mut rng = SplitMix64::new(0x_0001_5E10);
    let mut samples = Vec::new();
    let mut values = Vec::new();
    for k in 0..200 {
        let x = -5.0 + 0.05 * k as f64;
        // Noise ~ N(0, 0.3): small relative to the signal range (~20).
        let noise = 0.3 * rng.next_standard_normal();
        samples.push(vec![x]);
        values.push(3.0 + 2.0 * x + noise);
    }

    let surrogate = PolynomialSurrogate::fit(&samples, &values, 1).expect("fit ok");
    let r2 = surrogate.r_squared();
    assert!(
        r2 < 1.0 && r2 > 0.95,
        "noisy R² {r2} should be high but < 1"
    );

    // Slope/intercept recovered close to truth (coeffs: [const, x]).
    let c = surrogate.coefficients();
    assert!((c[0] - 3.0).abs() < 0.2, "intercept {} ~ 3", c[0]);
    assert!((c[1] - 2.0).abs() < 0.1, "slope {} ~ 2", c[1]);
}

// --- 5. Validation / degenerate handling ----------------------------------

#[test]
fn invalid_distributions_error() {
    assert!(matches!(
        Distribution::normal(0.0, 0.0),
        Err(UqError::InvalidDistribution(_))
    ));
    assert!(matches!(
        Distribution::normal(0.0, -1.0),
        Err(UqError::InvalidDistribution(_))
    ));
    assert!(matches!(
        Distribution::uniform(1.0, 1.0),
        Err(UqError::InvalidDistribution(_))
    ));
    assert!(matches!(
        Distribution::uniform(2.0, 1.0),
        Err(UqError::InvalidDistribution(_))
    ));
    // mode outside [lo, hi]
    assert!(matches!(
        Distribution::triangular(0.0, 5.0, 4.0),
        Err(UqError::InvalidDistribution(_))
    ));
    assert!(matches!(
        Distribution::triangular(0.0, -1.0, 4.0),
        Err(UqError::InvalidDistribution(_))
    ));
    // Valid ones succeed.
    assert!(Distribution::triangular(0.0, 2.0, 4.0).is_ok());
}

#[test]
fn empty_sample_statistics_do_not_panic_or_nan() {
    let empty: [f64; 0] = [];
    assert_eq!(statistics::mean(&empty), None);
    assert_eq!(statistics::variance(&empty), None);
    assert_eq!(statistics::std(&empty), None);
    assert!(matches!(
        statistics::percentile(&empty, 50.0),
        Err(UqError::EmptyInput(_))
    ));
    assert!(matches!(
        statistics::confidence_interval(&empty, 0.95),
        Err(UqError::EmptyInput(_))
    ));
    // Single-element variance is also undefined (needs n >= 2) — not NaN.
    assert_eq!(statistics::variance(&[1.0]), None);
}

#[test]
fn percentile_validates_p_range() {
    let data = [1.0, 2.0, 3.0, 4.0, 5.0];
    assert!(matches!(
        statistics::percentile(&data, -1.0),
        Err(UqError::OutOfRange(_))
    ));
    assert!(matches!(
        statistics::percentile(&data, 101.0),
        Err(UqError::OutOfRange(_))
    ));
    assert!(matches!(
        statistics::percentile(&data, f64::NAN),
        Err(UqError::OutOfRange(_))
    ));
    // Endpoints are valid and give min / max.
    assert_eq!(statistics::percentile(&data, 0.0).unwrap(), 1.0);
    assert_eq!(statistics::percentile(&data, 100.0).unwrap(), 5.0);
    // 50th percentile of 1..5 is 3.
    assert_eq!(statistics::percentile(&data, 50.0).unwrap(), 3.0);
    // Confidence level must be in (0, 1).
    assert!(matches!(
        statistics::confidence_interval(&data, 0.0),
        Err(UqError::OutOfRange(_))
    ));
    assert!(matches!(
        statistics::confidence_interval(&data, 1.0),
        Err(UqError::OutOfRange(_))
    ));
}

#[test]
fn surrogate_rejects_bad_inputs() {
    // Degree above the supported maximum.
    assert!(matches!(
        PolynomialSurrogate::fit(&[vec![0.0]], &[0.0], 3),
        Err(UqError::OutOfRange(_))
    ));
    // Empty samples.
    assert!(matches!(
        PolynomialSurrogate::fit(&[], &[], 1),
        Err(UqError::EmptyInput(_))
    ));
    // Mismatched sample/value counts.
    assert!(matches!(
        PolynomialSurrogate::fit(&[vec![0.0], vec![1.0]], &[0.0], 1),
        Err(UqError::DimensionMismatch(_))
    ));
    // Under-determined: a degree-2 fit in 2 inputs needs 6 terms, give 3.
    let few: Vec<Vec<f64>> = (0..3).map(|k| vec![k as f64, 0.0]).collect();
    assert!(matches!(
        PolynomialSurrogate::fit(&few, &[0.0, 1.0, 2.0], 2),
        Err(UqError::DimensionMismatch(_))
    ));
    // Predict with wrong input arity.
    let good: Vec<Vec<f64>> = (0..10).map(|k| vec![k as f64]).collect();
    let vals: Vec<f64> = (0..10).map(|k| k as f64).collect();
    let s = PolynomialSurrogate::fit(&good, &vals, 1).unwrap();
    assert!(matches!(
        s.predict(&[1.0, 2.0]),
        Err(UqError::DimensionMismatch(_))
    ));
}

#[test]
fn sobol_guards_ill_posed_inputs() {
    // Constant model ⇒ zero variance ⇒ None (no apportionable sensitivity).
    let constant = FnModel::new(2, 1, |_x| vec![42.0]);
    let dists = [
        Distribution::uniform(0.0, 1.0).unwrap(),
        Distribution::uniform(0.0, 1.0).unwrap(),
    ];
    let mut rng = SplitMix64::new(1);
    assert!(sobol_indices(&constant, &dists, 1000, 0, &mut rng).is_none());

    // n_base < 2, empty dists, and out-of-range output index all yield None.
    let m = FnModel::new(1, 1, |x| vec![x[0]]);
    let d1 = [Distribution::uniform(0.0, 1.0).unwrap()];
    assert!(sobol_indices(&m, &d1, 1, 0, &mut rng).is_none());
    assert!(sobol_indices(&m, &[], 100, 0, &mut rng).is_none());
    assert!(sobol_indices(&m, &d1, 100, 5, &mut rng).is_none());
}

// --- test-local helpers ----------------------------------------------------

/// Deterministic ground-truth mean of `f` over the box `[ax,bx]×[ay,by]`,
/// by averaging an `n×n` regular grid. No PRNG involved.
fn grid_mean(f: &impl Fn(&[f64]) -> f64, ax: f64, bx: f64, ay: f64, by: f64, n: usize) -> f64 {
    let mut acc = 0.0;
    for i in 0..n {
        for j in 0..n {
            let x = ax + (bx - ax) * (i as f64 + 0.5) / n as f64;
            let y = ay + (by - ay) * (j as f64 + 0.5) / n as f64;
            acc += f(&[x, y]);
        }
    }
    acc / (n * n) as f64
}
