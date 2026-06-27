//! Input-space sampling designs.
//!
//! Both routines turn a list of per-dimension [`Distribution`]s into `n`
//! sample points (each an input vector). They differ in *how* the points fill
//! the space:
//!
//! * [`monte_carlo`] draws every coordinate independently from its
//!   distribution — simple, unbiased, and with `O(1/√n)` error.
//! * [`latin_hypercube`] **stratifies** each dimension into `n` equal-
//!   probability bands and places exactly one sample in each band, then
//!   permutes the bands independently per dimension. This space-filling
//!   guarantee removes the clustering/gaps of plain Monte-Carlo and typically
//!   gives a lower-variance mean estimate for smooth functions at the same
//!   `n` (it improves the constant, not the asymptotic rate).

use crate::distribution::Distribution;
use crate::rng::SplitMix64;

/// Draw `n` independent Monte-Carlo samples.
///
/// Returns a `Vec` of `n` input vectors, each of length `dists.len()`, where
/// coordinate `j` of every sample is an independent draw from `dists[j]`.
///
/// An empty `dists` yields `n` empty vectors; `n == 0` yields an empty result.
#[must_use]
pub fn monte_carlo(n: usize, dists: &[Distribution], rng: &mut SplitMix64) -> Vec<Vec<f64>> {
    (0..n)
        .map(|_| dists.iter().map(|d| d.sample(rng)).collect())
        .collect()
}

/// Draw `n` Latin-hypercube samples.
///
/// For each of the `d = dists.len()` dimensions the unit interval is split into
/// `n` equal-probability strata `[k/n, (k+1)/n)`. One uniform draw is taken
/// inside each stratum and mapped through the distribution's inverse CDF; the
/// `n` resulting values are then assigned to sample rows according to an
/// independent random permutation of `0..n`. Each dimension therefore visits
/// every stratum exactly once (true stratification), while the per-dimension
/// permutation decorrelates the coordinates.
///
/// Returns a `Vec` of `n` input vectors of length `d`. `n == 0` yields an empty
/// result; an empty `dists` yields `n` empty vectors.
#[must_use]
pub fn latin_hypercube(n: usize, dists: &[Distribution], rng: &mut SplitMix64) -> Vec<Vec<f64>> {
    let d = dists.len();
    let mut samples = vec![vec![0.0_f64; d]; n];
    if n == 0 || d == 0 {
        return samples;
    }

    let inv_n = 1.0 / n as f64;
    for (j, dist) in dists.iter().enumerate() {
        // One jittered draw per stratum: u_k ∈ [k/n, (k+1)/n).
        let mut column: Vec<f64> = (0..n)
            .map(|k| {
                let u = (k as f64 + rng.next_f64()) * inv_n;
                dist.quantile(u)
            })
            .collect();
        // Permute this dimension independently (Fisher–Yates) so the
        // stratum order is randomised across rows.
        fisher_yates(&mut column, rng);
        for (i, value) in column.into_iter().enumerate() {
            samples[i][j] = value;
        }
    }
    samples
}

/// In-place Fisher–Yates shuffle driven by the crate PRNG.
fn fisher_yates(items: &mut [f64], rng: &mut SplitMix64) {
    let len = items.len();
    if len < 2 {
        return;
    }
    for i in (1..len).rev() {
        // Uniform index in 0..=i.
        let j = (rng.next_f64() * (i as f64 + 1.0)) as usize;
        // Guard the (vanishingly unlikely) next_f64() == 1.0 edge.
        let j = j.min(i);
        items.swap(i, j);
    }
}

/// Compute, for each dimension, the **stratum index** every sample falls into
/// under an `n`-way equal-probability split of that dimension's distribution.
///
/// Returns a `d × n` table: row `j` holds, for sample rows `0..n`, the stratum
/// (`0..n`) of that sample's `j`-th coordinate. For a correctly stratified
/// Latin-hypercube design each row is a **permutation of `0..n`** — exactly the
/// property the LHS test checks. Exposed for tests and for auditing a design.
///
/// Samples whose dimension is shorter than `d` contribute no stratum for the
/// missing coordinates (those entries are left at `0`); well-formed designs do
/// not hit this.
#[must_use]
pub fn stratum_indices(samples: &[Vec<f64>], dists: &[Distribution]) -> Vec<Vec<usize>> {
    let n = samples.len();
    let d = dists.len();
    let mut table = vec![vec![0_usize; n]; d];
    if n == 0 {
        return table;
    }
    for (j, dist) in dists.iter().enumerate() {
        for (i, sample) in samples.iter().enumerate() {
            let Some(&value) = sample.get(j) else {
                continue;
            };
            // Probability-integral transform → which of the n equal bands.
            let u = cdf(dist, value).clamp(0.0, 1.0);
            // Map u ∈ [0,1] to a stratum in 0..n, clamping the u == 1 edge.
            let s = ((u * n as f64) as usize).min(n - 1);
            table[j][i] = s;
        }
    }
    table
}

/// CDF of a [`Distribution`] at `x` — the inverse of
/// [`Distribution::quantile`]. Local to sampling because it is only needed to
/// recover stratum membership for auditing/tests.
fn cdf(dist: &Distribution, x: f64) -> f64 {
    match *dist {
        Distribution::Uniform { lo, hi } => ((x - lo) / (hi - lo)).clamp(0.0, 1.0),
        Distribution::Normal { mean, std } => standard_normal_cdf((x - mean) / std),
        Distribution::Triangular { lo, mode, hi } => {
            if x <= lo {
                0.0
            } else if x >= hi {
                1.0
            } else {
                let span = hi - lo;
                if x <= mode {
                    ((x - lo) * (x - lo)) / (span * (mode - lo))
                } else {
                    1.0 - ((hi - x) * (hi - x)) / (span * (hi - mode))
                }
            }
        }
    }
}

/// CDF of the standard normal via the error function, using Abramowitz &
/// Stegun 7.1.26 for `erf` (max abs error ≈ `1.5e-7`).
fn standard_normal_cdf(z: f64) -> f64 {
    0.5 * (1.0 + erf(z / std::f64::consts::SQRT_2))
}

/// Error function `erf(x)` via Abramowitz & Stegun rational approximation
/// 7.1.26. Accurate to about `1.5e-7` — well below any sampling tolerance.
fn erf(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    const A1: f64 = 0.254_829_592;
    const A2: f64 = -0.284_496_736;
    const A3: f64 = 1.421_413_741;
    const A4: f64 = -1.453_152_027;
    const A5: f64 = 1.061_405_429;
    const P: f64 = 0.327_591_1;
    let t = 1.0 / (1.0 + P * x);
    let y = 1.0 - (((((A5 * t + A4) * t) + A3) * t + A2) * t + A1) * t * (-x * x).exp();
    sign * y
}
