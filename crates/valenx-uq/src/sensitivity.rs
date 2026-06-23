//! Global sensitivity analysis — which inputs drive the output variance?
//!
//! Two complementary, widely-used global methods are provided:
//!
//! * [`sobol_indices`] — **variance-based** Sobol indices. The first-order
//!   index `S_i` is the fraction of the output variance explained by input `i`
//!   acting alone; the total index `S_Ti` adds in every interaction involving
//!   `i`. Estimated with the standard **Saltelli** A/B/AB cross-sampling
//!   design (`A`, `B`, and the `d` hybrid matrices `AB_i`), using the
//!   Jansen (1999) / Saltelli (2010) estimators. Cost: `N·(d + 2)` model
//!   evaluations.
//! * [`morris`] — the **Morris elementary-effects** screening method.
//!   Cheap one-at-a-time trajectory sampling yields, per input, `mu_star`
//!   (mean absolute elementary effect — overall importance) and `sigma`
//!   (their spread — a sign of non-linearity or interaction). Cost:
//!   `r·(d + 1)` model evaluations for `r` trajectories.
//!
//! Both scalarise a (possibly vector-valued) model to a single output index
//! chosen by the caller.
//!
//! ## Honesty note
//!
//! These are **estimates** from a finite sample. The Saltelli estimator is
//! unbiased in the limit but has finite-`N` variance: at small `N`, indices
//! can land slightly outside `[0, 1]` or fail to sum to 1. Morris is a
//! *screening* tool — it ranks and flags non-linearity, it does not partition
//! variance. Increase `N` / `r` (and re-seed to check stability) before
//! trusting tight numbers.

use crate::distribution::Distribution;
use crate::model::Model;
use crate::rng::SplitMix64;
use crate::statistics;

/// First-order and total Sobol indices for every input dimension.
#[derive(Debug, Clone, PartialEq)]
pub struct SobolIndices {
    /// First-order indices `S_i`, one per input. `S_i` is the share of output
    /// variance from input `i` alone (no interactions).
    pub first_order: Vec<f64>,
    /// Total-effect indices `S_Ti`, one per input. `S_Ti` is the share of
    /// output variance involving input `i`, including all interactions, so
    /// `S_Ti >= S_i` always (up to estimator noise).
    pub total: Vec<f64>,
}

/// Estimate Sobol first-order and total indices of `model`'s output
/// `output_index` with respect to its inputs, each distributed per `dists`.
///
/// `n_base` is the base sample size `N`; the total model-evaluation budget is
/// `N·(d + 2)`. The two independent base matrices `A` and `B` are drawn from
/// `dists` via the supplied deterministic PRNG.
///
/// Returns `None` if the analysis is ill-posed:
/// * `dists` is empty, or `n_base < 2`;
/// * `output_index >= model.n_outputs()`;
/// * the estimated output variance is zero (a constant response carries no
///   apportionable sensitivity).
#[must_use]
pub fn sobol_indices(
    model: &dyn Model,
    dists: &[Distribution],
    n_base: usize,
    output_index: usize,
    rng: &mut SplitMix64,
) -> Option<SobolIndices> {
    let d = dists.len();
    if d == 0 || n_base < 2 || output_index >= model.n_outputs() {
        return None;
    }

    // Two independent base sample matrices, N×d each.
    let a = sample_matrix(n_base, dists, rng);
    let b = sample_matrix(n_base, dists, rng);

    // Scalar model responses on A and B.
    let ya = eval_column(model, &a, output_index);
    let yb = eval_column(model, &b, output_index);

    // Total output variance, estimated from A and B pooled together so the
    // denominator is a stable, single estimate shared by every index.
    let mut pooled = ya.clone();
    pooled.extend_from_slice(&yb);
    let var_y = statistics::variance(&pooled)?;
    if var_y <= 0.0 {
        return None;
    }

    let mut first_order = Vec::with_capacity(d);
    let mut total = Vec::with_capacity(d);

    for i in 0..d {
        // AB_i = matrix A with its i-th column replaced by B's i-th column.
        let yab = eval_hybrid_column(model, &a, &b, i, output_index);

        // Jansen (1999) estimators:
        //   S_i   = 1 - (1/2N) Σ (yb - yab)²   / Var(Y)     (first order)
        // and Saltelli (2010):
        //   S_Ti  =      (1/2N) Σ (ya - yab)²   / Var(Y)     (total)
        let n = n_base as f64;

        let mut s_first_acc = 0.0;
        let mut s_total_acc = 0.0;
        for j in 0..n_base {
            let d_first = yb[j] - yab[j];
            s_first_acc += d_first * d_first;
            let d_total = ya[j] - yab[j];
            s_total_acc += d_total * d_total;
        }
        let s_i = 1.0 - (s_first_acc / (2.0 * n)) / var_y;
        let s_ti = (s_total_acc / (2.0 * n)) / var_y;

        first_order.push(s_i);
        total.push(s_ti);
    }

    Some(SobolIndices { first_order, total })
}

/// Per-input results of a Morris elementary-effects screening.
#[derive(Debug, Clone, PartialEq)]
pub struct MorrisResult {
    /// `mu_star_i` — the mean of the **absolute** elementary effects of input
    /// `i`. The primary importance ranking (Campolongo et al. 2007): large
    /// `mu_star` ⇒ influential input.
    pub mu_star: Vec<f64>,
    /// `sigma_i` — the standard deviation of input `i`'s elementary effects.
    /// Large `sigma` ⇒ the effect varies across the input space, i.e. the
    /// input is involved in non-linearities or interactions.
    pub sigma: Vec<f64>,
}

/// Run a Morris elementary-effects screening of `model`'s output
/// `output_index`.
///
/// Generates `n_trajectories` one-at-a-time trajectories on a `levels`-level
/// grid over the hyper-rectangle spanned by each input's effective range (for
/// an unbounded [`Distribution::Normal`] the range is taken as `mean ± 3σ`).
/// Each trajectory perturbs the inputs one at a time, giving one elementary
/// effect per input; the per-input mean-absolute (`mu_star`) and standard
/// deviation (`sigma`) of those effects are returned. Cost:
/// `n_trajectories·(d + 1)` evaluations.
///
/// Returns `None` if `dists` is empty, `n_trajectories == 0`, `levels < 2`, or
/// `output_index >= model.n_outputs()`.
#[must_use]
pub fn morris(
    model: &dyn Model,
    dists: &[Distribution],
    n_trajectories: usize,
    levels: usize,
    output_index: usize,
    rng: &mut SplitMix64,
) -> Option<MorrisResult> {
    let d = dists.len();
    if d == 0 || n_trajectories == 0 || levels < 2 || output_index >= model.n_outputs() {
        return None;
    }

    // The step Δ in the unit cube, the standard Morris choice p/(2(p-1)).
    let p = levels as f64;
    let delta_unit = p / (2.0 * (p - 1.0));
    let ranges: Vec<(f64, f64)> = dists.iter().map(effective_range).collect();

    // Collect elementary effects per input across all trajectories.
    let mut effects: Vec<Vec<f64>> = vec![Vec::with_capacity(n_trajectories); d];

    for _ in 0..n_trajectories {
        // Random starting point on the grid, in unit-cube coordinates.
        let mut base_unit: Vec<f64> = (0..d)
            .map(|_| {
                // Grid level in 0..=(levels-1-1) so a +Δ step stays in [0,1].
                let max_level = levels - 1;
                let lvl =
                    (rng.next_f64() * (max_level as f64 - delta_unit * (p - 1.0)).max(0.0)).floor();
                lvl / (p - 1.0)
            })
            .collect();

        // Random order in which the inputs are perturbed.
        let mut order: Vec<usize> = (0..d).collect();
        shuffle_indices(&mut order, rng);
        // Random sign per input (move forward or backward by Δ).
        let signs: Vec<f64> = (0..d)
            .map(|_| if rng.next_f64() < 0.5 { -1.0 } else { 1.0 })
            .collect();

        // Evaluate at the base point.
        let y_prev_point = unit_to_real(&base_unit, &ranges);
        let mut y_prev = model.evaluate(&y_prev_point)[output_index];

        for &i in &order {
            // Step input i by ±Δ, clamped into the unit cube.
            let mut moved = base_unit.clone();
            let step = signs[i] * delta_unit;
            let mut new_val = moved[i] + step;
            if !(0.0..=1.0).contains(&new_val) {
                // Reflect the step to stay inside the cube.
                new_val = moved[i] - step;
            }
            let actual_step = new_val - moved[i];
            moved[i] = new_val;

            let y_point = unit_to_real(&moved, &ranges);
            let y_new = model.evaluate(&y_point)[output_index];

            // Elementary effect: ΔY / Δx, where Δx is the change in the real
            // (unscaled) input. Guard a zero step.
            let (lo, hi) = ranges[i];
            let real_step = actual_step * (hi - lo);
            if real_step.abs() > f64::EPSILON {
                effects[i].push((y_new - y_prev) / real_step);
            }

            base_unit = moved;
            y_prev = y_new;
        }
    }

    let mut mu_star = Vec::with_capacity(d);
    let mut sigma = Vec::with_capacity(d);
    for ee in &effects {
        if ee.is_empty() {
            mu_star.push(0.0);
            sigma.push(0.0);
            continue;
        }
        let abs_mean = ee.iter().map(|e| e.abs()).sum::<f64>() / ee.len() as f64;
        mu_star.push(abs_mean);
        // Spread of the (signed) elementary effects.
        sigma.push(statistics::std(ee).unwrap_or(0.0));
    }

    Some(MorrisResult { mu_star, sigma })
}

// --- internal helpers ------------------------------------------------------

/// Draw an `n × d` sample matrix (rows are input vectors) from `dists`.
fn sample_matrix(n: usize, dists: &[Distribution], rng: &mut SplitMix64) -> Vec<Vec<f64>> {
    (0..n)
        .map(|_| dists.iter().map(|d| d.sample(rng)).collect())
        .collect()
}

/// Evaluate `model` on every row of `matrix`, keeping output `output_index`.
fn eval_column(model: &dyn Model, matrix: &[Vec<f64>], output_index: usize) -> Vec<f64> {
    matrix
        .iter()
        .map(|row| model.evaluate(row)[output_index])
        .collect()
}

/// Evaluate `model` on the hybrid matrix `AB_i` (rows of `a` with column `i`
/// taken from `b`), keeping output `output_index`.
fn eval_hybrid_column(
    model: &dyn Model,
    a: &[Vec<f64>],
    b: &[Vec<f64>],
    i: usize,
    output_index: usize,
) -> Vec<f64> {
    a.iter()
        .zip(b.iter())
        .map(|(ra, rb)| {
            let mut row = ra.clone();
            row[i] = rb[i];
            model.evaluate(&row)[output_index]
        })
        .collect()
}

/// The finite range used to embed a distribution into the Morris unit cube.
/// Bounded distributions use their support; the unbounded normal uses
/// `mean ± 3σ` (≈ 99.7 % of mass).
fn effective_range(dist: &Distribution) -> (f64, f64) {
    match *dist {
        Distribution::Uniform { lo, hi } => (lo, hi),
        Distribution::Triangular { lo, hi, .. } => (lo, hi),
        Distribution::Normal { mean, std } => (mean - 3.0 * std, mean + 3.0 * std),
    }
}

/// Map a unit-cube point to real input coordinates given per-input ranges.
fn unit_to_real(unit: &[f64], ranges: &[(f64, f64)]) -> Vec<f64> {
    unit.iter()
        .zip(ranges.iter())
        .map(|(&u, &(lo, hi))| lo + u * (hi - lo))
        .collect()
}

/// In-place Fisher–Yates shuffle of an index vector.
fn shuffle_indices(items: &mut [usize], rng: &mut SplitMix64) {
    let len = items.len();
    if len < 2 {
        return;
    }
    for i in (1..len).rev() {
        let j = ((rng.next_f64() * (i as f64 + 1.0)) as usize).min(i);
        items.swap(i, j);
    }
}
