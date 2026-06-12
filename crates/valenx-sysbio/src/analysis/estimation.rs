//! Parameter estimation - feature 36.
//!
//! Fit one or more model parameters to experimental time-course data
//! by minimising a residual `sum( (simulated[t_i] - observed[t_i])^2 )`.
//!
//! Two complementary algorithms are layered:
//!
//! - **Latin-hypercube + simulated-annealing pre-stage**. Parameter
//!   bounds rarely shape a convex landscape, so a deterministic
//!   `n_lhs` Latin-hypercube sample of the bounded box is evaluated
//!   and the best is taken as the warm start. An optional simulated-
//!   annealing refinement walks the best point further to escape
//!   shallow local minima.
//! - **Levenberg-Marquardt** refinement. Once the pre-stage has
//!   landed near a basin, LM does the precise nonlinear-least-squares
//!   fit. Each iteration builds a finite-difference Jacobian
//!   `J[i,j] = drr[i] / d theta[j]`, solves the damped normal
//!   equations `(J^T J + lambda I) delta = J^T r` for the step, and
//!   adjusts the trust radius via the standard up/down lambda rule.
//!
//! The driver returns the best-fit parameters, the final residual
//! sum-of-squares, per-parameter Hessian-based standard errors
//! (the inverse `J^T J` diagonals, scaled by the residual variance -
//! the classical asymptotic-normality formula used by COPASI's
//! parameter-estimation task), the number of LM iterations taken, and
//! the LHS / LM trace counts for diagnostics.
//!
//! ## v1 caveats
//!
//! The standard errors are computed from `J^T J` (a *Gauss-Newton*
//! Hessian approximation) rather than the full second derivatives of
//! the residual. This is the same approximation `scipy.optimize.curve_fit`
//! uses; it is exact for a linear model and asymptotically correct for
//! a well-conditioned nonlinear fit. A full Fisher-Information matrix
//! treatment (with profile likelihoods, identifiability analysis) is
//! out of scope for this v1 - the driver returns enough information
//! that a caller can layer that on top if it needs to.
//!
//! The nested-`for i in 0..n` patterns over the matrix scratch space
//! read more naturally than the `iter().enumerate()` adapters Clippy
//! would prefer - the `needless_range_loop` lint is suppressed locally
//! at the matrix-arithmetic sites.

use crate::analysis::param::ParamTarget;
use crate::error::{Result, SysbioError};
use crate::model::Model;
use crate::ode::linalg::solve_linear;
use crate::ode::{EventDrivenTimeCourse, OdeSystem, Rk45, Trajectory};
use crate::stochastic::rng::Rng;

/// One observation point in the experimental data set.
#[derive(Debug, Clone, PartialEq)]
pub struct ObservedPoint {
    /// The time at which the observation was made.
    pub time: f64,
    /// The species index whose amount was observed.
    pub species: usize,
    /// The observed value.
    pub value: f64,
    /// Optional 1-sigma weight - residuals are scaled by `1 / sigma`
    /// so each observation contributes its standard chi-squared term.
    /// `None` uses `1.0`.
    pub sigma: Option<f64>,
}

/// A bounded parameter to fit.
#[derive(Debug, Clone, PartialEq)]
pub struct EstimationTarget {
    /// Which knob of the model to vary.
    pub param: ParamTarget,
    /// Lower bound.
    pub lower: f64,
    /// Upper bound (must exceed `lower`).
    pub upper: f64,
    /// Optional initial guess. If `None` the LHS pre-stage picks one.
    pub initial: Option<f64>,
}

/// Knobs of the estimation driver.
#[derive(Debug, Clone)]
pub struct EstimationOptions {
    /// Number of Latin-hypercube samples in the pre-stage. Zero
    /// disables it (use the supplied `initial` guesses directly).
    pub lhs_samples: usize,
    /// Number of simulated-annealing refinement steps after the LHS
    /// pre-stage. Zero disables annealing.
    pub anneal_steps: usize,
    /// Maximum LM iterations.
    pub max_lm_iter: usize,
    /// LM convergence tolerance on the relative change in the
    /// residual sum of squares.
    pub tol: f64,
    /// PRNG seed for the LHS + SA stages (deterministic from the seed).
    pub seed: u64,
    /// Final integration time for each simulation; the driver picks
    /// the latest observation time if this is smaller.
    pub t_end: f64,
    /// Number of output points sampled on the simulation grid.
    /// Observations are linearly interpolated onto that grid.
    pub n_points: usize,
    /// Use the SBML-L3 event / rule-aware driver instead of the plain
    /// reaction-only ODE integrator. Default `true` so a model with
    /// events / rules is fit correctly out of the box.
    pub use_event_driver: bool,
}

impl Default for EstimationOptions {
    fn default() -> Self {
        EstimationOptions {
            lhs_samples: 16,
            anneal_steps: 0,
            max_lm_iter: 50,
            tol: 1e-8,
            seed: 2026,
            t_end: 0.0,
            n_points: 200,
            use_event_driver: true,
        }
    }
}

/// The output of a parameter-estimation run.
#[derive(Debug, Clone, PartialEq)]
pub struct EstimationReport {
    /// Best-fit parameter values, in the same order as the input
    /// `targets`.
    pub best: Vec<f64>,
    /// Final residual sum of squares.
    pub residual_ss: f64,
    /// Per-parameter Hessian-based standard errors (one entry per
    /// `targets`, `None` if the linear system was singular).
    pub std_errors: Vec<Option<f64>>,
    /// Number of LM iterations taken.
    pub lm_iterations: usize,
    /// Number of model evaluations (LHS + SA + LM).
    pub model_evals: usize,
    /// `true` if LM converged within `tol` before hitting
    /// `max_lm_iter`.
    pub converged: bool,
}

/// Run a parameter-estimation fit (feature 36).
///
/// Returns the best-fit parameters, residual, standard errors and
/// convergence diagnostics. Errors out only on malformed input
/// (an out-of-range bound, an empty target / observation list, a
/// model whose `validate` fails).
pub fn estimate_parameters(
    model: &Model,
    targets: &[EstimationTarget],
    observations: &[ObservedPoint],
    opts: &EstimationOptions,
) -> Result<EstimationReport> {
    if targets.is_empty() {
        return Err(SysbioError::invalid(
            "targets",
            "need at least one parameter to fit",
        ));
    }
    if observations.is_empty() {
        return Err(SysbioError::invalid(
            "observations",
            "need at least one data point",
        ));
    }
    for (i, t) in targets.iter().enumerate() {
        if !t.lower.is_finite() || !t.upper.is_finite() || t.upper <= t.lower {
            return Err(SysbioError::invalid(
                "targets",
                format!("target {i}: upper must exceed lower and both must be finite"),
            ));
        }
        // The supplied target must apply to the model (so we surface
        // a "wrong rate-law" error here rather than deep in LM).
        t.param.read(model)?;
    }
    model.validate()?;

    // Effective simulation horizon - cover every observation.
    let max_obs_t = observations.iter().map(|o| o.time).fold(0.0_f64, f64::max);
    let t_end = opts.t_end.max(max_obs_t * 1.05);
    if t_end <= 0.0 {
        return Err(SysbioError::invalid(
            "observations",
            "observation times must be positive",
        ));
    }

    let mut rng = Rng::new(opts.seed);
    let mut evals = 0usize;

    // Helper: clamp a candidate to the bounds.
    let clamp_to_bounds = |x: &mut [f64]| {
        for (i, v) in x.iter_mut().enumerate() {
            if *v < targets[i].lower {
                *v = targets[i].lower;
            }
            if *v > targets[i].upper {
                *v = targets[i].upper;
            }
        }
    };

    // Helper: evaluate the residual vector at `theta`. Each residual
    // is `(simulated - observed) / sigma`. Returns the residual vec
    // and the per-step trajectory used to compute it (for downstream
    // diagnostics, not re-exported).
    let residuals = |theta: &[f64], evals: &mut usize| -> Result<Vec<f64>> {
        *evals += 1;
        let mut m = model.clone();
        for (k, t) in targets.iter().enumerate() {
            m = t.param.apply(&m, theta[k])?;
        }
        let traj = simulate_for_fit(&m, t_end, opts)?;
        Ok(observation_residuals(&traj, observations))
    };

    // ----- Stage 1: Latin-hypercube pre-stage -----------------------
    let n_params = targets.len();
    let mut best_theta: Vec<f64> = targets
        .iter()
        .map(|t| t.initial.unwrap_or(0.5 * (t.lower + t.upper)))
        .collect();
    let mut best_ss = residual_ss(&residuals(&best_theta, &mut evals)?);

    if opts.lhs_samples > 0 {
        for sample in latin_hypercube_samples(&mut rng, n_params, opts.lhs_samples) {
            // Map each sample column from [0,1) to [lower, upper).
            let mut theta: Vec<f64> = sample
                .iter()
                .enumerate()
                .map(|(i, &u)| {
                    let lo = targets[i].lower;
                    let hi = targets[i].upper;
                    lo + (hi - lo) * u
                })
                .collect();
            clamp_to_bounds(&mut theta);
            let r = residuals(&theta, &mut evals)?;
            let ss = residual_ss(&r);
            if ss < best_ss {
                best_ss = ss;
                best_theta = theta;
            }
        }
    }

    // ----- Stage 2: simulated-annealing refinement -----------------
    if opts.anneal_steps > 0 {
        let mut current = best_theta.clone();
        let mut current_ss = best_ss;
        let t0 = 1.0_f64.max(current_ss * 0.1);
        for step in 0..opts.anneal_steps {
            // Cool linearly from t0 to t0/100.
            let temp = t0 * (1.0 - (step as f64) / (opts.anneal_steps as f64) * 0.99);
            // Propose a Gaussian step scaled to 5% of each bound width.
            let mut proposal = current.clone();
            for i in 0..n_params {
                let w = (targets[i].upper - targets[i].lower) * 0.05;
                proposal[i] = current[i] + w * rng.normal();
            }
            clamp_to_bounds(&mut proposal);
            let r = residuals(&proposal, &mut evals)?;
            let ss = residual_ss(&r);
            let delta = ss - current_ss;
            let accept = delta <= 0.0 || rng.uniform() < (-delta / temp.max(1e-12)).exp();
            if accept {
                current = proposal;
                current_ss = ss;
                if current_ss < best_ss {
                    best_ss = current_ss;
                    best_theta = current.clone();
                }
            }
        }
    }

    // ----- Stage 3: Levenberg-Marquardt ----------------------------
    let mut theta = best_theta;
    let mut r = residuals(&theta, &mut evals)?;
    let mut ss = residual_ss(&r);
    let mut lambda = 1e-3;
    let mut converged = false;
    let mut lm_iters = 0usize;

    for _ in 0..opts.max_lm_iter {
        lm_iters += 1;
        // Build the Jacobian by finite differences.
        let jac = jacobian(&residuals, &theta, &r, &mut evals)?;
        let n_res = r.len();
        // Form J^T J + lambda * diag(J^T J) and J^T r.
        let mut jtj = vec![vec![0.0; n_params]; n_params];
        let mut jtr = vec![0.0; n_params];
        #[allow(clippy::needless_range_loop)]
        for i in 0..n_params {
            for j in 0..n_params {
                let mut acc = 0.0;
                for k in 0..n_res {
                    acc += jac[k][i] * jac[k][j];
                }
                jtj[i][j] = acc;
            }
        }
        #[allow(clippy::needless_range_loop)]
        for i in 0..n_params {
            let mut acc = 0.0;
            for k in 0..n_res {
                acc += jac[k][i] * r[k];
            }
            jtr[i] = acc;
        }
        // Damped normal equations - Marquardt's diagonal scaling.
        let mut damped = jtj.clone();
        #[allow(clippy::needless_range_loop)]
        for i in 0..n_params {
            damped[i][i] *= 1.0 + lambda;
            // Floor the diagonal so a zero column does not collapse.
            if damped[i][i] < lambda * 1e-12 {
                damped[i][i] = lambda * 1e-12;
            }
        }
        let neg_jtr: Vec<f64> = jtr.iter().map(|x| -x).collect();
        let delta = match solve_linear(&damped, &neg_jtr) {
            Some(d) => d,
            None => {
                // Singular - try a larger damping; if still bad, exit.
                lambda *= 10.0;
                if lambda > 1e10 {
                    break;
                }
                continue;
            }
        };
        // Trial step.
        let mut trial: Vec<f64> = theta.iter().zip(&delta).map(|(a, d)| a + d).collect();
        clamp_to_bounds(&mut trial);
        let r_trial = residuals(&trial, &mut evals)?;
        let ss_trial = residual_ss(&r_trial);
        if ss_trial < ss {
            let rel = (ss - ss_trial) / ss.max(1e-300);
            theta = trial;
            r = r_trial;
            ss = ss_trial;
            lambda = (lambda / 3.0).max(1e-12);
            if rel < opts.tol {
                converged = true;
                break;
            }
        } else {
            lambda *= 3.0;
            if lambda > 1e12 {
                break;
            }
        }
    }

    // ----- Standard errors: (J^T J)^-1 diagonals * sigma^2 ---------
    let std_errors = compute_std_errors(&residuals, &theta, &r, n_params, &mut evals)?;

    Ok(EstimationReport {
        best: theta,
        residual_ss: ss,
        std_errors,
        lm_iterations: lm_iters,
        model_evals: evals,
        converged,
    })
}

/// Run the simulator at parameter set `theta`.
fn simulate_for_fit(model: &Model, t_end: f64, opts: &EstimationOptions) -> Result<Trajectory> {
    if opts.use_event_driver
        && (!model.events.is_empty()
            || !model.rules.assignments.is_empty()
            || !model.rules.rates.is_empty())
    {
        let driver = EventDrivenTimeCourse {
            t0: 0.0,
            t_end,
            n_points: opts.n_points,
            ..EventDrivenTimeCourse::new(t_end)
        };
        let evt = driver.run(model)?;
        Ok(evt.trajectory)
    } else {
        let sys = OdeSystem::from_model(model);
        let rk = Rk45::default();
        let traj = rk.integrate(&sys, &model.initial_state(), 0.0, t_end)?;
        Ok(traj)
    }
}

/// Linear interpolation of `traj` at `t` for species `i`.
fn interp_series_at(traj: &Trajectory, t: f64, species: usize) -> f64 {
    if traj.is_empty() {
        return 0.0;
    }
    if t <= traj.times[0] {
        return traj.states[0][species];
    }
    let last = traj.len() - 1;
    if t >= traj.times[last] {
        return traj.states[last][species];
    }
    let mut lo = 0usize;
    let mut hi = last;
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if traj.times[mid] <= t {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let (t0, t1) = (traj.times[lo], traj.times[hi]);
    let w = if t1 > t0 { (t - t0) / (t1 - t0) } else { 0.0 };
    traj.states[lo][species] * (1.0 - w) + traj.states[hi][species] * w
}

/// Build the residual vector for `traj` against the observation list.
fn observation_residuals(traj: &Trajectory, observations: &[ObservedPoint]) -> Vec<f64> {
    observations
        .iter()
        .map(|o| {
            let sim = interp_series_at(traj, o.time, o.species);
            let sigma = o.sigma.unwrap_or(1.0).max(1e-12);
            (sim - o.value) / sigma
        })
        .collect()
}

/// Residual sum of squares.
fn residual_ss(r: &[f64]) -> f64 {
    r.iter().map(|x| x * x).sum()
}

/// Finite-difference Jacobian of the residual vector with respect to
/// each parameter (central difference, relative step `sqrt(eps)`).
fn jacobian<F>(residuals: &F, theta: &[f64], r0: &[f64], evals: &mut usize) -> Result<Vec<Vec<f64>>>
where
    F: Fn(&[f64], &mut usize) -> Result<Vec<f64>>,
{
    let n_p = theta.len();
    let n_r = r0.len();
    let mut jac = vec![vec![0.0; n_p]; n_r];
    let mut th = theta.to_vec();
    let rel = 1e-5;
    for j in 0..n_p {
        let h = rel * theta[j].abs().max(1e-6);
        th[j] = theta[j] + h;
        let r_plus = residuals(&th, evals)?;
        th[j] = theta[j] - h;
        let r_minus = residuals(&th, evals)?;
        th[j] = theta[j];
        let inv = 1.0 / (2.0 * h);
        for k in 0..n_r {
            jac[k][j] = (r_plus[k] - r_minus[k]) * inv;
        }
    }
    Ok(jac)
}

/// Per-parameter Hessian-based standard errors.
///
/// `std_error[j] = sqrt( sigma_hat^2 * inv(J^T J)[j, j] )` where
/// `sigma_hat^2 = sum(r^2) / max(dof, 1)`. Returns `None` for a
/// parameter whose diagonal is non-positive (degenerate column).
fn compute_std_errors<F>(
    residuals: &F,
    theta: &[f64],
    r: &[f64],
    n_p: usize,
    evals: &mut usize,
) -> Result<Vec<Option<f64>>>
where
    F: Fn(&[f64], &mut usize) -> Result<Vec<f64>>,
{
    let jac = jacobian(residuals, theta, r, evals)?;
    let n_res = r.len();
    let mut jtj = vec![vec![0.0; n_p]; n_p];
    #[allow(clippy::needless_range_loop)]
    for i in 0..n_p {
        for j in 0..n_p {
            let mut acc = 0.0;
            for k in 0..n_res {
                acc += jac[k][i] * jac[k][j];
            }
            jtj[i][j] = acc;
        }
    }
    let dof = (n_res as isize - n_p as isize).max(1) as f64;
    let sigma2 = r.iter().map(|x| x * x).sum::<f64>() / dof;
    let inv = invert_matrix(&jtj);
    Ok((0..n_p)
        .map(|i| {
            inv.as_ref().and_then(|m| {
                let v = m[i][i] * sigma2;
                if v.is_finite() && v >= 0.0 {
                    Some(v.sqrt())
                } else {
                    None
                }
            })
        })
        .collect())
}

/// Inverse of a small dense matrix via Gauss-Jordan elimination.
fn invert_matrix(a: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    if n == 0 || a.iter().any(|r| r.len() != n) {
        return None;
    }
    let mut m: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut row = a[i].clone();
            for j in 0..n {
                row.push(if i == j { 1.0 } else { 0.0 });
            }
            row
        })
        .collect();
    for col in 0..n {
        // Pivot.
        let mut pivot = col;
        let mut best = m[col][col].abs();
        #[allow(clippy::needless_range_loop)]
        for r in (col + 1)..n {
            if m[r][col].abs() > best {
                best = m[r][col].abs();
                pivot = r;
            }
        }
        if best < 1e-300 {
            return None;
        }
        m.swap(col, pivot);
        let p = m[col][col];
        for j in 0..(2 * n) {
            m[col][j] /= p;
        }
        for r in 0..n {
            if r == col {
                continue;
            }
            let factor = m[r][col];
            if factor != 0.0 {
                for j in 0..(2 * n) {
                    m[r][j] -= factor * m[col][j];
                }
            }
        }
    }
    Some(m.into_iter().map(|row| row[n..(2 * n)].to_vec()).collect())
}

/// Latin-hypercube design generator over `[0, 1)^k`.
///
/// Each parameter axis is partitioned into `n` strata of equal width;
/// each stratum is sampled once with a uniform jitter; columns are
/// independently shuffled so the design is a true LHS (one sample per
/// row of every parameter's stratum grid). Returns a `Vec<Vec<f64>>`
/// shape `n x k`.
fn latin_hypercube_samples(rng: &mut Rng, k: usize, n: usize) -> Vec<Vec<f64>> {
    let mut samples: Vec<Vec<f64>> = (0..n).map(|_| vec![0.0; k]).collect();
    for col in 0..k {
        // Stratum sample per row.
        let mut col_samples: Vec<f64> = (0..n)
            .map(|row| (row as f64 + rng.uniform()) / n as f64)
            .collect();
        // Fisher-Yates shuffle so different columns are
        // independently re-ordered.
        for i in (1..n).rev() {
            let j = (rng.next_u64() as usize) % (i + 1);
            col_samples.swap(i, j);
        }
        for row in 0..n {
            samples[row][col] = col_samples[row];
        }
    }
    samples
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analysis::param::ParamTarget;
    use crate::model::{RateLaw, Reaction, Species};

    /// Build a model with a single exponential-decay rate `k` to fit.
    fn decay_model(k: f64, a0: f64) -> Model {
        let mut m = Model::new("decay_fit");
        let a = m.add_species(Species::new("A", a0));
        m.add_reaction(Reaction {
            id: "d".into(),
            reactants: vec![(a, 1.0)],
            products: vec![],
            rate_law: RateLaw::MassAction {
                k,
                reactants: vec![(a, 1.0)],
            },
            reversible: false,
        });
        m
    }

    /// Source-decay model, 2 params - source rate `s` and decay rate `k`.
    fn source_decay_model(s: f64, k: f64, a0: f64) -> Model {
        let mut m = Model::new("source_decay");
        let a = m.add_species(Species::new("A", a0));
        m.add_reaction(Reaction {
            id: "src".into(),
            reactants: vec![],
            products: vec![(a, 1.0)],
            rate_law: RateLaw::Constant { rate: s },
            reversible: false,
        });
        m.add_reaction(Reaction {
            id: "dec".into(),
            reactants: vec![(a, 1.0)],
            products: vec![],
            rate_law: RateLaw::MassAction {
                k,
                reactants: vec![(a, 1.0)],
            },
            reversible: false,
        });
        m
    }

    /// Generate synthetic observations from a model run at a known
    /// parameter value.
    fn synthetic_observations(
        model: &Model,
        species: usize,
        times: &[f64],
        t_end: f64,
    ) -> Vec<ObservedPoint> {
        let sys = OdeSystem::from_model(model);
        let rk = Rk45::default();
        let traj = rk
            .integrate(&sys, &model.initial_state(), 0.0, t_end)
            .unwrap();
        times
            .iter()
            .map(|&t| ObservedPoint {
                time: t,
                species,
                value: interp_series_at(&traj, t, species),
                sigma: None,
            })
            .collect()
    }

    #[test]
    fn fits_single_decay_rate_from_synthetic_data() {
        // Generate data at k_true = 1.7; start LM from k = 0.5; LHS
        // pre-stage should bracket the true value within a few
        // samples; LM polishes.
        let truth = decay_model(1.7, 100.0);
        let times = vec![0.1, 0.3, 0.6, 1.0, 1.5, 2.0, 3.0];
        let obs = synthetic_observations(&truth, 0, &times, 3.0);

        // Build a fresh model with the *wrong* k.
        let m_wrong = decay_model(0.5, 100.0);
        let targets = vec![EstimationTarget {
            param: ParamTarget::MassActionK { reaction: 0 },
            lower: 0.01,
            upper: 10.0,
            initial: Some(0.5),
        }];
        let report = estimate_parameters(
            &m_wrong,
            &targets,
            &obs,
            &EstimationOptions {
                lhs_samples: 24,
                anneal_steps: 0,
                max_lm_iter: 80,
                tol: 1e-10,
                seed: 7,
                t_end: 3.0,
                ..EstimationOptions::default()
            },
        )
        .unwrap();
        let k_est = report.best[0];
        assert!(
            (k_est - 1.7).abs() < 0.02,
            "k_est = {k_est:.4} (want 1.7), converged={}, ss={}, iters={}",
            report.converged,
            report.residual_ss,
            report.lm_iterations,
        );
        assert!(report.residual_ss < 1e-4);
        // Standard error should be a small, finite number.
        let se = report.std_errors[0].expect("se finite");
        assert!(se.is_finite() && se < 0.5, "se = {se}");
    }

    #[test]
    fn fits_two_parameter_source_decay() {
        // Fit both source `s` and decay `k` simultaneously.
        let truth = source_decay_model(4.5, 0.9, 0.0);
        let times = vec![0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0];
        let obs = synthetic_observations(&truth, 0, &times, 10.0);

        let m_wrong = source_decay_model(1.0, 0.1, 0.0);
        let targets = vec![
            EstimationTarget {
                param: ParamTarget::ConstantRate { reaction: 0 },
                lower: 0.1,
                upper: 20.0,
                initial: Some(1.0),
            },
            EstimationTarget {
                param: ParamTarget::MassActionK { reaction: 1 },
                lower: 0.01,
                upper: 5.0,
                initial: Some(0.1),
            },
        ];
        let report = estimate_parameters(
            &m_wrong,
            &targets,
            &obs,
            &EstimationOptions {
                lhs_samples: 32,
                max_lm_iter: 120,
                tol: 1e-10,
                seed: 99,
                t_end: 10.0,
                ..EstimationOptions::default()
            },
        )
        .unwrap();
        let s_est = report.best[0];
        let k_est = report.best[1];
        assert!(
            (s_est - 4.5).abs() < 0.1,
            "s_est = {s_est:.3} (want 4.5), converged={}, ss={}, iters={}",
            report.converged,
            report.residual_ss,
            report.lm_iterations,
        );
        assert!((k_est - 0.9).abs() < 0.02, "k_est = {k_est:.4} (want 0.9)");
    }

    #[test]
    fn rejects_bad_target_bounds() {
        let m = decay_model(1.0, 1.0);
        let targets = vec![EstimationTarget {
            param: ParamTarget::MassActionK { reaction: 0 },
            lower: 1.0,
            upper: 0.5, // upper < lower
            initial: None,
        }];
        let obs = vec![ObservedPoint {
            time: 1.0,
            species: 0,
            value: 0.5,
            sigma: None,
        }];
        assert!(estimate_parameters(&m, &targets, &obs, &EstimationOptions::default()).is_err());
    }

    #[test]
    fn empty_targets_is_an_error() {
        let m = decay_model(1.0, 1.0);
        let obs = vec![ObservedPoint {
            time: 1.0,
            species: 0,
            value: 0.5,
            sigma: None,
        }];
        assert!(estimate_parameters(&m, &[], &obs, &EstimationOptions::default()).is_err());
    }

    #[test]
    fn empty_observations_is_an_error() {
        let m = decay_model(1.0, 1.0);
        let targets = vec![EstimationTarget {
            param: ParamTarget::MassActionK { reaction: 0 },
            lower: 0.1,
            upper: 5.0,
            initial: None,
        }];
        assert!(estimate_parameters(&m, &targets, &[], &EstimationOptions::default()).is_err());
    }

    #[test]
    fn latin_hypercube_is_a_valid_design() {
        let mut rng = Rng::new(42);
        let n = 8;
        let k = 3;
        let samples = latin_hypercube_samples(&mut rng, k, n);
        assert_eq!(samples.len(), n);
        for s in &samples {
            assert_eq!(s.len(), k);
            for &v in s {
                assert!((0.0..1.0).contains(&v));
            }
        }
        // Every column should have one sample in each [k/n, (k+1)/n)
        // stratum - that's the defining LHS property.
        for col in 0..k {
            let mut strata = vec![false; n];
            for s in &samples {
                let stratum = (s[col] * n as f64).floor() as usize;
                let idx = stratum.min(n - 1);
                assert!(!strata[idx], "col {col}: stratum {idx} sampled twice");
                strata[idx] = true;
            }
            assert!(strata.iter().all(|&b| b), "col {col}: gap in design");
        }
    }

    #[test]
    fn invert_matrix_round_trip_identity() {
        // Identity inverts to itself.
        let i = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let inv = invert_matrix(&i).unwrap();
        #[allow(clippy::needless_range_loop)]
        for r in 0..3 {
            for c in 0..3 {
                let expect = if r == c { 1.0 } else { 0.0 };
                assert!((inv[r][c] - expect).abs() < 1e-12);
            }
        }
    }
}
