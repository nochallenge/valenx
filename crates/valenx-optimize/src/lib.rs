//! # valenx-optimize
//!
//! Parameter sweeps + optimization. First concrete chunk of
//! [RFC 0011](../../../rfcs/0011-parameter-sweep-optimization.md):
//!
//! - [`SweepConfig`] parser for the `[sweep]` block in `case.toml`.
//! - [`Optimizer`] trait describing the plan/step contract.
//! - [`GridOptimizer`] — full Cartesian product, in-scope reference
//!   implementation.
//!
//! The Latin-Hypercube and gradient-descent optimizers from the RFC
//! are sketched but not implemented; they land as follow-up commits.
//! Same for the wiring into the app's run pipeline (the optimizer
//! produces `DerivedCase` entries; turning each one into a real
//! adapter run is the app's job).

#![forbid(unsafe_code)]
#![allow(missing_docs)]

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Parsed `[sweep]` block from a `case.toml`.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SweepConfig {
    /// Which built-in optimizer drives this sweep.
    pub optimizer: OptimizerKind,
    /// Parameters to vary. The order matters — Cartesian-product
    /// optimizers iterate the last parameter fastest.
    pub parameters: Vec<ParameterSpec>,
    /// Optional objective declaration. Required for optimizers that
    /// need a fitness signal (gradient descent); optional for the
    /// grid sweep (which doesn't pick a "best" run).
    pub objective: Option<ObjectiveSpec>,
    /// Latin Hypercube settings — only consulted when
    /// `optimizer = "latin-hypercube"`.
    #[serde(default)]
    pub latin_hypercube: Option<LatinHypercubeConfig>,
    /// Gradient-descent settings — only consulted when
    /// `optimizer = "gradient-descent"`.
    #[serde(default)]
    pub gradient_descent: Option<GradientDescentConfig>,
}

/// Configuration for the Latin Hypercube optimizer. Each parameter
/// must declare exactly two `values` (interpreted as `[min, max]`);
/// the optimizer divides each parameter's range into `n_samples`
/// equal-width strata, picks the centre of each stratum, then
/// permutes per-parameter using the seeded RNG.
#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize)]
pub struct LatinHypercubeConfig {
    /// Total number of samples to draw.
    pub n_samples: usize,
    /// RNG seed for reproducibility. Two runs with the same seed
    /// produce identical sample sequences.
    #[serde(default)]
    pub seed: u64,
}

/// Configuration for the gradient-descent optimizer. Each parameter
/// declares `[min, max]` box bounds via its two `values`; the optimizer
/// starts at `initial`, then iteratively walks downhill by computing
/// a central-difference gradient and stepping `learning_rate * grad`.
///
/// `initial` keys the parameters by their JSON-Pointer path (same
/// strings as `ParameterSpec.path`).
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct GradientDescentConfig {
    /// Starting point. Must specify a value for every parameter in
    /// `SweepConfig.parameters` and lie within each one's [min, max]
    /// bounds.
    pub initial: BTreeMap<String, f64>,
    /// Central-difference step size used to estimate the gradient.
    /// Same units as the parameter; rule of thumb ≈ 1 % of the
    /// parameter's range.
    pub epsilon: f64,
    /// Multiplier applied to the negative gradient to take one step
    /// downhill. Larger = faster but risk of overshoot.
    pub learning_rate: f64,
    /// Hard cap on the number of descent iterations. The optimizer
    /// stops early if the gradient magnitude drops below numerical
    /// noise.
    pub max_iterations: usize,
}

/// Which optimizer drives a sweep.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum OptimizerKind {
    /// Full Cartesian product. Number of derived runs is the product
    /// of every parameter's value count. **In scope today.**
    Grid,
    /// Latin Hypercube Sampling — N samples maximally spread across
    /// the parameter space. Sketched in RFC 0011, not yet
    /// implemented.
    LatinHypercube,
    /// Finite-difference gradient + line search. Sketched in
    /// RFC 0011, not yet implemented.
    GradientDescent,
}

/// One parameter-value sweep declaration: a TOML-pointer path into
/// the base case + the values to substitute.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ParameterSpec {
    /// JSON-Pointer-style path into the base case.toml. Each segment
    /// after the leading `/` is a TOML key or numeric array index.
    /// Example: `/boundaries/inlet/velocity/0` to vary the X
    /// component of the inlet's velocity vector.
    pub path: String,
    /// Values to substitute at `path`. Type is JSON to allow
    /// scalars or strings (turbulence model names, etc.).
    pub values: Vec<serde_json::Value>,
}

/// Optimizer objective. The optimizer reads each completed run's
/// `Results.scalars` to extract `metric` and tries to minimise or
/// maximise it.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ObjectiveSpec {
    pub metric: String,
    #[serde(default)]
    pub direction: ObjectiveDirection,
}

/// Whether the optimizer should drive the objective metric down or
/// up.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ObjectiveDirection {
    /// Find runs with the smallest metric value (default).
    #[default]
    Minimize,
    /// Find runs with the largest metric value.
    Maximize,
}

/// One concrete run the optimizer wants the harness to execute.
#[derive(Clone, Debug)]
pub struct DerivedCase {
    /// Stable id within this sweep — used as the case directory name
    /// under the parent workdir. Format: `<base-name>-<seq>` where
    /// `seq` is a zero-padded index.
    pub id: String,
    /// Substitutions the harness must apply to the base case.toml
    /// before invoking the adapter. Values are by JSON-Pointer path.
    pub substitutions: BTreeMap<String, serde_json::Value>,
}

/// One completed run plus the optimizer-relevant scalars extracted
/// from it. Populated by the harness; passed back to
/// [`Optimizer::step`] so adaptive optimizers can decide what to
/// run next.
#[derive(Clone, Debug)]
pub struct CompletedRun {
    pub id: String,
    pub objective_value: Option<f64>,
    pub succeeded: bool,
}

/// Iterative optimizer contract. `plan()` returns the initial set of
/// runs; once those finish, `step()` decides whether more runs are
/// needed. Single-shot optimizers (Grid, LHS) return
/// [`OptimizerStep::Done`] from the first `step()` call.
///
/// `plan()` takes `&mut self` because adaptive optimizers
/// (gradient-descent etc.) need to capture the sweep config + initial
/// state for use in subsequent `step()` calls. Pure-batch optimizers
/// just don't mutate.
pub trait Optimizer {
    /// Optimizer id for UI display + audit logging.
    fn id(&self) -> &str;

    /// Initial set of derived runs. Adaptive optimizers also use this
    /// to seed their internal state for subsequent `step()` calls.
    fn plan(&mut self, sweep: &SweepConfig) -> Result<Vec<DerivedCase>, OptimizerError>;

    /// Process a completed batch and decide whether to schedule
    /// more.
    fn step(&mut self, completed: &[CompletedRun]) -> OptimizerStep;
}

/// Result of one [`Optimizer::step`] call.
#[derive(Clone, Debug)]
pub enum OptimizerStep {
    /// Sweep is complete. Harness shows the final summary.
    Done,
    /// Schedule these additional runs, then call `step()` again
    /// when they finish.
    More(Vec<DerivedCase>),
}

/// Errors raised by the sweep optimizers.
#[derive(Debug, Error)]
pub enum OptimizerError {
    /// The user-supplied [`SweepConfig`] is missing a required field
    /// or references an unknown metric.
    #[error("invalid sweep config: {0}")]
    InvalidConfig(String),
    /// The requested [`OptimizerKind`] is sketched but not yet wired
    /// up.
    #[error("optimizer `{0}` is not yet implemented")]
    NotImplemented(&'static str),
    /// Round-14 M6 (round-4 sister gap): the requested sample /
    /// iteration count exceeds the per-optimizer cap. Today the only
    /// site that raises this is LHS via [`MAX_LHS_SAMPLES`]; the
    /// variant is shaped to take any future optimizer's cap.
    #[error("{optimizer}: requested {requested} samples exceeds the {cap}-sample cap")]
    TooManySamples {
        /// Optimizer id (e.g. `latin-hypercube`).
        optimizer: &'static str,
        /// Sample count the user supplied.
        requested: usize,
        /// Cap that was exceeded.
        cap: usize,
    },
}

// ---------------------------------------------------------------------------
// Grid optimizer
// ---------------------------------------------------------------------------

/// Round-4 hardening: upper bound on the total number of cells in a
/// Cartesian-product grid sweep. 10 million is generous — production
/// sweeps run dozens to thousands of cases, not millions — and small
/// enough that the worst-case `Vec::with_capacity(MAX)` allocation
/// stays well inside any reasonable host's RAM budget.
pub const MAX_GRID_CELLS: usize = 10_000_000;

/// Round-14 M6 (round-4 sister gap): upper bound on the LHS
/// optimizer's `n_samples`. Pre-fix `cfg.n_samples` was honoured as
/// the user supplied it — a hostile / accidentally-typed
/// `n_samples = 10_000_000_000` would allocate the strata Vec (one
/// `Vec<f64>` per parameter, each `n_samples` long) before any
/// optimizer iteration ran, OOMing the host. 1 million is well past
/// any realistic LHS sweep (production runs top out at a few
/// thousand samples to keep wall-time tractable) and matches the
/// grid cap's order of magnitude.
pub const MAX_LHS_SAMPLES: usize = 1_000_000;

/// Full Cartesian product of every parameter's `values` array.
///
/// For a sweep declaring 3 × 5 × 2 parameters, `plan()` returns
/// 30 [`DerivedCase`] entries with substitutions covering every
/// combination. `step()` always returns [`OptimizerStep::Done`] —
/// the grid sweep schedules its full run set up front.
pub struct GridOptimizer;

impl GridOptimizer {
    /// New, stateless grid optimizer.
    pub fn new() -> Self {
        Self
    }
}

impl Default for GridOptimizer {
    fn default() -> Self {
        Self::new()
    }
}

impl Optimizer for GridOptimizer {
    fn id(&self) -> &str {
        "grid"
    }

    fn plan(&mut self, sweep: &SweepConfig) -> Result<Vec<DerivedCase>, OptimizerError> {
        if sweep.parameters.is_empty() {
            return Err(OptimizerError::InvalidConfig(
                "grid sweep needs at least one [[sweep.parameter]]".into(),
            ));
        }
        for p in &sweep.parameters {
            if p.values.is_empty() {
                return Err(OptimizerError::InvalidConfig(format!(
                    "parameter `{}` has empty values list",
                    p.path
                )));
            }
        }
        // Iterate the Cartesian product. Last parameter varies
        // fastest so the row order reads like nested loops in
        // declaration order.
        //
        // Round-4 hardening: a sweep with 10 parameters × 1000 values
        // each is 1e30 cells — the silent `.product()` would wrap
        // around `usize` and either underallocate or panic on
        // `Vec::with_capacity`. `checked_mul` over the iterator stops
        // the wrap, and a hard cap of 10M cells protects against the
        // OOM path even when the math doesn't overflow.
        let total: usize = sweep
            .parameters
            .iter()
            .map(|p| p.values.len())
            .try_fold(1usize, |acc, n| acc.checked_mul(n))
            .ok_or_else(|| {
                OptimizerError::InvalidConfig(format!(
                    "grid sweep cell count overflowed usize — \
                     reduce per-parameter value counts (limit: {MAX_GRID_CELLS} cells)"
                ))
            })?;
        if total > MAX_GRID_CELLS {
            return Err(OptimizerError::InvalidConfig(format!(
                "grid sweep would produce {total} cells, exceeding the \
                 {MAX_GRID_CELLS}-cell cap — reduce per-parameter value counts \
                 or split into multiple sweeps"
            )));
        }
        let pad_width = total.to_string().len();
        let mut out: Vec<DerivedCase> = Vec::with_capacity(total);
        let mut indices = vec![0usize; sweep.parameters.len()];
        for seq in 0..total {
            let mut substitutions: BTreeMap<String, serde_json::Value> = BTreeMap::new();
            for (param, &idx) in sweep.parameters.iter().zip(indices.iter()) {
                substitutions.insert(param.path.clone(), param.values[idx].clone());
            }
            out.push(DerivedCase {
                id: format!("sweep-{seq:0pad_width$}"),
                substitutions,
            });
            // Increment indices like an odometer — last position
            // varies fastest.
            for i in (0..indices.len()).rev() {
                indices[i] += 1;
                if indices[i] < sweep.parameters[i].values.len() {
                    break;
                }
                indices[i] = 0;
            }
        }
        Ok(out)
    }

    fn step(&mut self, _completed: &[CompletedRun]) -> OptimizerStep {
        // Grid sweeps schedule everything up front — once the harness
        // calls step(), there's nothing more to plan.
        OptimizerStep::Done
    }
}

/// Build the right [`Optimizer`] for the given [`OptimizerKind`].
/// Unknown / not-yet-implemented kinds return [`OptimizerError::NotImplemented`].
pub fn make_optimizer(kind: OptimizerKind) -> Result<Box<dyn Optimizer>, OptimizerError> {
    match kind {
        OptimizerKind::Grid => Ok(Box::new(GridOptimizer::new())),
        OptimizerKind::LatinHypercube => Ok(Box::new(LatinHypercubeOptimizer::new())),
        OptimizerKind::GradientDescent => Ok(Box::new(GradientDescentOptimizer::new())),
    }
}

// ---------------------------------------------------------------------------
// Gradient-descent optimizer
// ---------------------------------------------------------------------------

/// Finite-difference gradient + fixed-step descent. On each
/// iteration the optimizer schedules `1 + 2 * D` derived runs:
///
/// - 1 evaluation at the current point.
/// - 1 evaluation per parameter at `current + epsilon * e_j`.
/// - 1 evaluation per parameter at `current - epsilon * e_j`.
///
/// After all `1 + 2D` runs complete, [`Optimizer::step`] computes the
/// central-difference gradient, takes a step
/// `current ← current − learning_rate * grad` (clamped to the box
/// bounds), and either schedules the next iteration's batch or
/// returns [`OptimizerStep::Done`] when:
///
/// - the iteration count reaches `max_iterations`, or
/// - the gradient norm drops below 1e-12 (numerical floor — same
///   tolerance scipy uses by default).
///
/// Genuinely stateful — a single instance carries the current point,
/// iteration counter, and last-batch metadata. Don't share across
/// sweeps.
pub struct GradientDescentOptimizer {
    /// `None` until `plan()` runs; `Some(_)` while a sweep is active.
    state: Option<GdState>,
}

struct GdState {
    config: GradientDescentConfig,
    objective: ObjectiveSpec,
    parameters: Vec<ParameterSpec>,
    /// Current point, indexed in declaration order to match `parameters`.
    point: Vec<f64>,
    iteration: usize,
    /// One trace point per completed iteration. Useful for the UI's
    /// optimization history / convergence chart.
    trace: Vec<GdTracePoint>,
}

/// One row of a [`GradientDescentOptimizer`]'s convergence trace.
/// Records the iteration index, the parameter point that was
/// evaluated, the centre's objective value, and the gradient norm
/// the step was based on (for the convergence chart).
#[derive(Clone, Debug, PartialEq)]
pub struct GdTracePoint {
    pub iteration: usize,
    /// Parameter values in declaration order, matching the sweep's
    /// `parameters` array.
    pub point: Vec<f64>,
    /// Centre evaluation of the objective at this point. `None` when
    /// the centre's CompletedRun didn't carry a value.
    pub objective_value: Option<f64>,
    /// L2 norm of the gradient at this point (for "did the descent
    /// flatten out?" diagnostics). `None` for the seed iteration
    /// (no gradient has been computed yet).
    pub gradient_norm: Option<f64>,
}

impl Default for GradientDescentOptimizer {
    fn default() -> Self {
        Self::new()
    }
}

impl GradientDescentOptimizer {
    /// New optimizer with empty state — the first [`Optimizer::plan`]
    /// call seeds the descent from the sweep config.
    pub fn new() -> Self {
        Self { state: None }
    }

    /// Convergence trace produced so far. Empty until `plan()`
    /// runs; one entry per iteration after each `step()` call.
    /// Pair with the UI's convergence chart to surface "is this
    /// optimisation actually converging" without re-running.
    pub fn trace(&self) -> &[GdTracePoint] {
        match &self.state {
            Some(s) => &s.trace,
            None => &[],
        }
    }
}

impl Optimizer for GradientDescentOptimizer {
    fn id(&self) -> &str {
        "gradient-descent"
    }

    fn plan(&mut self, sweep: &SweepConfig) -> Result<Vec<DerivedCase>, OptimizerError> {
        let cfg = sweep.gradient_descent.as_ref().ok_or_else(|| {
            OptimizerError::InvalidConfig(
                "[sweep.gradient_descent] block missing — required for `optimizer = \"gradient-descent\"`".into(),
            )
        })?;
        if sweep.objective.is_none() {
            return Err(OptimizerError::InvalidConfig(
                "gradient-descent requires [sweep.objective] (which scalar to minimise/maximise)"
                    .into(),
            ));
        }
        if sweep.parameters.is_empty() {
            return Err(OptimizerError::InvalidConfig(
                "no parameters declared".into(),
            ));
        }
        if cfg.epsilon <= 0.0 {
            return Err(OptimizerError::InvalidConfig(format!(
                "epsilon must be > 0; got {}",
                cfg.epsilon
            )));
        }
        if cfg.max_iterations == 0 {
            return Err(OptimizerError::InvalidConfig(
                "max_iterations must be > 0".into(),
            ));
        }

        // Validate every parameter has [min, max] bounds and the
        // initial point lies within them. Pre-extract the starting
        // point in declaration order.
        let mut point: Vec<f64> = Vec::with_capacity(sweep.parameters.len());
        for p in &sweep.parameters {
            let (min, max) = parameter_bounds(p)?;
            let v = cfg.initial.get(&p.path).copied().ok_or_else(|| {
                OptimizerError::InvalidConfig(format!(
                    "[sweep.gradient_descent.initial] missing entry for `{}`",
                    p.path
                ))
            })?;
            if v < min || v > max {
                return Err(OptimizerError::InvalidConfig(format!(
                    "initial value {v} for `{}` is outside bounds [{min}, {max}]",
                    p.path
                )));
            }
            point.push(v);
        }

        let batch = build_gd_batch(0, &point, &sweep.parameters, cfg);
        // Capture the sweep + initial point so subsequent step()
        // calls can compute gradients + step the same point. Seed
        // the trace with the iteration-0 point (objective_value +
        // gradient_norm fill in once step() processes the first
        // batch).
        let seed_trace = GdTracePoint {
            iteration: 0,
            point: point.clone(),
            objective_value: None,
            gradient_norm: None,
        };
        self.state = Some(GdState {
            config: cfg.clone(),
            objective: sweep
                .objective
                .clone()
                .expect("objective presence checked above"),
            parameters: sweep.parameters.clone(),
            point,
            iteration: 0,
            trace: vec![seed_trace],
        });
        Ok(batch)
    }

    fn step(&mut self, completed: &[CompletedRun]) -> OptimizerStep {
        // Adaptive optimizers need persistent state across calls.
        // The harness is responsible for calling plan() once + step()
        // repeatedly on the *same* instance; we carry state on `self`.
        let state = match &mut self.state {
            Some(s) => s,
            None => {
                // Harness called step() before plan(), or the
                // optimizer was constructed but never planned. Treat
                // it as a "nothing to schedule" answer rather than
                // panicking — the harness should call plan() first.
                return OptimizerStep::Done;
            }
        };

        // Extract the central + perturbation values in the order our
        // batch was constructed. The harness guarantees `completed`
        // covers exactly the batch we last emitted.
        let cfg = &state.config;
        let d = state.parameters.len();
        let expected = 1 + 2 * d;
        if completed.len() != expected {
            // Partial batch — surface as Done so the harness
            // surfaces the failure rather than looping forever.
            return OptimizerStep::Done;
        }

        // Map batch ids back to gradient-FD indices.
        let center_id = format!("gd-iter{:0>2}-center", state.iteration);
        let center = completed
            .iter()
            .find(|c| c.id == center_id)
            .and_then(|c| c.objective_value);
        let mut grad: Vec<f64> = Vec::with_capacity(d);
        for j in 0..d {
            let plus_id = format!("gd-iter{:0>2}-p{:0>2}-plus", state.iteration, j);
            let minus_id = format!("gd-iter{:0>2}-p{:0>2}-minus", state.iteration, j);
            let plus = completed
                .iter()
                .find(|c| c.id == plus_id)
                .and_then(|c| c.objective_value);
            let minus = completed
                .iter()
                .find(|c| c.id == minus_id)
                .and_then(|c| c.objective_value);
            match (plus, minus) {
                (Some(p), Some(m)) => grad.push((p - m) / (2.0 * cfg.epsilon)),
                _ => {
                    // Missing FD evaluation = can't proceed.
                    return OptimizerStep::Done;
                }
            }
        }

        // Maximise = walk uphill (descend the negation).
        let direction_sign = match state.objective.direction {
            ObjectiveDirection::Minimize => -1.0,
            ObjectiveDirection::Maximize => 1.0,
        };

        // Step + clamp to bounds.
        let mut grad_norm_sq = 0.0;
        for (j, &g) in grad.iter().enumerate().take(d) {
            grad_norm_sq += g * g;
            let new_v = state.point[j] + direction_sign * cfg.learning_rate * g;
            let (min, max) =
                parameter_bounds(&state.parameters[j]).expect("validated at plan() time");
            state.point[j] = new_v.clamp(min, max);
        }
        let grad_norm = grad_norm_sq.sqrt();

        // Backfill the iteration's seed trace point with the centre
        // value + gradient norm computed this step. The next step's
        // record uses the new point captured below.
        if let Some(last) = state.trace.last_mut() {
            if last.iteration == state.iteration {
                last.objective_value = center;
                last.gradient_norm = Some(grad_norm);
            }
        }

        // Numerical-floor stop. scipy's L-BFGS uses 1e-12 by default
        // for the gradient norm; matching the convention.
        if grad_norm_sq < 1e-24 {
            return OptimizerStep::Done;
        }

        state.iteration += 1;
        // Record the new point as the seed for the next iteration.
        // objective + gradient_norm fill in next step().
        state.trace.push(GdTracePoint {
            iteration: state.iteration,
            point: state.point.clone(),
            objective_value: None,
            gradient_norm: None,
        });
        if state.iteration >= cfg.max_iterations {
            return OptimizerStep::Done;
        }

        let next = build_gd_batch(state.iteration, &state.point, &state.parameters, cfg);
        OptimizerStep::More(next)
    }
}

/// Build one iteration's batch of `1 + 2D` derived cases:
/// center + plus/minus perturbation per parameter. Pulled out so
/// `plan()` and `step()` build identical-shape batches.
fn build_gd_batch(
    iteration: usize,
    point: &[f64],
    parameters: &[ParameterSpec],
    cfg: &GradientDescentConfig,
) -> Vec<DerivedCase> {
    let d = parameters.len();
    let mut out: Vec<DerivedCase> = Vec::with_capacity(1 + 2 * d);
    // Center: parameters set to current point.
    out.push(DerivedCase {
        id: format!("gd-iter{iteration:0>2}-center"),
        substitutions: pack_substitutions(parameters, point),
    });
    // Per-parameter +epsilon and -epsilon perturbations.
    for j in 0..d {
        let mut plus = point.to_vec();
        let mut minus = point.to_vec();
        plus[j] += cfg.epsilon;
        minus[j] -= cfg.epsilon;
        out.push(DerivedCase {
            id: format!("gd-iter{iteration:0>2}-p{j:0>2}-plus"),
            substitutions: pack_substitutions(parameters, &plus),
        });
        out.push(DerivedCase {
            id: format!("gd-iter{iteration:0>2}-p{j:0>2}-minus"),
            substitutions: pack_substitutions(parameters, &minus),
        });
    }
    out
}

fn pack_substitutions(
    parameters: &[ParameterSpec],
    point: &[f64],
) -> BTreeMap<String, serde_json::Value> {
    parameters
        .iter()
        .zip(point.iter())
        .map(|(p, v)| (p.path.clone(), serde_json::Value::from(*v)))
        .collect()
}

fn parameter_bounds(p: &ParameterSpec) -> Result<(f64, f64), OptimizerError> {
    if p.values.len() != 2 {
        return Err(OptimizerError::InvalidConfig(format!(
            "parameter `{}`: gradient-descent needs exactly 2 values \
             ([min, max]); got {}",
            p.path,
            p.values.len()
        )));
    }
    let min = p.values[0].as_f64().ok_or_else(|| {
        OptimizerError::InvalidConfig(format!(
            "parameter `{}`: min value must be a number",
            p.path
        ))
    })?;
    let max = p.values[1].as_f64().ok_or_else(|| {
        OptimizerError::InvalidConfig(format!(
            "parameter `{}`: max value must be a number",
            p.path
        ))
    })?;
    // Reject NaN bounds and the empty range max == min as well —
    // wrap the comparison explicitly so clippy doesn't read the
    // negated `>` as a partially-ordered-comparison lint.
    if !(min.is_finite() && max.is_finite() && max > min) {
        return Err(OptimizerError::InvalidConfig(format!(
            "parameter `{}`: max ({max}) must be > min ({min})",
            p.path
        )));
    }
    Ok((min, max))
}

// ---------------------------------------------------------------------------
// Latin Hypercube optimizer
// ---------------------------------------------------------------------------

/// Latin Hypercube Sampling (LHS) — N samples maximally spread
/// across a continuous parameter space. Each parameter declares
/// `[min, max]` bounds via its two `values`; the optimizer divides
/// each parameter's range into N equal-width strata, picks the
/// centre of each stratum, then permutes per-parameter using a
/// seeded RNG so the combinations are diverse.
///
/// Why LHS over a grid: with D parameters and N samples per
/// dimension, a grid produces N^D runs; LHS gives N runs at the same
/// per-dimension resolution. For a typical 5-parameter sweep at 10
/// samples per dim that's 10 vs 100,000 runs.
#[derive(Default)]
pub struct LatinHypercubeOptimizer;

impl LatinHypercubeOptimizer {
    /// New, stateless LHS optimizer.
    pub fn new() -> Self {
        Self
    }
}

impl Optimizer for LatinHypercubeOptimizer {
    fn id(&self) -> &str {
        "latin-hypercube"
    }

    fn plan(&mut self, sweep: &SweepConfig) -> Result<Vec<DerivedCase>, OptimizerError> {
        let cfg = sweep.latin_hypercube.as_ref().ok_or_else(|| {
            OptimizerError::InvalidConfig(
                "[sweep.latin_hypercube] block missing — required for `optimizer = \"latin-hypercube\"`".into(),
            )
        })?;
        if cfg.n_samples == 0 {
            return Err(OptimizerError::InvalidConfig(
                "n_samples must be > 0".into(),
            ));
        }
        // Round-14 M6 (round-4 sister gap): refuse pathologically
        // large `n_samples` before the strata Vec allocation. Pre-fix
        // `n_samples = 10_000_000_000` flowed into the strata builder
        // and OOMed before any optimizer iteration ran.
        if cfg.n_samples > MAX_LHS_SAMPLES {
            return Err(OptimizerError::TooManySamples {
                optimizer: "latin-hypercube",
                requested: cfg.n_samples,
                cap: MAX_LHS_SAMPLES,
            });
        }
        if sweep.parameters.is_empty() {
            return Err(OptimizerError::InvalidConfig(
                "no parameters declared".into(),
            ));
        }

        // Each parameter's `values` is interpreted as `[min, max]`.
        // We pre-extract the bounds + reject anything that doesn't
        // fit so failures are caught before the loop body fires.
        let mut bounds: Vec<(f64, f64)> = Vec::with_capacity(sweep.parameters.len());
        for p in &sweep.parameters {
            if p.values.len() != 2 {
                return Err(OptimizerError::InvalidConfig(format!(
                    "parameter `{}`: latin-hypercube needs exactly 2 values \
                     ([min, max]); got {}",
                    p.path,
                    p.values.len()
                )));
            }
            let min = p.values[0].as_f64().ok_or_else(|| {
                OptimizerError::InvalidConfig(format!(
                    "parameter `{}`: min value must be a number",
                    p.path
                ))
            })?;
            let max = p.values[1].as_f64().ok_or_else(|| {
                OptimizerError::InvalidConfig(format!(
                    "parameter `{}`: max value must be a number",
                    p.path
                ))
            })?;
            // Same NaN-/infinite-rejection pattern as parameter_bounds
            // — explicit so clippy doesn't fire on the negated `>`.
            if !(min.is_finite() && max.is_finite() && max > min) {
                return Err(OptimizerError::InvalidConfig(format!(
                    "parameter `{}`: max ({max}) must be > min ({min})",
                    p.path
                )));
            }
            bounds.push((min, max));
        }

        // Build N equal-width stratum centres per parameter.
        let n = cfg.n_samples;
        let mut strata: Vec<Vec<f64>> = bounds
            .iter()
            .map(|(min, max)| {
                let step = (max - min) / n as f64;
                (0..n).map(|k| min + (k as f64 + 0.5) * step).collect()
            })
            .collect();

        // Per-parameter permutation: Fisher-Yates with the seeded
        // RNG. We mix the parameter index into the seed so each
        // dimension gets a distinct permutation even with seed=0.
        for (dim, values) in strata.iter_mut().enumerate() {
            let mut rng = SplitMix64::new(
                cfg.seed
                    .wrapping_add(dim as u64)
                    .wrapping_mul(0x9E3779B97F4A7C15),
            );
            for i in (1..values.len()).rev() {
                let j = (rng.next_u64() as usize) % (i + 1);
                values.swap(i, j);
            }
        }

        // Assemble N derived cases — sample i takes the i-th value
        // from each parameter's permuted strata.
        let id_width = id_width_for(n);
        let mut out: Vec<DerivedCase> = Vec::with_capacity(n);
        for sample in 0..n {
            let mut substitutions = std::collections::BTreeMap::new();
            for (param, values) in sweep.parameters.iter().zip(strata.iter()) {
                let v = values[sample];
                substitutions.insert(param.path.clone(), serde_json::Value::from(v));
            }
            out.push(DerivedCase {
                id: format!("sweep-{sample:0>id_width$}"),
                substitutions,
            });
        }
        Ok(out)
    }

    fn step(&mut self, _completed: &[CompletedRun]) -> OptimizerStep {
        // Single-shot like the grid optimizer.
        OptimizerStep::Done
    }
}

/// SplitMix64 — a tiny seedable RNG with great statistical
/// properties for what it is. Same algorithm Java's SplittableRandom
/// uses for its golden-section seed mixing. Inlined so we don't pull
/// in `rand` for one optimizer.
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
}

/// Width needed to zero-pad sample ids so they sort lexicographically.
fn id_width_for(n_samples: usize) -> usize {
    n_samples.saturating_sub(1).to_string().len().max(3)
}

// ---------------------------------------------------------------------------
// TOML-Pointer mutator + DerivedCase materialisation
// ---------------------------------------------------------------------------

/// Apply a [`DerivedCase`]'s substitutions to a base TOML document.
/// Returns the mutated TOML as a string ready to be written into the
/// derived run's `case.toml`.
///
/// Each substitution path is a JSON-Pointer-style string (RFC 6901
/// shape — leading `/`, segments separated by `/`). Numeric segments
/// index TOML arrays; alphanumeric segments index TOML tables. Paths
/// that don't resolve (missing key, out-of-range index, type mismatch)
/// produce a structured error rather than silently no-op'ing.
pub fn materialise_case(base_toml: &str, derived: &DerivedCase) -> Result<String, OptimizerError> {
    let mut doc: toml::Value = toml::from_str(base_toml)
        .map_err(|e| OptimizerError::InvalidConfig(format!("base case.toml parse: {e}")))?;
    for (path, json_value) in &derived.substitutions {
        let toml_value = json_to_toml(json_value)
            .map_err(|e| OptimizerError::InvalidConfig(format!("substitution at `{path}`: {e}")))?;
        apply_pointer(&mut doc, path, toml_value)
            .map_err(|e| OptimizerError::InvalidConfig(format!("substitute `{path}`: {e}")))?;
    }

    // Stamp a `[sweep.derived]` block recording the numeric inputs
    // for this derived run. The dataset assembler reads this back to
    // recover the per-sample input vector. We use short keys derived
    // from each pointer's last segment so the block is readable; the
    // full pointer is preserved as `<short-name>__path` if there's
    // any ambiguity (rare today, future-proof).
    let mut derived_table = toml::value::Table::new();
    let mut paths_table = toml::value::Table::new();
    for (path, json_value) in &derived.substitutions {
        let short = pointer_short_name(path);
        if let Some(num) = json_value.as_f64() {
            derived_table.insert(short.clone(), toml::Value::Float(num));
        } else if let Some(num) = json_value.as_i64() {
            derived_table.insert(short.clone(), toml::Value::Integer(num));
        }
        // String / bool inputs go into a sidecar map for traceability
        // but aren't picked up as numeric features.
        paths_table.insert(short.clone(), toml::Value::String(path.clone()));
    }
    if !derived_table.is_empty() || !paths_table.is_empty() {
        let mut sweep_block = match doc.get("sweep").and_then(|v| v.as_table()) {
            Some(t) => t.clone(),
            None => toml::value::Table::new(),
        };
        sweep_block.insert("derived".into(), toml::Value::Table(derived_table));
        sweep_block.insert("derived_paths".into(), toml::Value::Table(paths_table));
        if let toml::Value::Table(top) = &mut doc {
            top.insert("sweep".into(), toml::Value::Table(sweep_block));
        }
    }

    toml::to_string_pretty(&doc)
        .map_err(|e| OptimizerError::InvalidConfig(format!("re-serialise mutated case: {e}")))
}

/// Pull a short, file-friendly key out of a JSON-Pointer path.
/// `/boundaries/inlet/velocity/0` → `velocity_0`. Used by
/// `materialise_case` to label the `[sweep.derived]` block.
fn pointer_short_name(pointer: &str) -> String {
    let segs: Vec<&str> = pointer.trim_start_matches('/').split('/').collect();
    if segs.len() <= 1 {
        return segs.first().copied().unwrap_or("input").to_string();
    }
    // Last 2 segments joined by underscore — "velocity/0" -> "velocity_0".
    let last = segs[segs.len() - 1];
    let prev = segs[segs.len() - 2];
    if last.parse::<usize>().is_ok() {
        format!("{prev}_{last}")
    } else {
        last.to_string()
    }
}

/// Apply a single substitution at `pointer`. Walks down the document
/// segment-by-segment, replacing the leaf value with `new_value`.
fn apply_pointer(
    doc: &mut toml::Value,
    pointer: &str,
    new_value: toml::Value,
) -> Result<(), String> {
    let segments = parse_pointer(pointer)?;
    if segments.is_empty() {
        return Err("empty pointer (use `/key/...`)".to_string());
    }
    let (last, intermediate) = segments.split_last().expect("non-empty");
    let mut cur = doc;
    for seg in intermediate {
        cur = descend(cur, seg)?;
    }
    set_leaf(cur, last, new_value)
}

fn descend<'a>(cur: &'a mut toml::Value, seg: &str) -> Result<&'a mut toml::Value, String> {
    match cur {
        toml::Value::Table(t) => t
            .get_mut(seg)
            .ok_or_else(|| format!("missing key `{seg}` in table")),
        toml::Value::Array(a) => {
            let idx: usize = seg
                .parse()
                .map_err(|_| format!("array index `{seg}` is not a number"))?;
            let len = a.len();
            a.get_mut(idx)
                .ok_or_else(|| format!("array index `{idx}` out of range (len {len})"))
        }
        other => Err(format!(
            "cannot descend into segment `{seg}`: parent is {} not table/array",
            value_kind(other)
        )),
    }
}

fn set_leaf(cur: &mut toml::Value, last: &str, new_value: toml::Value) -> Result<(), String> {
    match cur {
        toml::Value::Table(t) => {
            t.insert(last.to_string(), new_value);
            Ok(())
        }
        toml::Value::Array(a) => {
            let idx: usize = last
                .parse()
                .map_err(|_| format!("array index `{last}` is not a number"))?;
            if idx >= a.len() {
                return Err(format!(
                    "array index `{idx}` out of range (len {})",
                    a.len()
                ));
            }
            a[idx] = new_value;
            Ok(())
        }
        other => Err(format!(
            "cannot set on `{last}`: parent is {} not table/array",
            value_kind(other)
        )),
    }
}

/// Split a JSON-Pointer-style path into its segments. Strips the
/// leading `/` and rejects empty paths.
fn parse_pointer(pointer: &str) -> Result<Vec<String>, String> {
    if !pointer.starts_with('/') {
        return Err(format!("pointer `{pointer}` must start with `/`"));
    }
    Ok(pointer[1..].split('/').map(|s| s.to_string()).collect())
}

/// Convert a serde_json::Value to a toml::Value. Carries the obvious
/// scalars + arrays + objects across; rejects nulls (TOML has no
/// null) with a clear error.
fn json_to_toml(v: &serde_json::Value) -> Result<toml::Value, String> {
    Ok(match v {
        serde_json::Value::Null => {
            return Err("TOML has no null type — pick a different sentinel".into());
        }
        serde_json::Value::Bool(b) => toml::Value::Boolean(*b),
        serde_json::Value::Number(n) => {
            if let Some(i) = n.as_i64() {
                toml::Value::Integer(i)
            } else if let Some(f) = n.as_f64() {
                toml::Value::Float(f)
            } else {
                return Err(format!("number `{n}` doesn't fit i64 or f64"));
            }
        }
        serde_json::Value::String(s) => toml::Value::String(s.clone()),
        serde_json::Value::Array(arr) => {
            let mut out = Vec::with_capacity(arr.len());
            for v in arr {
                out.push(json_to_toml(v)?);
            }
            toml::Value::Array(out)
        }
        serde_json::Value::Object(map) => {
            let mut t = toml::value::Table::new();
            for (k, v) in map {
                t.insert(k.clone(), json_to_toml(v)?);
            }
            toml::Value::Table(t)
        }
    })
}

fn value_kind(v: &toml::Value) -> &'static str {
    match v {
        toml::Value::String(_) => "string",
        toml::Value::Integer(_) => "integer",
        toml::Value::Float(_) => "float",
        toml::Value::Boolean(_) => "boolean",
        toml::Value::Datetime(_) => "datetime",
        toml::Value::Array(_) => "array",
        toml::Value::Table(_) => "table",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn two_param_sweep() -> SweepConfig {
        SweepConfig {
            optimizer: OptimizerKind::Grid,
            parameters: vec![
                ParameterSpec {
                    path: "/boundaries/inlet/velocity/0".into(),
                    values: vec![json!(10.0), json!(20.0), json!(30.0)],
                },
                ParameterSpec {
                    path: "/flow/turbulence".into(),
                    values: vec![json!("kEpsilon"), json!("kOmegaSST")],
                },
            ],
            objective: None,
            latin_hypercube: None,
            gradient_descent: None,
        }
    }

    #[test]
    fn grid_emits_full_cartesian_product() {
        let sweep = two_param_sweep();
        let plan = GridOptimizer::new().plan(&sweep).expect("plan");
        // 3 × 2 = 6 derived runs.
        assert_eq!(plan.len(), 6);
    }

    #[test]
    fn grid_substitution_carries_both_parameter_values() {
        let sweep = two_param_sweep();
        let plan = GridOptimizer::new().plan(&sweep).expect("plan");
        // Last parameter varies fastest: first two runs share the
        // same velocity (10) with two different turbulence models.
        let first = &plan[0];
        let second = &plan[1];
        assert_eq!(
            first.substitutions["/boundaries/inlet/velocity/0"],
            json!(10.0)
        );
        assert_eq!(first.substitutions["/flow/turbulence"], json!("kEpsilon"));
        assert_eq!(
            second.substitutions["/boundaries/inlet/velocity/0"],
            json!(10.0)
        );
        assert_eq!(second.substitutions["/flow/turbulence"], json!("kOmegaSST"));
        // Third run rolls over to the next velocity.
        let third = &plan[2];
        assert_eq!(
            third.substitutions["/boundaries/inlet/velocity/0"],
            json!(20.0)
        );
        assert_eq!(third.substitutions["/flow/turbulence"], json!("kEpsilon"));
    }

    #[test]
    fn grid_sweep_step_is_immediately_done() {
        let mut opt = GridOptimizer::new();
        match opt.step(&[]) {
            OptimizerStep::Done => {}
            other => panic!("expected Done, got {other:?}"),
        }
    }

    #[test]
    fn empty_parameter_list_errors_cleanly() {
        let sweep = SweepConfig {
            optimizer: OptimizerKind::Grid,
            parameters: vec![],
            objective: None,
            latin_hypercube: None,
            gradient_descent: None,
        };
        let err = GridOptimizer::new().plan(&sweep).unwrap_err();
        assert!(matches!(err, OptimizerError::InvalidConfig(_)));
    }

    #[test]
    fn empty_values_list_errors_cleanly() {
        let sweep = SweepConfig {
            optimizer: OptimizerKind::Grid,
            parameters: vec![ParameterSpec {
                path: "/x".into(),
                values: vec![],
            }],
            objective: None,
            latin_hypercube: None,
            gradient_descent: None,
        };
        let err = GridOptimizer::new().plan(&sweep).unwrap_err();
        match err {
            OptimizerError::InvalidConfig(reason) => {
                assert!(reason.contains("empty values list"));
            }
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn make_optimizer_returns_grid_for_grid_kind() {
        let opt = make_optimizer(OptimizerKind::Grid).expect("grid");
        assert_eq!(opt.id(), "grid");
    }

    #[test]
    fn make_optimizer_now_returns_gradient_descent_for_gd_kind() {
        // GradientDescentOptimizer landed alongside the other two —
        // every OptimizerKind now resolves.
        let opt = make_optimizer(OptimizerKind::GradientDescent).expect("ok");
        assert_eq!(opt.id(), "gradient-descent");
    }

    #[test]
    fn materialise_substitutes_into_table_value() {
        let base = r#"
[flow]
turbulence = "kEpsilon"
"#;
        let derived = DerivedCase {
            id: "sweep-0".into(),
            substitutions: {
                let mut m = BTreeMap::new();
                m.insert("/flow/turbulence".into(), json!("kOmegaSST"));
                m
            },
        };
        let out = materialise_case(base, &derived).expect("materialise");
        assert!(out.contains("turbulence = \"kOmegaSST\""), "got:\n{out}");
        // The original value is gone.
        assert!(!out.contains("kEpsilon"));
    }

    #[test]
    fn materialise_substitutes_into_array_index() {
        let base = r#"
[boundaries.inlet]
velocity = [10.0, 0.0, 0.0]
"#;
        let derived = DerivedCase {
            id: "sweep-0".into(),
            substitutions: {
                let mut m = BTreeMap::new();
                m.insert("/boundaries/inlet/velocity/0".into(), json!(50.0));
                m
            },
        };
        let out = materialise_case(base, &derived).expect("materialise");
        // The X-component changes from 10 to 50; Y/Z stay zero.
        assert!(out.contains("50"), "got:\n{out}");
        // We don't assert on the exact array spelling because TOML's
        // pretty printer chooses single- or multi-line layout based
        // on length.
    }

    #[test]
    fn materialise_rejects_missing_key() {
        let base = "[a]\nb = 1\n";
        let derived = DerivedCase {
            id: "sweep-0".into(),
            substitutions: {
                let mut m = BTreeMap::new();
                m.insert("/a/c".into(), json!(2));
                m
            },
        };
        // Setting a/c just inserts into table `a` — that should
        // succeed (TOML inserts are creative). Test the OTHER
        // failure mode: missing intermediate.
        let _ = materialise_case(base, &derived).expect("simple insert is fine");

        let derived = DerivedCase {
            id: "sweep-0".into(),
            substitutions: {
                let mut m = BTreeMap::new();
                m.insert("/missing/intermediate/key".into(), json!(2));
                m
            },
        };
        let err = materialise_case(base, &derived).unwrap_err();
        match err {
            OptimizerError::InvalidConfig(msg) => {
                assert!(msg.contains("missing key"), "got: {msg}");
            }
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn materialise_rejects_array_index_out_of_range() {
        let base = "v = [1.0, 2.0]\n";
        let derived = DerivedCase {
            id: "sweep-0".into(),
            substitutions: {
                let mut m = BTreeMap::new();
                m.insert("/v/5".into(), json!(99.0));
                m
            },
        };
        let err = materialise_case(base, &derived).unwrap_err();
        match err {
            OptimizerError::InvalidConfig(msg) => {
                assert!(msg.contains("out of range"), "got: {msg}");
            }
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn materialise_rejects_pointer_without_leading_slash() {
        let base = "[a]\nb = 1\n";
        let derived = DerivedCase {
            id: "sweep-0".into(),
            substitutions: {
                let mut m = BTreeMap::new();
                m.insert("a/b".into(), json!(99));
                m
            },
        };
        let err = materialise_case(base, &derived).unwrap_err();
        match err {
            OptimizerError::InvalidConfig(msg) => {
                assert!(msg.contains("must start with `/`"), "got: {msg}");
            }
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn materialise_handles_string_array_substitution() {
        // Parameter sweeping over discrete choices like turbulence
        // model names.
        let base = r#"
[flow]
turbulence = "laminar"
"#;
        let derived = DerivedCase {
            id: "sweep-0".into(),
            substitutions: {
                let mut m = BTreeMap::new();
                m.insert("/flow/turbulence".into(), json!("SpalartAllmaras"));
                m
            },
        };
        let out = materialise_case(base, &derived).expect("materialise");
        assert!(out.contains("SpalartAllmaras"));
    }

    #[test]
    fn materialise_stamps_a_sweep_derived_block_with_short_keys() {
        // The derived block lets the dataset assembler recover the
        // per-sample input vector without having to re-run the
        // optimizer.
        let base = r#"
[boundaries.inlet]
velocity = [10.0, 0.0, 0.0]

[flow]
turbulence = "kEpsilon"
"#;
        let derived = DerivedCase {
            id: "sweep-0".into(),
            substitutions: {
                let mut m = BTreeMap::new();
                m.insert("/boundaries/inlet/velocity/0".into(), json!(50.0));
                m.insert("/flow/turbulence".into(), json!("kOmegaSST"));
                m
            },
        };
        let out = materialise_case(base, &derived).expect("materialise");
        // [sweep.derived] picks up the numeric input only; the
        // string substitution lands in [sweep.derived_paths] as
        // free-form metadata.
        let value: toml::Value = toml::from_str(&out).expect("parse");
        let derived_block = value
            .get("sweep")
            .and_then(|s| s.get("derived"))
            .and_then(|d| d.as_table())
            .expect("[sweep.derived] block missing");
        assert!(derived_block.contains_key("velocity_0"));
        assert_eq!(
            derived_block.get("velocity_0").and_then(|v| v.as_float()),
            Some(50.0)
        );
        // String substitution is NOT in the numeric block.
        assert!(!derived_block.contains_key("turbulence"));
        // …but it IS in the paths sidecar.
        let paths_block = value
            .get("sweep")
            .and_then(|s| s.get("derived_paths"))
            .and_then(|d| d.as_table())
            .expect("[sweep.derived_paths] block missing");
        assert_eq!(
            paths_block.get("turbulence").and_then(|v| v.as_str()),
            Some("/flow/turbulence")
        );
    }

    #[test]
    fn pointer_short_name_handles_array_index_segments() {
        assert_eq!(
            super::pointer_short_name("/boundaries/inlet/velocity/0"),
            "velocity_0"
        );
        assert_eq!(super::pointer_short_name("/flow/turbulence"), "turbulence");
        assert_eq!(super::pointer_short_name("/x"), "x");
    }

    // -----------------------------------------------------------------
    // Latin Hypercube optimizer
    // -----------------------------------------------------------------

    fn lhs_sweep(n: usize, seed: u64) -> SweepConfig {
        SweepConfig {
            optimizer: OptimizerKind::LatinHypercube,
            parameters: vec![
                ParameterSpec {
                    path: "/aoa".into(),
                    values: vec![json!(-5.0), json!(15.0)],
                },
                ParameterSpec {
                    path: "/re".into(),
                    values: vec![json!(1e6), json!(1e7)],
                },
            ],
            objective: None,
            latin_hypercube: Some(LatinHypercubeConfig { n_samples: n, seed }),
            gradient_descent: None,
        }
    }

    #[test]
    fn latin_hypercube_emits_n_samples_distinct_runs() {
        let plan = LatinHypercubeOptimizer::new()
            .plan(&lhs_sweep(8, 42))
            .expect("plan");
        assert_eq!(plan.len(), 8);
        // Each sample id is unique.
        let mut ids: Vec<&str> = plan.iter().map(|p| p.id.as_str()).collect();
        ids.sort();
        ids.dedup();
        assert_eq!(ids.len(), 8);
    }

    #[test]
    fn latin_hypercube_each_parameter_covers_its_full_range_exactly_once() {
        // The defining property: every stratum is hit exactly once
        // per parameter. With 4 samples and aoa range [-5, 15], the
        // four centres are -2.5, 2.5, 7.5, 12.5. They appear in some
        // order across the 4 derived cases.
        let plan = LatinHypercubeOptimizer::new()
            .plan(&lhs_sweep(4, 0))
            .expect("plan");
        let mut aoa_values: Vec<f64> = plan
            .iter()
            .map(|c| {
                c.substitutions
                    .get("/aoa")
                    .and_then(|v| v.as_f64())
                    .unwrap()
            })
            .collect();
        aoa_values.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((aoa_values[0] - (-2.5)).abs() < 1e-9, "{aoa_values:?}");
        assert!((aoa_values[1] - 2.5).abs() < 1e-9, "{aoa_values:?}");
        assert!((aoa_values[2] - 7.5).abs() < 1e-9, "{aoa_values:?}");
        assert!((aoa_values[3] - 12.5).abs() < 1e-9, "{aoa_values:?}");
    }

    #[test]
    fn latin_hypercube_is_deterministic_for_a_given_seed() {
        let p1 = LatinHypercubeOptimizer::new()
            .plan(&lhs_sweep(8, 42))
            .expect("plan");
        let p2 = LatinHypercubeOptimizer::new()
            .plan(&lhs_sweep(8, 42))
            .expect("plan");
        let v1: Vec<(String, f64)> = p1
            .iter()
            .map(|c| {
                (
                    c.id.clone(),
                    c.substitutions
                        .get("/aoa")
                        .and_then(|v| v.as_f64())
                        .unwrap(),
                )
            })
            .collect();
        let v2: Vec<(String, f64)> = p2
            .iter()
            .map(|c| {
                (
                    c.id.clone(),
                    c.substitutions
                        .get("/aoa")
                        .and_then(|v| v.as_f64())
                        .unwrap(),
                )
            })
            .collect();
        assert_eq!(v1, v2);
    }

    #[test]
    fn latin_hypercube_different_seeds_produce_different_orderings() {
        let p1 = LatinHypercubeOptimizer::new()
            .plan(&lhs_sweep(8, 1))
            .expect("plan");
        let p2 = LatinHypercubeOptimizer::new()
            .plan(&lhs_sweep(8, 2))
            .expect("plan");
        let v1: Vec<f64> = p1
            .iter()
            .map(|c| {
                c.substitutions
                    .get("/aoa")
                    .and_then(|v| v.as_f64())
                    .unwrap()
            })
            .collect();
        let v2: Vec<f64> = p2
            .iter()
            .map(|c| {
                c.substitutions
                    .get("/aoa")
                    .and_then(|v| v.as_f64())
                    .unwrap()
            })
            .collect();
        assert_ne!(v1, v2);
    }

    #[test]
    fn latin_hypercube_rejects_missing_config_block() {
        let mut sweep = lhs_sweep(4, 0);
        sweep.latin_hypercube = None;
        let err = LatinHypercubeOptimizer::new()
            .plan(&sweep)
            .expect_err("fail");
        assert!(matches!(err, OptimizerError::InvalidConfig(_)));
    }

    #[test]
    fn latin_hypercube_rejects_parameters_without_two_bounds() {
        let mut sweep = lhs_sweep(4, 0);
        sweep.parameters[0].values = vec![json!(1.0)]; // only one value
        let err = LatinHypercubeOptimizer::new()
            .plan(&sweep)
            .expect_err("fail");
        let msg = format!("{err}");
        assert!(msg.contains("min, max"), "got {msg}");
    }

    #[test]
    fn latin_hypercube_rejects_inverted_bounds() {
        let mut sweep = lhs_sweep(4, 0);
        sweep.parameters[0].values = vec![json!(15.0), json!(-5.0)]; // max < min
        let err = LatinHypercubeOptimizer::new()
            .plan(&sweep)
            .expect_err("fail");
        let msg = format!("{err}");
        assert!(msg.contains("must be >"), "got {msg}");
    }

    #[test]
    fn make_optimizer_returns_lhs_for_kind_latin_hypercube() {
        let opt = make_optimizer(OptimizerKind::LatinHypercube).expect("ok");
        assert_eq!(opt.id(), "latin-hypercube");
    }

    /// Round-14 M6 RED→GREEN (round-4 sister gap): an LHS config
    /// with a pathologically large `n_samples` must be rejected
    /// before the strata Vec is allocated. Pre-fix
    /// `n_samples = 10_000_000_000` flowed straight into the strata
    /// builder and OOMed the host.
    #[test]
    fn latin_hypercube_rejects_oversized_n_samples() {
        // Build a sweep with one parameter and 10 billion samples.
        let mut sweep = lhs_sweep(0, 0);
        sweep.latin_hypercube = Some(LatinHypercubeConfig {
            n_samples: 10_000_000_000,
            seed: 0,
        });
        let err = LatinHypercubeOptimizer::new()
            .plan(&sweep)
            .expect_err("oversize n_samples must be rejected");
        match err {
            OptimizerError::TooManySamples {
                optimizer,
                requested,
                cap,
            } => {
                assert_eq!(optimizer, "latin-hypercube");
                assert_eq!(requested, 10_000_000_000);
                assert_eq!(cap, MAX_LHS_SAMPLES);
            }
            other => panic!("expected TooManySamples, got {other:?}"),
        }
    }

    // -----------------------------------------------------------------
    // Gradient descent optimizer
    // -----------------------------------------------------------------

    fn gd_sweep() -> SweepConfig {
        let mut initial = BTreeMap::new();
        initial.insert("/x".into(), 0.0);
        initial.insert("/y".into(), 0.0);
        SweepConfig {
            optimizer: OptimizerKind::GradientDescent,
            parameters: vec![
                ParameterSpec {
                    path: "/x".into(),
                    values: vec![json!(-10.0), json!(10.0)],
                },
                ParameterSpec {
                    path: "/y".into(),
                    values: vec![json!(-10.0), json!(10.0)],
                },
            ],
            objective: Some(ObjectiveSpec {
                metric: "f".into(),
                direction: ObjectiveDirection::Minimize,
            }),
            latin_hypercube: None,
            gradient_descent: Some(GradientDescentConfig {
                initial,
                epsilon: 0.1,
                learning_rate: 0.5,
                max_iterations: 3,
            }),
        }
    }

    #[test]
    fn gradient_descent_plan_emits_one_plus_two_d_runs() {
        // 2 params -> 1 center + 4 perturbations = 5 derived runs.
        let mut opt = GradientDescentOptimizer::new();
        let plan = opt.plan(&gd_sweep()).expect("plan");
        assert_eq!(plan.len(), 5);
        // Center first.
        assert_eq!(plan[0].id, "gd-iter00-center");
        // Then plus/minus pairs in declaration order.
        assert_eq!(plan[1].id, "gd-iter00-p00-plus");
        assert_eq!(plan[2].id, "gd-iter00-p00-minus");
        assert_eq!(plan[3].id, "gd-iter00-p01-plus");
        assert_eq!(plan[4].id, "gd-iter00-p01-minus");
    }

    #[test]
    fn gradient_descent_step_advances_the_point_downhill() {
        // Synthetic objective: f(x, y) = x^2 + y^2 with global min
        // at (0, 0). Starting at (3, 4), one descent step should
        // bring us closer to the origin.
        let mut sweep = gd_sweep();
        sweep
            .gradient_descent
            .as_mut()
            .unwrap()
            .initial
            .insert("/x".into(), 3.0);
        sweep
            .gradient_descent
            .as_mut()
            .unwrap()
            .initial
            .insert("/y".into(), 4.0);

        let mut opt = GradientDescentOptimizer::new();
        let plan = opt.plan(&sweep).expect("plan");
        // Build CompletedRuns with f = x^2 + y^2 evaluated at each
        // derived case's substitutions.
        let completed: Vec<CompletedRun> = plan
            .iter()
            .map(|d| {
                let x = d.substitutions["/x"].as_f64().unwrap();
                let y = d.substitutions["/y"].as_f64().unwrap();
                CompletedRun {
                    id: d.id.clone(),
                    objective_value: Some(x * x + y * y),
                    succeeded: true,
                }
            })
            .collect();

        let next = opt.step(&completed);
        let new_batch = match next {
            OptimizerStep::More(b) => b,
            OptimizerStep::Done => panic!("expected More, got Done"),
        };
        let center = new_batch
            .iter()
            .find(|d| d.id == "gd-iter01-center")
            .expect("center of iter 1");
        let new_x = center.substitutions["/x"].as_f64().unwrap();
        let new_y = center.substitutions["/y"].as_f64().unwrap();
        // grad(f) at (3, 4) = (6, 8). Step = -lr * grad = (-3, -4).
        // So new point = (0, 0) — exact for this convex quadratic.
        assert!((new_x - 0.0).abs() < 1e-9, "x: {new_x}");
        assert!((new_y - 0.0).abs() < 1e-9, "y: {new_y}");
    }

    #[test]
    fn gradient_descent_stops_at_max_iterations() {
        let mut opt = GradientDescentOptimizer::new();
        let sweep = gd_sweep(); // max_iterations = 3
        let mut current_batch = opt.plan(&sweep).expect("plan");
        // 3 iterations means: plan() = iter 0, then step three times
        // before Done. Hand the optimizer a non-flat objective so
        // it doesn't trip the gradient-norm early-stop.
        for _ in 0..3 {
            let completed: Vec<CompletedRun> = current_batch
                .iter()
                .map(|d| {
                    let x = d.substitutions["/x"].as_f64().unwrap();
                    let y = d.substitutions["/y"].as_f64().unwrap();
                    // Linear objective so the gradient never vanishes.
                    CompletedRun {
                        id: d.id.clone(),
                        objective_value: Some(x + 2.0 * y),
                        succeeded: true,
                    }
                })
                .collect();
            match opt.step(&completed) {
                OptimizerStep::More(b) => current_batch = b,
                OptimizerStep::Done => return, // Done before max — fine.
            }
        }
        // After max_iterations steps, a further step() must return Done.
        let completed: Vec<CompletedRun> = current_batch
            .iter()
            .map(|d| CompletedRun {
                id: d.id.clone(),
                objective_value: Some(0.0),
                succeeded: true,
            })
            .collect();
        assert!(matches!(opt.step(&completed), OptimizerStep::Done));
    }

    #[test]
    fn gradient_descent_clamps_step_to_box_bounds() {
        // Param bounds are [-10, 10]. With a ridiculous learning
        // rate (1e6) the unbounded step would shoot past 10.
        let mut sweep = gd_sweep();
        sweep.gradient_descent.as_mut().unwrap().learning_rate = 1e6;
        sweep
            .gradient_descent
            .as_mut()
            .unwrap()
            .initial
            .insert("/x".into(), 5.0);
        sweep
            .gradient_descent
            .as_mut()
            .unwrap()
            .initial
            .insert("/y".into(), 0.0);

        let mut opt = GradientDescentOptimizer::new();
        let plan = opt.plan(&sweep).expect("plan");
        // grad(f = -x) = -1, so descent direction (Minimize) takes a
        // huge positive step in x. Clamp must catch it at the upper
        // bound of 10.
        let completed: Vec<CompletedRun> = plan
            .iter()
            .map(|d| {
                let x = d.substitutions["/x"].as_f64().unwrap();
                CompletedRun {
                    id: d.id.clone(),
                    objective_value: Some(-x),
                    succeeded: true,
                }
            })
            .collect();
        let OptimizerStep::More(next) = opt.step(&completed) else {
            panic!("expected More");
        };
        let center = next
            .iter()
            .find(|d| d.id == "gd-iter01-center")
            .expect("iter 1 center");
        assert_eq!(center.substitutions["/x"].as_f64().unwrap(), 10.0);
    }

    #[test]
    fn gradient_descent_rejects_initial_outside_bounds() {
        let mut sweep = gd_sweep();
        sweep
            .gradient_descent
            .as_mut()
            .unwrap()
            .initial
            .insert("/x".into(), 999.0); // way outside [-10, 10]
        let mut opt = GradientDescentOptimizer::new();
        let err = opt.plan(&sweep).expect_err("must reject");
        let msg = format!("{err}");
        assert!(msg.contains("outside bounds"), "got {msg}");
    }

    #[test]
    fn gradient_descent_rejects_missing_initial_for_a_parameter() {
        let mut sweep = gd_sweep();
        sweep
            .gradient_descent
            .as_mut()
            .unwrap()
            .initial
            .remove("/y");
        let mut opt = GradientDescentOptimizer::new();
        let err = opt.plan(&sweep).expect_err("must reject");
        let msg = format!("{err}");
        assert!(msg.contains("missing entry for `/y`"), "got {msg}");
    }

    #[test]
    fn gradient_descent_requires_objective() {
        let mut sweep = gd_sweep();
        sweep.objective = None;
        let mut opt = GradientDescentOptimizer::new();
        let err = opt.plan(&sweep).expect_err("must reject");
        let msg = format!("{err}");
        assert!(msg.contains("objective"), "got {msg}");
    }

    #[test]
    fn gradient_descent_maximize_walks_uphill_not_downhill() {
        // Direction sign should flip for Maximize.
        let mut sweep = gd_sweep();
        sweep.objective = Some(ObjectiveSpec {
            metric: "f".into(),
            direction: ObjectiveDirection::Maximize,
        });
        sweep
            .gradient_descent
            .as_mut()
            .unwrap()
            .initial
            .insert("/x".into(), 1.0);
        sweep
            .gradient_descent
            .as_mut()
            .unwrap()
            .initial
            .insert("/y".into(), 0.0);

        let mut opt = GradientDescentOptimizer::new();
        let plan = opt.plan(&sweep).expect("plan");
        // Objective f(x, y) = x. grad = (1, 0). Maximize -> step
        // in +grad direction -> x grows.
        let completed: Vec<CompletedRun> = plan
            .iter()
            .map(|d| {
                let x = d.substitutions["/x"].as_f64().unwrap();
                CompletedRun {
                    id: d.id.clone(),
                    objective_value: Some(x),
                    succeeded: true,
                }
            })
            .collect();
        let OptimizerStep::More(next) = opt.step(&completed) else {
            panic!("expected More");
        };
        let center = next.iter().find(|d| d.id == "gd-iter01-center").unwrap();
        let new_x = center.substitutions["/x"].as_f64().unwrap();
        assert!(new_x > 1.0, "Maximize must increase x; got {new_x}");
    }

    #[test]
    fn gradient_descent_trace_empty_before_plan() {
        let opt = GradientDescentOptimizer::new();
        assert!(opt.trace().is_empty());
    }

    #[test]
    fn gradient_descent_trace_seeded_after_plan() {
        let mut opt = GradientDescentOptimizer::new();
        let _ = opt.plan(&gd_sweep()).expect("plan");
        // After plan() there's one seed entry — iteration 0, the
        // initial point, no objective / gradient yet.
        let trace = opt.trace();
        assert_eq!(trace.len(), 1);
        assert_eq!(trace[0].iteration, 0);
        assert_eq!(trace[0].point, vec![0.0, 0.0]);
        assert!(trace[0].objective_value.is_none());
        assert!(trace[0].gradient_norm.is_none());
    }

    #[test]
    fn gradient_descent_trace_records_per_iteration_data() {
        // Same fixture as the canonical "advances downhill" test.
        // Three step calls -> trace has 4 entries (iter 0 through 3),
        // with iter 0..=2 backfilled with center + gradient_norm.
        let mut sweep = gd_sweep();
        sweep
            .gradient_descent
            .as_mut()
            .unwrap()
            .initial
            .insert("/x".into(), 3.0);
        sweep
            .gradient_descent
            .as_mut()
            .unwrap()
            .initial
            .insert("/y".into(), 4.0);

        let mut opt = GradientDescentOptimizer::new();
        let mut batch = opt.plan(&sweep).expect("plan");
        for _ in 0..3 {
            let completed: Vec<CompletedRun> = batch
                .iter()
                .map(|d| {
                    let x = d.substitutions["/x"].as_f64().unwrap();
                    let y = d.substitutions["/y"].as_f64().unwrap();
                    CompletedRun {
                        id: d.id.clone(),
                        // x^2 + y^2; gradient is (2x, 2y).
                        objective_value: Some(x * x + y * y),
                        succeeded: true,
                    }
                })
                .collect();
            match opt.step(&completed) {
                OptimizerStep::More(b) => batch = b,
                OptimizerStep::Done => break,
            }
        }
        let trace = opt.trace();
        // After plan() + 3 steps: iter 0..=N entries; the most-
        // recent entry's objective hasn't been recorded yet (next
        // step would do it). The seed iteration's center value
        // gets filled in on the FIRST step.
        assert!(trace.len() >= 2);
        assert!(
            trace[0].objective_value.is_some(),
            "iter-0 trace should have objective backfilled after first step"
        );
        assert!(
            trace[0].gradient_norm.is_some(),
            "iter-0 trace should have gradient norm backfilled"
        );
        // iter 0 starts at (3, 4); its objective is 25.
        assert!((trace[0].objective_value.unwrap() - 25.0).abs() < 1e-9);
    }

    #[test]
    fn gradient_descent_trace_carries_gradient_norm_for_diagnostics() {
        let mut sweep = gd_sweep();
        sweep
            .gradient_descent
            .as_mut()
            .unwrap()
            .initial
            .insert("/x".into(), 3.0);
        sweep
            .gradient_descent
            .as_mut()
            .unwrap()
            .initial
            .insert("/y".into(), 4.0);
        let mut opt = GradientDescentOptimizer::new();
        let plan = opt.plan(&sweep).expect("plan");
        let completed: Vec<CompletedRun> = plan
            .iter()
            .map(|d| {
                let x = d.substitutions["/x"].as_f64().unwrap();
                let y = d.substitutions["/y"].as_f64().unwrap();
                CompletedRun {
                    id: d.id.clone(),
                    objective_value: Some(x * x + y * y),
                    succeeded: true,
                }
            })
            .collect();
        let _ = opt.step(&completed);
        let trace = opt.trace();
        // Gradient at (3, 4) for f = x^2 + y^2 is (6, 8); ||grad|| = 10.
        let g = trace[0].gradient_norm.expect("gradient norm filled");
        assert!((g - 10.0).abs() < 1e-9, "got grad norm {g}");
    }

    #[test]
    fn ids_are_zero_padded_to_total_width() {
        // A 100-run sweep should give ids sweep-000 .. sweep-099 so
        // they sort lexicographically in a directory listing.
        let sweep = SweepConfig {
            optimizer: OptimizerKind::Grid,
            parameters: vec![ParameterSpec {
                path: "/x".into(),
                values: (0..100).map(|i| json!(i as f64)).collect(),
            }],
            objective: None,
            latin_hypercube: None,
            gradient_descent: None,
        };
        let plan = GridOptimizer::new().plan(&sweep).expect("plan");
        assert_eq!(plan[0].id, "sweep-000");
        assert_eq!(plan[99].id, "sweep-099");
    }
}
