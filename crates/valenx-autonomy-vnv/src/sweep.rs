//! [`run_suite`] / [`SweepResult`] — run a whole suite and aggregate.
//!
//! A *sweep* runs every [`Scenario`] in a [`ScenarioSuite`], evaluates the same
//! [`RequirementSet`] against each resulting [`crate::Trace`], and aggregates the
//! verdicts into a [`SweepResult`]: the overall pass-rate, the worst-case
//! (minimum) margin per requirement and overall, and — when a [`ParamGrid`] is
//! supplied — parameter-space coverage. It also exposes builders that *generate*
//! a suite, either over the full grid ([`grid_suite`]) or by seeded Monte-Carlo
//! sampling of axis ranges ([`monte_carlo_suite`], reusing the deterministic
//! `valenx-sensors` [`SplitMix64`] so a seed reproduces the same suite exactly).

use std::collections::BTreeMap;

use valenx_sensors::SplitMix64;

use crate::coverage::{parameter_coverage, ParamCoverage, ParamGrid};
use crate::error::VnvError;
use crate::report::{evaluate, VnvReport};
use crate::requirement::RequirementSet;
use crate::scenario::{Scenario, ScenarioSuite};
use crate::trace::run_scenario;

/// The aggregated result of running a whole suite against a requirement set.
#[derive(Debug, Clone, PartialEq)]
pub struct SweepResult {
    /// The suite name.
    pub suite: String,
    /// One report per scenario, in suite order.
    pub reports: Vec<VnvReport>,
    /// Number of scenarios run.
    pub total_runs: usize,
    /// Number of scenarios whose every requirement passed.
    pub passed_runs: usize,
    /// The worst (minimum) margin observed for each requirement label across all
    /// runs, in deterministic label order. The minimum is meaningful because all
    /// margins share the sign convention (negative = violated), so this is the
    /// nearest-miss / deepest-violation per requirement over the whole sweep.
    pub worst_margin_by_requirement: BTreeMap<String, f64>,
    /// Parameter-space coverage, if a [`ParamGrid`] was supplied to the sweep.
    pub parameter_coverage: Option<ParamCoverage>,
}

impl SweepResult {
    /// The pass-rate in `[0, 1]` (`passed / total`), or `0.0` for no runs.
    #[must_use]
    pub fn pass_rate(&self) -> f64 {
        if self.total_runs == 0 {
            0.0
        } else {
            self.passed_runs as f64 / self.total_runs as f64
        }
    }

    /// The single worst margin across every requirement and every run, or `None`
    /// if nothing was evaluated. This is the global nearest-miss / deepest
    /// violation for the whole sweep.
    #[must_use]
    pub fn worst_margin(&self) -> Option<f64> {
        self.worst_margin_by_requirement
            .values()
            .copied()
            .fold(None, |acc, m| Some(acc.map_or(m, |a: f64| a.min(m))))
    }

    /// Whether every run passed.
    #[must_use]
    pub fn all_passed(&self) -> bool {
        self.passed_runs == self.total_runs && self.total_runs > 0
    }
}

/// Run `suite` against `requirements`, aggregating into a [`SweepResult`].
///
/// Validates the suite and requirement set first (fail loud on an empty suite, a
/// bad scenario, or a bad requirement). Each scenario is run to a [`crate::Trace`] and
/// scored; the per-requirement worst margins and the pass-rate are accumulated.
/// If `grid` is `Some`, parameter-space coverage of the suite against it is also
/// computed and attached.
///
/// # Errors
/// Any [`ScenarioSuite::validate`], [`RequirementSet::validate`],
/// [`run_scenario`], or [`evaluate`] error — surfaced rather than swallowed, so
/// a misconfigured sweep does not silently report a bogus pass-rate.
pub fn run_suite(
    suite: &ScenarioSuite,
    requirements: &RequirementSet,
    grid: Option<&ParamGrid>,
) -> Result<SweepResult, VnvError> {
    suite.validate()?;
    requirements.validate()?;

    let mut reports = Vec::with_capacity(suite.len());
    let mut worst: BTreeMap<String, f64> = BTreeMap::new();
    let mut passed_runs = 0;

    for scenario in &suite.scenarios {
        let trace = run_scenario(scenario)?;
        let report = evaluate(requirements, &trace)?;
        if report.overall_pass {
            passed_runs += 1;
        }
        for o in &report.outcomes {
            worst
                .entry(o.label.clone())
                .and_modify(|m| *m = m.min(o.margin))
                .or_insert(o.margin);
        }
        reports.push(report);
    }

    let parameter_coverage = match grid {
        Some(g) => Some(parameter_coverage(suite, g)?),
        None => None,
    };

    Ok(SweepResult {
        suite: suite.name.clone(),
        total_runs: reports.len(),
        passed_runs,
        worst_margin_by_requirement: worst,
        parameter_coverage,
        reports,
    })
}

/// Build a [`ScenarioSuite`] spanning every cell of `grid`.
///
/// `build` is called once per grid cell with the cell's `{axis → value}` map (in
/// deterministic, axis-name-sorted order) and the cell index; it returns the
/// [`Scenario`] for that cell. The returned scenarios are expected to tag their
/// [`crate::Scenario::params`] from the cell (the helper does **not** auto-tag,
/// so the builder stays in control — but see [`grid_suite_auto`] for the common
/// case that just copies the cell into the params).
///
/// # Errors
/// [`ParamGrid::validate`] errors, or any error the `build` closure returns.
pub fn grid_suite<F>(
    name: impl Into<String>,
    grid: &ParamGrid,
    mut build: F,
) -> Result<ScenarioSuite, VnvError>
where
    F: FnMut(&BTreeMap<String, f64>, usize) -> Result<Scenario, VnvError>,
{
    grid.validate()?;
    let cells = grid.cells();
    let mut scenarios = Vec::with_capacity(cells.len());
    for (i, cell) in cells.iter().enumerate() {
        scenarios.push(build(cell, i)?);
    }
    Ok(ScenarioSuite::new(name, scenarios))
}

/// Like [`grid_suite`], but the cell's `{axis → value}` map is automatically
/// copied into each scenario's [`crate::Scenario::params`] *after* `build`
/// returns, guaranteeing the suite has full parameter coverage against `grid`.
/// The `build` closure only has to turn a cell into a scenario; it may set
/// additional params of its own (they are preserved; the grid axes overwrite any
/// same-named entry so coverage matching is exact).
///
/// # Errors
/// As [`grid_suite`].
pub fn grid_suite_auto<F>(
    name: impl Into<String>,
    grid: &ParamGrid,
    mut build: F,
) -> Result<ScenarioSuite, VnvError>
where
    F: FnMut(&BTreeMap<String, f64>, usize) -> Result<Scenario, VnvError>,
{
    grid.validate()?;
    let cells = grid.cells();
    let mut scenarios = Vec::with_capacity(cells.len());
    for (i, cell) in cells.iter().enumerate() {
        let mut scenario = build(cell, i)?;
        for (k, v) in cell {
            scenario.set_param(k.clone(), *v);
        }
        scenarios.push(scenario);
    }
    Ok(ScenarioSuite::new(name, scenarios))
}

/// One continuous axis to Monte-Carlo sample: a name and an inclusive
/// `[lo, hi]` range. Sampling draws uniformly in `[lo, hi]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SampleAxis {
    /// Inclusive lower bound.
    pub lo: f64,
    /// Inclusive upper bound (`>= lo`).
    pub hi: f64,
}

impl SampleAxis {
    /// Build a sample axis.
    ///
    /// # Errors
    /// [`VnvError::InvalidConfig`] if `lo`/`hi` are non-finite or `hi < lo`.
    pub fn new(lo: f64, hi: f64) -> Result<Self, VnvError> {
        if !(lo.is_finite() && hi.is_finite() && hi >= lo) {
            return Err(VnvError::InvalidConfig(format!(
                "sample axis must have finite lo ≤ hi, got [{lo}, {hi}]"
            )));
        }
        Ok(Self { lo, hi })
    }

    /// Map a uniform `u ∈ [0, 1)` to a value in `[lo, hi]`.
    fn sample(&self, u: f64) -> f64 {
        self.lo + (self.hi - self.lo) * u
    }
}

/// Build a [`ScenarioSuite`] of `count` scenarios by seeded Monte-Carlo
/// sampling of continuous parameter `axes`.
///
/// For each of the `count` samples, every axis is drawn uniformly in its range
/// from a single seeded [`SplitMix64`] (so the whole suite is reproducible from
/// `seed`), and `build` turns the `{axis → value}` sample into a [`Scenario`].
/// The sampled values are automatically copied into each scenario's
/// [`crate::Scenario::params`].
///
/// Note: continuous MC samples are (almost surely) distinct, so they do **not**
/// align to a discrete [`ParamGrid`]'s cells; use the grid builders for cell
/// coverage. MC is for stress-testing the *pass-rate* over a continuous
/// envelope.
///
/// # Errors
/// [`VnvError::InvalidConfig`] for a zero `count`; any axis error; or any error
/// the `build` closure returns.
pub fn monte_carlo_suite<F>(
    name: impl Into<String>,
    axes: &BTreeMap<String, SampleAxis>,
    count: usize,
    seed: u64,
    mut build: F,
) -> Result<ScenarioSuite, VnvError>
where
    F: FnMut(&BTreeMap<String, f64>, usize) -> Result<Scenario, VnvError>,
{
    if count == 0 {
        return Err(VnvError::InvalidConfig(
            "monte_carlo_suite count must be ≥ 1".into(),
        ));
    }
    for (axis_name, axis) in axes {
        if !(axis.lo.is_finite() && axis.hi.is_finite() && axis.hi >= axis.lo) {
            return Err(VnvError::InvalidConfig(format!(
                "sample axis '{axis_name}' must have finite lo ≤ hi, got [{}, {}]",
                axis.lo, axis.hi
            )));
        }
    }

    let mut rng = SplitMix64::new(seed);
    let mut scenarios = Vec::with_capacity(count);
    for i in 0..count {
        let mut sample: BTreeMap<String, f64> = BTreeMap::new();
        // Draw in axis-name order so the stream → sample mapping is
        // deterministic and stable.
        for (axis_name, axis) in axes {
            let u = rng.next_f64();
            sample.insert(axis_name.clone(), axis.sample(u));
        }
        let mut scenario = build(&sample, i)?;
        for (k, v) in &sample {
            scenario.set_param(k.clone(), *v);
        }
        scenarios.push(scenario);
    }
    Ok(ScenarioSuite::new(name, scenarios))
}
