//! Digital-engineering / MBSE spine: the **digital thread** that ties design
//! parameters → requirements → trade studies into one versionable record.
//!
//! This is the systems-engineering layer of [`crate`]: pure
//! systems-engineering software with **zero dual-use concern**. It composes
//! transparent, auditable book-keeping — it adds no predictive power of its own.
//! The numbers it checks come from valenx's physics solvers (or, in a trade
//! study, from any evaluator you hand it).
//!
//! ## What it provides
//!
//! - **Parameter & requirement model** — a [`Parameter`] is a named, united
//!   scalar; a [`Requirement`] is a named acceptance criterion over one *metric*
//!   with a [`Comparator`] (`Min` / `Max` / `Range`) and a target.
//!   [`Requirement::evaluate`] returns a [`Verdict`] with a **signed margin**:
//!   how far *inside* (positive) or *outside* (negative) the target a metric
//!   value lands.
//! - **Design point + compliance** — a [`DesignPoint`] bundles parameters with a
//!   map of computed metric values; [`DesignPoint::check`] produces a
//!   [`ComplianceReport`] (per-requirement verdicts, `overall_pass`, and the
//!   `worst_margin`).
//! - **Trade study / DOE driver** — [`TradeStudy`] runs a **full-factorial**
//!   sweep over discrete parameter levels, calling an evaluator at each grid
//!   point to get the metrics, then attaches compliance and computes the
//!   **Pareto-non-dominated front** over chosen [`Objective`]s.
//! - **Traceability** — [`coverage`] maps requirements ↔ the metrics that
//!   satisfy them and flags any metric left **unconstrained** by a requirement.
//!
//! ## Honest scope & cost
//!
//! - A full-factorial sweep evaluates **the product of the sweep sizes** of every
//!   parameter (a 3×4 sweep is 12 points; adding a third 5-level parameter makes
//!   it 60). This grows combinatorially — it is intended for small, deliberate
//!   design-of-experiments grids, not high-dimensional exploration.
//! - The [`pareto_front`] helper is an `O(n²)` pairwise dominance scan, fine for
//!   the small candidate sets a trade study produces.
//! - **Future hook (not a dependency now):** sampling-based DOE (Monte-Carlo /
//!   Latin-hypercube), surrogates and global sensitivity belong in the planned
//!   `valenx-uq` track and would plug in *behind* the [`Evaluator`] seam — a
//!   sampler choosing parameter vectors instead of the full grid. This module
//!   takes **no** dependency on `valenx-uq`; the seam is the only contract.
//!
//! ## Example
//!
//! ```
//! use valenx_orchestrator::digeng::{
//!     Comparator, DesignPoint, Parameter, Requirement,
//! };
//!
//! // One requirement: mass must stay at or below 5.0 kg.
//! let reqs = vec![Requirement::new(
//!     "R1",
//!     "mass budget",
//!     "mass_kg",
//!     Comparator::Max,
//!     5.0,
//! )];
//!
//! // A design point whose computed mass is 4.2 kg.
//! let dp = DesignPoint::new("baseline")
//!     .with_parameter(Parameter::new("wall_mm", 2.0, "mm"))
//!     .with_metric("mass_kg", 4.2);
//!
//! let report = dp.check(&reqs);
//! assert!(report.overall_pass);
//! assert!((report.worst_margin - 0.8).abs() < 1e-12); // 5.0 - 4.2, inside
//! ```

use std::collections::{BTreeMap, BTreeSet};

use serde::{Deserialize, Serialize};

/// A named, united scalar design parameter.
///
/// `unit` is carried for provenance/reporting only — this module does **not**
/// perform unit conversion; the evaluator and requirements are responsible for
/// using consistent units.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Parameter {
    /// Parameter name (e.g. `"wall_thickness"`).
    pub name: String,
    /// Parameter value, in `unit`.
    pub value: f64,
    /// Unit string for provenance (e.g. `"mm"`); not interpreted.
    pub unit: String,
}

impl Parameter {
    /// Construct a parameter from a name, value and unit.
    pub fn new(name: impl Into<String>, value: f64, unit: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            value,
            unit: unit.into(),
        }
    }
}

/// How a [`Requirement`] compares a metric value against its target.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum Comparator {
    /// The metric must be **at least** the target (`value >= target`). A
    /// floor; e.g. a minimum factor of safety.
    Min,
    /// The metric must be **at most** the target (`value <= target`). A
    /// ceiling; e.g. a mass or cost budget.
    Max,
    /// The metric must lie **within** an inclusive band `[lo, hi]`.
    Range {
        /// Inclusive lower bound.
        lo: f64,
        /// Inclusive upper bound.
        hi: f64,
    },
}

/// The outcome of checking one metric value against one [`Requirement`].
///
/// `margin` is **signed**: it is `>= 0` exactly when [`Verdict::pass`] is `true`
/// (the value satisfies the comparator), and its magnitude is the distance to
/// the nearest violated bound — i.e. how much head-room remains (when passing)
/// or how far the value must move to comply (when failing).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Verdict {
    /// Whether the metric value satisfies the comparator.
    pub pass: bool,
    /// Signed margin: `>= 0` when inside the target, `< 0` when outside.
    pub margin: f64,
}

/// A single acceptance criterion over one named metric.
///
/// The `metric` field is the *key* linking this requirement to a value in a
/// [`DesignPoint`]'s metric map (and to a trade study's evaluator output via
/// [`TradeStudy::metric_names`]).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Requirement {
    /// Stable requirement identifier (e.g. `"R-001"`).
    pub id: String,
    /// Human-readable requirement name.
    pub name: String,
    /// The metric key this requirement constrains.
    pub metric: String,
    /// How the metric is compared to `target`.
    pub comparator: Comparator,
    /// The target / bound value (the band bounds live in `comparator` for a
    /// [`Comparator::Range`]; `target` is then unused and conventionally `0.0`).
    pub target: f64,
}

impl Requirement {
    /// Construct a requirement.
    ///
    /// For a [`Comparator::Range`] the band lives in the comparator itself; pass
    /// any `target` (`0.0` by convention) — it is ignored for ranges.
    pub fn new(
        id: impl Into<String>,
        name: impl Into<String>,
        metric: impl Into<String>,
        comparator: Comparator,
        target: f64,
    ) -> Self {
        Self {
            id: id.into(),
            name: name.into(),
            metric: metric.into(),
            comparator,
            target,
        }
    }

    /// Evaluate `metric_value` against this requirement.
    ///
    /// Returns a [`Verdict`] whose `margin` is signed (positive = inside the
    /// target with that much head-room; negative = outside by that much):
    ///
    /// - [`Comparator::Min`]: `margin = value - target` (pass when `>= 0`).
    /// - [`Comparator::Max`]: `margin = target - value` (pass when `>= 0`).
    /// - [`Comparator::Range`]: `margin = min(value - lo, hi - value)` — the
    ///   distance to the nearer band edge (pass when `>= 0`, i.e. inside).
    pub fn evaluate(&self, metric_value: f64) -> Verdict {
        let margin = match self.comparator {
            Comparator::Min => metric_value - self.target,
            Comparator::Max => self.target - metric_value,
            Comparator::Range { lo, hi } => (metric_value - lo).min(hi - metric_value),
        };
        Verdict {
            pass: margin >= 0.0,
            margin,
        }
    }
}

/// A named design point: a set of [`Parameter`]s plus the metric values computed
/// for them.
///
/// Metrics are keyed by name so they line up with each [`Requirement::metric`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DesignPoint {
    /// Name / label for this design point.
    pub name: String,
    /// The parameters defining this point.
    pub parameters: Vec<Parameter>,
    /// Computed metric values, keyed by metric name.
    pub metrics: BTreeMap<String, f64>,
}

impl DesignPoint {
    /// Create an empty design point with the given name.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            parameters: Vec::new(),
            metrics: BTreeMap::new(),
        }
    }

    /// Builder: add a parameter.
    pub fn with_parameter(mut self, p: Parameter) -> Self {
        self.parameters.push(p);
        self
    }

    /// Builder: set a computed metric value.
    pub fn with_metric(mut self, name: impl Into<String>, value: f64) -> Self {
        self.metrics.insert(name.into(), value);
        self
    }

    /// Look up a metric value by name.
    pub fn metric(&self, name: &str) -> Option<f64> {
        self.metrics.get(name).copied()
    }

    /// Check this design point against `requirements`, producing a
    /// [`ComplianceReport`].
    ///
    /// Each requirement is evaluated against the matching metric. A requirement
    /// whose metric is **missing** from this point is recorded as a failing
    /// verdict with `margin = NEG_INFINITY` (an unmet requirement is never
    /// silently treated as a pass — mirroring the funnel's "absence is not
    /// evidence of safety" stance).
    ///
    /// With **no** requirements the report is vacuously passing
    /// (`overall_pass = true`) and its `worst_margin` is
    /// [`f64::INFINITY`] (no constraint binds).
    pub fn check(&self, requirements: &[Requirement]) -> ComplianceReport {
        let mut per_requirement = Vec::with_capacity(requirements.len());
        let mut overall_pass = true;
        let mut worst_margin = f64::INFINITY;

        for req in requirements {
            let verdict = match self.metric(&req.metric) {
                Some(v) => req.evaluate(v),
                None => Verdict {
                    pass: false,
                    margin: f64::NEG_INFINITY,
                },
            };
            overall_pass &= verdict.pass;
            if verdict.margin < worst_margin {
                worst_margin = verdict.margin;
            }
            per_requirement.push((req.id.clone(), verdict));
        }

        ComplianceReport {
            design_point: self.name.clone(),
            per_requirement,
            overall_pass,
            worst_margin,
        }
    }
}

/// The result of checking a [`DesignPoint`] against a set of [`Requirement`]s.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ComplianceReport {
    /// Name of the design point this report covers.
    pub design_point: String,
    /// One `(requirement id, verdict)` pair per requirement checked, in order.
    pub per_requirement: Vec<(String, Verdict)>,
    /// `true` iff **every** requirement passed.
    pub overall_pass: bool,
    /// The smallest (most binding) signed margin across all requirements;
    /// [`f64::INFINITY`] when there are no requirements.
    pub worst_margin: f64,
}

impl ComplianceReport {
    /// Render the report as a short human-readable block.
    pub fn render(&self) -> String {
        let mut s = String::new();
        s.push_str(&format!(
            "[compliance] {} — {}\n",
            self.design_point,
            if self.overall_pass { "PASS" } else { "FAIL" }
        ));
        for (id, v) in &self.per_requirement {
            s.push_str(&format!(
                "  {id}: {} (margin {:.4})\n",
                if v.pass { "pass" } else { "FAIL" },
                v.margin
            ));
        }
        s.push_str(&format!("  worst margin: {:.4}\n", self.worst_margin));
        s
    }
}

/// Whether a trade-study [`Objective`] is to be minimized or maximized when
/// computing the Pareto front.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Direction {
    /// Smaller is better (e.g. mass, cost).
    Minimize,
    /// Larger is better (e.g. payload, factor of safety).
    Maximize,
}

/// One objective for Pareto analysis: a metric and the direction that improves it.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Objective {
    /// The metric key to optimize (must be produced by the evaluator).
    pub metric: String,
    /// Whether to minimize or maximize the metric.
    pub direction: Direction,
}

impl Objective {
    /// An objective to minimize `metric`.
    pub fn minimize(metric: impl Into<String>) -> Self {
        Self {
            metric: metric.into(),
            direction: Direction::Minimize,
        }
    }

    /// An objective to maximize `metric`.
    pub fn maximize(metric: impl Into<String>) -> Self {
        Self {
            metric: metric.into(),
            direction: Direction::Maximize,
        }
    }
}

/// One discrete parameter sweep: a name, a unit, and the levels to enumerate.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ParameterSweep {
    /// Parameter name.
    pub name: String,
    /// Unit for provenance (carried onto each generated [`Parameter`]).
    pub unit: String,
    /// The discrete levels to sweep over.
    pub levels: Vec<f64>,
}

impl ParameterSweep {
    /// Construct a sweep from a name, unit and explicit levels.
    pub fn new(name: impl Into<String>, unit: impl Into<String>, levels: Vec<f64>) -> Self {
        Self {
            name: name.into(),
            unit: unit.into(),
            levels,
        }
    }
}

/// An evaluator maps a parameter vector (one value per swept parameter, in sweep
/// order) to a metric vector (one value per [`TradeStudy::metric_names`], in
/// order).
///
/// This is the seam a future sampling/UQ driver would sit behind: the same
/// `Fn` contract, fed sampled vectors instead of a full grid.
pub type Evaluator<'a> = dyn Fn(&[f64]) -> Vec<f64> + 'a;

/// A full-factorial trade study (design of experiments).
///
/// Holds the parameter sweeps, the ordered metric names the evaluator returns,
/// and the requirements to check each design point against. Run it with
/// [`TradeStudy::run`].
pub struct TradeStudy {
    /// The parameters to sweep (full-factorial across their levels).
    pub sweeps: Vec<ParameterSweep>,
    /// Names of the metrics the evaluator returns, **in evaluator output order**.
    pub metric_names: Vec<String>,
    /// Requirements every generated design point is checked against.
    pub requirements: Vec<Requirement>,
}

impl TradeStudy {
    /// Construct a trade study from sweeps, the evaluator's metric names (in
    /// output order) and the requirements to check.
    pub fn new(
        sweeps: Vec<ParameterSweep>,
        metric_names: Vec<String>,
        requirements: Vec<Requirement>,
    ) -> Self {
        Self {
            sweeps,
            metric_names,
            requirements,
        }
    }

    /// The number of design points a full-factorial run will produce: the
    /// product of every sweep's level count.
    ///
    /// An empty sweep set yields `0` (there is nothing to vary, so no grid). A
    /// sweep with zero levels collapses the whole product to `0`.
    pub fn point_count(&self) -> usize {
        if self.sweeps.is_empty() {
            return 0;
        }
        self.sweeps.iter().map(|s| s.levels.len()).product()
    }

    /// Run the full-factorial sweep, calling `evaluator` at every grid point.
    ///
    /// For each combination of one level per parameter (Cartesian product, in
    /// sweep order with the **last** parameter varying fastest), the evaluator is
    /// called with that parameter vector. Its returned metric vector is paired
    /// positionally with [`TradeStudy::metric_names`] to build the design point's
    /// metric map, the requirements are checked, and a [`TradeResult`] is
    /// recorded.
    ///
    /// Returns an empty [`TradeStudyOutcome`] (no panic) when there is nothing to
    /// sweep — no sweeps, or any sweep with zero levels.
    ///
    /// # Panics
    ///
    /// Never. A short metric vector from the evaluator is tolerated: only the
    /// metric names it covers are populated (a name with no corresponding output
    /// value is simply absent from that point's metric map, and any requirement
    /// on it then records the missing-metric failing verdict).
    pub fn run(&self, evaluator: &Evaluator<'_>) -> TradeStudyOutcome {
        let mut results = Vec::new();

        for combo in cartesian_product(&self.sweeps) {
            let metric_values = evaluator(&combo);

            let mut dp = DesignPoint::new(point_label(&self.sweeps, &combo));
            for (sweep, &value) in self.sweeps.iter().zip(combo.iter()) {
                dp = dp.with_parameter(Parameter::new(
                    sweep.name.clone(),
                    value,
                    sweep.unit.clone(),
                ));
            }
            for (name, &value) in self.metric_names.iter().zip(metric_values.iter()) {
                dp = dp.with_metric(name.clone(), value);
            }

            let compliance = dp.check(&self.requirements);
            results.push(TradeResult {
                design_point: dp,
                compliance,
            });
        }

        TradeStudyOutcome { results }
    }
}

/// One row of a [`TradeStudyOutcome`]: a generated design point and its
/// compliance against the study's requirements.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TradeResult {
    /// The generated design point (parameters + evaluated metrics).
    pub design_point: DesignPoint,
    /// Its compliance report against the study's requirements.
    pub compliance: ComplianceReport,
}

/// The full output of a [`TradeStudy::run`]: every design point with its
/// compliance.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TradeStudyOutcome {
    /// One result per full-factorial grid point.
    pub results: Vec<TradeResult>,
}

impl TradeStudyOutcome {
    /// The indices of the [`TradeResult`]s on the Pareto-non-dominated front for
    /// `objectives`.
    ///
    /// A result is on the front when **no other** result is at least as good in
    /// every objective and strictly better in at least one (standard Pareto
    /// dominance). A result missing any objective's metric cannot dominate and is
    /// itself never dominated on the missing axis — it is conservatively kept off
    /// no front purely by an absent value (see [`pareto_front`]).
    ///
    /// With no objectives, **all** indices are returned (nothing can dominate).
    pub fn pareto_indices(&self, objectives: &[Objective]) -> Vec<usize> {
        let points: Vec<&DesignPoint> = self.results.iter().map(|r| &r.design_point).collect();
        pareto_front(&points, objectives)
    }

    /// The [`TradeResult`]s on the Pareto-non-dominated front for `objectives`.
    pub fn pareto_front(&self, objectives: &[Objective]) -> Vec<&TradeResult> {
        self.pareto_indices(objectives)
            .into_iter()
            .map(|i| &self.results[i])
            .collect()
    }

    /// The results that satisfy **all** requirements (`overall_pass`).
    pub fn compliant(&self) -> Vec<&TradeResult> {
        self.results
            .iter()
            .filter(|r| r.compliance.overall_pass)
            .collect()
    }
}

/// Compute the Pareto-non-dominated front over `points` for `objectives`.
///
/// Returns the indices into `points` that are not dominated by any other point.
/// Point `a` *dominates* `b` when, across every objective, `a` is no worse than
/// `b` (in that objective's [`Direction`]) and strictly better in at least one.
///
/// A point missing an objective's metric is treated as having no value on that
/// axis: it can be neither better nor worse there, so it can never strictly
/// dominate (it lacks the "strictly better somewhere" it might have needed) and
/// is harder to dominate. This keeps incomplete points from silently swallowing
/// the front.
///
/// Complexity is `O(n² · m)` for `n` points and `m` objectives — intended for
/// the small candidate sets a trade study produces.
///
/// With no objectives every index is returned (the front is everything).
pub fn pareto_front(points: &[&DesignPoint], objectives: &[Objective]) -> Vec<usize> {
    if objectives.is_empty() {
        return (0..points.len()).collect();
    }
    (0..points.len())
        .filter(|&i| {
            !(0..points.len()).any(|j| j != i && dominates(points[j], points[i], objectives))
        })
        .collect()
}

/// Does `a` Pareto-dominate `b` over `objectives`?
fn dominates(a: &DesignPoint, b: &DesignPoint, objectives: &[Objective]) -> bool {
    let mut strictly_better_somewhere = false;
    for obj in objectives {
        let (av, bv) = match (a.metric(&obj.metric), b.metric(&obj.metric)) {
            (Some(av), Some(bv)) => (av, bv),
            // `a` lacks this metric: it cannot be "no worse" here, so it cannot
            // dominate `b` at all.
            (None, _) => return false,
            // `b` lacks this metric: `a` is not worse here, but gains no strict
            // edge from a missing value; treat this axis as a tie.
            (Some(_), None) => continue,
        };
        // Normalize to "smaller is better" so one comparison covers both.
        let (a_cost, b_cost) = match obj.direction {
            Direction::Minimize => (av, bv),
            Direction::Maximize => (-av, -bv),
        };
        if a_cost > b_cost {
            return false; // strictly worse in this objective -> cannot dominate
        }
        if a_cost < b_cost {
            strictly_better_somewhere = true;
        }
    }
    strictly_better_somewhere
}

/// A requirements-↔-metrics traceability report.
///
/// Built by [`coverage`]; answers "which requirements exercise which metrics,
/// and which produced metrics are left unconstrained?".
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CoverageReport {
    /// One `(requirement id, metric)` link per requirement, in input order.
    pub requirement_metrics: Vec<(String, String)>,
    /// Metric names that are produced but constrained by **no** requirement,
    /// sorted and de-duplicated.
    pub unconstrained_metrics: Vec<String>,
    /// Requirement ids whose metric is **not** among the produced metrics
    /// (a requirement with no data to check it), sorted and de-duplicated.
    pub uncovered_requirements: Vec<String>,
}

impl CoverageReport {
    /// `true` when every produced metric is constrained and every requirement
    /// has a metric to check it against.
    pub fn is_fully_covered(&self) -> bool {
        self.unconstrained_metrics.is_empty() && self.uncovered_requirements.is_empty()
    }

    /// Render the coverage report as a short human-readable block.
    pub fn render(&self) -> String {
        let mut s = String::from("[traceability]\n");
        for (id, metric) in &self.requirement_metrics {
            s.push_str(&format!("  {id} -> {metric}\n"));
        }
        if self.unconstrained_metrics.is_empty() {
            s.push_str("  unconstrained metrics: (none)\n");
        } else {
            s.push_str(&format!(
                "  unconstrained metrics: {}\n",
                self.unconstrained_metrics.join(", ")
            ));
        }
        if !self.uncovered_requirements.is_empty() {
            s.push_str(&format!(
                "  requirements with no metric data: {}\n",
                self.uncovered_requirements.join(", ")
            ));
        }
        s
    }
}

/// Build a [`CoverageReport`] linking `requirements` to the set of
/// `produced_metrics` (the metric names an evaluator / design point actually
/// produces).
///
/// - `requirement_metrics` lists each requirement's `(id, metric)` link.
/// - `unconstrained_metrics` are produced metrics no requirement constrains.
/// - `uncovered_requirements` are requirements whose metric is absent from the
///   produced set.
///
/// Both empty inputs are handled gracefully (no panic): with no requirements,
/// every produced metric is reported unconstrained; with no produced metrics,
/// every requirement is reported uncovered.
pub fn coverage(requirements: &[Requirement], produced_metrics: &[String]) -> CoverageReport {
    let constrained: BTreeSet<&str> = requirements.iter().map(|r| r.metric.as_str()).collect();
    let produced: BTreeSet<&str> = produced_metrics.iter().map(|s| s.as_str()).collect();

    let requirement_metrics = requirements
        .iter()
        .map(|r| (r.id.clone(), r.metric.clone()))
        .collect();

    let unconstrained_metrics = produced
        .iter()
        .filter(|m| !constrained.contains(*m))
        .map(|m| m.to_string())
        .collect();

    let uncovered_requirements = requirements
        .iter()
        .filter(|r| !produced.contains(r.metric.as_str()))
        .map(|r| r.id.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();

    CoverageReport {
        requirement_metrics,
        unconstrained_metrics,
        uncovered_requirements,
    }
}

// --- internal: full-factorial enumeration --------------------------------

/// Enumerate the Cartesian product of the sweeps' levels, in sweep order with
/// the last sweep varying fastest. Empty when there are no sweeps or any sweep
/// has zero levels.
fn cartesian_product(sweeps: &[ParameterSweep]) -> Vec<Vec<f64>> {
    if sweeps.is_empty() || sweeps.iter().any(|s| s.levels.is_empty()) {
        return Vec::new();
    }
    let mut out: Vec<Vec<f64>> = vec![Vec::new()];
    for sweep in sweeps {
        let mut next = Vec::with_capacity(out.len() * sweep.levels.len());
        for prefix in &out {
            for &level in &sweep.levels {
                let mut row = prefix.clone();
                row.push(level);
                next.push(row);
            }
        }
        out = next;
    }
    out
}

/// A compact label for a grid point, e.g. `"wall=2.0; rib=3.0"`.
fn point_label(sweeps: &[ParameterSweep], combo: &[f64]) -> String {
    sweeps
        .iter()
        .zip(combo.iter())
        .map(|(s, v)| format!("{}={}", s.name, v))
        .collect::<Vec<_>>()
        .join("; ")
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- Test 1: Requirement::evaluate for each comparator -----------------

    #[test]
    fn evaluate_min_comparator() {
        // FoS must be >= 1.5.
        let req = Requirement::new("R", "factor of safety", "fos", Comparator::Min, 1.5);

        let pass = req.evaluate(2.0);
        assert!(pass.pass);
        assert!((pass.margin - 0.5).abs() < 1e-12); // 2.0 - 1.5

        let fail = req.evaluate(1.0);
        assert!(!fail.pass);
        assert!((fail.margin - (-0.5)).abs() < 1e-12); // 1.0 - 1.5

        // Exactly on the bound passes with zero margin.
        let edge = req.evaluate(1.5);
        assert!(edge.pass);
        assert!(edge.margin.abs() < 1e-12);
    }

    #[test]
    fn evaluate_max_comparator() {
        // Mass must be <= 5.0.
        let req = Requirement::new("R", "mass budget", "mass", Comparator::Max, 5.0);

        let pass = req.evaluate(4.2);
        assert!(pass.pass);
        assert!((pass.margin - 0.8).abs() < 1e-12); // 5.0 - 4.2

        let fail = req.evaluate(6.0);
        assert!(!fail.pass);
        assert!((fail.margin - (-1.0)).abs() < 1e-12); // 5.0 - 6.0
    }

    #[test]
    fn evaluate_range_comparator() {
        // Operating temp must be within [20, 80].
        let req = Requirement::new(
            "R",
            "operating temperature",
            "temp_c",
            Comparator::Range { lo: 20.0, hi: 80.0 },
            0.0,
        );

        // Inside, nearer the high edge: margin = min(50-20, 80-50)=30 -> wait,
        // value 50 -> min(30, 30) = 30.
        let mid = req.evaluate(50.0);
        assert!(mid.pass);
        assert!((mid.margin - 30.0).abs() < 1e-12);

        // Inside, nearer the low edge: value 25 -> min(5, 55) = 5.
        let near_lo = req.evaluate(25.0);
        assert!(near_lo.pass);
        assert!((near_lo.margin - 5.0).abs() < 1e-12);

        // Below the band: value 10 -> min(-10, 70) = -10 (fail).
        let below = req.evaluate(10.0);
        assert!(!below.pass);
        assert!((below.margin - (-10.0)).abs() < 1e-12);

        // Above the band: value 100 -> min(80, -20) = -20 (fail).
        let above = req.evaluate(100.0);
        assert!(!above.pass);
        assert!((above.margin - (-20.0)).abs() < 1e-12);
    }

    // --- Test 2: ComplianceReport over a mixed pass/fail design point ------

    #[test]
    fn compliance_report_mixed() {
        let reqs = vec![
            Requirement::new("R1", "mass", "mass", Comparator::Max, 5.0), // 4.0 -> +1.0 pass
            Requirement::new("R2", "fos", "fos", Comparator::Min, 1.5),   // 1.2 -> -0.3 FAIL
            Requirement::new("R3", "cost", "cost", Comparator::Max, 100.0), // 90 -> +10 pass
        ];
        let dp = DesignPoint::new("mixed")
            .with_metric("mass", 4.0)
            .with_metric("fos", 1.2)
            .with_metric("cost", 90.0);

        let report = dp.check(&reqs);
        assert_eq!(report.per_requirement.len(), 3);
        assert!(!report.overall_pass); // R2 fails
                                       // Worst margin is R2's -0.3.
        assert!((report.worst_margin - (-0.3)).abs() < 1e-12);

        // Per-requirement verdicts line up by id and pass-state.
        assert!(report.per_requirement[0].1.pass);
        assert!(!report.per_requirement[1].1.pass);
        assert!(report.per_requirement[2].1.pass);

        // An all-pass point reports overall_pass and the *smallest* head-room.
        let dp_ok = DesignPoint::new("ok")
            .with_metric("mass", 4.0) // +1.0
            .with_metric("fos", 2.0) // +0.5  <- worst (most binding)
            .with_metric("cost", 90.0); // +10
        let ok = dp_ok.check(&reqs);
        assert!(ok.overall_pass);
        assert!((ok.worst_margin - 0.5).abs() < 1e-12);
    }

    #[test]
    fn compliance_missing_metric_fails_not_passes() {
        // A requirement whose metric the design point never computed must FAIL,
        // never silently pass.
        let reqs = vec![Requirement::new("R1", "mass", "mass", Comparator::Max, 5.0)];
        let dp = DesignPoint::new("no-mass").with_metric("cost", 10.0);
        let report = dp.check(&reqs);
        assert!(!report.overall_pass);
        assert!(report.worst_margin.is_infinite() && report.worst_margin < 0.0);
    }

    // --- Test 3: full-factorial DOE produces every grid point --------------

    #[test]
    fn full_factorial_3x4_yields_12_points() {
        let sweeps = vec![
            ParameterSweep::new("a", "mm", vec![1.0, 2.0, 3.0]), // 3 levels
            ParameterSweep::new("b", "mm", vec![10.0, 20.0, 30.0, 40.0]), // 4 levels
        ];
        // Evaluator returns two metrics: sum and product of the two params.
        let study = TradeStudy::new(
            sweeps,
            vec!["sum".to_string(), "product".to_string()],
            Vec::new(),
        );
        assert_eq!(study.point_count(), 12);

        let outcome = study.run(&|p: &[f64]| vec![p[0] + p[1], p[0] * p[1]]);
        assert_eq!(outcome.results.len(), 12);

        // Spot-check: the (a=2, b=30) point has sum=32, product=60, and both
        // parameters recorded with their units.
        let r = outcome
            .results
            .iter()
            .find(|r| {
                r.design_point
                    .parameters
                    .iter()
                    .any(|p| p.name == "a" && p.value == 2.0)
                    && r.design_point
                        .parameters
                        .iter()
                        .any(|p| p.name == "b" && p.value == 30.0)
            })
            .expect("grid must contain the (2,30) point");
        assert_eq!(r.design_point.metric("sum"), Some(32.0));
        assert_eq!(r.design_point.metric("product"), Some(60.0));
        assert_eq!(r.design_point.parameters.len(), 2);

        // Every grid point carries both evaluated metrics.
        for res in &outcome.results {
            assert!(res.design_point.metric("sum").is_some());
            assert!(res.design_point.metric("product").is_some());
        }

        // The full Cartesian product is covered (all 12 unique (a,b) pairs).
        let mut pairs: Vec<(f64, f64)> = outcome
            .results
            .iter()
            .map(|r| {
                let a = r
                    .design_point
                    .parameters
                    .iter()
                    .find(|p| p.name == "a")
                    .unwrap()
                    .value;
                let b = r
                    .design_point
                    .parameters
                    .iter()
                    .find(|p| p.name == "b")
                    .unwrap()
                    .value;
                (a, b)
            })
            .collect();
        pairs.sort_by(|x, y| x.partial_cmp(y).unwrap());
        pairs.dedup();
        assert_eq!(pairs.len(), 12);
    }

    // --- Test 4: Pareto front on a known small set -------------------------

    #[test]
    fn pareto_front_known_set() {
        // Two objectives: minimize mass, maximize payload.
        // Points (mass, payload):
        //   A (2, 10)  -- best payload, not lightest
        //   B (1,  5)  -- lightest, modest payload
        //   C (3,  4)  -- dominated by both A and B (heavier AND less payload than A; heavier+less than... B? B=1,5 vs C=3,4: B lighter and more payload -> B dominates C)
        //   D (2,  8)  -- dominated by A (same mass, less payload)
        // Non-dominated front = {A, B}.
        let mk = |name: &str, mass: f64, payload: f64| {
            DesignPoint::new(name)
                .with_metric("mass", mass)
                .with_metric("payload", payload)
        };
        let a = mk("A", 2.0, 10.0);
        let b = mk("B", 1.0, 5.0);
        let c = mk("C", 3.0, 4.0);
        let d = mk("D", 2.0, 8.0);
        let pts = vec![&a, &b, &c, &d];

        let objs = vec![Objective::minimize("mass"), Objective::maximize("payload")];
        let front = pareto_front(&pts, &objs);
        // Indices 0 (A) and 1 (B) only.
        assert_eq!(front, vec![0, 1]);

        // Sanity: A does not dominate B (A heavier) and B does not dominate A
        // (A more payload) -> both on the front.
        assert!(!dominates(&a, &b, &objs));
        assert!(!dominates(&b, &a, &objs));
        // A dominates D (same mass, more payload).
        assert!(dominates(&a, &d, &objs));
        // B dominates C (lighter and more payload).
        assert!(dominates(&b, &c, &objs));

        // Same front when driven through a TradeStudyOutcome.
        let outcome = TradeStudyOutcome {
            results: pts
                .iter()
                .map(|dp| TradeResult {
                    design_point: (*dp).clone(),
                    compliance: dp.check(&[]),
                })
                .collect(),
        };
        let front_results = outcome.pareto_front(&objs);
        let names: Vec<&str> = front_results
            .iter()
            .map(|r| r.design_point.name.as_str())
            .collect();
        assert_eq!(names, vec!["A", "B"]);
    }

    #[test]
    fn pareto_no_objectives_keeps_all() {
        let a = DesignPoint::new("A").with_metric("m", 1.0);
        let b = DesignPoint::new("B").with_metric("m", 2.0);
        let pts = vec![&a, &b];
        assert_eq!(pareto_front(&pts, &[]), vec![0, 1]);
    }

    // --- Test 5: degenerate inputs are graceful ----------------------------

    #[test]
    fn empty_requirements_is_vacuously_passing() {
        let dp = DesignPoint::new("d").with_metric("m", 1.0);
        let report = dp.check(&[]);
        assert!(report.overall_pass);
        assert!(report.worst_margin.is_infinite() && report.worst_margin > 0.0);
        assert!(report.per_requirement.is_empty());
    }

    #[test]
    fn empty_sweep_runs_without_panic() {
        // No sweeps at all -> zero design points, no panic.
        let study = TradeStudy::new(Vec::new(), vec!["m".to_string()], Vec::new());
        assert_eq!(study.point_count(), 0);
        let outcome = study.run(&|_p: &[f64]| vec![0.0]);
        assert!(outcome.results.is_empty());

        // A sweep declared with zero levels also collapses to nothing.
        let study2 = TradeStudy::new(
            vec![ParameterSweep::new("a", "mm", Vec::new())],
            vec!["m".to_string()],
            Vec::new(),
        );
        assert_eq!(study2.point_count(), 0);
        let outcome2 = study2.run(&|_p: &[f64]| vec![0.0]);
        assert!(outcome2.results.is_empty());

        // Pareto over an empty outcome is empty, not a panic.
        assert!(outcome2
            .pareto_front(&[Objective::minimize("m")])
            .is_empty());
    }

    #[test]
    fn coverage_flags_unconstrained_metric() {
        // Two requirements constrain "mass" and "fos"; the evaluator also
        // produces "cost", which nothing constrains.
        let reqs = vec![
            Requirement::new("R1", "mass", "mass", Comparator::Max, 5.0),
            Requirement::new("R2", "fos", "fos", Comparator::Min, 1.5),
        ];
        let produced = vec!["mass".to_string(), "fos".to_string(), "cost".to_string()];
        let cov = coverage(&reqs, &produced);

        assert_eq!(cov.unconstrained_metrics, vec!["cost".to_string()]);
        assert!(cov.uncovered_requirements.is_empty());
        assert!(!cov.is_fully_covered());
        assert_eq!(cov.requirement_metrics.len(), 2);

        // Degenerate: no requirements -> every produced metric is unconstrained.
        let cov_none = coverage(&[], &produced);
        assert_eq!(cov_none.unconstrained_metrics.len(), 3);
        assert!(cov_none.uncovered_requirements.is_empty());

        // Degenerate: no produced metrics -> every requirement is uncovered.
        let cov_empty = coverage(&reqs, &[]);
        assert!(cov_empty.unconstrained_metrics.is_empty());
        assert_eq!(cov_empty.uncovered_requirements.len(), 2);
        assert!(!cov_empty.is_fully_covered());
    }
}
