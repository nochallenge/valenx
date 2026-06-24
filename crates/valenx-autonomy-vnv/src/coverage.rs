//! Coverage metrics over a swept suite.
//!
//! Two complementary notions of coverage, both reported as exact fractions:
//!
//! * **Parameter-space coverage** ([`ParamCoverage`]) — given a [`ParamGrid`]
//!   (named axes, each a finite list of values, defining a discrete cell grid),
//!   which cells did the suite's scenarios actually land in? The fraction is
//!   `cells_exercised / total_cells`. This answers "did we test the whole
//!   parameter envelope we said we would?".
//!
//! * **Requirement coverage** ([`RequirementCoverage`]) — across a set of
//!   evaluated runs, which requirements were *exercised* (actually evaluated to a
//!   real pass/fail), and which were *triggered* (driven to a failure at least
//!   once — i.e. the test was strong enough to make the requirement bite)? This
//!   answers "did our scenarios actually stress each requirement?". A requirement
//!   that never fails anywhere may simply be safe — or may never have been
//!   exercised hard enough; coverage makes that distinction visible.
//!
//! Coverage is computed from the *labels* a scenario carries
//! ([`crate::Scenario::params`]) matched against the grid; it does not re-run
//! anything.

use std::collections::{BTreeMap, BTreeSet};

use crate::error::VnvError;
use crate::report::VnvReport;
use crate::scenario::ScenarioSuite;

/// A discrete parameter grid: a set of named axes, each with a finite list of
/// allowed values. The total number of cells is the product of the axis
/// lengths. Used as the *denominator* for parameter-space coverage.
#[derive(Debug, Clone, PartialEq)]
pub struct ParamGrid {
    /// Axis name → its ordered list of values. A `BTreeMap` so the axis order
    /// (and any derived cell key) is deterministic.
    pub axes: BTreeMap<String, Vec<f64>>,
}

impl ParamGrid {
    /// An empty grid (no axes; one trivial cell — the empty tuple).
    #[must_use]
    pub fn new() -> Self {
        Self {
            axes: BTreeMap::new(),
        }
    }

    /// Add an axis (builder style). Replaces any existing axis of the same name.
    #[must_use]
    pub fn with_axis(mut self, name: impl Into<String>, values: Vec<f64>) -> Self {
        self.axes.insert(name.into(), values);
        self
    }

    /// The total number of cells = product of axis lengths. An empty grid has
    /// exactly one (trivial) cell.
    #[must_use]
    pub fn total_cells(&self) -> usize {
        self.axes.values().map(Vec::len).product::<usize>().max(1)
    }

    /// Validate the grid: every axis non-empty and every value finite.
    ///
    /// # Errors
    /// [`VnvError::InvalidConfig`] for a zero-length axis;
    /// [`VnvError::NonFinite`] for a non-finite value.
    pub fn validate(&self) -> Result<(), VnvError> {
        for (name, values) in &self.axes {
            if values.is_empty() {
                return Err(VnvError::InvalidConfig(format!(
                    "parameter grid axis '{name}' has no values"
                )));
            }
            for v in values {
                if !v.is_finite() {
                    return Err(VnvError::NonFinite(format!(
                        "parameter grid axis '{name}' has a non-finite value ({v})"
                    )));
                }
            }
        }
        Ok(())
    }

    /// Enumerate every cell as an ordered `BTreeMap` of axis→value, in a
    /// deterministic (axis-name-sorted, row-major) order. Used to generate a
    /// full-coverage suite. An empty grid yields a single empty cell.
    #[must_use]
    pub fn cells(&self) -> Vec<BTreeMap<String, f64>> {
        let mut cells = vec![BTreeMap::new()];
        for (name, values) in &self.axes {
            let mut next = Vec::with_capacity(cells.len() * values.len());
            for cell in &cells {
                for &v in values {
                    let mut c = cell.clone();
                    c.insert(name.clone(), v);
                    next.push(c);
                }
            }
            cells = next;
        }
        cells
    }

    /// The cell key a set of scenario parameters maps to: for each grid axis,
    /// the index of the value that (exactly) matches the scenario's value on
    /// that axis. Returns `None` if the scenario is missing an axis or its value
    /// is not one of the axis's listed values (i.e. it falls *outside* the grid
    /// and so covers no cell).
    ///
    /// Matching is exact equality on the `f64` (the sweep generates scenario
    /// params *from* the grid values, so they are bit-identical; this is
    /// deliberate — coverage should not silently snap an off-grid value to a
    /// cell).
    #[must_use]
    pub fn cell_index_of(&self, params: &BTreeMap<String, f64>) -> Option<Vec<usize>> {
        let mut key = Vec::with_capacity(self.axes.len());
        for (name, values) in &self.axes {
            let v = params.get(name)?;
            let idx = values.iter().position(|x| x == v)?;
            key.push(idx);
        }
        Some(key)
    }
}

impl Default for ParamGrid {
    fn default() -> Self {
        Self::new()
    }
}

/// Parameter-space coverage of a suite against a [`ParamGrid`].
#[derive(Debug, Clone, PartialEq)]
pub struct ParamCoverage {
    /// Total number of cells in the grid (the denominator).
    pub total_cells: usize,
    /// Number of distinct cells the suite exercised (the numerator).
    pub covered_cells: usize,
}

impl ParamCoverage {
    /// The coverage fraction in `[0, 1]` (`covered / total`). A grid with one
    /// trivial cell that is covered reads `1.0`.
    #[must_use]
    pub fn fraction(&self) -> f64 {
        if self.total_cells == 0 {
            0.0
        } else {
            self.covered_cells as f64 / self.total_cells as f64
        }
    }

    /// Whether every cell was covered.
    #[must_use]
    pub fn is_complete(&self) -> bool {
        self.covered_cells >= self.total_cells
    }
}

/// Compute parameter-space coverage of `suite` against `grid`.
///
/// Each scenario's [`crate::Scenario::params`] is mapped to a grid cell via
/// [`ParamGrid::cell_index_of`]; scenarios that fall outside the grid are
/// ignored (they cover no cell). The covered count is the number of *distinct*
/// in-grid cells hit.
///
/// # Errors
/// [`ParamGrid::validate`] / [`ScenarioSuite::validate`] errors.
pub fn parameter_coverage(
    suite: &ScenarioSuite,
    grid: &ParamGrid,
) -> Result<ParamCoverage, VnvError> {
    grid.validate()?;
    suite.validate()?;
    let mut covered: BTreeSet<Vec<usize>> = BTreeSet::new();
    for s in &suite.scenarios {
        if let Some(key) = grid.cell_index_of(&s.params) {
            covered.insert(key);
        }
    }
    Ok(ParamCoverage {
        total_cells: grid.total_cells(),
        covered_cells: covered.len(),
    })
}

/// Requirement coverage across a batch of evaluated runs.
#[derive(Debug, Clone, PartialEq)]
pub struct RequirementCoverage {
    /// For each requirement label seen, whether it was *exercised* (evaluated to
    /// a real pass/fail) anywhere, in deterministic label order.
    pub exercised: BTreeMap<String, bool>,
    /// For each requirement label seen, whether it was *triggered* — driven to a
    /// failure (`pass == false`) in at least one run.
    pub triggered: BTreeMap<String, bool>,
}

impl RequirementCoverage {
    /// The number of distinct requirements seen.
    #[must_use]
    pub fn num_requirements(&self) -> usize {
        self.exercised.len()
    }

    /// Fraction of distinct requirements that were exercised somewhere.
    #[must_use]
    pub fn exercised_fraction(&self) -> f64 {
        if self.exercised.is_empty() {
            0.0
        } else {
            let n = self.exercised.values().filter(|&&b| b).count();
            n as f64 / self.exercised.len() as f64
        }
    }

    /// Fraction of distinct requirements that were triggered (failed at least
    /// once). A high pass-rate suite has a *low* triggered fraction — useful for
    /// spotting requirements no scenario ever stressed.
    #[must_use]
    pub fn triggered_fraction(&self) -> f64 {
        if self.triggered.is_empty() {
            0.0
        } else {
            let n = self.triggered.values().filter(|&&b| b).count();
            n as f64 / self.triggered.len() as f64
        }
    }
}

/// Compute requirement coverage from a batch of [`VnvReport`]s (one per run).
///
/// A requirement (keyed by its label) is *exercised* if any outcome for it has
/// `exercised == true`, and *triggered* if any outcome for it failed
/// (`pass == false`).
#[must_use]
pub fn requirement_coverage(reports: &[VnvReport]) -> RequirementCoverage {
    let mut exercised: BTreeMap<String, bool> = BTreeMap::new();
    let mut triggered: BTreeMap<String, bool> = BTreeMap::new();
    for report in reports {
        for o in &report.outcomes {
            *exercised.entry(o.label.clone()).or_insert(false) |= o.exercised;
            *triggered.entry(o.label.clone()).or_insert(false) |= !o.pass;
        }
    }
    RequirementCoverage {
        exercised,
        triggered,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_grid_has_one_trivial_cell() {
        let grid = ParamGrid::new();
        assert_eq!(grid.total_cells(), 1);
        assert_eq!(grid.cells().len(), 1);
        assert!(grid.cells()[0].is_empty());
    }

    #[test]
    fn grid_cell_count_is_product_of_axes() {
        let grid = ParamGrid::new()
            .with_axis("a", vec![1.0, 2.0, 3.0])
            .with_axis("b", vec![10.0, 20.0]);
        assert_eq!(grid.total_cells(), 6);
        let cells = grid.cells();
        assert_eq!(cells.len(), 6);
        // Every cell carries both axes.
        assert!(cells
            .iter()
            .all(|c| c.contains_key("a") && c.contains_key("b")));
        // Cells are distinct.
        let mut keys: Vec<_> = cells
            .iter()
            .map(|c| grid.cell_index_of(c).unwrap())
            .collect();
        keys.sort();
        keys.dedup();
        assert_eq!(keys.len(), 6, "all six cells map to distinct indices");
    }

    #[test]
    fn cell_index_of_rejects_off_grid_and_missing() {
        let grid = ParamGrid::new().with_axis("a", vec![1.0, 2.0]);
        let mut p = BTreeMap::new();
        p.insert("a".to_string(), 1.0);
        assert_eq!(grid.cell_index_of(&p), Some(vec![0]));
        // Off-grid value ⇒ None.
        p.insert("a".to_string(), 1.5);
        assert_eq!(grid.cell_index_of(&p), None);
        // Missing axis ⇒ None.
        assert_eq!(grid.cell_index_of(&BTreeMap::new()), None);
    }
}
