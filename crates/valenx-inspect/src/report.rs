//! Inspection report aggregation.
//!
//! An [`InspectReport`] is a flat list of [`ReportRow`] — one row per
//! (measurement, tolerance, result) triple. Rows are appended in
//! order; the report can be serialised to RON via [`crate::InspectFile`]
//! or exported to CSV via [`InspectReport::to_csv`] for hand-off to a
//! metrology workflow.

use serde::{Deserialize, Serialize};

use crate::measurement::{compute, Measurement};
use crate::tolerance::Tolerance;

/// Result of evaluating one measurement against a tolerance band.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum CheckResult {
    /// Value inside the inner Pass band.
    Pass,
    /// Value inside the outer band but outside the inner Pass band.
    Warning,
    /// Value outside the outer band.
    Fail,
}

impl CheckResult {
    /// Single-character status — `'P'` / `'W'` / `'F'`.
    pub fn glyph(self) -> char {
        match self {
            Self::Pass => 'P',
            Self::Warning => 'W',
            Self::Fail => 'F',
        }
    }

    /// Long-form label.
    pub fn label(self) -> &'static str {
        match self {
            Self::Pass => "Pass",
            Self::Warning => "Warning",
            Self::Fail => "Fail",
        }
    }
}

/// One (measurement, tolerance, result) triple.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ReportRow {
    /// Optional human-readable label (e.g. "Bore A diameter").
    pub label: Option<String>,
    /// What was measured.
    pub measurement: Measurement,
    /// Tolerance band the measurement was checked against.
    pub tolerance: Tolerance,
    /// Pass / Warning / Fail.
    pub result: CheckResult,
}

/// Ordered list of rows. The full inspection-run record.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct InspectReport {
    /// Optional title (e.g. "Final inspection — part 1234").
    pub title: String,
    /// Rows in run order.
    pub rows: Vec<ReportRow>,
}

impl InspectReport {
    /// Empty report (no title, no rows).
    pub fn new() -> Self {
        Self::default()
    }

    /// Empty report with title.
    pub fn with_title(title: impl Into<String>) -> Self {
        Self {
            title: title.into(),
            rows: Vec::new(),
        }
    }

    /// Append a pre-computed row.
    pub fn add_row(&mut self, measurement: Measurement, tolerance: Tolerance, result: CheckResult) {
        self.rows.push(ReportRow {
            label: None,
            measurement,
            tolerance,
            result,
        });
    }

    /// Append a row with a custom label.
    pub fn add_labelled(
        &mut self,
        label: impl Into<String>,
        measurement: Measurement,
        tolerance: Tolerance,
        result: CheckResult,
    ) {
        self.rows.push(ReportRow {
            label: Some(label.into()),
            measurement,
            tolerance,
            result,
        });
    }

    /// Convenience: compute + evaluate a measurement against a band
    /// and append. Returns the appended row's index, or
    /// `Err(InspectError)` if measurement computation fails.
    pub fn measure_and_check(
        &mut self,
        label: Option<String>,
        measurement: Measurement,
        tolerance: Tolerance,
    ) -> Result<usize, crate::error::InspectError> {
        let actual = compute(&measurement)?;
        let result = tolerance.evaluate(actual);
        let row = ReportRow {
            label,
            measurement,
            tolerance,
            result,
        };
        self.rows.push(row);
        Ok(self.rows.len() - 1)
    }

    /// Count of Pass rows.
    pub fn pass_count(&self) -> usize {
        self.rows
            .iter()
            .filter(|r| r.result == CheckResult::Pass)
            .count()
    }

    /// Count of Warning rows.
    pub fn warning_count(&self) -> usize {
        self.rows
            .iter()
            .filter(|r| r.result == CheckResult::Warning)
            .count()
    }

    /// Count of Fail rows.
    pub fn fail_count(&self) -> usize {
        self.rows
            .iter()
            .filter(|r| r.result == CheckResult::Fail)
            .count()
    }

    /// `true` when every row is a Pass.
    pub fn all_pass(&self) -> bool {
        self.rows.iter().all(|r| r.result == CheckResult::Pass)
    }

    /// CSV with header row. Columns: label, kind, nominal, upper_dev,
    /// lower_dev, actual, result.
    pub fn to_csv(&self) -> String {
        let mut out = String::new();
        out.push_str("label,kind,nominal,upper_dev,lower_dev,actual,result\n");
        for r in &self.rows {
            let label = r.label.clone().unwrap_or_default();
            let kind = r.measurement.label();
            // Best-effort actual; falls back to "<err>" if recomputation fails.
            let actual = compute(&r.measurement)
                .map(|v| format!("{v}"))
                .unwrap_or_else(|_| "<err>".into());
            out.push_str(&format!(
                "{label},{kind},{},{},{},{actual},{}\n",
                r.tolerance.nominal,
                r.tolerance.upper_dev,
                r.tolerance.lower_dev,
                r.result.label(),
            ));
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tolerance::Tolerance;
    use nalgebra::Vector3;

    #[test]
    fn measure_and_check_pass() {
        let mut r = InspectReport::with_title("Test");
        let idx = r
            .measure_and_check(
                Some("AB".into()),
                Measurement::Distance {
                    from: Vector3::zeros(),
                    to: Vector3::new(10.0, 0.0, 0.0),
                },
                Tolerance::symmetric(10.0, 0.1),
            )
            .unwrap();
        assert_eq!(idx, 0);
        assert_eq!(r.rows[0].result, CheckResult::Pass);
        assert_eq!(r.pass_count(), 1);
        assert!(r.all_pass());
    }

    #[test]
    fn csv_has_header_and_row() {
        let mut r = InspectReport::new();
        r.add_row(
            Measurement::Distance {
                from: Vector3::zeros(),
                to: Vector3::new(1.0, 0.0, 0.0),
            },
            Tolerance::symmetric(1.0, 0.05),
            CheckResult::Pass,
        );
        let s = r.to_csv();
        assert!(s.starts_with("label,kind,"));
        assert!(s.contains("Distance"));
        assert!(s.contains("Pass"));
    }

    #[test]
    fn glyph_round_trip() {
        assert_eq!(CheckResult::Pass.glyph(), 'P');
        assert_eq!(CheckResult::Warning.glyph(), 'W');
        assert_eq!(CheckResult::Fail.glyph(), 'F');
    }
}
