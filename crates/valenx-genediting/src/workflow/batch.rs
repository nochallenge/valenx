//! Feature 29 — the batch editing-design driver.
//!
//! Editing-design work is rarely a single target — a screen designs
//! reagents for dozens or hundreds of loci at once. This module runs
//! the top-level [`crate::workflow::driver`] over many
//! [`EditingRequest`]s, **captures per-target failures rather than
//! aborting the run**, and assembles a result table.
//!
//! ## v1 scope
//!
//! The batch driver runs each target sequentially (no thread pool — a
//! deliberate dependency-minimal choice consistent with the rest of
//! the crate). A target whose design fails is recorded in
//! [`BatchResult::failures`] with its error; the run always completes.

use crate::error::Result;
use crate::workflow::advisor::EditApproach;
use crate::workflow::driver::{run_editing_design, EditingReport, EditingRequest};
use serde::{Deserialize, Serialize};

/// One labelled target in a batch editing-design run.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BatchTarget {
    /// A human-readable label (gene name, variant id, …).
    pub label: String,
    /// The editing request for this target.
    pub request: EditingRequest,
}

/// One successfully designed target.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BatchEntry {
    /// The target label.
    pub label: String,
    /// The editing report for this target.
    pub report: EditingReport,
}

/// One failed target in a batch run.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BatchFailure {
    /// The target label.
    pub label: String,
    /// The stable error code of the failure.
    pub error_code: String,
    /// A human-readable failure message.
    pub error: String,
}

/// The result of a batch editing-design run (feature 29).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BatchResult {
    /// Successfully designed targets.
    pub entries: Vec<BatchEntry>,
    /// Targets whose design failed (the run continued past them).
    pub failures: Vec<BatchFailure>,
}

impl BatchResult {
    /// Total number of targets attempted.
    pub fn total(&self) -> usize {
        self.entries.len() + self.failures.len()
    }

    /// The success fraction in `[0, 1]` (`1.0` for an empty batch).
    pub fn success_rate(&self) -> f64 {
        if self.total() == 0 {
            return 1.0;
        }
        self.entries.len() as f64 / self.total() as f64
    }

    /// How many successful designs used `approach`.
    pub fn count_by_approach(&self, approach: EditApproach) -> usize {
        self.entries
            .iter()
            .filter(|e| e.report.approach == approach)
            .count()
    }

    /// A compact text summary table of the run — one line per target.
    pub fn summary_table(&self) -> String {
        let mut out = String::new();
        out.push_str("label\tstatus\tapproach\n");
        for e in &self.entries {
            out.push_str(&format!("{}\tok\t{}\n", e.label, e.report.approach.name()));
        }
        for f in &self.failures {
            out.push_str(&format!("{}\tfailed\t{}\n", f.label, f.error_code));
        }
        out
    }
}

/// Runs a batch editing design over many targets (feature 29).
///
/// Each [`BatchTarget`] is dispatched through [`run_editing_design`];
/// successes land in [`BatchResult::entries`], failures (captured, not
/// propagated) in [`BatchResult::failures`]. The run always completes
/// and the result is `Ok` — an individual design failure is data, not
/// a run-level error.
///
/// # Errors
/// This function does not itself fail — it returns `Result` only for
/// signature symmetry with the rest of the crate; the result is always
/// `Ok`.
pub fn run_batch_design(targets: &[BatchTarget]) -> Result<BatchResult> {
    let mut entries = Vec::new();
    let mut failures = Vec::new();
    for t in targets {
        match run_editing_design(&t.request) {
            Ok(report) => entries.push(BatchEntry {
                label: t.label.clone(),
                report,
            }),
            Err(e) => failures.push(BatchFailure {
                label: t.label.clone(),
                error_code: e.code().to_string(),
                error: e.to_string(),
            }),
        }
    }
    Ok(BatchResult { entries, failures })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crispr::nuclease::NucleaseId;
    use crate::workflow::driver::EditingTask;

    fn knockout_cds() -> Vec<u8> {
        let mut s = Vec::new();
        for _ in 0..6 {
            s.extend_from_slice(b"ACGTACGTACGTACGTACGTAGG");
        }
        s
    }

    fn good_target(label: &str) -> BatchTarget {
        BatchTarget {
            label: label.to_string(),
            request: EditingRequest {
                task: EditingTask::KnockoutGene {
                    cds: knockout_cds(),
                    exon_boundaries: vec![0, 69],
                    nuclease: NucleaseId::SpCas9,
                },
            },
        }
    }

    fn bad_target(label: &str) -> BatchTarget {
        BatchTarget {
            label: label.to_string(),
            request: EditingRequest {
                task: EditingTask::KnockoutGene {
                    cds: knockout_cds(),
                    exon_boundaries: vec![7], // must start at 0 → fails
                    nuclease: NucleaseId::SpCas9,
                },
            },
        }
    }

    #[test]
    fn batch_designs_every_good_target() {
        let targets = vec![good_target("A"), good_target("B"), good_target("C")];
        let result = run_batch_design(&targets).unwrap();
        assert_eq!(result.entries.len(), 3);
        assert!(result.failures.is_empty());
        assert_eq!(result.total(), 3);
        assert!((result.success_rate() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn batch_captures_failures_without_aborting() {
        let targets = vec![good_target("A"), bad_target("B"), good_target("C")];
        let result = run_batch_design(&targets).unwrap();
        // The bad target is captured; the good ones still designed.
        assert_eq!(result.entries.len(), 2);
        assert_eq!(result.failures.len(), 1);
        assert_eq!(result.failures[0].label, "B");
        assert!(result.failures[0].error_code.starts_with("genediting."));
    }

    #[test]
    fn success_rate_reflects_failures() {
        let targets = vec![good_target("A"), bad_target("B")];
        let result = run_batch_design(&targets).unwrap();
        assert!((result.success_rate() - 0.5).abs() < 1e-9);
    }

    #[test]
    fn empty_batch_is_ok() {
        let result = run_batch_design(&[]).unwrap();
        assert_eq!(result.total(), 0);
        assert!((result.success_rate() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn count_by_approach_tallies() {
        let targets = vec![good_target("A"), good_target("B")];
        let result = run_batch_design(&targets).unwrap();
        // Both are knockouts.
        assert_eq!(result.count_by_approach(EditApproach::NucleaseNhej), 2);
        assert_eq!(result.count_by_approach(EditApproach::BaseEditing), 0);
    }

    #[test]
    fn summary_table_has_a_row_per_target() {
        let targets = vec![good_target("A"), bad_target("B")];
        let result = run_batch_design(&targets).unwrap();
        let table = result.summary_table();
        assert!(table.contains("A\tok"));
        assert!(table.contains("B\tfailed"));
        // header + 2 rows = 3 lines.
        assert_eq!(table.lines().count(), 3);
    }
}
