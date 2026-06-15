//! The calibration status of a run — calibrated, or honestly blocked.

use serde::{Deserialize, Serialize};

/// Whether the run's confidences are calibrated against a labelled benchmark,
/// or blocked because no such benchmark was available.
///
/// The [`CalibrationStatus::Blocked`] variant is the honest output when a
/// calibration dataset (e.g. SKEMPI for binding ΔΔG) is absent: the dossier
/// records the missing dependency instead of presenting an uncalibrated score
/// as if it were a probability.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum CalibrationStatus {
    /// Confidences were calibrated by `method` against `calibration_set`.
    Calibrated {
        /// Calibration method (e.g. `"isotonic"`, `"Platt"`).
        method: String,
        /// The labelled set used (e.g. `"SKEMPI v2"`).
        calibration_set: String,
    },
    /// No calibration was possible; `dependency` names what is missing.
    Blocked {
        /// The missing dependency (e.g. `"SKEMPI dataset"`).
        dependency: String,
    },
}

impl CalibrationStatus {
    /// A calibrated status.
    pub fn calibrated(method: impl Into<String>, calibration_set: impl Into<String>) -> Self {
        CalibrationStatus::Calibrated {
            method: method.into(),
            calibration_set: calibration_set.into(),
        }
    }

    /// A blocked status naming the missing `dependency`.
    pub fn blocked(dependency: impl Into<String>) -> Self {
        CalibrationStatus::Blocked {
            dependency: dependency.into(),
        }
    }

    /// Whether confidences in this run are calibrated.
    pub fn is_calibrated(&self) -> bool {
        matches!(self, CalibrationStatus::Calibrated { .. })
    }

    /// A one-line summary for the dossier record.
    pub fn summary(&self) -> String {
        match self {
            CalibrationStatus::Calibrated {
                method,
                calibration_set,
            } => format!("calibrated via {method} on {calibration_set}"),
            CalibrationStatus::Blocked { dependency } => {
                format!("BLOCKED: {dependency} (no calibrated confidence available)")
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn blocked_reports_dependency_and_is_not_calibrated() {
        let c = CalibrationStatus::blocked("SKEMPI dataset");
        assert!(!c.is_calibrated());
        assert!(c.summary().contains("BLOCKED: SKEMPI dataset"));
    }

    #[test]
    fn calibrated_records_method_and_set() {
        let c = CalibrationStatus::calibrated("isotonic", "SKEMPI v2");
        assert!(c.is_calibrated());
        assert!(c.summary().contains("isotonic"));
        assert!(c.summary().contains("SKEMPI v2"));
    }
}
