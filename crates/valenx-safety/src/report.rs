//! The consolidated per-candidate risk report.

use std::fmt::Write as _;

use serde::{Deserialize, Serialize};

use crate::error::SafetyError;
use crate::flag::RiskFlag;
use crate::severity::Severity;

/// All of a candidate's safety flags in one place, plus an optional calibrated
/// confidence. **Never carries a "safe" verdict** — see
/// [`RiskReport::requires_human_review`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RiskReport {
    /// The candidate this report is about.
    pub candidate_id: String,
    /// Every flag raised by the candidate's screens.
    pub flags: Vec<RiskFlag>,
    /// Optional calibrated confidence in the overall assessment, in `[0, 1]`.
    pub calibrated_confidence: Option<f64>,
}

impl RiskReport {
    /// A new, empty report for `candidate_id` (non-empty).
    pub fn new(candidate_id: impl Into<String>) -> Result<Self, SafetyError> {
        let candidate_id = candidate_id.into();
        if candidate_id.trim().is_empty() {
            return Err(SafetyError::Empty {
                what: "candidate id",
            });
        }
        Ok(Self {
            candidate_id,
            flags: Vec::new(),
            calibrated_confidence: None,
        })
    }

    /// Add a flag.
    pub fn with_flag(mut self, flag: RiskFlag) -> Self {
        self.flags.push(flag);
        self
    }

    /// Attach a calibrated confidence (`[0, 1]`).
    pub fn with_confidence(mut self, confidence: f64) -> Result<Self, SafetyError> {
        if !confidence.is_finite() || !(0.0..=1.0).contains(&confidence) {
            return Err(SafetyError::ConfidenceOutOfRange { value: confidence });
        }
        self.calibrated_confidence = Some(confidence);
        Ok(self)
    }

    /// The worst flag's severity, or `None` when no flag was raised.
    ///
    /// `None` means *no screen flagged this candidate* — it does **not** mean the
    /// candidate is safe.
    pub fn aggregate_severity(&self) -> Option<Severity> {
        self.flags.iter().map(|f| f.severity).max()
    }

    /// Whether any flag was raised.
    pub fn has_flags(&self) -> bool {
        !self.flags.is_empty()
    }

    /// Always `true`. This report organises evidence; it never approves a
    /// candidate. Promotion to wet-lab testing is a human decision.
    pub fn requires_human_review(&self) -> bool {
        true
    }

    /// A human-readable report — flags worst-first, ending with the mandatory
    /// no-auto-safe notice.
    pub fn render(&self) -> String {
        let mut s = String::new();
        let _ = writeln!(s, "Risk report: candidate {}", self.candidate_id);
        match self.aggregate_severity() {
            Some(sev) => {
                let _ = writeln!(s, "  aggregate severity: {}", sev.as_str());
            }
            None => {
                let _ = writeln!(
                    s,
                    "  aggregate severity: no flags raised (NOT an assertion of safety)"
                );
            }
        }
        if let Some(c) = self.calibrated_confidence {
            let _ = writeln!(s, "  calibrated confidence: {:.1}%", c * 100.0);
        }
        let mut flags: Vec<&RiskFlag> = self.flags.iter().collect();
        flags.sort_by(|a, b| b.severity.cmp(&a.severity));
        for f in flags {
            let ev = f.evidence.map(|v| format!(" [{v:.3}]")).unwrap_or_default();
            let _ = writeln!(
                s,
                "  [{}] {}: {}{ev}",
                f.severity.as_str(),
                f.source,
                f.detail
            );
        }
        let _ = writeln!(s);
        let _ = writeln!(
            s,
            "NOTICE: flags + confidence only. No candidate is marked safe or \
             approved; promotion to wet-lab testing requires explicit human review."
        );
        s
    }
}

/// Consolidate a candidate's flags into one report.
pub fn consolidate(
    candidate_id: impl Into<String>,
    flags: Vec<RiskFlag>,
) -> Result<RiskReport, SafetyError> {
    let mut r = RiskReport::new(candidate_id)?;
    r.flags = flags;
    Ok(r)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::flag::{crispr_offtarget_flag, immunogenicity_flag, offtarget_flag};

    #[test]
    fn aggregate_is_worst_flag() {
        let flags = vec![
            immunogenicity_flag(0.3, 0.1).unwrap().unwrap(), // Moderate
            offtarget_flag("GDF-11", 0.899, 0.8).unwrap().unwrap(), // High
        ];
        let r = consolidate("design_A", flags).unwrap();
        assert_eq!(r.aggregate_severity(), Some(Severity::High));
        assert!(r.has_flags());
    }

    #[test]
    fn no_flags_is_not_safe() {
        let r = RiskReport::new("clean_design").unwrap();
        assert!(r.aggregate_severity().is_none());
        assert!(!r.has_flags());
        assert!(r.requires_human_review()); // still requires review
        let report = r.render();
        assert!(report.contains("NOT an assertion of safety"));
        assert!(report.contains("human review"));
    }

    #[test]
    fn always_requires_human_review() {
        let r = RiskReport::new("x")
            .unwrap()
            .with_flag(crispr_offtarget_flag(12, 1).unwrap().unwrap());
        assert!(r.requires_human_review());
        assert!(!r.render().to_lowercase().contains("approved for"));
    }

    #[test]
    fn rejects_bad_input() {
        assert_eq!(RiskReport::new("  ").unwrap_err().code(), "empty");
        assert_eq!(
            RiskReport::new("x")
                .unwrap()
                .with_confidence(2.0)
                .unwrap_err()
                .code(),
            "confidence_out_of_range"
        );
    }

    #[test]
    fn serde_round_trips() {
        let r = consolidate(
            "design_A",
            vec![offtarget_flag("GDF-11", 0.899, 0.8).unwrap().unwrap()],
        )
        .unwrap();
        let j = serde_json::to_string(&r).unwrap();
        let back: RiskReport = serde_json::from_str(&j).unwrap();
        assert_eq!(r, back);
    }
}
