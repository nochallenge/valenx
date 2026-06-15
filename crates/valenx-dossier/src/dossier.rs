//! The run dossier: assemble, hash, and render.

use std::fmt::Write as _;

use serde::{Deserialize, Serialize};
use valenx_repro::{Artifact, ArtifactRole, Parameter, ProvenanceStep, ReproBundle, SoftwareRef};

use crate::calibration::CalibrationStatus;
use crate::candidate::ScoredCandidate;
use crate::error::DossierError;

/// A complete biologic-design run record: the goal, its ranked candidates, the
/// calibration status, and the software manifest. Build it, then read off
/// [`RunDossier::fingerprint`] (reproducible hash) and [`RunDossier::render`]
/// (human report).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RunDossier {
    /// The design goal (e.g. "Inhibit myostatin (GDF-8)").
    pub goal: String,
    /// Candidates (any order; [`RunDossier::ranked`] sorts them).
    pub candidates: Vec<ScoredCandidate>,
    /// Whether confidences are calibrated, or blocked.
    pub calibration: CalibrationStatus,
    /// Software manifest `(name, version)`.
    pub software: Vec<(String, String)>,
}

impl RunDossier {
    /// Start a dossier for `goal` (non-empty) with a calibration status.
    pub fn new(
        goal: impl Into<String>,
        calibration: CalibrationStatus,
    ) -> Result<Self, DossierError> {
        let goal = goal.into();
        if goal.trim().is_empty() {
            return Err(DossierError::Empty { what: "goal" });
        }
        Ok(Self {
            goal,
            candidates: Vec::new(),
            calibration,
            software: Vec::new(),
        })
    }

    /// Add a candidate.
    pub fn with_candidate(mut self, candidate: ScoredCandidate) -> Self {
        self.candidates.push(candidate);
        self
    }

    /// Record a piece of software used in the run.
    pub fn with_software(mut self, name: impl Into<String>, version: impl Into<String>) -> Self {
        self.software.push((name.into(), version.into()));
        self
    }

    /// Candidates sorted best-first by comparable score (ties broken by id).
    pub fn ranked(&self) -> Vec<&ScoredCandidate> {
        let mut v: Vec<&ScoredCandidate> = self.candidates.iter().collect();
        v.sort_by(|a, b| {
            b.comparable_score
                .partial_cmp(&a.comparable_score)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.id.cmp(&b.id))
        });
        v
    }

    /// Always `true`: a dossier ranks and flags but never approves. Promotion of
    /// any candidate to wet-lab testing is a human decision.
    pub fn requires_human_signoff(&self) -> bool {
        true
    }

    /// Build the content-hashed [`ReproBundle`] backing this dossier: candidates
    /// as ordered output artifacts, calibration + sign-off as parameters, the
    /// software manifest, and the assembly step.
    pub fn to_repro_bundle(&self) -> Result<ReproBundle, DossierError> {
        let mut b = ReproBundle::new(
            format!("Biologic-design dossier: {}", self.goal),
            format!(
                "Ranked, safety-flagged candidate shortlist. Calibration: {}. \
                 No candidate is marked safe; wet-lab promotion requires human sign-off.",
                self.calibration.summary()
            ),
        )?
        .with_parameter(Parameter::new("calibration", self.calibration.summary()))
        .with_parameter(Parameter::new("requires_human_signoff", "true"))
        .with_parameter(Parameter::new(
            "n_candidates",
            self.candidates.len().to_string(),
        ));

        for (name, version) in &self.software {
            b = b.with_software(SoftwareRef::new(name.clone(), version.clone()));
        }

        for (rank, c) in self.ranked().iter().enumerate() {
            let conf = c
                .calibrated_confidence
                .map(|x| format!("{x:.3}"))
                .unwrap_or_else(|| "uncalibrated".to_string());
            let flags = if c.safety_flags.is_empty() {
                "none (NOT a safety assertion)".to_string()
            } else {
                c.safety_flags.join(",")
            };
            let line = format!(
                "rank={rank} id={} score={:.4} confidence={conf} flags={flags}",
                c.id, c.comparable_score
            );
            b = b.with_artifact(Artifact::from_bytes(
                format!("candidate_{rank}"),
                ArtifactRole::Output,
                line.as_bytes(),
            ));
        }

        b = b.with_step(ProvenanceStep::new(
            1,
            "valenx-dossier",
            env!("CARGO_PKG_VERSION"),
            self.goal.clone(),
        ))?;
        Ok(b)
    }

    /// The reproducible SHA-256 fingerprint of the dossier (via [`ReproBundle`]).
    pub fn fingerprint(&self) -> Result<String, DossierError> {
        Ok(self.to_repro_bundle()?.fingerprint())
    }

    /// A human-readable dossier report — the thing an operator reads to cut a
    /// shortlist. Ends with the mandatory no-auto-safe / human-sign-off notice.
    pub fn render(&self) -> String {
        let mut s = String::new();
        let _ = writeln!(s, "Biologic-design dossier");
        let _ = writeln!(s, "  goal       : {}", self.goal);
        let _ = writeln!(s, "  calibration: {}", self.calibration.summary());
        let _ = writeln!(s, "  candidates : {}", self.candidates.len());
        let _ = writeln!(s);
        for (rank, c) in self.ranked().iter().enumerate() {
            let conf = c
                .calibrated_confidence
                .map(|x| format!("{:.1}%", x * 100.0))
                .unwrap_or_else(|| "uncalibrated".to_string());
            let _ = writeln!(
                s,
                "  #{rank}  {:<16} score {:.4}  confidence {conf}",
                c.id, c.comparable_score
            );
            if !c.components.is_empty() {
                let comps: Vec<String> = c
                    .components
                    .iter()
                    .map(|(n, v)| format!("{n}={v:.3}"))
                    .collect();
                let _ = writeln!(s, "        components: {}", comps.join(", "));
            }
            let flags = if c.safety_flags.is_empty() {
                "none raised (NOT a safety assertion)".to_string()
            } else {
                c.safety_flags.join(", ")
            };
            let _ = writeln!(s, "        safety flags: {flags}");
        }
        let _ = writeln!(s);
        let _ = writeln!(
            s,
            "  fingerprint (SHA-256): {}",
            self.fingerprint()
                .unwrap_or_else(|_| "<unavailable>".to_string())
        );
        let _ = writeln!(s);
        let _ = writeln!(
            s,
            "NOTICE: This is a ranked, flagged shortlist. No candidate is marked \
             safe or approved. Promotion of any candidate to wet-lab testing \
             requires explicit human operator sign-off."
        );
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> RunDossier {
        RunDossier::new(
            "Inhibit myostatin (GDF-8)",
            CalibrationStatus::blocked("SKEMPI dataset"),
        )
        .unwrap()
        .with_software("valenx-score", "0.1.0")
        .with_candidate(ScoredCandidate::new("design_A", 0.82).unwrap())
        .with_candidate(
            ScoredCandidate::new("design_B", 0.91)
                .unwrap()
                .with_flag("GDF11_crossreactivity"),
        )
    }

    #[test]
    fn ranked_is_best_first() {
        let d = sample();
        let r = d.ranked();
        assert_eq!(r[0].id, "design_B"); // 0.91
        assert_eq!(r[1].id, "design_A"); // 0.82
    }

    #[test]
    fn rejects_empty_goal() {
        assert_eq!(
            RunDossier::new("  ", CalibrationStatus::blocked("x"))
                .unwrap_err()
                .code(),
            "empty"
        );
    }

    #[test]
    fn fingerprint_is_deterministic_64_hex() {
        let d = sample();
        assert_eq!(d.fingerprint().unwrap(), d.fingerprint().unwrap());
        assert_eq!(d.fingerprint().unwrap().len(), 64);
        // changing a score changes the fingerprint (tamper-evidence)
        let mut d2 = sample();
        d2.candidates[0].comparable_score = 0.5;
        assert_ne!(d.fingerprint().unwrap(), d2.fingerprint().unwrap());
    }

    #[test]
    fn always_requires_human_signoff_and_never_says_safe() {
        let d = sample();
        assert!(d.requires_human_signoff());
        let report = d.render();
        assert!(report.contains("human operator sign-off"));
        assert!(!report.to_lowercase().contains("approved for"));
        // a flagless candidate is explicitly NOT asserted safe
        assert!(report.contains("NOT a safety assertion"));
    }

    #[test]
    fn render_surfaces_blocked_calibration() {
        let report = sample().render();
        assert!(report.contains("BLOCKED: SKEMPI dataset"));
    }

    #[test]
    fn serde_round_trips() {
        let d = sample();
        let j = serde_json::to_string(&d).unwrap();
        let back: RunDossier = serde_json::from_str(&j).unwrap();
        assert_eq!(d, back);
    }
}
