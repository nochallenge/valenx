//! The funnel: chain selection → safety → dossier into one run.

use serde::{Deserialize, Serialize};

use valenx_dossier::{CalibrationStatus, RunDossier, ScoredCandidate};
use valenx_safety::flag::{crispr_offtarget_flag, immunogenicity_flag, offtarget_flag};
use valenx_safety::{consolidate, RiskFlag, RiskReport, Severity};
use valenx_select::{select_shortlist, Candidate as SelectCandidate, Shortlist};

use crate::error::OrchestratorError;

/// A sequence-identity hit against a known off-target, the evidence an
/// off-target screen produces. `identity` is a fraction in `[0, 1]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OfftargetEvidence {
    /// The off-target the candidate resembles (e.g. `"GDF-11"`).
    pub reference: String,
    /// Fractional sequence identity to `reference`, in `[0, 1]`.
    pub identity: f64,
}

/// One candidate entering the funnel.
///
/// `method_scores` and `features` drive selection; the screen inputs drive
/// safety. **An absent optional screen input (`None`) means the screen was not
/// run — it is recorded as such, never treated as a clean pass.** (For
/// sequence-driven runs, [`crate::run_funnel_seqs`] fills these from the actual
/// screens.)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FunnelCandidate {
    /// Stable candidate identifier.
    pub id: String,
    /// One score per orthogonal scoring method (e.g. docking, MM-GBSA,
    /// interface confidence). Every candidate must carry the same number.
    pub method_scores: Vec<f64>,
    /// Descriptor vector used for diversity selection. Every candidate must
    /// carry the same dimension.
    pub features: Vec<f64>,
    /// Optional calibrated confidence in this candidate, in `[0, 1]`.
    pub calibrated_confidence: Option<f64>,
    /// Off-target screen result, or `None` if the screen was not run.
    pub offtarget: Option<OfftargetEvidence>,
    /// Predicted T-cell epitope density, or `None` if the screen was not run.
    pub immunogenicity: Option<f64>,
    /// Predicted CRISPR off-target edit-site count, or `None` if not run.
    pub crispr_offtarget_sites: Option<u32>,
    /// Developability liability flags (e.g. aggregation-prone region, extreme
    /// pI). Empty means none raised; it does **not** mean the screen was skipped.
    pub developability_flags: Vec<String>,
    /// Count of predicted linear B-cell epitope regions, or `None` if the screen
    /// was not run.
    pub bcell_epitope_regions: Option<usize>,
}

/// An upstream pipeline stage that needs a gated resource (a trained model, a
/// docking engine, an experimental structure, a GPU or a license).
///
/// The orchestrator runs the selection → safety → dossier core with none of
/// these. Any stage listed in [`FunnelConfig::blocked_stages`] is recorded as a
/// [`BlockedStage`] (`BLOCKED: <dep>`) and skipped — it is **never faked**.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GatedStage {
    /// De-novo candidate design via a generative model.
    Generate,
    /// Structure-based docking against a target structure.
    Dock,
    /// Physics-based binding-affinity scoring (FEP / MM-GBSA).
    Score,
}

impl GatedStage {
    /// The stage's short name.
    pub fn name(self) -> &'static str {
        match self {
            GatedStage::Generate => "generate",
            GatedStage::Dock => "dock",
            GatedStage::Score => "score",
        }
    }

    /// The gated dependency this stage needs to run for real.
    pub fn dependency(self) -> &'static str {
        match self {
            GatedStage::Generate => "generative model (GPU + trained weights)",
            GatedStage::Dock => "docking engine + experimental target structure",
            GatedStage::Score => "physics-based affinity (FEP/MM-GBSA; GPU or licensed)",
        }
    }
}

/// A gated stage that was skipped because its dependency is unavailable.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct BlockedStage {
    /// The stage that could not run.
    pub stage: GatedStage,
    /// The missing dependency.
    pub dependency: String,
    /// The honest, human-readable `BLOCKED: <dep>` message.
    pub message: String,
}

impl BlockedStage {
    fn from_stage(stage: GatedStage) -> Self {
        let dependency = stage.dependency().to_string();
        Self {
            message: format!("BLOCKED: {} — needs {dependency}", stage.name()),
            dependency,
            stage,
        }
    }
}

/// Tunables for one funnel run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FunnelConfig {
    /// Where the operator cuts the shortlist (the diverse top-`N`).
    pub top_n: usize,
    /// Sphere-exclusion radius in feature space (must be `> 0`); larger = more
    /// diverse / fewer.
    pub diversity_radius: f64,
    /// Off-target identity at/above which a flag is raised.
    pub offtarget_threshold: f64,
    /// Epitope density at/above which an immunogenicity flag is raised.
    pub immunogenicity_threshold: f64,
    /// CRISPR off-target site count at/above which a flag is raised.
    pub crispr_threshold: u32,
    /// Calibration provenance recorded in the dossier (calibrated, or an honest
    /// `BLOCKED` note when no ground truth exists).
    pub calibration: CalibrationStatus,
    /// Upstream gated stages declared unavailable for this run.
    pub blocked_stages: Vec<GatedStage>,
}

/// The result of one funnel run: the assembled dossier plus the intermediate
/// products and the honest record of which gated stages were skipped.
#[derive(Debug, Clone, PartialEq)]
pub struct FunnelOutcome {
    /// The run goal.
    pub goal: String,
    /// The signed, fingerprintable dossier (ranked candidates + flags).
    pub dossier: RunDossier,
    /// The pre-safety diversity cut from the selection stage.
    pub shortlist: Shortlist,
    /// One consolidated safety report per shortlisted candidate.
    pub reports: Vec<RiskReport>,
    /// Gated upstream stages that were skipped, each with its `BLOCKED` note.
    pub blocked: Vec<BlockedStage>,
}

impl FunnelOutcome {
    /// Always `true`: the funnel ranks and flags, it never approves. Promotion
    /// of any candidate to wet-lab testing is a human decision.
    pub fn requires_human_signoff(&self) -> bool {
        true
    }

    /// Content-hash fingerprint of the assembled dossier.
    pub fn fingerprint(&self) -> Result<String, OrchestratorError> {
        Ok(self.dossier.fingerprint()?)
    }

    /// A full human-readable report: goal, gated-stage status, the selection
    /// cut, the dossier, and every safety report.
    pub fn render(&self) -> String {
        let mut s = String::new();
        s.push_str("=== valenx-orchestrator: design funnel ===\n");
        s.push_str(&format!("goal: {}\n\n", self.goal));

        s.push_str("[upstream gated stages]\n");
        if self.blocked.is_empty() {
            s.push_str("  (none declared unavailable)\n");
        } else {
            for b in &self.blocked {
                s.push_str(&format!("  {}\n", b.message));
            }
        }
        s.push('\n');

        s.push_str(&format!(
            "[selection] {} shortlisted (requested top-{})\n",
            self.shortlist.entries.len(),
            self.shortlist.requested_n
        ));
        for e in &self.shortlist.entries {
            s.push_str(&format!(
                "  {}  consensus={:.4}  disagreement={:.4}\n",
                e.id, e.consensus_score, e.disagreement
            ));
        }
        s.push('\n');

        s.push_str(&self.dossier.render());
        s.push('\n');

        s.push_str("[safety reports]\n");
        for r in &self.reports {
            s.push_str(&r.render());
            s.push('\n');
        }

        s.push_str(
            "NOTE: This funnel ranks and flags candidates; it NEVER approves one. \
             No candidate here is asserted safe. Promotion to wet-lab testing is a \
             human decision requiring explicit operator sign-off.\n",
        );
        s
    }
}

/// Run the full design funnel for `goal` over `candidates` with `config`.
///
/// Stages, in order:
/// 1. **selection** ([`valenx_select`]): consensus-rank across the candidates'
///    method scores, then take a feature-diverse top-`N` by sphere exclusion.
/// 2. **safety** ([`valenx_safety`]): for each shortlisted candidate, consolidate
///    its off-target / immunogenicity / CRISPR screens into one risk report.
///    A screen with no input is recorded as `Info: not run`, **never** as a pass.
/// 3. **dossier** ([`valenx_dossier`]): assemble the goal, ranked candidates
///    (with broken-out scores, confidence and flags), calibration provenance and
///    software list into one fingerprintable [`RunDossier`].
///
/// Gated upstream stages (generate / dock / score) declared in
/// [`FunnelConfig::blocked_stages`] are recorded as [`BlockedStage`]s and skipped;
/// the core runs on whatever real candidate scores were supplied. Nothing is
/// fabricated and no candidate is ever marked safe.
///
/// # Errors
///
/// Returns [`OrchestratorError`] if `candidates` is empty, if any stage's
/// underlying crate errors (inconsistent score/feature lengths, out-of-range
/// confidence, …), or if an internal invariant is violated.
pub fn run_funnel(
    goal: &str,
    candidates: &[FunnelCandidate],
    config: &FunnelConfig,
) -> Result<FunnelOutcome, OrchestratorError> {
    if candidates.is_empty() {
        return Err(OrchestratorError::Empty { what: "candidates" });
    }

    // --- Gated upstream stages: record BLOCKED, never fabricate. ---
    let blocked: Vec<BlockedStage> = config
        .blocked_stages
        .iter()
        .copied()
        .map(BlockedStage::from_stage)
        .collect();

    // --- Stage 1: selection (consensus -> diversify). No gated deps. ---
    let select_candidates: Vec<SelectCandidate> = candidates
        .iter()
        .map(|c| SelectCandidate {
            id: c.id.clone(),
            method_scores: c.method_scores.clone(),
            features: c.features.clone(),
            calibrated_confidence: c.calibrated_confidence,
            // Safety is consolidated downstream on the shortlist only.
            safety_flags: Vec::new(),
        })
        .collect();
    let shortlist = select_shortlist(&select_candidates, config.top_n, config.diversity_radius)?;

    // --- Stage 2: safety consolidation on the shortlist. ---
    let by_id: std::collections::HashMap<&str, &FunnelCandidate> =
        candidates.iter().map(|c| (c.id.as_str(), c)).collect();
    let mut reports = Vec::with_capacity(shortlist.entries.len());
    for entry in &shortlist.entries {
        let cand = by_id
            .get(entry.id.as_str())
            .copied()
            .ok_or(OrchestratorError::Internal {
                what: "shortlisted id not found among input candidates",
            })?;
        reports.push(consolidate_candidate(cand, config)?);
    }

    // --- Stage 3: dossier assembly (fingerprintable). ---
    let mut dossier = RunDossier::new(goal, config.calibration.clone())?
        .with_software("valenx-orchestrator", env!("CARGO_PKG_VERSION"))
        .with_software("valenx-select", "in-tree")
        .with_software("valenx-safety", "in-tree")
        .with_software("valenx-dossier", "in-tree");

    for (entry, report) in shortlist.entries.iter().zip(reports.iter()) {
        let mut sc = ScoredCandidate::new(entry.id.clone(), entry.consensus_score)?
            .with_component("consensus", entry.consensus_score)?
            .with_component("disagreement", entry.disagreement)?;
        if let Some(conf) = entry.calibrated_confidence {
            sc = sc.with_confidence(conf)?;
        }
        for flag in &report.flags {
            sc = sc.with_flag(render_flag(flag));
        }
        if let Some(sev) = report.aggregate_severity() {
            sc = sc.with_flag(format!("aggregate severity: {}", sev.as_str()));
        }
        dossier = dossier.with_candidate(sc);
    }

    Ok(FunnelOutcome {
        goal: goal.to_string(),
        dossier,
        shortlist,
        reports,
        blocked,
    })
}

/// Consolidate one candidate's three screens into a [`RiskReport`]. A screen
/// with no input becomes an `Info: not run` flag — absence is never a pass.
fn consolidate_candidate(
    cand: &FunnelCandidate,
    config: &FunnelConfig,
) -> Result<RiskReport, OrchestratorError> {
    let mut flags: Vec<RiskFlag> = Vec::new();

    match &cand.offtarget {
        Some(ev) => {
            if let Some(f) = offtarget_flag(&ev.reference, ev.identity, config.offtarget_threshold)?
            {
                flags.push(f);
            }
        }
        None => flags.push(not_run_flag("off-target")),
    }

    match cand.immunogenicity {
        Some(density) => {
            if let Some(f) = immunogenicity_flag(density, config.immunogenicity_threshold)? {
                flags.push(f);
            }
        }
        None => flags.push(not_run_flag("immunogenicity")),
    }

    match cand.crispr_offtarget_sites {
        Some(sites) => {
            if let Some(f) = crispr_offtarget_flag(sites, config.crispr_threshold)? {
                flags.push(f);
            }
        }
        None => flags.push(not_run_flag("crispr-off-target")),
    }

    // Developability liabilities (manufacturability concerns, Low severity).
    for detail in &cand.developability_flags {
        flags.push(RiskFlag::new(
            "developability",
            Severity::Low,
            detail.clone(),
        ));
    }

    // Predicted linear B-cell epitope regions (antibody-binding surface).
    match cand.bcell_epitope_regions {
        Some(n) if n >= 3 => flags.push(
            RiskFlag::new(
                "bcell-epitope",
                Severity::Moderate,
                format!("{n} predicted linear B-cell epitope region(s)"),
            )
            .with_evidence(n as f64),
        ),
        Some(n) if n >= 1 => flags.push(
            RiskFlag::new(
                "bcell-epitope",
                Severity::Low,
                format!("{n} predicted linear B-cell epitope region(s)"),
            )
            .with_evidence(n as f64),
        ),
        Some(_) => {} // zero regions: ran, nothing to flag
        None => flags.push(not_run_flag("bcell-epitope")),
    }

    Ok(consolidate(cand.id.clone(), flags)?)
}

/// An `Info` flag recording that a screen was not run.
fn not_run_flag(screen: &str) -> RiskFlag {
    RiskFlag::new(
        screen,
        Severity::Info,
        "screen not run — absence of a flag is NOT evidence of safety",
    )
}

/// Render one risk flag as a single line.
fn render_flag(flag: &RiskFlag) -> String {
    match flag.evidence {
        Some(ev) => format!(
            "[{}] {}: {} (evidence {ev:.4})",
            flag.severity.as_str(),
            flag.source,
            flag.detail
        ),
        None => format!(
            "[{}] {}: {}",
            flag.severity.as_str(),
            flag.source,
            flag.detail
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cand(id: &str, scores: Vec<f64>, feats: Vec<f64>) -> FunnelCandidate {
        FunnelCandidate {
            id: id.to_string(),
            method_scores: scores,
            features: feats,
            calibrated_confidence: None,
            offtarget: None,
            immunogenicity: None,
            crispr_offtarget_sites: None,
            developability_flags: Vec::new(),
            bcell_epitope_regions: None,
        }
    }

    fn config(top_n: usize) -> FunnelConfig {
        FunnelConfig {
            top_n,
            diversity_radius: 0.5,
            offtarget_threshold: 0.8,
            immunogenicity_threshold: 0.1,
            crispr_threshold: 1,
            calibration: CalibrationStatus::blocked("no held-out ground truth"),
            blocked_stages: vec![GatedStage::Generate, GatedStage::Dock, GatedStage::Score],
        }
    }

    #[test]
    fn empty_candidates_is_error() {
        let err = run_funnel("goal", &[], &config(3)).unwrap_err();
        assert_eq!(err.code(), "empty");
    }

    #[test]
    fn funnel_ranks_and_shortlists() {
        let cands = vec![
            cand("a", vec![3.0, 9.0], vec![0.0, 0.0]),
            cand("b", vec![1.0, 1.0], vec![1.0, 0.0]),
            cand("c", vec![2.0, 5.0], vec![0.0, 1.0]),
        ];
        let out = run_funnel("rank test", &cands, &config(3)).unwrap();
        // 'a' is best in both methods -> top of the consensus.
        assert_eq!(out.shortlist.entries[0].id, "a");
        assert_eq!(out.dossier.ranked()[0].id, "a");
        assert_eq!(out.shortlist.entries.len(), 3);
    }

    #[test]
    fn gated_stages_are_blocked_not_faked() {
        let cands = vec![cand("a", vec![1.0], vec![0.0])];
        let out = run_funnel("blocked test", &cands, &config(1)).unwrap();
        assert_eq!(out.blocked.len(), 3);
        assert!(out
            .blocked
            .iter()
            .all(|b| b.message.starts_with("BLOCKED:")));
        let names: Vec<&str> = out.blocked.iter().map(|b| b.stage.name()).collect();
        assert!(names.contains(&"generate") && names.contains(&"dock") && names.contains(&"score"));
    }

    #[test]
    fn unrun_screen_is_recorded_not_treated_as_safe() {
        let cands = vec![cand("a", vec![1.0], vec![0.0])];
        let out = run_funnel("unrun test", &cands, &config(1)).unwrap();
        let report = &out.reports[0];
        // All four optional screens were None/empty -> four Info "not run" flags
        // (off-target, immunogenicity, CRISPR off-target, B-cell epitope).
        assert_eq!(report.flags.len(), 4);
        assert!(report.flags.iter().all(|f| f.severity == Severity::Info));
        assert!(report.flags.iter().all(|f| f.detail.contains("not run")));
        // A report with only "not run" flags still demands human review.
        assert!(report.requires_human_review());
    }

    #[test]
    fn real_screen_flags_propagate_to_dossier() {
        let mut c = cand("risky", vec![5.0], vec![0.0]);
        c.offtarget = Some(OfftargetEvidence {
            reference: "GDF-11".to_string(),
            identity: 0.899,
        });
        c.immunogenicity = Some(0.05); // below threshold 0.1 -> no flag
        c.crispr_offtarget_sites = Some(0); // below threshold 1 -> no flag
        let out = run_funnel("flag test", &[c], &config(1)).unwrap();
        let report = &out.reports[0];
        // off-target fires High (0.899 in [0.85,0.95)); immuno/crispr ran but
        // were below threshold -> no "not run" flags for them.
        assert_eq!(report.aggregate_severity(), Some(Severity::High));
        let dossier_flags = &out.dossier.candidates[0].safety_flags;
        assert!(dossier_flags.iter().any(|f| f.contains("off-target")));
        assert!(dossier_flags
            .iter()
            .any(|f| f.contains("aggregate severity: high")));
    }

    #[test]
    fn never_approves_and_fingerprints() {
        let cands = vec![cand("a", vec![1.0], vec![0.0])];
        let out = run_funnel("signoff test", &cands, &config(1)).unwrap();
        assert!(out.requires_human_signoff());
        assert!(out.dossier.requires_human_signoff());
        let fp = out.fingerprint().unwrap();
        assert!(!fp.is_empty());
        // Rendered report carries the no-approval banner.
        assert!(out.render().contains("NEVER approves"));
    }
}
