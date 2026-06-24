//! Multi-run **campaign**: sweep a set of parameter configs through the funnel.
//!
//! A [`Campaign`] runs the existing [`run_funnel`] once per named
//! [`FunnelConfig`] in a sweep, over a single shared candidate set, and collects
//! each run's [`FunnelOutcome`] together with its **provenance**: the dossier's
//! content-hash fingerprint, the count of gated stages that were
//! [`BLOCKED`](crate::BlockedStage), and whether the run succeeded. From those it
//! produces a [`CampaignSummary`] (best / worst run by a chosen metric, plus
//! success / blocked tallies).
//!
//! This is *plumbing*: it composes the funnel, it does not add scoring power.
//! Crucially it preserves the funnel's honesty contract — a config that declares
//! gated stages still records them as `BLOCKED: <dep>` (never fabricated), and no
//! candidate is ever marked safe. A run whose config errors (e.g. inconsistent
//! score lengths) is recorded as a **failed** [`CampaignRun`] carrying the stable
//! error [`code`](crate::OrchestratorError::code) — the campaign does not abort on
//! the first failure, so one bad config does not lose the rest of the sweep.
//!
//! ## What "best" means
//!
//! Ranking is by a [`CampaignMetric`] read off each successful run:
//! - [`CampaignMetric::TopConsensus`] — the top shortlisted candidate's consensus
//!   score (higher is better);
//! - [`CampaignMetric::ShortlistSize`] — how many candidates survived the
//!   diversity cut (higher is better);
//! - [`CampaignMetric::FewestBlockedStages`] — how few gated stages were blocked
//!   (fewer is better) — useful when sweeping which resources to make available.
//!
//! A *failed* run is never "best"; it sorts after every successful run.
//!
//! ```
//! use valenx_dossier::CalibrationStatus;
//! use valenx_orchestrator::{
//!     Campaign, CampaignMetric, FunnelCandidate, FunnelConfig, GatedStage,
//! };
//!
//! let candidates = vec![
//!     FunnelCandidate {
//!         id: "A".into(),
//!         method_scores: vec![9.0, 8.0],
//!         features: vec![0.0, 0.0],
//!         calibrated_confidence: None,
//!         offtarget: None,
//!         immunogenicity: None,
//!         crispr_offtarget_sites: None,
//!         developability_flags: Vec::new(),
//!         bcell_epitope_regions: None,
//!     },
//!     FunnelCandidate {
//!         id: "B".into(),
//!         method_scores: vec![1.0, 1.0],
//!         features: vec![1.0, 1.0],
//!         calibrated_confidence: None,
//!         offtarget: None,
//!         immunogenicity: None,
//!         crispr_offtarget_sites: None,
//!         developability_flags: Vec::new(),
//!         bcell_epitope_regions: None,
//!     },
//! ];
//!
//! let mk = |top_n, radius| FunnelConfig {
//!     top_n,
//!     diversity_radius: radius,
//!     offtarget_threshold: 0.8,
//!     immunogenicity_threshold: 0.1,
//!     crispr_threshold: 1,
//!     calibration: CalibrationStatus::blocked("no held-out ground truth"),
//!     blocked_stages: vec![GatedStage::Generate, GatedStage::Dock, GatedStage::Score],
//! };
//!
//! let campaign = Campaign::new("inhibit target X")
//!     .with_config("loose", mk(2, 0.5))
//!     .with_config("tight", mk(1, 2.0));
//!
//! let report = campaign.run(&candidates, CampaignMetric::ShortlistSize).unwrap();
//! assert_eq!(report.runs.len(), 2);
//! assert_eq!(report.summary.total, 2);
//! assert_eq!(report.summary.succeeded, 2);
//! // Every gated stage stayed BLOCKED, never faked.
//! assert!(report.runs.iter().all(|r| r.blocked_stages == 3));
//! ```

use serde::{Deserialize, Serialize};

use crate::error::OrchestratorError;
use crate::funnel::{run_funnel, FunnelCandidate, FunnelConfig, FunnelOutcome};

/// Which metric a [`Campaign`] ranks its runs by.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CampaignMetric {
    /// The top shortlisted candidate's consensus score; **higher is better**.
    TopConsensus,
    /// The number of candidates that survived the diversity cut; **higher is
    /// better**.
    ShortlistSize,
    /// The number of gated stages recorded as `BLOCKED`; **fewer is better**.
    FewestBlockedStages,
}

impl CampaignMetric {
    /// The metric's short name.
    pub fn name(self) -> &'static str {
        match self {
            CampaignMetric::TopConsensus => "top-consensus",
            CampaignMetric::ShortlistSize => "shortlist-size",
            CampaignMetric::FewestBlockedStages => "fewest-blocked-stages",
        }
    }

    /// Whether a larger metric value is the better one.
    pub fn higher_is_better(self) -> bool {
        match self {
            CampaignMetric::TopConsensus | CampaignMetric::ShortlistSize => true,
            CampaignMetric::FewestBlockedStages => false,
        }
    }
}

/// The provenance of one campaign run that **succeeded**: the dossier
/// fingerprint plus the headline numbers the campaign ranks and tallies on.
///
/// This is the per-run content-addressed lineage — the dossier `fingerprint`
/// ties the run back to its exact inputs and is reproducible (see
/// [`FunnelOutcome::fingerprint`]).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RunProvenance {
    /// The assembled dossier's reproducible SHA-256 fingerprint.
    pub fingerprint: String,
    /// The top shortlisted candidate's consensus score, if any candidate
    /// survived the diversity cut.
    pub top_consensus: Option<f64>,
    /// How many candidates survived the diversity cut.
    pub shortlist_size: usize,
}

/// One run in a campaign: its config label and either the run's products +
/// provenance (success), or the stable error code (failure).
///
/// The campaign never aborts on a failing config — it records the failure here
/// (with the underlying [`OrchestratorError::code`]) and carries on, so a single
/// bad config cannot lose the rest of the sweep.
#[derive(Debug, Clone)]
pub struct CampaignRun {
    /// The label given to this run's config in the sweep.
    pub label: String,
    /// The full funnel outcome, or `None` if the run failed.
    pub outcome: Option<FunnelOutcome>,
    /// The run's provenance (fingerprint + headline numbers), or `None` if it
    /// failed or its fingerprint could not be computed.
    pub provenance: Option<RunProvenance>,
    /// How many gated stages this run's config declared `BLOCKED` (recorded even
    /// for a failed run, straight from the config — these are never faked).
    pub blocked_stages: usize,
    /// The stable error code if the run failed, else `None`.
    pub error_code: Option<&'static str>,
}

impl CampaignRun {
    /// Whether this run completed successfully.
    pub fn succeeded(&self) -> bool {
        self.outcome.is_some() && self.error_code.is_none()
    }

    /// The value of `metric` for this run, or `None` if the run failed (a failed
    /// run has no metric and is never "best").
    pub fn metric_value(&self, metric: CampaignMetric) -> Option<f64> {
        let prov = self.provenance.as_ref()?;
        match metric {
            CampaignMetric::TopConsensus => prov.top_consensus,
            CampaignMetric::ShortlistSize => Some(prov.shortlist_size as f64),
            CampaignMetric::FewestBlockedStages => Some(self.blocked_stages as f64),
        }
    }
}

/// The roll-up over a whole campaign: tallies plus the best / worst run labels.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CampaignSummary {
    /// Total runs attempted (one per config in the sweep).
    pub total: usize,
    /// How many runs succeeded.
    pub succeeded: usize,
    /// How many runs failed (config or stage error).
    pub failed: usize,
    /// Total gated stages blocked across all runs (honesty tally — every one of
    /// these was recorded as `BLOCKED`, never fabricated).
    pub total_blocked_stages: usize,
    /// The metric the runs were ranked by.
    pub metric: CampaignMetric,
    /// Label of the best successful run by `metric`, or `None` if none succeeded.
    pub best_label: Option<String>,
    /// Label of the worst successful run by `metric`, or `None` if none
    /// succeeded.
    pub worst_label: Option<String>,
}

/// The full result of running a campaign: every run plus the summary.
#[derive(Debug, Clone)]
pub struct CampaignReport {
    /// The campaign goal.
    pub goal: String,
    /// One entry per config in the sweep, in sweep order.
    pub runs: Vec<CampaignRun>,
    /// The roll-up.
    pub summary: CampaignSummary,
}

impl CampaignReport {
    /// A human-readable report: the summary line, then one line per run.
    pub fn render(&self) -> String {
        let mut s = String::new();
        s.push_str("=== valenx-orchestrator: campaign ===\n");
        s.push_str(&format!("goal: {}\n", self.goal));
        s.push_str(&format!(
            "runs: {} ({} succeeded, {} failed)  metric: {}\n",
            self.summary.total,
            self.summary.succeeded,
            self.summary.failed,
            self.summary.metric.name(),
        ));
        s.push_str(&format!(
            "best: {}   worst: {}\n",
            self.summary.best_label.as_deref().unwrap_or("(none)"),
            self.summary.worst_label.as_deref().unwrap_or("(none)"),
        ));
        s.push_str(&format!(
            "gated stages blocked across all runs: {} (never faked)\n\n",
            self.summary.total_blocked_stages,
        ));
        for run in &self.runs {
            if run.succeeded() {
                let metric = run
                    .metric_value(self.summary.metric)
                    .map(|m| format!("{m:.4}"))
                    .unwrap_or_else(|| "n/a".to_string());
                let fp = run
                    .provenance
                    .as_ref()
                    .map(|p| p.fingerprint.as_str())
                    .unwrap_or("<none>");
                s.push_str(&format!(
                    "  [ok]   {:<16} {}={metric}  blocked={}  fp={fp}\n",
                    run.label,
                    self.summary.metric.name(),
                    run.blocked_stages,
                ));
            } else {
                s.push_str(&format!(
                    "  [FAIL] {:<16} error={}\n",
                    run.label,
                    run.error_code.unwrap_or("unknown"),
                ));
            }
        }
        s.push_str(
            "\nNOTE: A campaign sweeps configs and ranks runs; it NEVER approves a \
             candidate. No run here asserts any candidate is safe. Promotion to \
             wet-lab testing is a human decision requiring explicit sign-off.\n",
        );
        s
    }
}

/// A multi-run sweep: one [`run_funnel`] per labelled [`FunnelConfig`], over a
/// shared candidate set, rolled up into a [`CampaignReport`].
///
/// Build it with [`Campaign::new`] and [`Campaign::with_config`], then
/// [`Campaign::run`].
#[derive(Debug, Clone)]
pub struct Campaign {
    goal: String,
    configs: Vec<(String, FunnelConfig)>,
}

impl Campaign {
    /// Start an empty campaign for `goal`.
    pub fn new(goal: impl Into<String>) -> Self {
        Self {
            goal: goal.into(),
            configs: Vec::new(),
        }
    }

    /// Add a labelled config to the sweep. Labels are free-form (e.g.
    /// `"loose"`, `"radius=0.5"`); they are echoed back in the report and do not
    /// need to be unique, though unique labels make a report easier to read.
    pub fn with_config(mut self, label: impl Into<String>, config: FunnelConfig) -> Self {
        self.configs.push((label.into(), config));
        self
    }

    /// The number of configs in the sweep.
    pub fn len(&self) -> usize {
        self.configs.len()
    }

    /// Whether the sweep is empty.
    pub fn is_empty(&self) -> bool {
        self.configs.is_empty()
    }

    /// Run every config in the sweep over `candidates`, ranking by `metric`.
    ///
    /// Each config is run via [`run_funnel`]. A config that errors is recorded as
    /// a failed [`CampaignRun`] (carrying the stable error code) and the sweep
    /// continues — the campaign as a whole only errors if it has **nothing** to
    /// run.
    ///
    /// # Errors
    ///
    /// Returns [`OrchestratorError::Empty`] if the sweep has no configs, or if
    /// `candidates` is empty (each run would fail identically, so this is
    /// surfaced once up front rather than as N identical failed runs).
    pub fn run(
        &self,
        candidates: &[FunnelCandidate],
        metric: CampaignMetric,
    ) -> Result<CampaignReport, OrchestratorError> {
        if self.configs.is_empty() {
            return Err(OrchestratorError::Empty {
                what: "campaign configs",
            });
        }
        if candidates.is_empty() {
            return Err(OrchestratorError::Empty { what: "candidates" });
        }

        let mut runs = Vec::with_capacity(self.configs.len());
        for (label, config) in &self.configs {
            // A config's blocked stages are known up front and recorded even if
            // the run errors — they are read from the config, never fabricated.
            let blocked_stages = config.blocked_stages.len();
            match run_funnel(&self.goal, candidates, config) {
                Ok(outcome) => {
                    let provenance = outcome.fingerprint().ok().map(|fingerprint| {
                        let top_consensus =
                            outcome.shortlist.entries.first().map(|e| e.consensus_score);
                        RunProvenance {
                            fingerprint,
                            top_consensus,
                            shortlist_size: outcome.shortlist.entries.len(),
                        }
                    });
                    runs.push(CampaignRun {
                        label: label.clone(),
                        outcome: Some(outcome),
                        provenance,
                        blocked_stages,
                        error_code: None,
                    });
                }
                Err(e) => runs.push(CampaignRun {
                    label: label.clone(),
                    outcome: None,
                    provenance: None,
                    blocked_stages,
                    error_code: Some(e.code()),
                }),
            }
        }

        let summary = summarize(&runs, metric);
        Ok(CampaignReport {
            goal: self.goal.clone(),
            runs,
            summary,
        })
    }
}

/// Roll the per-run results into a [`CampaignSummary`].
fn summarize(runs: &[CampaignRun], metric: CampaignMetric) -> CampaignSummary {
    let total = runs.len();
    let succeeded = runs.iter().filter(|r| r.succeeded()).count();
    let failed = total - succeeded;
    let total_blocked_stages = runs.iter().map(|r| r.blocked_stages).sum();

    // Best / worst over successful runs only. A run with no metric value (e.g. an
    // empty shortlist under `TopConsensus`) is excluded from ranking — it cannot
    // be "best" or "worst" on a metric it does not have.
    let mut best: Option<(&str, f64)> = None;
    let mut worst: Option<(&str, f64)> = None;
    for run in runs.iter().filter(|r| r.succeeded()) {
        let Some(value) = run.metric_value(metric) else {
            continue;
        };
        let better_than = |challenger: f64, incumbent: f64| {
            if metric.higher_is_better() {
                challenger > incumbent
            } else {
                challenger < incumbent
            }
        };
        match best {
            Some((_, b)) if !better_than(value, b) => {}
            _ => best = Some((run.label.as_str(), value)),
        }
        match worst {
            Some((_, w)) if better_than(value, w) => {}
            _ => worst = Some((run.label.as_str(), value)),
        }
    }

    CampaignSummary {
        total,
        succeeded,
        failed,
        total_blocked_stages,
        metric,
        best_label: best.map(|(l, _)| l.to_string()),
        worst_label: worst.map(|(l, _)| l.to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::funnel::GatedStage;
    use valenx_dossier::CalibrationStatus;

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

    fn config(top_n: usize, radius: f64) -> FunnelConfig {
        FunnelConfig {
            top_n,
            diversity_radius: radius,
            offtarget_threshold: 0.8,
            immunogenicity_threshold: 0.1,
            crispr_threshold: 1,
            calibration: CalibrationStatus::blocked("no held-out ground truth"),
            blocked_stages: vec![GatedStage::Generate, GatedStage::Dock, GatedStage::Score],
        }
    }

    fn three_candidates() -> Vec<FunnelCandidate> {
        vec![
            cand("a", vec![9.0, 8.0], vec![0.0, 0.0]),
            cand("b", vec![5.0, 5.0], vec![1.0, 0.0]),
            cand("c", vec![1.0, 1.0], vec![0.0, 1.0]),
        ]
    }

    #[test]
    fn empty_sweep_is_error() {
        let err = Campaign::new("g")
            .run(&three_candidates(), CampaignMetric::ShortlistSize)
            .unwrap_err();
        assert_eq!(err.code(), "empty");
    }

    #[test]
    fn empty_candidates_is_error() {
        let err = Campaign::new("g")
            .with_config("x", config(2, 0.5))
            .run(&[], CampaignMetric::ShortlistSize)
            .unwrap_err();
        assert_eq!(err.code(), "empty");
    }

    #[test]
    fn sweeps_all_configs_and_counts() {
        let report = Campaign::new("inhibit X")
            .with_config("loose", config(3, 0.5))
            .with_config("tight", config(1, 5.0))
            .run(&three_candidates(), CampaignMetric::ShortlistSize)
            .unwrap();
        assert_eq!(report.runs.len(), 2);
        assert_eq!(report.summary.total, 2);
        assert_eq!(report.summary.succeeded, 2);
        assert_eq!(report.summary.failed, 0);
    }

    #[test]
    fn shortlist_size_picks_the_looser_radius_as_best() {
        // A small radius keeps more candidates; a huge radius keeps just one.
        let report = Campaign::new("g")
            .with_config("loose", config(3, 0.1))
            .with_config("tight", config(3, 100.0))
            .run(&three_candidates(), CampaignMetric::ShortlistSize)
            .unwrap();
        assert_eq!(report.summary.best_label.as_deref(), Some("loose"));
        assert_eq!(report.summary.worst_label.as_deref(), Some("tight"));
        // The looser run kept more than the tighter run.
        let loose = report.runs.iter().find(|r| r.label == "loose").unwrap();
        let tight = report.runs.iter().find(|r| r.label == "tight").unwrap();
        assert!(
            loose.metric_value(CampaignMetric::ShortlistSize).unwrap()
                > tight.metric_value(CampaignMetric::ShortlistSize).unwrap()
        );
    }

    #[test]
    fn gated_stages_stay_blocked_never_faked() {
        let report = Campaign::new("g")
            .with_config("a", config(2, 0.5))
            .with_config("b", config(2, 0.5))
            .run(&three_candidates(), CampaignMetric::TopConsensus)
            .unwrap();
        // Each run declared 3 gated stages -> each records 3 blocked, never run.
        assert!(report.runs.iter().all(|r| r.blocked_stages == 3));
        assert_eq!(report.summary.total_blocked_stages, 6);
        for run in &report.runs {
            let blocked = &run.outcome.as_ref().unwrap().blocked;
            assert_eq!(blocked.len(), 3);
            assert!(blocked.iter().all(|b| b.message.starts_with("BLOCKED:")));
        }
    }

    #[test]
    fn each_successful_run_has_a_reproducible_fingerprint() {
        let report = Campaign::new("g")
            .with_config("only", config(2, 0.5))
            .run(&three_candidates(), CampaignMetric::TopConsensus)
            .unwrap();
        let prov = report.runs[0].provenance.as_ref().unwrap();
        assert_eq!(prov.fingerprint.len(), 64);
        assert!(prov.fingerprint.chars().all(|c| c.is_ascii_hexdigit()));
        // top_consensus matches the dossier's top-ranked candidate.
        assert!(prov.top_consensus.is_some());
    }

    #[test]
    fn a_failing_config_is_recorded_not_aborting() {
        // Candidates with mismatched method-score lengths make the selection
        // stage error; the other config still runs.
        let bad = vec![
            cand("a", vec![1.0, 2.0], vec![0.0, 0.0]),
            cand("b", vec![1.0], vec![1.0, 0.0]), // shorter scores -> select error
        ];
        let report = Campaign::new("g")
            .with_config("will-fail", config(2, 0.5))
            .run(&bad, CampaignMetric::ShortlistSize)
            .unwrap();
        assert_eq!(report.summary.total, 1);
        assert_eq!(report.summary.failed, 1);
        assert_eq!(report.summary.succeeded, 0);
        assert!(!report.runs[0].succeeded());
        // The stable error code is carried; it is the selection stage.
        assert_eq!(report.runs[0].error_code, Some("select"));
        // With no successful run there is no best/worst.
        assert!(report.summary.best_label.is_none());
        assert!(report.summary.worst_label.is_none());
    }

    #[test]
    fn fewest_blocked_metric_prefers_the_less_gated_config() {
        let mut few = config(2, 0.5);
        few.blocked_stages = vec![GatedStage::Generate]; // only one gated stage
        let many = config(2, 0.5); // three gated stages
        let report = Campaign::new("g")
            .with_config("few-blocked", few)
            .with_config("many-blocked", many)
            .run(&three_candidates(), CampaignMetric::FewestBlockedStages)
            .unwrap();
        // Fewer blocked is better.
        assert_eq!(report.summary.best_label.as_deref(), Some("few-blocked"));
        assert_eq!(report.summary.worst_label.as_deref(), Some("many-blocked"));
    }

    #[test]
    fn render_is_honest_and_lists_runs() {
        let report = Campaign::new("inhibit X")
            .with_config("loose", config(3, 0.1))
            .with_config("tight", config(1, 100.0))
            .run(&three_candidates(), CampaignMetric::ShortlistSize)
            .unwrap();
        let text = report.render();
        assert!(text.contains("campaign"));
        assert!(text.contains("loose"));
        assert!(text.contains("tight"));
        // The no-approval banner must survive.
        assert!(text.contains("NEVER approves"));
        assert!(text.contains("never faked"));
    }
}
