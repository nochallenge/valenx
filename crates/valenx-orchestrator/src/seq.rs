//! Sequence-driven funnel: run the screens from amino-acid sequences.
//!
//! [`run_funnel_seqs`] takes candidate **sequences** (plus the caller's upstream
//! method scores) and a reference panel, *runs the actual screens* to compute the
//! safety evidence — off-target identity, T-cell epitope density, developability
//! and linear B-cell epitopes — derives a diversity descriptor from amino-acid
//! composition, then funnels through the same selection → safety → dossier path
//! as [`crate::run_funnel`]. Nothing is fabricated; a screen that cannot run on a
//! given sequence (e.g. shorter than the epitope matrix) is recorded as not run.

use serde::{Deserialize, Serialize};

use valenx_developability::assess;
use valenx_epitope_map::{hydrophilicity_kd, linear_epitope_regions};
use valenx_immuno::{epitope_density, library::illustrative_hla_a0201};
use valenx_offtarget::best_ungapped_identity;

use crate::error::OrchestratorError;
use crate::funnel::{run_funnel, FunnelCandidate, FunnelConfig, FunnelOutcome, OfftargetEvidence};

/// The 20 standard amino acids, in the order used for the composition feature.
const AA20: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

/// One candidate identified by its amino-acid sequence.
///
/// The funnel computes the safety evidence and the diversity descriptor from the
/// sequence; the caller still supplies `method_scores` (these come from upstream
/// scoring methods — docking, affinity, interface confidence — which the funnel
/// will not invent).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SeqCandidate {
    /// Stable candidate identifier.
    pub id: String,
    /// The candidate's amino-acid sequence (one-letter codes).
    pub sequence: String,
    /// One score per orthogonal scoring method; every candidate must carry the
    /// same number (this is what consensus ranking aggregates).
    pub method_scores: Vec<f64>,
}

/// Reference panel and tunables for sequence-driven screening.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScreenConfig {
    /// `(id, sequence)` off-target reference panel — known proteins the candidate
    /// must not cross-react with. Empty means the off-target screen is not run.
    pub offtarget_references: Vec<(String, String)>,
    /// Per-peptide score at/above which a 9-mer window counts as a predicted
    /// T-cell epitope (illustrative; depends on the matrix scale).
    pub immuno_score_threshold: f64,
    /// Sliding-window length for linear B-cell epitope propensity.
    pub bcell_window: usize,
    /// Smoothed-propensity threshold above which a window is called an epitope.
    pub bcell_threshold: f64,
}

impl ScreenConfig {
    /// A config for `offtarget_references` with illustrative default tunables
    /// (`immuno_score_threshold = 1.0`, `bcell_window = 7`, `bcell_threshold =
    /// 0.0`). Override the fields directly to tune.
    pub fn new(offtarget_references: Vec<(String, String)>) -> Self {
        Self {
            offtarget_references,
            immuno_score_threshold: 1.0,
            bcell_window: 7,
            bcell_threshold: 0.0,
        }
    }
}

/// Amino-acid composition: the fraction of each of the 20 standard residues, in
/// [`AA20`] order. Non-standard characters count toward the length but not toward
/// any bucket, so the vector sums to (standard residues / total length).
fn aa_composition(seq: &str) -> Vec<f64> {
    let mut counts = [0u32; 20];
    let mut total = 0u32;
    for b in seq.bytes() {
        let up = b.to_ascii_uppercase();
        if let Some(i) = AA20.iter().position(|&a| a == up) {
            counts[i] += 1;
        }
        total += 1;
    }
    let denom = if total == 0 { 1.0 } else { f64::from(total) };
    counts.iter().map(|&c| f64::from(c) / denom).collect()
}

/// Run the design funnel directly from candidate **sequences**.
///
/// For each candidate this computes:
/// - **off-target**: the worst ungapped identity over `screens.offtarget_references`
///   ([`valenx_offtarget`]); skipped (recorded as not run) if the panel is empty;
/// - **immunogenicity**: predicted T-cell epitope density ([`valenx_immuno`]),
///   skipped if the sequence is shorter than the epitope matrix;
/// - **developability**: liability flags ([`valenx_developability`]);
/// - **B-cell epitopes**: count of predicted linear epitope regions
///   ([`valenx_epitope_map`]);
///
/// then derives a diversity descriptor from amino-acid composition and funnels
/// through [`run_funnel`]. CRISPR off-target is not applicable to a protein
/// binder, so it is recorded as not run.
///
/// # Errors
///
/// Returns [`OrchestratorError`] if `candidates` is empty, a candidate sequence
/// is empty, any screen errors, or the downstream funnel errors.
pub fn run_funnel_seqs(
    goal: &str,
    candidates: &[SeqCandidate],
    screens: &ScreenConfig,
    config: &FunnelConfig,
) -> Result<FunnelOutcome, OrchestratorError> {
    if candidates.is_empty() {
        return Err(OrchestratorError::Empty { what: "candidates" });
    }

    let pssm = illustrative_hla_a0201();
    let scale = hydrophilicity_kd();

    let mut funnel_candidates = Vec::with_capacity(candidates.len());
    for c in candidates {
        if c.sequence.trim().is_empty() {
            return Err(OrchestratorError::Empty {
                what: "candidate sequence",
            });
        }
        let seq = c.sequence.as_str();

        // Off-target: worst ungapped identity over the reference panel.
        let offtarget = if screens.offtarget_references.is_empty() {
            None
        } else {
            let mut best_id = -1.0_f64;
            let mut best_ref = "";
            for (rid, rseq) in &screens.offtarget_references {
                let id = best_ungapped_identity(seq, rseq)?;
                if id > best_id {
                    best_id = id;
                    best_ref = rid.as_str();
                }
            }
            Some(OfftargetEvidence {
                reference: best_ref.to_string(),
                identity: best_id,
            })
        };

        // T-cell epitope density (needs a sequence at least as long as the matrix).
        let immunogenicity = if seq.len() >= pssm.length() {
            Some(epitope_density(&pssm, seq, screens.immuno_score_threshold)?)
        } else {
            None
        };

        // Developability liabilities.
        let developability_flags = assess(seq)?.flags;

        // Linear B-cell epitope regions (window clamped to the sequence length).
        let window = screens.bcell_window.clamp(1, seq.len());
        let bcell_epitope_regions =
            Some(linear_epitope_regions(seq, &scale, window, screens.bcell_threshold)?.len());

        funnel_candidates.push(FunnelCandidate {
            id: c.id.clone(),
            method_scores: c.method_scores.clone(),
            features: aa_composition(seq),
            calibrated_confidence: None,
            offtarget,
            immunogenicity,
            crispr_offtarget_sites: None, // not applicable to a protein binder
            developability_flags,
            bcell_epitope_regions,
        });
    }

    run_funnel(goal, &funnel_candidates, config)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::GatedStage;
    use valenx_dossier::CalibrationStatus;

    fn funnel_config(top_n: usize) -> FunnelConfig {
        FunnelConfig {
            top_n,
            diversity_radius: 0.05,
            offtarget_threshold: 0.8,
            immunogenicity_threshold: 0.25,
            crispr_threshold: 3,
            calibration: CalibrationStatus::blocked("no held-out ground truth"),
            blocked_stages: vec![GatedStage::Generate, GatedStage::Dock, GatedStage::Score],
        }
    }

    #[test]
    fn composition_is_20d_and_sums_to_one() {
        let f = aa_composition("ACDEFGHIKL");
        assert_eq!(f.len(), 20);
        let s: f64 = f.iter().sum();
        assert!(
            (s - 1.0).abs() < 1e-9,
            "composition should sum to 1, got {s}"
        );
    }

    #[test]
    fn empty_candidates_is_error() {
        let err =
            run_funnel_seqs("g", &[], &ScreenConfig::new(vec![]), &funnel_config(2)).unwrap_err();
        assert_eq!(err.code(), "empty");
    }

    #[test]
    fn empty_sequence_is_error() {
        let cands = vec![SeqCandidate {
            id: "x".into(),
            sequence: "   ".into(),
            method_scores: vec![1.0],
        }];
        let err = run_funnel_seqs("g", &cands, &ScreenConfig::new(vec![]), &funnel_config(1))
            .unwrap_err();
        assert_eq!(err.code(), "empty");
    }

    #[test]
    fn identical_to_reference_flags_critical_offtarget() {
        let refs = vec![("self-target".to_string(), "ACDEFGHIKLMNPQRST".to_string())];
        let cands = vec![
            SeqCandidate {
                id: "exact".into(),
                sequence: "ACDEFGHIKLMNPQRST".into(), // identical to the reference
                method_scores: vec![5.0],
            },
            SeqCandidate {
                id: "distinct".into(),
                sequence: "WYWYWYWYWYWYWYWYW".into(), // shares almost nothing
                method_scores: vec![3.0],
            },
        ];
        let out = run_funnel_seqs(
            "inhibit X",
            &cands,
            &ScreenConfig::new(refs),
            &funnel_config(2),
        )
        .unwrap();
        // 'exact' is 100% identical to the reference -> Critical off-target.
        let exact = out
            .reports
            .iter()
            .find(|r| r.candidate_id == "exact")
            .unwrap();
        assert_eq!(
            exact.aggregate_severity(),
            Some(valenx_safety::Severity::Critical)
        );
        // The funnel never approves.
        assert!(out.requires_human_signoff());
    }

    #[test]
    fn runs_all_screens_ranks_and_fingerprints() {
        let refs = vec![("ref".to_string(), "ACDEFGHIKL".to_string())];
        let cands = vec![
            SeqCandidate {
                id: "a".into(),
                sequence: "ACDEFGHIKLMNPQRSTVWY".into(),
                method_scores: vec![9.0],
            },
            SeqCandidate {
                id: "b".into(),
                sequence: "MNPQRSTVWYMNPQRSTVWY".into(),
                method_scores: vec![1.0],
            },
        ];
        let out =
            run_funnel_seqs("g", &cands, &ScreenConfig::new(refs), &funnel_config(2)).unwrap();
        assert_eq!(out.reports.len(), 2);
        // Every report carries flags (at minimum the CRISPR 'not run' record).
        assert!(out.reports.iter().all(|r| r.has_flags()));
        // The CRISPR screen is recorded as not run, never as a clean pass.
        let a = out.reports.iter().find(|r| r.candidate_id == "a").unwrap();
        assert!(a
            .flags
            .iter()
            .any(|f| f.source == "crispr-off-target" && f.detail.contains("not run")));
        // The assembled dossier is fingerprintable.
        assert!(!out.fingerprint().unwrap().is_empty());
    }
}
