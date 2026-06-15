//! Demo: run the design funnel end-to-end on **synthetic** candidates.
//!
//! The scores, features and screen results below are illustrative inputs — they
//! are NOT real predictions. The point is to show the selection → safety →
//! dossier pipeline composes into one command, blocks gated stages honestly, and
//! never approves a candidate.
//!
//! Run with: `cargo run -p valenx-orchestrator`

use valenx_dossier::CalibrationStatus;
use valenx_orchestrator::{
    run_funnel, FunnelCandidate, FunnelConfig, GatedStage, OfftargetEvidence,
};

fn main() {
    // Six synthetic candidates with two orthogonal method scores and a 2-D
    // feature vector each. Two carry off-target evidence; the rest have screens
    // that were "not run" (recorded as such, never as a clean pass).
    let candidates = vec![
        FunnelCandidate {
            id: "design_01".into(),
            method_scores: vec![9.1, 8.7],
            features: vec![0.10, 0.90],
            calibrated_confidence: None,
            offtarget: Some(OfftargetEvidence {
                reference: "GDF-11".into(),
                identity: 0.899, // High band
            }),
            immunogenicity: Some(0.05),
            crispr_offtarget_sites: Some(0),
        },
        FunnelCandidate {
            id: "design_02".into(),
            method_scores: vec![8.4, 9.0],
            features: vec![0.85, 0.20],
            calibrated_confidence: None,
            offtarget: Some(OfftargetEvidence {
                reference: "activin-A".into(),
                identity: 0.183, // below threshold -> no flag
            }),
            immunogenicity: Some(0.32), // Moderate band
            crispr_offtarget_sites: Some(2),
        },
        FunnelCandidate {
            id: "design_03".into(),
            method_scores: vec![7.0, 7.4],
            features: vec![0.50, 0.55],
            calibrated_confidence: None,
            offtarget: None, // screen not run
            immunogenicity: None,
            crispr_offtarget_sites: None,
        },
        FunnelCandidate {
            id: "design_04".into(),
            method_scores: vec![6.2, 5.9],
            features: vec![0.20, 0.15],
            calibrated_confidence: None,
            offtarget: None,
            immunogenicity: None,
            crispr_offtarget_sites: None,
        },
        FunnelCandidate {
            id: "design_05".into(),
            method_scores: vec![5.5, 6.8],
            features: vec![0.95, 0.90],
            calibrated_confidence: None,
            offtarget: None,
            immunogenicity: None,
            crispr_offtarget_sites: None,
        },
        FunnelCandidate {
            id: "design_06".into(),
            method_scores: vec![4.0, 4.2],
            features: vec![0.45, 0.50],
            calibrated_confidence: None,
            offtarget: None,
            immunogenicity: None,
            crispr_offtarget_sites: None,
        },
    ];

    let config = FunnelConfig {
        top_n: 4,
        diversity_radius: 0.25,
        offtarget_threshold: 0.80,
        immunogenicity_threshold: 0.25,
        crispr_threshold: 3,
        // No held-out ground truth here -> calibration is honestly BLOCKED.
        calibration: CalibrationStatus::blocked("no held-out experimental ground truth"),
        // De-novo design, docking and physics affinity all need gated resources.
        blocked_stages: vec![GatedStage::Generate, GatedStage::Dock, GatedStage::Score],
    };

    let outcome = match run_funnel(
        "inhibit myostatin (GDF-8) — synthetic demo",
        &candidates,
        &config,
    ) {
        Ok(o) => o,
        Err(e) => {
            eprintln!("funnel failed [{}]: {e}", e.code());
            std::process::exit(1);
        }
    };

    println!("{}", outcome.render());

    match outcome.fingerprint() {
        Ok(fp) => println!("dossier SHA-256 fingerprint: {fp}"),
        Err(e) => eprintln!("fingerprint failed [{}]: {e}", e.code()),
    }
}
