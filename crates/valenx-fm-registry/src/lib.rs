//! # valenx-fm-registry
//!
//! A registry of the **gated frontier foundation models** the design pipeline
//! would call when the resources exist — and an honest readiness probe that
//! reports `BLOCKED` when they do not.
//!
//! The valenx funnel ([`valenx-orchestrator`]) runs its selection → safety →
//! dossier core entirely in-house. The high-end *generate* / *fold* / *perturb*
//! stages, by contrast, depend on large neural models that need a **GPU and
//! downloaded weights** (and sometimes a licensed runtime). This crate records,
//! for each such model, **verified public metadata** — task, license, upstream —
//! plus the exact gated dependency it needs, and a [`probe`] that checks whether
//! local weights have been staged.
//!
//! **This crate never runs inference and never fabricates a result.** In any
//! environment without the gated resources (CI, this repo) every probe returns
//! [`ProbeStatus::Blocked`]. Even when weights are staged, [`probe`] only reports
//! readiness — actually running a model is the upstream tool's job, out of scope
//! here. This is the connective-tissue placeholder for the pipeline's
//! resource-gated frontier, kept honest.
//!
//! ## Models
//!
//! - **Evo 2** (Arc Institute) — genome-scale DNA/RNA/protein design;
//!   NVIDIA Open Model License; <https://arcinstitute.org/tools/evo>
//! - **Boltz-2** (MIT / Recursion) — joint structure + binding-affinity;
//!   MIT license; <https://boltz.bio/boltz2>
//! - **scGPT** (CZI Virtual Cells) — single-cell perturbation prediction;
//!   MIT license; <https://virtualcellmodels.cziscience.com/model/scgpt>
//!
//! ## Example
//!
//! ```
//! use valenx_fm_registry::{registry, probe, ProbeStatus};
//!
//! // With no weights staged (the default), every model is BLOCKED — never faked.
//! for m in registry() {
//!     let status = probe(&m);
//!     assert!(status.is_blocked(), "{} should be blocked without weights", m.id);
//!     println!("{}", status.message(&m)); // "BLOCKED: … — needs GPU + weights …"
//! }
//! ```
//!
//! [`valenx-orchestrator`]: https://docs.rs/valenx-orchestrator

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use std::path::{Path, PathBuf};

use serde::Serialize;

/// What a foundation model does in the pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum FmTask {
    /// Genome-scale DNA / RNA / protein modelling and design.
    GenomeDesign,
    /// Joint 3-D structure and binding-affinity prediction.
    StructureAffinity,
    /// Single-cell transcriptomics and genetic-perturbation prediction.
    SingleCellPerturbation,
}

impl FmTask {
    /// A short label for the task.
    pub fn as_str(self) -> &'static str {
        match self {
            FmTask::GenomeDesign => "genome-design",
            FmTask::StructureAffinity => "structure-affinity",
            FmTask::SingleCellPerturbation => "single-cell-perturbation",
        }
    }
}

/// A gated foundation model the pipeline could call. Every field is **verified
/// public metadata** (see `upstream_url`); this struct carries no weights and
/// runs no inference.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub struct FmModel {
    /// Stable short identifier (e.g. `"evo2"`).
    pub id: &'static str,
    /// Human-readable name and source.
    pub display_name: &'static str,
    /// What the model does.
    pub task: FmTask,
    /// The model's published license.
    pub license: &'static str,
    /// Canonical upstream URL.
    pub upstream_url: &'static str,
    /// The gated dependency required to actually run the model.
    pub needs: &'static str,
    /// Environment variable that, if set to a readable path, points at
    /// locally-staged model weights.
    pub weights_env: &'static str,
}

/// The registry of known gated foundation models, with verified metadata.
pub fn registry() -> Vec<FmModel> {
    vec![
        FmModel {
            id: "evo2",
            display_name: "Evo 2 (Arc Institute)",
            task: FmTask::GenomeDesign,
            license: "NVIDIA Open Model License",
            upstream_url: "https://arcinstitute.org/tools/evo",
            needs: "GPU (multi-H100 class) + downloaded weights (7B/40B) via NVIDIA BioNeMo",
            weights_env: "VALENX_EVO2_WEIGHTS",
        },
        FmModel {
            id: "boltz2",
            display_name: "Boltz-2 (MIT / Recursion)",
            task: FmTask::StructureAffinity,
            license: "MIT",
            upstream_url: "https://boltz.bio/boltz2",
            needs: "GPU + downloaded weights + the Boltz Python runtime",
            weights_env: "VALENX_BOLTZ2_WEIGHTS",
        },
        FmModel {
            id: "scgpt",
            display_name: "scGPT (CZI Virtual Cells)",
            task: FmTask::SingleCellPerturbation,
            license: "MIT",
            upstream_url: "https://virtualcellmodels.cziscience.com/model/scgpt",
            needs: "GPU + downloaded weights + a single-cell expression matrix (AnnData)",
            weights_env: "VALENX_SCGPT_WEIGHTS",
        },
    ]
}

/// Look up a model by its `id`, or `None` if unknown.
pub fn model(id: &str) -> Option<FmModel> {
    registry().into_iter().find(|m| m.id == id)
}

/// The readiness of a model in the current environment.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub enum ProbeStatus {
    /// Weights are staged at `weights`. Running still needs the model's GPU
    /// runtime, so this crate will not run inference even now — it only reports
    /// that the weights are present.
    WeightsStaged {
        /// Path to the staged weights.
        weights: PathBuf,
    },
    /// **BLOCKED**: the gated dependency is absent. Carries the dependency text.
    Blocked {
        /// What is missing (GPU, weights, licensed runtime, …).
        dependency: String,
    },
}

impl ProbeStatus {
    /// Whether the model is blocked (cannot run here).
    pub fn is_blocked(&self) -> bool {
        matches!(self, ProbeStatus::Blocked { .. })
    }

    /// A human-readable status line for `model`.
    pub fn message(&self, model: &FmModel) -> String {
        match self {
            ProbeStatus::Blocked { dependency } => {
                format!("BLOCKED: {} — needs {dependency}", model.display_name)
            }
            ProbeStatus::WeightsStaged { weights } => format!(
                "weights staged for {} at {} (running still needs the upstream GPU runtime; \
                 not run here)",
                model.display_name,
                weights.display()
            ),
        }
    }
}

/// Probe a model's readiness. If its `weights_env` variable is set to a path
/// that exists, the result is [`ProbeStatus::WeightsStaged`]; otherwise — the
/// default in any environment without the gated resources — it is
/// [`ProbeStatus::Blocked`]. This function **never runs the model and never
/// returns a prediction.**
pub fn probe(model: &FmModel) -> ProbeStatus {
    match std::env::var_os(model.weights_env) {
        Some(p) if !p.is_empty() && Path::new(&p).exists() => ProbeStatus::WeightsStaged {
            weights: PathBuf::from(p),
        },
        _ => ProbeStatus::Blocked {
            dependency: model.needs.to_string(),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn registry_has_three_unique_models() {
        let reg = registry();
        assert_eq!(reg.len(), 3);
        let mut ids: Vec<&str> = reg.iter().map(|m| m.id).collect();
        ids.sort_unstable();
        ids.dedup();
        assert_eq!(ids.len(), 3, "model ids must be unique");
    }

    #[test]
    fn lookup_by_id() {
        assert_eq!(model("evo2").unwrap().task, FmTask::GenomeDesign);
        assert_eq!(model("boltz2").unwrap().task, FmTask::StructureAffinity);
        assert_eq!(model("scgpt").unwrap().task, FmTask::SingleCellPerturbation);
        assert!(model("does-not-exist").is_none());
    }

    #[test]
    fn every_model_is_blocked_without_staged_weights() {
        // The honest default: no weights env set anywhere -> BLOCKED, never faked.
        for m in registry() {
            let status = probe(&m);
            assert!(status.is_blocked(), "{} should be blocked", m.id);
            let msg = status.message(&m);
            assert!(msg.starts_with("BLOCKED:"), "got: {msg}");
            assert!(msg.contains("GPU"), "dependency should mention GPU: {msg}");
        }
    }

    #[test]
    fn weights_staged_status_is_not_blocked() {
        // Construct the staged variant directly (no global env mutation).
        let m = model("boltz2").unwrap();
        let status = ProbeStatus::WeightsStaged {
            weights: PathBuf::from("/tmp/boltz2.ckpt"),
        };
        assert!(!status.is_blocked());
        let msg = status.message(&m);
        assert!(msg.contains("weights staged"));
        assert!(msg.contains("not run here"));
    }
}
