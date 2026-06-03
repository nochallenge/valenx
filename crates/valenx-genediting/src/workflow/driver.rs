//! Feature 30 — top-level drivers, [`EditingReport`], and the
//! LLM / MCP-controllable surface.
//!
//! This module is the crate's front door. It provides:
//!
//! - a single typed [`EditingRequest`] → [`EditingReport`] driver
//!   ([`run_editing_design`]) that dispatches to the right Group-A/B/C
//!   designer for a desired change;
//! - a unified [`GeneditingRequest`] / [`GeneditingResponse`] envelope
//!   — a **typed request / response surface an external LLM (over an
//!   MCP tool) can drive** without touching the internal modules. It
//!   covers all three workflows (the strategy advisor, an editing
//!   design, and a full mRNA-therapeutic design via the
//!   [`crate::mrna::design`] driver). Every variant is
//!   `serde`-serialisable, so the whole crate is reachable as
//!   structured JSON-shaped data.
//!
//! ## v1 scope
//!
//! The driver chooses *one* design per request (the advisor's
//! recommendation, with a prime-editing fallback). It is a dispatcher,
//! not an exhaustive search; callers that want every guide call the
//! per-module designers directly. All scores are the transparent
//! heuristics of the source modules.

use crate::base_edit::design::BaseEditDesign;
use crate::crispr::knockout::{design_knockout, GeneModel, KnockoutRequest, KnockoutStrategy};
use crate::crispr::nuclease::NucleaseId;
use crate::error::{GeneditingError, Result};
use crate::mrna::design::{design_mrna, MrnaDesignReport, MrnaDesignRequest};
use crate::prime_edit::pegrna::PegRna;
use crate::prime_edit::strategy::{model_strategy, StrategyModel};
use crate::workflow::advisor::{advise_strategy, DesiredChange, EditApproach, StrategyAdvice};
use crate::workflow::variant::{
    plan_variant_correction, CorrectionReagent, PathogenicVariant, VariantCorrectionPlan,
    VariantCorrectionRequest,
};
use serde::{Deserialize, Serialize};

/// What an [`EditingRequest`] asks the driver to do.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum EditingTask {
    /// Correct a pathogenic variant — design the reverting reagent.
    CorrectVariant {
        /// The wild-type reference window around the variant.
        reference: Vec<u8>,
        /// 0-based variant position in the reference.
        variant_pos: usize,
        /// The pathogenic variant.
        variant: PathogenicVariant,
    },
    /// Knock out a gene — design an NHEJ frameshift strategy.
    KnockoutGene {
        /// The coding sequence of the gene.
        cds: Vec<u8>,
        /// Exon-boundary offsets (0-based, the first must be `0`).
        exon_boundaries: Vec<usize>,
        /// Which nuclease to use.
        nuclease: NucleaseId,
    },
}

/// A typed request to the top-level editing driver.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EditingRequest {
    /// The task to run.
    pub task: EditingTask,
}

/// The designed reagent an [`EditingReport`] carries.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum DesignedReagent {
    /// A base-editing guide.
    BaseEdit(BaseEditDesign),
    /// A prime-editing pegRNA plus its strategy model.
    Prime {
        /// The designed pegRNA.
        pegrna: PegRna,
        /// The strategy model (efficiency heuristic + rationale).
        strategy: StrategyModel,
    },
    /// An NHEJ knockout strategy — a ranked set of knockout guides.
    Knockout(KnockoutStrategy),
}

impl DesignedReagent {
    /// The best single guide / pegRNA spacer sequence, for a quick
    /// summary. Returns an empty string for an empty knockout strategy.
    pub fn primary_spacer(&self) -> String {
        match self {
            DesignedReagent::BaseEdit(d) => d.protospacer.clone(),
            DesignedReagent::Prime { pegrna, .. } => {
                String::from_utf8_lossy(&pegrna.spacer).into_owned()
            }
            DesignedReagent::Knockout(k) => k
                .best()
                .map(|g| g.guide.protospacer.clone())
                .unwrap_or_default(),
        }
    }
}

/// The bundled result of a top-level editing-design run (feature 30).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EditingReport {
    /// The editing modality the driver chose.
    pub approach: EditApproach,
    /// The designed reagent.
    pub reagent: DesignedReagent,
    /// A human-readable explanation of the design.
    pub rationale: String,
    /// Advisory notes (one line per consideration).
    pub notes: Vec<String>,
}

impl EditingReport {
    /// `true` when the chosen approach makes a double-strand break.
    pub fn makes_dsb(&self) -> bool {
        self.approach.makes_dsb()
    }
}

/// Runs a top-level editing design (feature 30).
///
/// Dispatches on [`EditingTask`]: a `CorrectVariant` task is routed
/// through the [`crate::workflow::variant`] planner; a `KnockoutGene`
/// task through the [`crate::crispr::knockout`] designer. The result is
/// a uniform [`EditingReport`] regardless of which modality was used.
///
/// # Errors
/// Propagates [`GeneditingError`] from the underlying designer — most
/// often [`GeneditingError::NoValidDesign`] when no reagent reaches the
/// target.
pub fn run_editing_design(req: &EditingRequest) -> Result<EditingReport> {
    match &req.task {
        EditingTask::CorrectVariant {
            reference,
            variant_pos,
            variant,
        } => {
            let plan = plan_variant_correction(&VariantCorrectionRequest {
                reference: reference.clone(),
                variant_pos: *variant_pos,
                variant: variant.clone(),
            })?;
            Ok(variant_plan_to_report(plan))
        }
        EditingTask::KnockoutGene {
            cds,
            exon_boundaries,
            nuclease,
        } => {
            let gene = GeneModel::from_cds(cds.clone(), exon_boundaries)?;
            let strategy = design_knockout(&KnockoutRequest::new(gene, *nuclease))?;
            let mut notes = vec![strategy.rationale.clone()];
            notes.push(
                "An NHEJ knockout relies on a frameshift; verify the indel \
                 spectrum experimentally (a CRISPResso-class analysis)."
                    .to_string(),
            );
            let rationale = strategy.rationale.clone();
            Ok(EditingReport {
                approach: EditApproach::NucleaseNhej,
                reagent: DesignedReagent::Knockout(strategy),
                rationale,
                notes,
            })
        }
    }
}

/// Converts a [`VariantCorrectionPlan`] into an [`EditingReport`].
fn variant_plan_to_report(plan: VariantCorrectionPlan) -> EditingReport {
    let mut notes = Vec::new();
    let reagent = match plan.reagent {
        CorrectionReagent::BaseEdit { design, .. } => {
            notes.push(format!(
                "Base-editing product purity (heuristic): {:.2}; {} bystander \
                 editable base(s) in the window.",
                design.product_purity,
                design.window.bystander_count(),
            ));
            DesignedReagent::BaseEdit(design)
        }
        CorrectionReagent::Prime { pegrna, editor } => {
            let strategy = model_strategy(&pegrna, editor);
            notes.push(strategy.rationale.clone());
            DesignedReagent::Prime { pegrna, strategy }
        }
    };
    EditingReport {
        approach: plan.approach,
        reagent,
        rationale: plan.rationale,
        notes,
    }
}

// --- The LLM / MCP-controllable surface -------------------------------

/// A unified, `serde`-serialisable request envelope for an external
/// LLM driving the crate over an MCP tool (feature 30).
///
/// Every editing and mRNA workflow is reachable through this one type
/// — an LLM emits a [`GeneditingRequest`] as structured data and reads
/// back a [`GeneditingResponse`], with no access to the internal
/// module APIs.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum GeneditingRequest {
    /// Recommend an editing strategy for a desired change (the
    /// advisor) — no sequence design, just the modality recommendation.
    AdviseStrategy {
        /// The desired genomic change.
        change: DesiredChange,
        /// Whether a PAM-adjacent guide is known to exist at the locus.
        pam_available: bool,
    },
    /// Run a full editing design (variant correction or knockout).
    DesignEdit(EditingRequest),
    /// Run a full mRNA-therapeutic design.
    DesignMrna(MrnaDesignRequest),
}

/// The matching `serde`-serialisable response envelope (feature 30).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum GeneditingResponse {
    /// A strategy recommendation.
    Strategy(StrategyAdvice),
    /// A completed editing design.
    Edit(EditingReport),
    /// A completed mRNA design.
    Mrna(MrnaDesignReport),
    /// The request failed; carries the stable error code and a
    /// human-readable message.
    Error {
        /// The [`GeneditingError::code`] of the failure.
        code: String,
        /// A human-readable failure message.
        message: String,
    },
}

impl GeneditingResponse {
    /// `true` when the response is *not* an [`GeneditingResponse::Error`].
    pub fn is_ok(&self) -> bool {
        !matches!(self, GeneditingResponse::Error { .. })
    }
}

/// The single entry point an external LLM / MCP tool calls (feature
/// 30).
///
/// Dispatches a [`GeneditingRequest`] to the matching workflow and
/// always returns a [`GeneditingResponse`] — errors are *captured into*
/// [`GeneditingResponse::Error`] rather than propagated, so the caller
/// (the LLM) always receives a well-formed, serialisable answer.
pub fn handle_request(req: &GeneditingRequest) -> GeneditingResponse {
    match req {
        GeneditingRequest::AdviseStrategy {
            change,
            pam_available,
        } => GeneditingResponse::Strategy(advise_strategy(change, *pam_available)),
        GeneditingRequest::DesignEdit(edit_req) => match run_editing_design(edit_req) {
            Ok(report) => GeneditingResponse::Edit(report),
            Err(e) => error_response(&e),
        },
        GeneditingRequest::DesignMrna(mrna_req) => match design_mrna(mrna_req) {
            Ok(report) => GeneditingResponse::Mrna(report),
            Err(e) => error_response(&e),
        },
    }
}

/// Builds a [`GeneditingResponse::Error`] from a [`GeneditingError`].
fn error_response(e: &GeneditingError) -> GeneditingResponse {
    GeneditingResponse::Error {
        code: e.code().to_string(),
        message: e.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mrna::tailcap::MrnaUseCase;

    fn knockout_cds() -> Vec<u8> {
        let mut s = Vec::new();
        for _ in 0..6 {
            s.extend_from_slice(b"ACGTACGTACGTACGTACGTAGG");
        }
        s
    }

    fn variant_reference() -> Vec<u8> {
        b"ACGTACGTACGTACGTACGTAGGCATGCATGCATGCATGCATGCATGC".to_vec()
    }

    #[test]
    fn drives_a_knockout_design() {
        let req = EditingRequest {
            task: EditingTask::KnockoutGene {
                cds: knockout_cds(),
                exon_boundaries: vec![0, 69],
                nuclease: NucleaseId::SpCas9,
            },
        };
        let report = run_editing_design(&req).unwrap();
        assert_eq!(report.approach, EditApproach::NucleaseNhej);
        assert!(report.makes_dsb());
        assert!(matches!(report.reagent, DesignedReagent::Knockout(_)));
        assert!(!report.reagent.primary_spacer().is_empty());
    }

    #[test]
    fn drives_a_variant_correction() {
        // variant_reference() index 25 is `T` (…A G G C A[24] T[25] G…),
        // so the pathogenic substitution declares ref_base = T.
        let req = EditingRequest {
            task: EditingTask::CorrectVariant {
                reference: variant_reference(),
                variant_pos: 25,
                variant: PathogenicVariant::Substitution {
                    ref_base: b'T',
                    variant_base: b'A',
                },
            },
        };
        let report = run_editing_design(&req).unwrap();
        // A point correction → base editing or prime editing.
        assert!(matches!(
            report.approach,
            EditApproach::BaseEditing | EditApproach::PrimeEditing
        ));
        assert!(!report.notes.is_empty());
    }

    #[test]
    fn knockout_with_bad_boundaries_errors() {
        let req = EditingRequest {
            task: EditingTask::KnockoutGene {
                cds: knockout_cds(),
                exon_boundaries: vec![5, 10], // must start at 0
                nuclease: NucleaseId::SpCas9,
            },
        };
        assert!(run_editing_design(&req).is_err());
    }

    #[test]
    fn llm_surface_advises_a_strategy() {
        let req = GeneditingRequest::AdviseStrategy {
            change: DesiredChange::GeneKnockout,
            pam_available: true,
        };
        let resp = handle_request(&req);
        assert!(resp.is_ok());
        match resp {
            GeneditingResponse::Strategy(a) => {
                assert_eq!(a.recommended(), EditApproach::NucleaseNhej);
            }
            _ => panic!("expected a Strategy response"),
        }
    }

    #[test]
    fn llm_surface_runs_an_edit() {
        let req = GeneditingRequest::DesignEdit(EditingRequest {
            task: EditingTask::KnockoutGene {
                cds: knockout_cds(),
                exon_boundaries: vec![0, 69],
                nuclease: NucleaseId::SpCas9,
            },
        });
        let resp = handle_request(&req);
        assert!(resp.is_ok());
        assert!(matches!(resp, GeneditingResponse::Edit(_)));
    }

    #[test]
    fn llm_surface_runs_an_mrna_design() {
        let req = GeneditingRequest::DesignMrna(MrnaDesignRequest::new(
            b"ATGGCCCTGCTGGAAGAATAA".to_vec(),
            MrnaUseCase::Vaccine,
        ));
        let resp = handle_request(&req);
        assert!(resp.is_ok());
        assert!(matches!(resp, GeneditingResponse::Mrna(_)));
    }

    #[test]
    fn llm_surface_captures_errors() {
        // A knockout with bad boundaries → an Error response, not a
        // panic or a propagated Err.
        let req = GeneditingRequest::DesignEdit(EditingRequest {
            task: EditingTask::KnockoutGene {
                cds: knockout_cds(),
                exon_boundaries: vec![3], // must start at 0
                nuclease: NucleaseId::SpCas9,
            },
        });
        let resp = handle_request(&req);
        assert!(!resp.is_ok());
        match resp {
            GeneditingResponse::Error { code, .. } => {
                assert!(code.starts_with("genediting."));
            }
            _ => panic!("expected an Error response"),
        }
    }

    #[test]
    fn request_response_types_are_clone_and_eq() {
        // The LLM surface types must be `Clone`/`PartialEq` (and, via
        // the derives, `serde`-serialisable) so an MCP layer can carry
        // them as structured data.
        let req = GeneditingRequest::AdviseStrategy {
            change: DesiredChange::PointMutation { from: b'A', to: b'G' },
            pam_available: true,
        };
        assert_eq!(req.clone(), req);
        let resp = handle_request(&req);
        assert_eq!(resp.clone(), resp);
        assert!(resp.is_ok());
    }

    /// Compile-time proof the LLM-surface envelopes implement `Serialize`
    /// + `Deserialize` (so an MCP tool can ferry them as JSON).
    #[test]
    fn llm_surface_is_serde_capable() {
        fn assert_serde<T: serde::Serialize + serde::de::DeserializeOwned>() {}
        assert_serde::<GeneditingRequest>();
        assert_serde::<GeneditingResponse>();
        assert_serde::<EditingReport>();
    }
}
