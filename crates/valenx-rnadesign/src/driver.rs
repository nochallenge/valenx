//! Group F — the top-level driver (features 23–25).
//!
//! This module is the crate's front door. It provides:
//!
//! - [`design_rna`] (feature 23) — the single top-level driver: it
//!   takes a [`DesignGoal`] and a [`DesignConstraints`] set, runs the
//!   whole Goal → Design → Optimize → Validate → Export workflow, and
//!   returns a [`RnaDesignReport`] (the design, the validation, the
//!   synthesis plan).
//! - [`design_rna_batch`] (feature 24) — a batch mode: design several
//!   candidates (from different random seeds) and rank them by quality.
//! - [`RnaDesignRequest`] / [`RnaDesignResponse`] and
//!   [`handle_request`] (feature 25) — a typed, `serde`-serialisable
//!   request / response surface an external LLM can drive over an MCP
//!   tool, with no access to the internal modules.
//!
//! ## v1 scope — honest framing
//!
//! [`design_rna`] runs one design per goal and one optimisation /
//! validation pass; [`design_rna_batch`] is independent restarts, not a
//! population search. Every [`RnaDesignReport`] carries the same honest
//! disclaimer the [`crate::validate::ValidationReport`] does: the output
//! is a strong in-silico candidate, not a guarantee — see the crate
//! root.

use crate::design::{
    design_coding, design_riboswitch, design_structural, motif_scaffold, CodingDesignParams,
    DesignKind, RnaDesign, StructuralDesignParams, TwoStateParams,
};
use crate::error::{Result, RnaDesignError};
use crate::export::export_design_report;
use crate::goal::{DesignConstraints, DesignGoal};
use crate::optimize::{optimize_design, OptimizationResult, OptimizeParams};
use crate::synthesis::{plan_synthesis, Promoter, SynthesisPackage};
use crate::validate::{validate_design, ValidationReport, ValidationVerdict};
use crate::workflow::{DesignSession, WorkflowStage};
use serde::{Deserialize, Serialize};

/// The bundled result of a top-level RNA-design run (feature 23).
///
/// # Honest framing
///
/// An [`RnaDesignReport`] is a *predicted* design — a strong candidate
/// validated against a nearest-neighbor energy model and rule-based
/// heuristics, NOT a guarantee of correct in-vivo behaviour. The
/// physical RNA must be synthesised and lab-validated. See the crate
/// root and [`crate::validate::ValidationReport`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RnaDesignReport {
    /// The final (optimised) RNA design candidate.
    pub design: RnaDesign,
    /// The optimisation result.
    pub optimization: OptimizationResult,
    /// The in-silico validation report.
    pub validation: ValidationReport,
    /// The synthesis-ready output package.
    pub synthesis: SynthesisPackage,
    /// A human-readable plain-text design report (the
    /// [`crate::export::export_design_report`] output).
    pub text_report: String,
    /// Top-level workflow notes (one line per stage).
    pub notes: Vec<String>,
}

impl RnaDesignReport {
    /// `true` when the design cleared in-silico validation
    /// ([`ValidationVerdict::Pass`]).
    ///
    /// Note: a `true` here means the *in-silico* checks passed — it is
    /// not a wet-lab result.
    pub fn validation_passed(&self) -> bool {
        self.validation.verdict == ValidationVerdict::Pass
    }

    /// A scalar **quality score** in `[0, 1]` for ranking candidates —
    /// a transparent blend of the validation verdict, the
    /// structure-match / ensemble metrics and the optimisation
    /// improvement. Higher is better.
    pub fn quality_score(&self) -> f64 {
        // Verdict component.
        let verdict_score = match self.validation.verdict {
            ValidationVerdict::Pass => 1.0,
            ValidationVerdict::Warn => 0.6,
            ValidationVerdict::Fail => 0.2,
        };
        // Structure component (if the design has a target).
        let structure_score = match (&self.validation.fold_back, &self.validation.ensemble) {
            (Some(fb), Some(en)) => {
                0.5 * (fb.structure_match_percent / 100.0)
                    + 0.5 * en.target_probability.clamp(0.0, 1.0)
            }
            _ => verdict_score, // a non-structural design leans on the verdict
        };
        // Pass fraction of the constraint checks.
        let total = self.validation.constraint_checks.len().max(1);
        let pass_fraction =
            self.validation.count_with(ValidationVerdict::Pass) as f64 / total as f64;

        (0.45 * verdict_score + 0.35 * structure_score + 0.20 * pass_fraction).clamp(0.0, 1.0)
    }
}

/// Designs an RNA candidate from a goal, running the whole workflow
/// (feature 23).
///
/// Drives a [`DesignSession`] through Goal → Design → Optimize →
/// Validate → Export and returns the bundled [`RnaDesignReport`]. The
/// designer is chosen from the [`DesignGoal`] shape; the synthesis plan
/// uses a T7 promoter.
///
/// # Errors
/// - [`RnaDesignError::Goal`] if the goal or constraints are invalid.
/// - [`RnaDesignError::NoDesign`] if no candidate can be produced.
/// - [`RnaDesignError::Upstream`] if a building-block call fails.
pub fn design_rna(goal: &DesignGoal, constraints: &DesignConstraints) -> Result<RnaDesignReport> {
    design_rna_seeded(goal, constraints, 0)
}

/// [`design_rna`] with an explicit base random seed — used by the batch
/// driver to produce distinct candidates.
///
/// # Errors
/// As [`design_rna`].
pub fn design_rna_seeded(
    goal: &DesignGoal,
    constraints: &DesignConstraints,
    seed: u64,
) -> Result<RnaDesignReport> {
    // --- Goal stage --------------------------------------------------
    let mut session = DesignSession::new(goal.clone(), constraints.clone())?;
    let mut notes = vec![format!(
        "Goal: design a {} under {} constraint(s).",
        goal.kind_name(),
        describe_constraints(constraints),
    )];

    // --- Design stage ------------------------------------------------
    let design = run_design(goal, constraints, seed)?;
    session.advance_design(design)?;
    notes.push(format!(
        "Design: produced a {}-nt candidate.",
        session.design().map(|d| d.len()).unwrap_or(0),
    ));

    // --- Optimize stage ----------------------------------------------
    let design = session.design().expect("design stage completed").clone();
    let opt_params = OptimizeParams {
        seed: seed.wrapping_add(0x0_9714),
        ..OptimizeParams::default()
    };
    let optimization = optimize_design(&design, constraints, opt_params)?;
    session.advance_optimize(optimization.clone())?;
    notes.push(format!(
        "Optimize: objective score {:.3} -> {:.3}.",
        optimization.score_before, optimization.score_after,
    ));

    // --- Validate stage ----------------------------------------------
    let optimized = optimization.design.clone();
    let validation = validate_design(&optimized, constraints)?;
    session.advance_validate(validation.clone())?;
    notes.push(format!(
        "Validate: in-silico verdict `{}`.",
        validation.verdict.name(),
    ));

    // --- Export stage ------------------------------------------------
    let synthesis = plan_synthesis(&optimized, constraints, Promoter::T7)?;
    session.advance_export(synthesis.clone())?;
    notes.push(format!(
        "Export: synthesis package ready — a {}-bp DNA template.",
        synthesis.template_len(),
    ));
    debug_assert_eq!(session.stage(), WorkflowStage::Export);

    let text_report = export_design_report(&optimized, &validation, &synthesis, goal.kind_name());

    notes.push(
        "Reminder: this report is an in-silico prediction — a strong validated candidate, \
         not a guarantee. Synthesise and lab-validate the design."
            .to_string(),
    );

    Ok(RnaDesignReport {
        design: optimized,
        optimization,
        validation,
        synthesis,
        text_report,
        notes,
    })
}

/// Runs the Design stage — dispatches on the goal shape to the right
/// designer.
fn run_design(goal: &DesignGoal, constraints: &DesignConstraints, seed: u64) -> Result<RnaDesign> {
    match goal {
        DesignGoal::Structural { target, .. } => {
            let params = StructuralDesignParams {
                seed: seed.wrapping_add(0x5EED),
                ..StructuralDesignParams::default()
            };
            design_structural(target, params)
        }
        DesignGoal::Coding { protein, use_case } => {
            let mut params = CodingDesignParams::new(*use_case);
            params.host = constraints.host;
            params.nucleoside = constraints.nucleoside;
            design_coding(protein, params)
        }
        DesignGoal::Hybrid {
            protein, use_case, ..
        } => {
            // v1: a hybrid design produces the coding mRNA; the required
            // structural elements are reported but the construct is the
            // coding construct. (A full hybrid designer that seats the
            // elements into the UTRs is future work — see the crate
            // note.)
            let mut params = CodingDesignParams::new(*use_case);
            params.host = constraints.host;
            params.nucleoside = constraints.nucleoside;
            let mut design = design_coding(protein, params)?;
            design.notes.push(
                "Hybrid goal: the coding mRNA was designed; the required structural \
                 elements are recorded but not yet auto-seated into the UTRs (v1 scope)."
                    .to_string(),
            );
            Ok(design)
        }
    }
}

/// A short human-readable summary of which constraints are active.
fn describe_constraints(c: &DesignConstraints) -> String {
    let mut parts = Vec::new();
    parts.push(format!(
        "GC {:.0}-{:.0}%",
        c.gc_min * 100.0,
        c.gc_max * 100.0
    ));
    if c.length_min != 0 || c.length_max != 0 {
        parts.push("length-bounded".to_string());
    }
    if !c.forbidden_motifs.is_empty() {
        parts.push(format!("{} forbidden motif(s)", c.forbidden_motifs.len()));
    }
    if !c.forbidden_restriction_sites.is_empty() {
        parts.push(format!(
            "{} forbidden site(s)",
            c.forbidden_restriction_sites.len()
        ));
    }
    parts.join(", ")
}

// ---------------------------------------------------------------------
// Feature 24 — batch mode
// ---------------------------------------------------------------------

/// The result of a batch RNA-design run (feature 24).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BatchDesignResult {
    /// The designed candidates, **ranked best-first** by
    /// [`RnaDesignReport::quality_score`].
    pub candidates: Vec<RnaDesignReport>,
}

impl BatchDesignResult {
    /// The best (highest-quality) candidate, or `None` if the batch is
    /// empty.
    pub fn best(&self) -> Option<&RnaDesignReport> {
        self.candidates.first()
    }

    /// The number of candidates in the batch.
    pub fn len(&self) -> usize {
        self.candidates.len()
    }

    /// `true` when the batch produced no candidate.
    pub fn is_empty(&self) -> bool {
        self.candidates.is_empty()
    }
}

/// Designs several RNA candidates and ranks them (feature 24).
///
/// Runs [`design_rna_seeded`] `count` times from distinct seeds and
/// returns the candidates sorted best-first by
/// [`RnaDesignReport::quality_score`]. A candidate whose design step
/// fails is skipped; the batch fails only if *every* candidate fails.
///
/// # Errors
/// - [`RnaDesignError::Invalid`] if `count == 0`.
/// - [`RnaDesignError::NoDesign`] if no candidate could be designed at
///   all.
pub fn design_rna_batch(
    goal: &DesignGoal,
    constraints: &DesignConstraints,
    count: usize,
) -> Result<BatchDesignResult> {
    if count == 0 {
        return Err(RnaDesignError::invalid(
            "count",
            "a batch must request at least one candidate",
        ));
    }
    // Validate the goal / constraints once up front so an invalid goal
    // fails fast rather than `count` times.
    goal.validate()?;
    constraints.validate()?;

    let mut candidates = Vec::new();
    let mut last_err: Option<RnaDesignError> = None;
    for k in 0..count {
        // Spread the seeds widely so the candidates genuinely differ.
        let seed = (k as u64).wrapping_mul(0x9E37_79B9).wrapping_add(1);
        match design_rna_seeded(goal, constraints, seed) {
            Ok(report) => candidates.push(report),
            Err(e) => last_err = Some(e),
        }
    }

    if candidates.is_empty() {
        return Err(last_err.unwrap_or_else(|| {
            RnaDesignError::no_design("batch", "no candidate could be designed")
        }));
    }

    // Rank best-first; a stable sort keeps earlier seeds ahead on a tie.
    candidates.sort_by(|a, b| {
        b.quality_score()
            .partial_cmp(&a.quality_score())
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    Ok(BatchDesignResult { candidates })
}

// ---------------------------------------------------------------------
// Feature 25 — the MCP / LLM-controllable surface
// ---------------------------------------------------------------------

/// A typed, `serde`-serialisable request envelope for an external LLM
/// driving the crate over an MCP tool (feature 25).
///
/// Every workflow this crate offers is reachable through this one type
/// — an LLM emits an [`RnaDesignRequest`] as structured data and reads
/// back an [`RnaDesignResponse`], with no access to the internal module
/// APIs.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum RnaDesignRequest {
    /// Design a single RNA candidate for a goal.
    Design {
        /// The design goal.
        goal: DesignGoal,
        /// The design constraints.
        constraints: DesignConstraints,
    },
    /// Design a batch of candidates and rank them.
    DesignBatch {
        /// The design goal.
        goal: DesignGoal,
        /// The design constraints.
        constraints: DesignConstraints,
        /// How many candidates to design.
        count: usize,
    },
    /// Return a built-in functional-motif scaffold as a starting design
    /// (no optimisation / validation — just the scaffold).
    GetScaffold {
        /// The scaffold id (see [`crate::design::motif`]).
        id: String,
    },
}

/// The matching `serde`-serialisable response envelope (feature 25).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum RnaDesignResponse {
    /// A completed single design.
    Design(Box<RnaDesignReport>),
    /// A completed batch design.
    Batch(BatchDesignResult),
    /// A functional-motif scaffold.
    Scaffold(RnaDesign),
    /// The request failed; carries the stable error code and a
    /// human-readable message.
    Error {
        /// The [`RnaDesignError::code`] of the failure.
        code: String,
        /// A human-readable failure message.
        message: String,
    },
}

impl RnaDesignResponse {
    /// `true` when the response is *not* an [`RnaDesignResponse::Error`].
    pub fn is_ok(&self) -> bool {
        !matches!(self, RnaDesignResponse::Error { .. })
    }
}

/// The single entry point an external LLM / MCP tool calls (feature 25).
///
/// Dispatches an [`RnaDesignRequest`] to the matching workflow and
/// always returns an [`RnaDesignResponse`] — errors are *captured into*
/// [`RnaDesignResponse::Error`] rather than propagated, so the caller
/// (the LLM) always receives a well-formed, serialisable answer.
pub fn handle_request(req: &RnaDesignRequest) -> RnaDesignResponse {
    match req {
        RnaDesignRequest::Design { goal, constraints } => match design_rna(goal, constraints) {
            Ok(report) => RnaDesignResponse::Design(Box::new(report)),
            Err(e) => error_response(&e),
        },
        RnaDesignRequest::DesignBatch {
            goal,
            constraints,
            count,
        } => match design_rna_batch(goal, constraints, *count) {
            Ok(batch) => RnaDesignResponse::Batch(batch),
            Err(e) => error_response(&e),
        },
        RnaDesignRequest::GetScaffold { id } => match motif_scaffold(id) {
            Ok(design) => RnaDesignResponse::Scaffold(design),
            Err(e) => error_response(&e),
        },
    }
}

/// Builds an [`RnaDesignResponse::Error`] from an [`RnaDesignError`].
fn error_response(e: &RnaDesignError) -> RnaDesignResponse {
    RnaDesignResponse::Error {
        code: e.code().to_string(),
        message: e.to_string(),
    }
}

// ---------------------------------------------------------------------
// A small two-state convenience driver (riboswitch goals)
// ---------------------------------------------------------------------

/// Designs a riboswitch end-to-end and returns a full
/// [`RnaDesignReport`] (a convenience wrapper that runs the two-state
/// designer, then the optimise / validate / export stages).
///
/// `free_dot_bracket` and `bound_dot_bracket` are the resting and
/// ligand-bound target structures.
///
/// # Errors
/// As [`design_rna`], plus [`RnaDesignError::Goal`] for unparseable
/// dot-bracket strings.
pub fn design_riboswitch_workflow(
    free_dot_bracket: &str,
    bound_dot_bracket: &str,
    constraints: &DesignConstraints,
) -> Result<RnaDesignReport> {
    use valenx_rnastruct::Structure;
    let free = Structure::from_dot_bracket(free_dot_bracket)
        .map_err(|e| RnaDesignError::goal("target", e.to_string()))?;
    let bound = Structure::from_dot_bracket(bound_dot_bracket)
        .map_err(|e| RnaDesignError::goal("target", e.to_string()))?;
    constraints.validate()?;

    // Design stage — the two-state designer.
    let rs = design_riboswitch(&free, &bound, TwoStateParams::default())?;
    let design = rs.design;

    // Optimize.
    let optimization = optimize_design(&design, constraints, OptimizeParams::default())?;
    let optimized = optimization.design.clone();

    // Validate + export.
    let validation = validate_design(&optimized, constraints)?;
    let synthesis = plan_synthesis(&optimized, constraints, Promoter::T7)?;
    let text_report = export_design_report(
        &optimized,
        &validation,
        &synthesis,
        "riboswitch (two-state)",
    );

    let notes = vec![
        "Riboswitch workflow: two-state design, optimised, validated, synthesis-planned."
            .to_string(),
        format!(
            "Free-state base-pair distance {}, bound-state energy gap {:.2} kcal/mol.",
            rs.free_state_distance, rs.energy_gap,
        ),
        "In-silico prediction — the ligand is not modelled; lab-validate the switch.".to_string(),
    ];

    Ok(RnaDesignReport {
        design: optimized,
        optimization,
        validation,
        synthesis,
        text_report,
        notes,
    })
}

/// `true` when a design report's design is a coding mRNA — a tiny
/// helper used by callers inspecting a batch.
pub fn report_is_coding(report: &RnaDesignReport) -> bool {
    matches!(report.design.kind, DesignKind::Coding)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::goal::StructuralClass;
    use valenx_genediting::mrna::tailcap::MrnaUseCase;

    #[test]
    fn design_rna_structural_end_to_end() {
        let goal = DesignGoal::structural("((((((....))))))", StructuralClass::Hairpin).unwrap();
        let report = design_rna(&goal, &DesignConstraints::default()).unwrap();
        assert_eq!(report.design.len(), 16);
        assert!(!report.text_report.is_empty());
        // The workflow ran every stage.
        assert!(report.notes.iter().any(|n| n.starts_with("Goal:")));
        assert!(report.notes.iter().any(|n| n.starts_with("Export:")));
        // The quality score is in range.
        assert!((0.0..=1.0).contains(&report.quality_score()));
    }

    #[test]
    fn design_rna_coding_end_to_end() {
        let goal = DesignGoal::coding(b"MKVLAGD".to_vec(), MrnaUseCase::Vaccine);
        let report = design_rna(&goal, &DesignConstraints::default()).unwrap();
        assert!(report.design.is_coding());
        assert!(report.synthesis.template_len() > report.design.len());
    }

    #[test]
    fn design_rna_rejects_bad_goal() {
        let bad = DesignGoal::coding(Vec::new(), MrnaUseCase::Vaccine);
        assert!(design_rna(&bad, &DesignConstraints::default()).is_err());
    }

    #[test]
    fn design_rna_is_deterministic() {
        let goal = DesignGoal::structural("(((....)))", StructuralClass::Hairpin).unwrap();
        let a = design_rna(&goal, &DesignConstraints::default()).unwrap();
        let b = design_rna(&goal, &DesignConstraints::default()).unwrap();
        assert_eq!(a.design.sequence, b.design.sequence);
    }

    #[test]
    fn batch_designs_and_ranks() {
        let goal = DesignGoal::structural("((((....))))", StructuralClass::Hairpin).unwrap();
        let batch = design_rna_batch(&goal, &DesignConstraints::default(), 3).unwrap();
        assert_eq!(batch.len(), 3);
        // The batch is ranked best-first.
        for w in batch.candidates.windows(2) {
            assert!(
                w[0].quality_score() >= w[1].quality_score() - 1e-9,
                "batch is not sorted best-first"
            );
        }
        assert!(batch.best().is_some());
    }

    #[test]
    fn batch_rejects_zero_count() {
        let goal = DesignGoal::structural("((((....))))", StructuralClass::Hairpin).unwrap();
        assert!(design_rna_batch(&goal, &DesignConstraints::default(), 0).is_err());
    }

    #[test]
    fn llm_surface_runs_a_design() {
        let req = RnaDesignRequest::Design {
            goal: DesignGoal::structural("(((....)))", StructuralClass::Hairpin).unwrap(),
            constraints: DesignConstraints::default(),
        };
        let resp = handle_request(&req);
        assert!(resp.is_ok());
        assert!(matches!(resp, RnaDesignResponse::Design(_)));
    }

    #[test]
    fn llm_surface_runs_a_batch() {
        let req = RnaDesignRequest::DesignBatch {
            goal: DesignGoal::structural("(((....)))", StructuralClass::Hairpin).unwrap(),
            constraints: DesignConstraints::default(),
            count: 2,
        };
        let resp = handle_request(&req);
        assert!(resp.is_ok());
        match resp {
            RnaDesignResponse::Batch(b) => assert_eq!(b.len(), 2),
            _ => panic!("expected a Batch response"),
        }
    }

    #[test]
    fn llm_surface_gets_a_scaffold() {
        let req = RnaDesignRequest::GetScaffold {
            id: "gnra_tetraloop".to_string(),
        };
        let resp = handle_request(&req);
        assert!(resp.is_ok());
        assert!(matches!(resp, RnaDesignResponse::Scaffold(_)));
    }

    #[test]
    fn llm_surface_captures_errors() {
        // An unknown scaffold id → an Error response, not a panic.
        let req = RnaDesignRequest::GetScaffold {
            id: "no_such_scaffold".to_string(),
        };
        let resp = handle_request(&req);
        assert!(!resp.is_ok());
        match resp {
            RnaDesignResponse::Error { code, .. } => {
                assert!(code.starts_with("rnadesign."));
            }
            _ => panic!("expected an Error response"),
        }
    }

    #[test]
    fn request_response_types_are_clone_and_eq() {
        let req = RnaDesignRequest::GetScaffold {
            id: "gnra_tetraloop".to_string(),
        };
        assert_eq!(req.clone(), req);
        let resp = handle_request(&req);
        assert_eq!(resp.clone(), resp);
        assert!(resp.is_ok());
    }

    /// Compile-time proof the LLM-surface envelopes implement
    /// `Serialize` + `Deserialize` (so an MCP tool can ferry them as
    /// JSON).
    #[test]
    fn llm_surface_is_serde_capable() {
        fn assert_serde<T: serde::Serialize + serde::de::DeserializeOwned>() {}
        assert_serde::<RnaDesignRequest>();
        assert_serde::<RnaDesignResponse>();
        assert_serde::<RnaDesignReport>();
        assert_serde::<BatchDesignResult>();
    }

    #[test]
    fn riboswitch_workflow_runs() {
        let report = design_riboswitch_workflow(
            "((((....))))....",
            "....((((....))))",
            &DesignConstraints::default(),
        )
        .unwrap();
        assert_eq!(report.design.len(), 16);
        assert!(matches!(report.design.kind, DesignKind::Riboswitch { .. }));
    }
}
