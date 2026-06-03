//! Features 3–4 — the design session and the workflow state machine.
//!
//! A synthetic-RNA design is not a single function call — it is a
//! **workflow** with distinct stages, each consuming the previous
//! stage's output:
//!
//! ```text
//!   Goal  ->  Design  ->  Optimize  ->  Validate  ->  Export
//! ```
//!
//! [`DesignSession`] is the state machine that drives that workflow. It
//! holds the immutable design intent (the [`DesignGoal`] and the
//! [`DesignConstraints`]), the current [`WorkflowStage`], and the
//! intermediate results produced so far. Each `advance_*` method runs
//! one stage, checks the stage order, and stores the result.
//!
//! A session is the *stateful* front end; callers that just want the
//! whole pipeline in one call use [`crate::driver::design_rna`], which
//! drives a session internally.
//!
//! ## v1 scope
//!
//! The state machine enforces stage *order* (you cannot validate before
//! you design) but allows a stage to be re-run (re-optimising discards
//! the later results). It does not persist a session to disk or support
//! branching / undo — it is a single linear workflow object.

use crate::design::RnaDesign;
use crate::error::{Result, RnaDesignError};
use crate::goal::{DesignConstraints, DesignGoal};
use crate::optimize::OptimizationResult;
use crate::synthesis::SynthesisPackage;
use crate::validate::ValidationReport;

/// A stage of the design workflow (feature 4).
///
/// The stages are totally ordered; [`WorkflowStage::index`] gives the
/// position so a session can check that work happens in order.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum WorkflowStage {
    /// The goal and constraints are set; no sequence designed yet.
    Goal,
    /// A candidate sequence has been designed.
    Design,
    /// The candidate has been multi-objective-optimised.
    Optimize,
    /// The optimised candidate has been validated in silico.
    Validate,
    /// The synthesis-ready output package has been produced.
    Export,
}

impl WorkflowStage {
    /// The stage's 0-based position in the workflow order.
    pub fn index(self) -> usize {
        match self {
            WorkflowStage::Goal => 0,
            WorkflowStage::Design => 1,
            WorkflowStage::Optimize => 2,
            WorkflowStage::Validate => 3,
            WorkflowStage::Export => 4,
        }
    }

    /// A short human-readable name.
    pub fn name(self) -> &'static str {
        match self {
            WorkflowStage::Goal => "goal",
            WorkflowStage::Design => "design",
            WorkflowStage::Optimize => "optimize",
            WorkflowStage::Validate => "validate",
            WorkflowStage::Export => "export",
        }
    }

    /// The next stage in the workflow, or `None` if this is the last.
    pub fn next(self) -> Option<WorkflowStage> {
        match self {
            WorkflowStage::Goal => Some(WorkflowStage::Design),
            WorkflowStage::Design => Some(WorkflowStage::Optimize),
            WorkflowStage::Optimize => Some(WorkflowStage::Validate),
            WorkflowStage::Validate => Some(WorkflowStage::Export),
            WorkflowStage::Export => None,
        }
    }
}

/// A progress snapshot of a [`DesignSession`] (feature 4).
#[derive(Clone, Debug, PartialEq)]
pub struct StageStatus {
    /// The stage the session has currently completed.
    pub current: WorkflowStage,
    /// Fraction of the workflow complete, in `[0, 1]`.
    pub progress: f64,
    /// `true` when the workflow has reached its final stage.
    pub done: bool,
    /// A one-line human-readable status message.
    pub message: String,
}

/// The design-workflow state machine (feature 3).
///
/// Holds the design intent and every intermediate result. Drive it with
/// [`advance_design`](Self::advance_design),
/// [`advance_optimize`](Self::advance_optimize),
/// [`advance_validate`](Self::advance_validate) and
/// [`advance_export`](Self::advance_export), or use
/// [`crate::driver::design_rna`] to run the whole pipeline at once.
#[derive(Clone, Debug)]
pub struct DesignSession {
    goal: DesignGoal,
    constraints: DesignConstraints,
    stage: WorkflowStage,
    design: Option<RnaDesign>,
    optimization: Option<OptimizationResult>,
    validation: Option<ValidationReport>,
    synthesis: Option<SynthesisPackage>,
}

impl DesignSession {
    /// Opens a new design session for a goal and a constraint set.
    ///
    /// Both the goal and the constraints are validated up front, so a
    /// malformed design intent is rejected before any work begins.
    ///
    /// # Errors
    /// [`RnaDesignError::Goal`] if the goal or the constraints fail
    /// validation.
    pub fn new(goal: DesignGoal, constraints: DesignConstraints) -> Result<Self> {
        goal.validate()?;
        constraints.validate()?;
        Ok(DesignSession {
            goal,
            constraints,
            stage: WorkflowStage::Goal,
            design: None,
            optimization: None,
            validation: None,
            synthesis: None,
        })
    }

    /// The design goal.
    pub fn goal(&self) -> &DesignGoal {
        &self.goal
    }

    /// The design constraints.
    pub fn constraints(&self) -> &DesignConstraints {
        &self.constraints
    }

    /// The stage the session has currently completed.
    pub fn stage(&self) -> WorkflowStage {
        self.stage
    }

    /// The designed candidate, once the Design stage has run.
    pub fn design(&self) -> Option<&RnaDesign> {
        self.design.as_ref()
    }

    /// The optimisation result, once the Optimize stage has run.
    pub fn optimization(&self) -> Option<&OptimizationResult> {
        self.optimization.as_ref()
    }

    /// The validation report, once the Validate stage has run.
    pub fn validation(&self) -> Option<&ValidationReport> {
        self.validation.as_ref()
    }

    /// The synthesis package, once the Export stage has run.
    pub fn synthesis(&self) -> Option<&SynthesisPackage> {
        self.synthesis.as_ref()
    }

    /// A progress snapshot of the session.
    pub fn status(&self) -> StageStatus {
        let progress = self.stage.index() as f64 / WorkflowStage::Export.index() as f64;
        let done = self.stage == WorkflowStage::Export;
        let message = if done {
            "Design workflow complete — a synthesis-ready candidate package is available."
                .to_string()
        } else {
            format!(
                "Completed the `{}` stage; next is `{}`.",
                self.stage.name(),
                self.stage.next().map(|s| s.name()).unwrap_or("(none)"),
            )
        };
        StageStatus {
            current: self.stage,
            progress,
            done,
            message,
        }
    }

    /// Records the result of the **Design** stage and advances to
    /// `Design`.
    ///
    /// The design is supplied by the caller (the driver runs the
    /// appropriate `valenx-rnadesign::design` function); the session
    /// just stores it and moves the stage forward. Re-running this
    /// stage discards any later results.
    ///
    /// # Errors
    /// [`RnaDesignError::Invalid`] if the session has already advanced
    /// past `Design` — call [`reset_to`](Self::reset_to) first to
    /// re-design.
    pub fn advance_design(&mut self, design: RnaDesign) -> Result<()> {
        // Designing is allowed from Goal, or as a re-run from Design.
        if self.stage > WorkflowStage::Design {
            return Err(RnaDesignError::invalid(
                "stage",
                format!(
                    "cannot design from the `{}` stage — reset the session first",
                    self.stage.name()
                ),
            ));
        }
        self.design = Some(design);
        self.optimization = None;
        self.validation = None;
        self.synthesis = None;
        self.stage = WorkflowStage::Design;
        Ok(())
    }

    /// Records the result of the **Optimize** stage and advances to
    /// `Optimize`.
    ///
    /// # Errors
    /// [`RnaDesignError::Invalid`] if the Design stage has not run, or
    /// the session has advanced past `Optimize`.
    pub fn advance_optimize(&mut self, result: OptimizationResult) -> Result<()> {
        if self.stage < WorkflowStage::Design {
            return Err(RnaDesignError::invalid(
                "stage",
                "cannot optimise before a candidate has been designed",
            ));
        }
        if self.stage > WorkflowStage::Optimize {
            return Err(RnaDesignError::invalid(
                "stage",
                "cannot optimise from a later stage — reset the session first",
            ));
        }
        // The optimised design becomes the working design.
        self.design = Some(result.design.clone());
        self.optimization = Some(result);
        self.validation = None;
        self.synthesis = None;
        self.stage = WorkflowStage::Optimize;
        Ok(())
    }

    /// Records the result of the **Validate** stage and advances to
    /// `Validate`.
    ///
    /// # Errors
    /// [`RnaDesignError::Invalid`] if the Optimize stage has not run, or
    /// the session has advanced past `Validate`.
    pub fn advance_validate(&mut self, report: ValidationReport) -> Result<()> {
        if self.stage < WorkflowStage::Optimize {
            return Err(RnaDesignError::invalid(
                "stage",
                "cannot validate before the candidate has been optimised",
            ));
        }
        if self.stage > WorkflowStage::Validate {
            return Err(RnaDesignError::invalid(
                "stage",
                "cannot validate from a later stage — reset the session first",
            ));
        }
        self.validation = Some(report);
        self.synthesis = None;
        self.stage = WorkflowStage::Validate;
        Ok(())
    }

    /// Records the result of the **Export** stage and advances to
    /// `Export` — the terminal stage.
    ///
    /// # Errors
    /// [`RnaDesignError::Invalid`] if the Validate stage has not run.
    pub fn advance_export(&mut self, package: SynthesisPackage) -> Result<()> {
        if self.stage < WorkflowStage::Validate {
            return Err(RnaDesignError::invalid(
                "stage",
                "cannot export before the candidate has been validated",
            ));
        }
        self.synthesis = Some(package);
        self.stage = WorkflowStage::Export;
        Ok(())
    }

    /// Rewinds the session to an earlier stage, discarding every result
    /// produced at or after `stage`.
    ///
    /// Use this to re-run part of the workflow (e.g. reset to `Goal` to
    /// design a fresh candidate).
    ///
    /// # Errors
    /// [`RnaDesignError::Invalid`] if `stage` is not earlier than the
    /// current stage (there is nothing to rewind).
    pub fn reset_to(&mut self, stage: WorkflowStage) -> Result<()> {
        if stage >= self.stage {
            return Err(RnaDesignError::invalid(
                "stage",
                format!(
                    "cannot reset to `{}` — the session is already at `{}`",
                    stage.name(),
                    self.stage.name()
                ),
            ));
        }
        if stage < WorkflowStage::Design {
            self.design = None;
        }
        if stage < WorkflowStage::Optimize {
            self.optimization = None;
        }
        if stage < WorkflowStage::Validate {
            self.validation = None;
        }
        if stage < WorkflowStage::Export {
            self.synthesis = None;
        }
        self.stage = stage;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::design::DesignKind;
    use crate::goal::StructuralClass;
    use valenx_rnastruct::Structure;

    fn demo_goal() -> DesignGoal {
        DesignGoal::structural("((((....))))", StructuralClass::Hairpin).unwrap()
    }

    fn demo_design() -> RnaDesign {
        RnaDesign {
            sequence: b"GGGGAAAACCCC".to_vec(),
            kind: DesignKind::Structural {
                target: Structure::from_dot_bracket("((((....))))").unwrap(),
            },
            cds_span: None,
            construct: None,
            notes: Vec::new(),
        }
    }

    #[test]
    fn stage_ordering() {
        assert!(WorkflowStage::Goal < WorkflowStage::Design);
        assert!(WorkflowStage::Design < WorkflowStage::Export);
        assert_eq!(WorkflowStage::Goal.index(), 0);
        assert_eq!(WorkflowStage::Export.index(), 4);
        assert_eq!(WorkflowStage::Goal.next(), Some(WorkflowStage::Design));
        assert_eq!(WorkflowStage::Export.next(), None);
    }

    #[test]
    fn new_session_starts_at_goal() {
        let s = DesignSession::new(demo_goal(), DesignConstraints::default()).unwrap();
        assert_eq!(s.stage(), WorkflowStage::Goal);
        assert!(s.design().is_none());
        let status = s.status();
        assert!((status.progress - 0.0).abs() < 1e-9);
        assert!(!status.done);
    }

    #[test]
    fn new_session_rejects_bad_goal() {
        let bad = DesignGoal::coding(Vec::new(), valenx_genediting::mrna::tailcap::MrnaUseCase::Vaccine);
        assert!(DesignSession::new(bad, DesignConstraints::default()).is_err());
    }

    #[test]
    fn advance_design_moves_stage() {
        let mut s = DesignSession::new(demo_goal(), DesignConstraints::default()).unwrap();
        s.advance_design(demo_design()).unwrap();
        assert_eq!(s.stage(), WorkflowStage::Design);
        assert!(s.design().is_some());
    }

    #[test]
    fn stages_must_run_in_order() {
        let mut s = DesignSession::new(demo_goal(), DesignConstraints::default()).unwrap();
        // Cannot validate before designing+optimising.
        let dummy_report = crate::validate::ValidationReport {
            verdict: crate::validate::ValidationVerdict::Fail,
            constraint_checks: Vec::new(),
            fold_back: None,
            ensemble: None,
            robustness: None,
            notes: Vec::new(),
        };
        assert!(s.advance_validate(dummy_report).is_err());
    }

    #[test]
    fn reset_rewinds_and_clears() {
        let mut s = DesignSession::new(demo_goal(), DesignConstraints::default()).unwrap();
        s.advance_design(demo_design()).unwrap();
        assert!(s.design().is_some());
        // Reset to Goal clears the design.
        s.reset_to(WorkflowStage::Goal).unwrap();
        assert_eq!(s.stage(), WorkflowStage::Goal);
        assert!(s.design().is_none());
        // Cannot reset forward.
        assert!(s.reset_to(WorkflowStage::Design).is_err());
    }

    #[test]
    fn status_progress_advances() {
        let mut s = DesignSession::new(demo_goal(), DesignConstraints::default()).unwrap();
        let p0 = s.status().progress;
        s.advance_design(demo_design()).unwrap();
        let p1 = s.status().progress;
        assert!(p1 > p0);
    }
}
