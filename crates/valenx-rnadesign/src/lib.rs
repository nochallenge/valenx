//! # valenx-rnadesign — unified synthetic-RNA design workflow
//!
//! Round 6 Block 16 of the Valenx roadmap — a Round 6 *extension* that
//! turns the scattered RNA-design capability of the earlier Round 6
//! crates into one guided, start-to-finish workflow.
//!
//! **This crate orchestrates; it does not reimplement folding.** Round
//! 6 already shipped the RNA building blocks — this crate is the
//! integration layer on top of them:
//!
//! - [`valenx_rnastruct`] (Block 6.5) — RNA secondary structure: Zuker
//!   minimum-free-energy folding, the McCaskill partition function and
//!   base-pair probabilities, inverse folding, structure distance.
//! - [`valenx_genediting`] (Block 6.13) — mRNA-therapeutic design: the
//!   five-part construct model, codon optimisation, start-region
//!   structure minimisation, uridine / pseudouridine planning.
//! - [`valenx_bioseq`] (Block 6.1) — the [`Seq`](valenx_bioseq::Seq)
//!   type, transcription, the restriction-enzyme database, primer
//!   design, and FASTA / GenBank I/O.
//!
//! It is the same software *category* as the unified synthetic-RNA
//! design environments — the design front-ends of Benchling / Geneious,
//! the NUPACK / ViennaRNA design suites, the IDT / GenScript synthesis
//! design portals — a single pipeline from a design goal to a
//! synthesis-ready candidate package.
//!
//! ## What it does
//!
//! - **Goal & workflow** ([`goal`], [`workflow`]) — a [`DesignGoal`]
//!   (a structural RNA targeting a dot-bracket; a coding mRNA encoding a
//!   protein; or a hybrid), a [`DesignConstraints`] set, and a
//!   [`DesignSession`] state machine that walks Goal → Design →
//!   Optimize → Validate → Export.
//! - **Sequence design** ([`design`]) — structural inverse folding,
//!   coding-mRNA construct assembly, riboswitch / two-state design,
//!   functional-motif scaffold templates, regulatory-element design.
//! - **Multi-objective optimisation** ([`optimize`]) — a weighted
//!   simulated-annealing optimiser balancing target-structure match,
//!   ensemble defect, GC content, repeats, restriction sites, codon
//!   adaptation, uridine content and off-target hairpins; plus the
//!   ensemble-defect metric, repeat / low-complexity scan,
//!   forbidden-motif removal and synthesizability scan.
//! - **In-silico validation** ([`validate`]) — fold-back validation,
//!   partition-function ensemble validation, robustness / melting /
//!   co-transcriptional sanity checks, and a [`ValidationReport`].
//! - **Output & synthesis plan** ([`synthesis`], [`export`]) — DNA
//!   template generation with T7 / SP6 promoters, an in-vitro-
//!   transcription plan, a synthesis-order package, and FASTA /
//!   GenBank / text-report export.
//! - **Driver** ([`driver`]) — the top-level [`design_rna`] entry
//!   point, a batch / ranking mode, and a typed MCP/LLM-controllable
//!   request / response surface.
//!
//! ## Errors
//!
//! Every fallible function returns
//! [`Result<_, RnaDesignError>`](error::RnaDesignError). The error type
//! carries stable [`code`](error::RnaDesignError::code) and
//! [`category`](error::RnaDesignError::category) accessors for
//! telemetry, and an [`Upstream`](error::RnaDesignError::Upstream)
//! variant that names the building-block crate when a folding /
//! codon-optimisation call fails underneath.
//!
//! ## Honest framing — read this
//!
//! **The output of this workflow is a strong, validated-*in-silico*
//! candidate — a *prediction* from an energy model, NOT a guarantee of
//! correct in-vivo behaviour.**
//!
//! Every structural claim this crate makes comes from a
//! nearest-neighbor secondary-structure energy model (the Turner-2004
//! parameters, as implemented by `valenx-rnastruct` — themselves a
//! faithful representative subset, see that crate's note). Every
//! translation / expression score comes from a transparent rule-based
//! heuristic in `valenx-genediting`. None of it is a measurement.
//!
//! Physical RNA is made by **chemical synthesis or in-vitro
//! transcription in a wet lab**; a design that this crate reports as
//! folding to its target, clearing every constraint and passing
//! validation is a *well-supported hypothesis* that must still be
//! synthesised and **lab-validated** — folded, assayed, and tested for
//! function. Nothing in this crate is named "verified" or "guaranteed
//! correct", and the [`ValidationReport`] carries this disclaimer in
//! its own type documentation. Treat the synthesis package as the
//! starting point for wet-lab work, not its conclusion.
//!
//! ## v1 scope
//!
//! This is a real working v1, not production parity with a commercial
//! design suite. Each module documents its own simplifications inline;
//! the load-bearing ones are:
//!
//! - Structural design is a multi-start of the `valenx-rnastruct`
//!   inverse-folding adaptive walk — not an exhaustive or
//!   constraint-programming designer.
//! - The riboswitch / two-state designer is a thermodynamic heuristic:
//!   it does not model the ligand, its binding energy, or switching
//!   kinetics.
//! - The multi-objective optimiser is a weighted-sum simulated
//!   annealing over synonymous (CDS) or structure-preserving mutations
//!   — a local search, not a global optimum.
//! - The ensemble-defect, off-target-hairpin and synthesizability
//!   metrics are transparent computed scores, not trained models.
//! - Co-transcriptional folding is a v1 sanity check (a 5′-window fold
//!   progression), not a full kinetic folding simulation.
//! - No `cargo test` is ever run here per the project lockdown — tests
//!   are compile-checked by `cargo clippy --all-targets` only, never
//!   executed.

#![forbid(unsafe_code)]

pub mod aptamer;
pub mod constraints;
pub mod design;
pub mod driver;
pub mod error;
pub mod export;
pub mod goal;
pub mod inverse;
pub mod lineardesign;
pub mod multistate;
pub mod optimize;
pub mod riboswitch_ed;
pub mod synthesis;
pub mod tube;
pub mod validate;
pub mod workflow;

// --- Convenience re-exports of the most-used types --------------------

pub use error::{ErrorCategory, Result, RnaDesignError};

pub use goal::{
    DesignConstraints, DesignGoal, ElementPlacement, StructuralClass, StructuralElement,
};
pub use workflow::{DesignSession, StageStatus, WorkflowStage};

pub use design::{DesignKind, RnaDesign};

pub use lineardesign::{
    linear_design, pareto_sweep, LinearDesignRequest, LinearDesignResult, ParetoPoint,
};
pub use inverse::{
    inverse_fold_constrained, inverse_fold_ensemble_defect, EnsembleDefectDesign,
    EnsembleDefectParams,
};
pub use constraints::{lock_entry, parse_locked, DesignConstraintSet};
pub use multistate::{design_multistate, MultiStateDesign, MultiStateParams, StateSpec};

pub use aptamer::{
    design_aptamer, extract_pockets, pharmacophore_pocket_score, AptamerDesign,
    AptamerDesignParams, BaseEdgeFeatures, FeatureKind, Pharmacophore, PharmacophoreFeature,
    Pocket, PocketKind,
};
pub use riboswitch_ed::{
    design_riboswitch_ed, LigandBindingSite, LigandConstraint, RiboswitchEdDesign,
    RiboswitchEdParams,
};
pub use tube::{
    design_tube, fold_all_complexes, solve_tube_equilibrium, ComplexEnergies, ComplexKind,
    TargetDistribution, TargetFraction, TubeDesign, TubeDesignParams, TubeEquilibrium,
    TubeStrand,
};

pub use optimize::{optimize_design, ObjectiveWeights, OptimizationResult};
pub use validate::{validate_design, ConstraintCheck, ValidationReport, ValidationVerdict};

pub use synthesis::{plan_synthesis, Promoter, SynthesisPackage};
pub use export::{export_design_report, export_fasta, export_genbank};

pub use driver::{
    design_rna, design_rna_batch, handle_request, RnaDesignReport, RnaDesignRequest,
    RnaDesignResponse,
};

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_genediting::mrna::tailcap::MrnaUseCase;

    /// End-to-end: a structural-RNA design goal runs the whole
    /// Goal → Export workflow and produces a coherent report.
    #[test]
    fn structural_design_end_to_end() {
        let goal =
            DesignGoal::structural("((((((....))))))", StructuralClass::Hairpin).unwrap();
        let report = design_rna(&goal, &DesignConstraints::default()).unwrap();
        // The report carries a design, a validation, and a synthesis plan.
        assert!(!report.design.is_empty());
        assert_eq!(report.design.len(), 16);
        // The synthesis package has a DNA template longer than the RNA
        // (it carries a promoter).
        assert!(report.synthesis.dna_template.len() > report.design.len());
    }

    /// End-to-end: a coding-mRNA design goal.
    #[test]
    fn coding_design_end_to_end() {
        let goal = DesignGoal::coding(b"MKVLAGD".to_vec(), MrnaUseCase::Vaccine);
        let report = design_rna(&goal, &DesignConstraints::default()).unwrap();
        assert!(report.design.is_coding());
        assert!(report.design.construct.is_some());
    }

    #[test]
    fn re_exports_are_wired() {
        // Touch a sampling of the convenience re-exports.
        let e = RnaDesignError::not_yet("x");
        assert_eq!(e.category_enum(), ErrorCategory::Capability);
        let c = DesignConstraints::default();
        assert!(c.validate().is_ok());
        let w = ObjectiveWeights::default();
        assert!(w.structure_match >= 0.0);
    }
}
