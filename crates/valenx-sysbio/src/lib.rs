//! # valenx-sysbio — systems & synthetic biology
//!
//! A native-Rust replacement
//! for the daily-driver core of the systems- and synthetic-biology
//! tool landscape: COPASI, Tellurium / libRoadRunner, BioNetGen,
//! CellDesigner, libSBML, PySB, iBioSim, SBOLCanvas / pySBOL, j5 and
//! Synbiopython — pure deterministic algorithms with a hand-rolled
//! PRNG, no neural-network weights and no external processes.
//!
//! It builds on [`valenx_bioseq`] (Block 6.1) for the
//! [`Seq`](valenx_bioseq::Seq) type used by the SBOL part model and
//! the DNA-assembly planners, and on `nalgebra` for the dense linear
//! algebra behind the ODE Jacobians, the BDF Newton solve, the
//! conserved-moiety null-space computation and the FBA layer.
//!
//! ## What it does
//!
//! - **Model** ([`mod@model`]) — a reaction-network [`Model`] of
//!   species, reactions, compartments and parameters; an SBML-subset
//!   reader and writer (now with SBML L3 events, assignment rules and
//!   rate rules); mass-action / Michaelis-Menten / Hill kinetic
//!   rate laws; a BioNetGen-class rule-based model expander; and a
//!   small expression AST ([`Expr`]) for triggers and rule formulas.
//! - **Deterministic simulation** ([`ode`]) — stoichiometry-matrix ODE
//!   assembly, RK4 / adaptive RK45 / implicit BDF integrators, a
//!   damped-Newton steady-state solver, an event-aware time-course
//!   driver, and an [`EventDrivenTimeCourse`] that detects SBML L3
//!   trigger crossings between integrator steps with bisection.
//! - **Stochastic simulation** ([`stochastic`]) — the Gillespie
//!   direct-method SSA, explicit tau-leaping, the Gibson-Bruck
//!   next-reaction method, and ensemble runs with mean / variance /
//!   percentile statistics.
//! - **Analysis** ([`analysis`]) — 1-D / 2-D parameter scans, local
//!   and Morris-global sensitivity analysis, conserved-moiety
//!   detection, a steady-state continuation bifurcation scan, and a
//!   Levenberg-Marquardt parameter-estimation driver with
//!   Latin-hypercube + simulated-annealing pre-stages and Hessian-
//!   based standard errors ([`estimate_parameters`]).
//! - **Constraint-based modelling** ([`fba`]) — a from-scratch simplex
//!   LP solver, flux-balance analysis, flux variability analysis and
//!   parsimonious FBA.
//! - **Synthetic biology** ([`synbio`]) — an SBOL-class genetic-design
//!   data model, a Cello-class genetic-circuit logic simulator, a
//!   gene-regulatory-network ODE model, Gibson / Golden-Gate /
//!   BioBrick DNA-assembly planners, and a standard part library.
//! - **Pipeline** ([`pipeline`]) — a bundled [`pipeline::SysbioReport`]
//!   and an SBML-model round-trip helper.
//!
//! ## Errors
//!
//! Every fallible function returns
//! [`Result<_, SysbioError>`](error::SysbioError). The error type
//! carries stable [`code`](error::SysbioError::code) and
//! [`category`](error::SysbioError::category) accessors for telemetry.
//!
//! ## v1 scope
//!
//! This is a real working v1, not production parity with the reference
//! tools. Each module documents its own simplifications inline; the
//! notable ones are: the SBML support is the *core subset* in its
//! text encoding only (no Level-3 packages, no full MathML evaluator
//! — kinetic laws, event triggers and rule formulas round-trip through
//! a compact annotation); SBML L3 algebraic rules `0 = f(...)` are
//! supported only in the explicit-substitution case (the implicit-DAE
//! form is an honest omission); the BioNetGen-class rule expander does
//! combinatorial *site-state* expansion but has no bonds or complexes;
//! the BDF integrator is order 1-2 (A-stable but not the order-5 of
//! CVODE / LSODA); the FBA LP solver is a dense two-phase simplex —
//! correct but sized for v1 toy metabolic models, not a 10 000-reaction
//! genome-scale reconstruction; the bifurcation scan follows a single
//! steady-state branch (no pseudo-arclength fold tracing); the
//! parameter-estimation standard errors come from the Gauss-Newton
//! `J^T J` approximation (correct for a well-conditioned fit, not a
//! profile-likelihood / Fisher-information treatment); and the
//! Cello-class circuit layer *evaluates and scores* a supplied gate
//! assignment rather than searching for one.

#![forbid(unsafe_code)]

pub mod analysis;
pub mod error;
pub mod fba;
pub mod model;
pub mod ode;
pub mod pipeline;
pub mod stochastic;
pub mod synbio;

// --- Convenience re-exports of the most-used types --------------------

pub use analysis::estimation::{
    estimate_parameters, EstimationOptions, EstimationReport, EstimationTarget,
    ObservedPoint,
};
pub use error::{ErrorCategory, Result, SysbioError};
pub use fba::{FbaProblem, FbaSolution};
pub use model::events::{
    AssignmentRule, EventAssignment, RateRule, SbmlEvent, SbmlRules, VarRef,
};
pub use model::expr::Expr;
pub use model::{read_sbml, write_sbml, Model, RateLaw, Reaction, Species};
pub use ode::{
    steady_state, EventDrivenTimeCourse, EventTrajectory, OdeSystem, TimeCourse, Trajectory,
};
pub use pipeline::{analyze_model, sbml_round_trip, SysbioReport};
pub use stochastic::{run_ensemble, StochasticModel};
pub use synbio::{Circuit, Device, GeneNetwork, Part, PartLibrary};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crate_error_reexported() {
        let e = SysbioError::not_yet("x");
        assert_eq!(e.category_enum(), ErrorCategory::Capability);
    }

    #[test]
    fn end_to_end_model_build_and_analyze() {
        // Build a tiny model through the public API and run the
        // top-level driver — the crate's headline workflow.
        let mut m = Model::new("smoke");
        let a = m.add_species(Species::new("A", 0.0));
        m.add_reaction(Reaction {
            id: "src".into(),
            reactants: vec![],
            products: vec![(a, 1.0)],
            rate_law: RateLaw::Constant { rate: 3.0 },
            reversible: false,
        });
        m.add_reaction(Reaction {
            id: "dec".into(),
            reactants: vec![(a, 1.0)],
            products: vec![],
            rate_law: RateLaw::MassAction {
                k: 1.0,
                reactants: vec![(a, 1.0)],
            },
            reversible: false,
        });
        let report = analyze_model(&m, 30.0, None).expect("driver runs");
        // Steady state A* = 3/1 = 3.
        assert!((report.steady_state.unwrap()[0] - 3.0).abs() < 1e-3);
    }
}
