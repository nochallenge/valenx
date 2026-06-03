//! Top-level driver and bundled report — feature 30.
//!
//! This module ties the crate together. It offers:
//!
//! - [`SysbioReport`] — a single bundled analysis of a reaction
//!   network: its steady state, its metabolic-flux distribution and a
//!   summary of a deterministic time course. One call gets a caller
//!   the headline numbers without orchestrating the ODE, FBA and
//!   steady-state layers by hand.
//! - [`analyze_model`] — the driver that produces a [`SysbioReport`].
//! - [`sbml_round_trip`] — an SBML-model round-trip helper: write a
//!   model to SBML and read it back, asserting structural equality —
//!   the standard interchange-fidelity check.
//!
//! The report is intentionally a *bundle of v1 summaries*, not a
//! COPASI-grade analysis suite — each underlying method documents its
//! own scope in its own module.

use crate::analysis::conservation::conservation_laws;
use crate::error::{Result, SysbioError};
use crate::fba::{FbaProblem, FbaSolution};
use crate::model::{read_sbml, write_sbml, Model};
use crate::ode::steady::steady_state;
use crate::ode::{OdeSystem, TimeCourse, Trajectory};

/// A compact summary of a deterministic time course.
#[derive(Debug, Clone, PartialEq)]
pub struct TimeCourseSummary {
    /// Simulated end time.
    pub t_end: f64,
    /// Number of output samples.
    pub n_samples: usize,
    /// Final amount of each species.
    pub final_state: Vec<f64>,
    /// Peak (maximum) amount reached by each species over the course.
    pub peak_state: Vec<f64>,
}

impl TimeCourseSummary {
    /// Summarise a [`Trajectory`].
    pub fn from_trajectory(traj: &Trajectory) -> Self {
        let dim = traj.states.first().map(|s| s.len()).unwrap_or(0);
        let mut peak = vec![f64::NEG_INFINITY; dim];
        for st in &traj.states {
            for (i, &v) in st.iter().enumerate() {
                if v > peak[i] {
                    peak[i] = v;
                }
            }
        }
        if peak.iter().any(|v| !v.is_finite()) {
            peak = vec![0.0; dim];
        }
        TimeCourseSummary {
            t_end: traj.times.last().copied().unwrap_or(0.0),
            n_samples: traj.len(),
            final_state: traj.final_state().map(|s| s.to_vec()).unwrap_or_default(),
            peak_state: peak,
        }
    }
}

/// A bundled systems-biology analysis of one reaction-network model.
#[derive(Debug, Clone, PartialEq)]
pub struct SysbioReport {
    /// The analysed model's id.
    pub model_id: String,
    /// Species count.
    pub n_species: usize,
    /// Reaction count.
    pub n_reactions: usize,
    /// Number of conserved moieties detected.
    pub n_conserved_moieties: usize,
    /// The located steady state, if the Newton solve converged.
    pub steady_state: Option<Vec<f64>>,
    /// Steady-state residual `‖f‖` (`None` if it did not converge).
    pub steady_residual: Option<f64>,
    /// The FBA flux distribution, if an objective was supplied and the
    /// LP was feasible.
    pub fba: Option<FbaSolution>,
    /// A summary of the deterministic time course.
    pub time_course: TimeCourseSummary,
}

impl SysbioReport {
    /// `true` if every headline analysis succeeded (steady state
    /// converged and, when requested, FBA was feasible).
    pub fn fully_resolved(&self) -> bool {
        self.steady_state.is_some()
            && self.fba.as_ref().map(|f| f.feasible).unwrap_or(true)
    }
}

/// Run the bundled analysis driver on `model` (feature 30).
///
/// Produces a [`SysbioReport`] containing:
/// 1. conserved-moiety count (stoichiometry null space),
/// 2. a deterministic time course to `t_end`,
/// 3. a steady-state solve warm-started from the time course's final
///    state (the robust choice — the trajectory has already flowed
///    toward the attractor), and
/// 4. optionally an FBA solve, if `fba_objective` names a reaction.
///
/// A non-converged steady-state solve or an infeasible FBA is recorded
/// in the report rather than raised — the report is a best-effort
/// bundle, and a partial result is still useful.
pub fn analyze_model(
    model: &Model,
    t_end: f64,
    fba_objective: Option<&str>,
) -> Result<SysbioReport> {
    model.validate()?;
    if t_end <= 0.0 {
        return Err(SysbioError::invalid("t_end", "t_end must be positive"));
    }

    // (1) Conservation analysis.
    let moieties = conservation_laws(model, 1e-9)?.len();

    // (2) Deterministic time course.
    let tc = TimeCourse::new(t_end);
    let traj = tc.run(model)?;
    let summary = TimeCourseSummary::from_trajectory(&traj);

    // (3) Steady state, warm-started from the time-course endpoint.
    let sys = OdeSystem::from_model(model);
    let warm = traj
        .final_state()
        .map(|s| s.to_vec())
        .unwrap_or_else(|| model.initial_state());
    let (steady_state_vec, residual) = match steady_state(&sys, &warm, 1e-7, 300) {
        Ok(ss) => (Some(ss.state), Some(ss.residual)),
        Err(_) => (None, None),
    };

    // (4) Optional FBA.
    let fba = match fba_objective {
        Some(obj) => {
            let mut prob = FbaProblem::from_model(model)?;
            prob.set_objective(obj)?;
            Some(prob.optimize()?)
        }
        None => None,
    };

    Ok(SysbioReport {
        model_id: model.id.clone(),
        n_species: model.species.len(),
        n_reactions: model.reactions.len(),
        n_conserved_moieties: moieties,
        steady_state: steady_state_vec,
        steady_residual: residual,
        fba,
        time_course: summary,
    })
}

/// SBML round-trip helper (feature 30).
///
/// Serialises `model` to SBML, parses it back, and verifies that the
/// reconstructed model has the same species, reactions, compartments
/// and parameters. Returns the re-parsed model on success, or a
/// [`SysbioError::Parse`] describing the first structural divergence.
///
/// This is the standard interchange-fidelity check — a design tool
/// runs it to be sure a model survives an export / import cycle.
pub fn sbml_round_trip(model: &Model) -> Result<Model> {
    let xml = write_sbml(model);
    let (back, report) = read_sbml(&xml)?;
    if back.species.len() != model.species.len() {
        return Err(SysbioError::parse(
            "sbml",
            "round-trip changed the species count",
        ));
    }
    if back.reactions.len() != model.reactions.len() {
        return Err(SysbioError::parse(
            "sbml",
            "round-trip changed the reaction count",
        ));
    }
    if back.compartments.len() != model.compartments.len() {
        return Err(SysbioError::parse(
            "sbml",
            "round-trip changed the compartment count",
        ));
    }
    if back.parameters.len() != model.parameters.len() {
        return Err(SysbioError::parse(
            "sbml",
            "round-trip changed the parameter count",
        ));
    }
    if report.placeholder_laws > 0 {
        return Err(SysbioError::parse(
            "sbml",
            "round-trip lost a kinetic law (placeholder substituted)",
        ));
    }
    if back.events.len() != model.events.len() {
        return Err(SysbioError::parse(
            "sbml",
            "round-trip changed the event count",
        ));
    }
    if back.rules.assignments.len() != model.rules.assignments.len()
        || back.rules.rates.len() != model.rules.rates.len()
    {
        return Err(SysbioError::parse(
            "sbml",
            "round-trip changed the rule count",
        ));
    }
    Ok(back)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{RateLaw, Reaction, Species};

    /// Source/decay model: 0 ->(rate 4) A ->(2 A) 0. Steady A* = 2.
    fn source_decay() -> Model {
        let mut m = Model::new("report_demo");
        let a = m.add_species(Species::new("A", 0.0));
        m.add_reaction(Reaction {
            id: "src".into(),
            reactants: vec![],
            products: vec![(a, 1.0)],
            rate_law: RateLaw::Constant { rate: 4.0 },
            reversible: false,
        });
        m.add_reaction(Reaction {
            id: "dec".into(),
            reactants: vec![(a, 1.0)],
            products: vec![],
            rate_law: RateLaw::MassAction {
                k: 2.0,
                reactants: vec![(a, 1.0)],
            },
            reversible: false,
        });
        m
    }

    #[test]
    fn report_bundles_steady_state_and_time_course() {
        let m = source_decay();
        let report = analyze_model(&m, 20.0, None).unwrap();
        assert_eq!(report.model_id, "report_demo");
        assert_eq!(report.n_species, 1);
        assert_eq!(report.n_reactions, 2);
        // Steady state A* = source/decay = 4/2 = 2.
        let ss = report.steady_state.as_ref().expect("converged");
        assert!((ss[0] - 2.0).abs() < 1e-4, "A* = {}", ss[0]);
        // Time course ends near the same value.
        assert!((report.time_course.final_state[0] - 2.0).abs() < 0.1);
        assert!(report.fully_resolved());
    }

    #[test]
    fn report_includes_fba_when_objective_given() {
        // A linear pathway model so FBA has something to maximise.
        let mut m = Model::new("path");
        let a = m.add_species(Species::new("A", 0.0));
        let b = m.add_species(Species::new("B", 0.0));
        m.add_reaction(Reaction {
            id: "imp".into(),
            reactants: vec![],
            products: vec![(a, 1.0)],
            rate_law: RateLaw::Constant { rate: 0.0 },
            reversible: false,
        });
        m.add_reaction(Reaction {
            id: "conv".into(),
            reactants: vec![(a, 1.0)],
            products: vec![(b, 1.0)],
            rate_law: RateLaw::Constant { rate: 0.0 },
            reversible: false,
        });
        m.add_reaction(Reaction {
            id: "exp".into(),
            reactants: vec![(b, 1.0)],
            products: vec![],
            rate_law: RateLaw::Constant { rate: 0.0 },
            reversible: false,
        });
        let report = analyze_model(&m, 5.0, Some("exp")).unwrap();
        let fba = report.fba.expect("fba ran");
        assert!(fba.feasible);
        // Default bounds cap flux at 1000.
        assert!((fba.objective_value - 1000.0).abs() < 1e-3);
    }

    #[test]
    fn report_records_conserved_moiety() {
        // A <-> B isomerisation conserves A + B.
        let mut m = Model::new("iso");
        let a = m.add_species(Species::new("A", 5.0));
        let b = m.add_species(Species::new("B", 5.0));
        m.add_reaction(Reaction {
            id: "f".into(),
            reactants: vec![(a, 1.0)],
            products: vec![(b, 1.0)],
            rate_law: RateLaw::MassAction {
                k: 1.0,
                reactants: vec![(a, 1.0)],
            },
            reversible: true,
        });
        m.add_reaction(Reaction {
            id: "r".into(),
            reactants: vec![(b, 1.0)],
            products: vec![(a, 1.0)],
            rate_law: RateLaw::MassAction {
                k: 1.0,
                reactants: vec![(b, 1.0)],
            },
            reversible: true,
        });
        let report = analyze_model(&m, 10.0, None).unwrap();
        assert_eq!(report.n_conserved_moieties, 1);
    }

    #[test]
    fn sbml_round_trip_preserves_the_model() {
        let m = source_decay();
        let back = sbml_round_trip(&m).unwrap();
        assert_eq!(back.species.len(), m.species.len());
        assert_eq!(back.reactions.len(), m.reactions.len());
        // The reconstructed model still simulates to the same steady
        // state.
        let report = analyze_model(&back, 20.0, None).unwrap();
        assert!((report.steady_state.unwrap()[0] - 2.0).abs() < 1e-3);
    }

    #[test]
    fn analyze_rejects_bad_horizon() {
        let m = source_decay();
        assert!(analyze_model(&m, -1.0, None).is_err());
    }

    #[test]
    fn analyze_rejects_unknown_fba_objective() {
        let m = source_decay();
        assert!(analyze_model(&m, 5.0, Some("ghost")).is_err());
    }

    #[test]
    fn sbml_round_trip_preserves_events_and_rules() {
        use crate::model::events::{
            AssignmentRule, EventAssignment, RateRule, SbmlEvent, VarRef,
        };
        use crate::model::expr::Expr;
        let mut m = Model::new("evt_rules");
        let _a = m.add_species(Species::new("A", 0.0));
        let _b = m.add_species(Species::new("B", 0.0));
        m.add_parameter(crate::model::Parameter::new("k", 1.0));
        m.add_reaction(Reaction {
            id: "noop".into(),
            reactants: vec![],
            products: vec![],
            rate_law: RateLaw::Constant { rate: 0.0 },
            reversible: false,
        });
        m.add_rate_rule(RateRule {
            target: VarRef::Species(0),
            formula: Expr::param(0),
        });
        m.add_assignment_rule(AssignmentRule {
            target: VarRef::Species(1),
            formula: Expr::mul(Expr::var(0), Expr::k(2.0)),
        });
        m.add_event(SbmlEvent::new(
            "halt",
            Expr::ge(Expr::Time, Expr::k(3.0)),
            vec![EventAssignment {
                target: VarRef::Parameter(0),
                formula: Expr::k(0.0),
            }],
        ));
        let back = sbml_round_trip(&m).unwrap();
        assert_eq!(back.events.len(), 1);
        assert_eq!(back.rules.rates.len(), 1);
        assert_eq!(back.rules.assignments.len(), 1);
    }

    #[test]
    fn time_course_summary_tracks_peak() {
        let traj = Trajectory {
            times: vec![0.0, 1.0, 2.0],
            states: vec![vec![1.0], vec![5.0], vec![3.0]],
        };
        let s = TimeCourseSummary::from_trajectory(&traj);
        assert_eq!(s.peak_state, vec![5.0]);
        assert_eq!(s.final_state, vec![3.0]);
        assert_eq!(s.n_samples, 3);
    }
}
