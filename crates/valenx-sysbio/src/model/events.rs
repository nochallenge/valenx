//! SBML L3 events, assignment rules and rate rules - features 32 - 34.
//!
//! Three pieces of commercial-depth SBML L3 machinery live here, each
//! using the [`Expr`] AST from `model::expr`:
//!
//! - [`SbmlEvent`] - a *trigger* expression (boolean-valued) plus zero
//!   or more [`EventAssignment`]s that fire when the trigger
//!   transitions from `false` to `true`. Triggers may be time-based
//!   (`t >= 10`), state-based (`s0 >= 100`) or arbitrary boolean
//!   combinations of both, and may carry an optional `delay` and
//!   `priority` (simultaneous-event ordering, biggest-first), plus the
//!   `use_values_from_trigger_time` flag that decides whether the
//!   assignment right-hand sides are evaluated at the trigger moment
//!   or at the delayed execution moment.
//! - [`AssignmentRule`] - `var := f(...)`, an algebraic constraint
//!   enforced every integrator output and every rule-graph
//!   re-evaluation. The set of rules is topologically sorted so they
//!   can be applied in one sweep with no recurrence.
//! - [`RateRule`] - `dvar/dt = f(...)`, folded into the ODE RHS.
//!   Targets a species or a global parameter; for the species case the
//!   contribution is added to whatever the stoichiometric law already
//!   produces (the COPASI convention).
//!
//! A general algebraic rule of the form `0 = f(...)` is handled by the
//! v1 path *only* in the explicit-substitution case where the
//! algebraic-rule expression can be rewritten as `var := g(...)` by
//! solving for the dependent variable (i.e. it is structurally an
//! assignment rule). Implicit DAEs are an honestly documented
//! omission in [`SbmlRules::topo_sort`].
//!
//! ## Where these get used
//!
//! The simulation surface is the existing time-course driver: it picks
//! up events and rules out of the [`crate::model::Model`] (no separate
//! API), checks
//! rules every integrator output, and watches triggers between RK / BDF
//! steps with a bisection root-find. The SBML reader / writer extends
//! the same `sbml:*`-attribute encoding the rate-law writer uses, so
//! event-carrying models round-trip just like the rest.

use serde::{Deserialize, Serialize};

use crate::error::{Result, SysbioError};
use crate::model::expr::Expr;

/// What kind of model variable an event assignment or a rule targets.
///
/// Species are addressed by their position in the species table;
/// parameters by their position in the global parameter table. The
/// model index is the same scheme the rate laws use.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum VarRef {
    /// A species.
    Species(usize),
    /// A global parameter.
    Parameter(usize),
}

/// A single SBML L3 event assignment: a target plus the right-hand
/// expression to assign at fire-time.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EventAssignment {
    /// The variable to overwrite.
    pub target: VarRef,
    /// The numeric expression whose value becomes the new amount.
    pub formula: Expr,
}

/// An SBML L3 event: a trigger plus zero or more assignments, with an
/// optional delay and priority.
///
/// The driver evaluates `trigger.value(...)` each integrator step.
/// The trigger is *fired* on a rising-edge crossing: it was non-
/// positive at the previous step (`<= 0`) and is now positive (`> 0`).
/// Bisection refines the crossing time to within `1e-6` of the step.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SbmlEvent {
    /// Stable identifier (kept for SBML round-trip; the simulator uses
    /// the event's index in the model).
    pub id: String,
    /// Trigger expression. Treated as boolean (positive = true), so
    /// `Gt`, `Ge`, `Lt`, `Le`, `And`, `Or` and `Not` all work as the
    /// user expects, and a raw `s0 - 5.0` also works as a real-valued
    /// crossing function.
    pub trigger: Expr,
    /// Assignments fired when the trigger goes from false to true.
    pub assignments: Vec<EventAssignment>,
    /// Optional delay (in simulation time units) between the trigger
    /// firing and the assignments being applied. `None` means
    /// instantaneous.
    pub delay: Option<f64>,
    /// Priority among simultaneously firing events. Higher fires
    /// first; ties broken by event index. Default `0.0`.
    pub priority: f64,
    /// If `true`, the right-hand sides are computed at the trigger
    /// moment (the SBML default) and frozen until the delayed
    /// execution applies them. If `false`, they are computed at the
    /// execution moment - the value can therefore depend on whatever
    /// other events fired in between.
    pub use_values_from_trigger_time: bool,
    /// Optional `t0` initial-state value of the trigger. Set to `false`
    /// on construction; the driver overwrites it once after the
    /// initial state is established so the very first crossing
    /// detection has a proper "previous" baseline.
    pub initial_value: bool,
}

impl SbmlEvent {
    /// Constructor with sensible defaults: no delay, priority zero,
    /// `use_values_from_trigger_time` true, `initial_value` false.
    pub fn new(id: impl Into<String>, trigger: Expr, assignments: Vec<EventAssignment>) -> Self {
        SbmlEvent {
            id: id.into(),
            trigger,
            assignments,
            delay: None,
            priority: 0.0,
            use_values_from_trigger_time: true,
            initial_value: false,
        }
    }

    /// Builder: set the delay.
    pub fn with_delay(mut self, d: f64) -> Self {
        self.delay = Some(d);
        self
    }
    /// Builder: set the priority.
    pub fn with_priority(mut self, p: f64) -> Self {
        self.priority = p;
        self
    }
    /// Builder: flip the "use values from execution time" flag.
    pub fn evaluate_at_execution_time(mut self) -> Self {
        self.use_values_from_trigger_time = false;
        self
    }
}

/// An SBML L3 assignment rule: `target := formula`.
///
/// Rules are applied every integrator output and every rule-graph
/// sweep. They are not part of the dynamics in the integrator's eyes -
/// they project the integrator's raw state onto the algebraic
/// constraint surface before any other code reads `y`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AssignmentRule {
    /// The targeted species or parameter.
    pub target: VarRef,
    /// The numeric right-hand side.
    pub formula: Expr,
}

/// An SBML L3 rate rule: `d target / dt = formula`.
///
/// The driver folds `formula(t, y, p)` into the species or parameter
/// derivative. For a species this is *additive* on top of the
/// stoichiometric law - so a species with both reaction kinetics and
/// a rate rule sees both contributions, the COPASI semantic. For a
/// parameter the rate rule is the *only* source of dynamics; the
/// parameter is treated as a state component.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RateRule {
    /// The targeted species or parameter.
    pub target: VarRef,
    /// The numeric derivative expression.
    pub formula: Expr,
}

/// The full SBML L3 rule set attached to a [`crate::model::Model`].
///
/// Stored on the model and consumed by the time-course driver. The
/// assignment-rule list is topologically sorted on demand
/// ([`SbmlRules::topo_sort`]) so the driver can apply them in one
/// pass.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct SbmlRules {
    /// Assignment rules.
    pub assignments: Vec<AssignmentRule>,
    /// Rate rules.
    pub rates: Vec<RateRule>,
}

impl SbmlRules {
    /// Topologically sort the assignment rules so each rule's
    /// right-hand side reads only species / parameters that have
    /// either been assigned by an earlier rule or are independent
    /// (driven by the integrator).
    ///
    /// Returns a permutation - an index into [`SbmlRules::assignments`]
    /// in execution order - or [`SysbioError::InvalidModel`] if the
    /// dependency graph has a cycle.
    pub fn topo_sort(&self) -> Result<Vec<usize>> {
        let n = self.assignments.len();
        if n == 0 {
            return Ok(Vec::new());
        }
        // assigned[i] = the var assigned by rule i.
        let assigned: Vec<&VarRef> = self.assignments.iter().map(|r| &r.target).collect();

        // For each rule, the set of *other* rules whose target appears
        // in its formula. That dependency edge means "the other rule
        // must run first".
        let mut deps: Vec<Vec<usize>> = vec![Vec::new(); n];
        for (i, rule) in self.assignments.iter().enumerate() {
            let vd = rule.formula.var_deps();
            let pd = rule.formula.param_deps();
            for (j, other_target) in assigned.iter().enumerate() {
                if i == j {
                    continue;
                }
                let depends_on = match other_target {
                    VarRef::Species(s) => vd.contains(s),
                    VarRef::Parameter(p) => pd.contains(p),
                };
                if depends_on {
                    deps[i].push(j);
                }
            }
        }

        // Kahn's algorithm.
        let mut in_deg: Vec<usize> = deps.iter().map(|d| d.len()).collect();
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
        for (i, d) in deps.iter().enumerate() {
            for &j in d {
                adj[j].push(i);
            }
        }
        let mut queue: Vec<usize> = (0..n).filter(|&i| in_deg[i] == 0).collect();
        let mut out = Vec::with_capacity(n);
        while let Some(node) = queue.pop() {
            out.push(node);
            for &child in &adj[node] {
                in_deg[child] -= 1;
                if in_deg[child] == 0 {
                    queue.push(child);
                }
            }
        }
        if out.len() != n {
            return Err(SysbioError::invalid_model(
                "reaction_network",
                "assignment-rule dependency graph is cyclic",
            ));
        }
        Ok(out)
    }

    /// Apply the assignment rules in topological order to `y` and `p`.
    /// `t` is the current simulation time. Returns the rule
    /// permutation used (for diagnostics).
    pub fn apply_assignments(
        &self,
        y: &mut [f64],
        p: &mut [f64],
        t: f64,
    ) -> Result<Vec<usize>> {
        let order = self.topo_sort()?;
        for &idx in &order {
            let rule = &self.assignments[idx];
            let v = rule.formula.value(y, p, t);
            match &rule.target {
                VarRef::Species(i) => {
                    if let Some(slot) = y.get_mut(*i) {
                        *slot = v;
                    }
                }
                VarRef::Parameter(i) => {
                    if let Some(slot) = p.get_mut(*i) {
                        *slot = v;
                    }
                }
            }
        }
        Ok(order)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn topo_sort_chains_in_dependency_order() {
        // s0 := 1, s1 := s0 + 1, s2 := s1 * 2.
        // Storage order is shuffled to make the test interesting.
        let rules = SbmlRules {
            assignments: vec![
                AssignmentRule {
                    target: VarRef::Species(2),
                    formula: Expr::mul(Expr::var(1), Expr::k(2.0)),
                },
                AssignmentRule {
                    target: VarRef::Species(0),
                    formula: Expr::k(1.0),
                },
                AssignmentRule {
                    target: VarRef::Species(1),
                    formula: Expr::add(Expr::var(0), Expr::k(1.0)),
                },
            ],
            rates: vec![],
        };
        let order = rules.topo_sort().unwrap();
        // rule 1 (s0) must come before rule 2 (s1) which must come
        // before rule 0 (s2).
        let pos = |k: usize| order.iter().position(|&x| x == k).unwrap();
        assert!(pos(1) < pos(2));
        assert!(pos(2) < pos(0));

        let mut y = vec![0.0_f64; 3];
        let mut p: Vec<f64> = vec![];
        rules.apply_assignments(&mut y, &mut p, 0.0).unwrap();
        assert_eq!(y, vec![1.0, 2.0, 4.0]);
    }

    #[test]
    fn cyclic_assignment_rules_are_rejected() {
        // s0 := s1, s1 := s0 + 1.
        let rules = SbmlRules {
            assignments: vec![
                AssignmentRule {
                    target: VarRef::Species(0),
                    formula: Expr::var(1),
                },
                AssignmentRule {
                    target: VarRef::Species(1),
                    formula: Expr::add(Expr::var(0), Expr::k(1.0)),
                },
            ],
            rates: vec![],
        };
        assert!(rules.topo_sort().is_err());
    }

    #[test]
    fn parameter_rule_can_assign_parameter() {
        let rules = SbmlRules {
            assignments: vec![AssignmentRule {
                target: VarRef::Parameter(0),
                formula: Expr::mul(Expr::var(0), Expr::k(3.0)),
            }],
            rates: vec![],
        };
        let mut y = vec![5.0];
        let mut p = vec![0.0];
        rules.apply_assignments(&mut y, &mut p, 0.0).unwrap();
        assert_eq!(p[0], 15.0);
    }

    #[test]
    fn event_builder_sets_defaults() {
        let ev = SbmlEvent::new(
            "e1",
            Expr::ge(Expr::Time, Expr::k(5.0)),
            vec![EventAssignment {
                target: VarRef::Species(0),
                formula: Expr::k(100.0),
            }],
        );
        assert_eq!(ev.id, "e1");
        assert!(ev.delay.is_none());
        assert_eq!(ev.priority, 0.0);
        assert!(ev.use_values_from_trigger_time);
        assert!(!ev.initial_value);
    }

    #[test]
    fn event_builder_chaining_works() {
        let ev = SbmlEvent::new(
            "e1",
            Expr::Const(0.0),
            vec![],
        )
        .with_delay(2.0)
        .with_priority(5.0)
        .evaluate_at_execution_time();
        assert_eq!(ev.delay, Some(2.0));
        assert_eq!(ev.priority, 5.0);
        assert!(!ev.use_values_from_trigger_time);
    }
}
