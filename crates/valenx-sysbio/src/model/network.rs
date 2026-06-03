//! Reaction-network model — feature 1.
//!
//! The [`Model`] is the central data structure of the crate: a list of
//! [`Species`] (each in a [`Compartment`]), a list of [`Reaction`]s
//! (each carrying a stoichiometry and a [`RateLaw`]),
//! and a table of named global [`Parameter`]s. Every downstream
//! engine — the ODE assembler, the three stochastic simulators, the
//! FBA solver — consumes a `Model`.
//!
//! Design notes:
//!
//! - Species, reactions and parameters are stored in `Vec`s; their
//!   array index *is* their identity in the numerical layer. A
//!   `name → index` map is kept alongside for the SBML reader and for
//!   human-facing queries.
//! - Stoichiometry is stored sparsely on each reaction as
//!   `(species index, signed coefficient)` pairs. The dense
//!   stoichiometry matrix `S` is materialised on demand by
//!   [`Model::stoichiometry_matrix`].
//! - A species may be flagged `constant` (a boundary / buffered
//!   species, e.g. an SBML boundary condition): its amount never
//!   changes, so the ODE assembler drops its row.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::error::{Result, SysbioError};
use crate::model::events::{SbmlEvent, SbmlRules, VarRef};
use crate::model::RateLaw;

/// A reaction compartment — a named volume.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Compartment {
    /// Unique compartment identifier.
    pub id: String,
    /// Compartment size (volume). Amounts are interpreted as molecule
    /// counts or concentrations depending on the caller's convention;
    /// the size matters for the SBML round-trip and for converting
    /// concentrations to counts in the SSA.
    pub size: f64,
}

impl Compartment {
    /// A unit-volume compartment with the given id.
    pub fn new(id: impl Into<String>) -> Self {
        Compartment {
            id: id.into(),
            size: 1.0,
        }
    }
}

/// A chemical species — one state variable of the model.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Species {
    /// Unique species identifier.
    pub id: String,
    /// Index of the owning compartment in [`Model::compartments`].
    pub compartment: usize,
    /// Initial amount (count or concentration).
    pub initial: f64,
    /// If `true` the amount is held fixed (an SBML boundary species or
    /// a buffered metabolite). The ODE assembler omits its row.
    pub constant: bool,
}

impl Species {
    /// A non-constant species in compartment 0 with the given initial
    /// amount.
    pub fn new(id: impl Into<String>, initial: f64) -> Self {
        Species {
            id: id.into(),
            compartment: 0,
            initial,
            constant: false,
        }
    }

    /// Builder: mark this species as a fixed boundary species.
    pub fn constant(mut self) -> Self {
        self.constant = true;
        self
    }

    /// Builder: place this species in compartment `c`.
    pub fn in_compartment(mut self, c: usize) -> Self {
        self.compartment = c;
        self
    }
}

/// A named global parameter (a rate constant, a Km, …).
///
/// Parameters give the analysis layer something to *scan* and to
/// *differentiate against* — the sensitivity and bifurcation modules
/// perturb parameter values, not the resolved constants baked into a
/// [`RateLaw`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Parameter {
    /// Unique parameter identifier.
    pub id: String,
    /// Current numeric value.
    pub value: f64,
}

impl Parameter {
    /// A parameter with the given id and value.
    pub fn new(id: impl Into<String>, value: f64) -> Self {
        Parameter {
            id: id.into(),
            value,
        }
    }
}

/// A reaction: a stoichiometry plus a kinetic law.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Reaction {
    /// Unique reaction identifier.
    pub id: String,
    /// Reactant `(species index, positive multiplicity)` pairs.
    pub reactants: Vec<(usize, f64)>,
    /// Product `(species index, positive multiplicity)` pairs.
    pub products: Vec<(usize, f64)>,
    /// Kinetic rate law (the *forward* flux for an irreversible
    /// reaction; FBA treats `reversible` separately).
    pub rate_law: RateLaw,
    /// Whether the reaction may carry negative flux (used by FBA to
    /// pick a default lower bound).
    pub reversible: bool,
}

impl Reaction {
    /// The signed stoichiometric change of each species in this
    /// reaction: `+products − reactants`, summed over duplicates.
    pub fn net_stoichiometry(&self) -> BTreeMap<usize, f64> {
        let mut net: BTreeMap<usize, f64> = BTreeMap::new();
        for &(i, c) in &self.reactants {
            *net.entry(i).or_insert(0.0) -= c;
        }
        for &(i, c) in &self.products {
            *net.entry(i).or_insert(0.0) += c;
        }
        net
    }
}

/// A complete reaction-network model.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct Model {
    /// Optional model identifier.
    pub id: String,
    /// All compartments.
    pub compartments: Vec<Compartment>,
    /// All species (state variables). Index == identity.
    pub species: Vec<Species>,
    /// All reactions. Index == identity.
    pub reactions: Vec<Reaction>,
    /// All global parameters.
    pub parameters: Vec<Parameter>,
    /// SBML L3 discrete events fired on rising-edge trigger
    /// crossings. The time-course driver evaluates each event's
    /// trigger between integrator steps. Empty by default - a model
    /// with no events behaves exactly as before this field was added.
    #[serde(default)]
    pub events: Vec<SbmlEvent>,
    /// SBML L3 assignment + rate rules. Applied every integrator
    /// output (assignment rules) and folded into the ODE RHS (rate
    /// rules). Empty by default - same back-compat note as `events`.
    #[serde(default)]
    pub rules: SbmlRules,
}

impl Model {
    /// An empty model with the given id and one unit compartment.
    pub fn new(id: impl Into<String>) -> Self {
        Model {
            id: id.into(),
            compartments: vec![Compartment::new("default")],
            species: Vec::new(),
            reactions: Vec::new(),
            parameters: Vec::new(),
            events: Vec::new(),
            rules: SbmlRules::default(),
        }
    }

    /// Append an event, returning its index.
    pub fn add_event(&mut self, e: SbmlEvent) -> usize {
        self.events.push(e);
        self.events.len() - 1
    }

    /// Append an assignment rule. The driver will topologically sort
    /// the assignment-rule set every time it sweeps them, so order of
    /// insertion is not significant.
    pub fn add_assignment_rule(&mut self, rule: crate::model::events::AssignmentRule) {
        self.rules.assignments.push(rule);
    }

    /// Append a rate rule. Multiple rate rules on the same target are
    /// summed by the driver.
    pub fn add_rate_rule(&mut self, rule: crate::model::events::RateRule) {
        self.rules.rates.push(rule);
    }

    /// Snapshot of the global parameter values as a parallel `f64`
    /// vector, one entry per [`Parameter`]. Used by the rule and
    /// event evaluators - they read parameter slots by index, so the
    /// driver can update parameter values without re-cloning the
    /// model.
    pub fn parameter_values(&self) -> Vec<f64> {
        self.parameters.iter().map(|p| p.value).collect()
    }

    /// Append a species, returning its index.
    pub fn add_species(&mut self, s: Species) -> usize {
        self.species.push(s);
        self.species.len() - 1
    }

    /// Append a reaction, returning its index.
    pub fn add_reaction(&mut self, r: Reaction) -> usize {
        self.reactions.push(r);
        self.reactions.len() - 1
    }

    /// Append a parameter, returning its index.
    pub fn add_parameter(&mut self, p: Parameter) -> usize {
        self.parameters.push(p);
        self.parameters.len() - 1
    }

    /// Index of the species with id `id`, if any.
    pub fn species_index(&self, id: &str) -> Option<usize> {
        self.species.iter().position(|s| s.id == id)
    }

    /// Index of the parameter with id `id`, if any.
    pub fn parameter_index(&self, id: &str) -> Option<usize> {
        self.parameters.iter().position(|p| p.id == id)
    }

    /// The initial-amount vector, one entry per species, in index
    /// order. This is the `y0` handed to every integrator.
    pub fn initial_state(&self) -> Vec<f64> {
        self.species.iter().map(|s| s.initial).collect()
    }

    /// Structural validation. Checks that every stoichiometry index is
    /// in range and that every rate-law species reference is in range.
    /// Returns [`SysbioError::InvalidModel`] on the first problem.
    pub fn validate(&self) -> Result<()> {
        let ns = self.species.len();
        let nc = self.compartments.len();
        if ns == 0 {
            return Err(SysbioError::invalid_model(
                "reaction_network",
                "model has no species",
            ));
        }
        for (i, s) in self.species.iter().enumerate() {
            if s.compartment >= nc {
                return Err(SysbioError::invalid_model(
                    "reaction_network",
                    format!("species {i} ({}) references missing compartment", s.id),
                ));
            }
        }
        for (ri, r) in self.reactions.iter().enumerate() {
            for &(i, c) in r.reactants.iter().chain(r.products.iter()) {
                if i >= ns {
                    return Err(SysbioError::invalid_model(
                        "reaction_network",
                        format!("reaction {ri} ({}) references species index {i}", r.id),
                    ));
                }
                if c < 0.0 {
                    return Err(SysbioError::invalid_model(
                        "reaction_network",
                        format!("reaction {ri} ({}) has negative stoichiometry", r.id),
                    ));
                }
            }
            for dep in r.rate_law.dependencies() {
                if dep >= ns {
                    return Err(SysbioError::invalid_model(
                        "reaction_network",
                        format!("reaction {ri} ({}) rate law reads species {dep}", r.id),
                    ));
                }
            }
            // Biochemical-constant validity: unphysical rate-law constants
            // (a negative Michaelis constant, a non-positive Hill
            // dissociation constant, or a negative Hill coefficient) would
            // silently produce wrong rates, so reject them at validation
            // rather than mask them at evaluation time.
            match &r.rate_law {
                RateLaw::MichaelisMenten { km, .. } => {
                    if !km.is_finite() || *km < 0.0 {
                        return Err(SysbioError::invalid_model(
                            "reaction_network",
                            format!(
                                "reaction {ri} ({}) Michaelis constant Km must be finite and >= 0",
                                r.id
                            ),
                        ));
                    }
                }
                RateLaw::Hill { kd, n, .. } => {
                    if !kd.is_finite() || *kd <= 0.0 {
                        return Err(SysbioError::invalid_model(
                            "reaction_network",
                            format!(
                                "reaction {ri} ({}) Hill dissociation constant Kd must be finite and > 0",
                                r.id
                            ),
                        ));
                    }
                    if !n.is_finite() || *n < 0.0 {
                        return Err(SysbioError::invalid_model(
                            "reaction_network",
                            format!(
                                "reaction {ri} ({}) Hill coefficient n must be finite and >= 0",
                                r.id
                            ),
                        ));
                    }
                }
                _ => {}
            }
        }
        // Events: trigger + assignment formulas must index valid
        // species / parameters; assignment targets must be in range.
        let np = self.parameters.len();
        let check_ref = |r: &VarRef, ctx: &str| -> Result<()> {
            match r {
                VarRef::Species(i) => {
                    if *i >= ns {
                        return Err(SysbioError::invalid_model(
                            "reaction_network",
                            format!("{ctx} references missing species index {i}"),
                        ));
                    }
                }
                VarRef::Parameter(i) => {
                    if *i >= np {
                        return Err(SysbioError::invalid_model(
                            "reaction_network",
                            format!("{ctx} references missing parameter index {i}"),
                        ));
                    }
                }
            }
            Ok(())
        };
        let check_expr_indices = |e: &crate::model::expr::Expr, ctx: &str| -> Result<()> {
            for v in e.var_deps() {
                if v >= ns {
                    return Err(SysbioError::invalid_model(
                        "reaction_network",
                        format!("{ctx} reads missing species index {v}"),
                    ));
                }
            }
            for p in e.param_deps() {
                if p >= np {
                    return Err(SysbioError::invalid_model(
                        "reaction_network",
                        format!("{ctx} reads missing parameter index {p}"),
                    ));
                }
            }
            Ok(())
        };
        for (ei, ev) in self.events.iter().enumerate() {
            let ctx = format!("event {ei} ({})", ev.id);
            check_expr_indices(&ev.trigger, &format!("{ctx} trigger"))?;
            for (ai, a) in ev.assignments.iter().enumerate() {
                check_ref(&a.target, &format!("{ctx} assignment {ai}"))?;
                check_expr_indices(&a.formula, &format!("{ctx} assignment {ai}"))?;
            }
            if let Some(d) = ev.delay {
                if d < 0.0 {
                    return Err(SysbioError::invalid_model(
                        "reaction_network",
                        format!("{ctx} has a negative delay"),
                    ));
                }
            }
        }
        for (i, rule) in self.rules.assignments.iter().enumerate() {
            let ctx = format!("assignment rule {i}");
            check_ref(&rule.target, &ctx)?;
            check_expr_indices(&rule.formula, &ctx)?;
        }
        for (i, rule) in self.rules.rates.iter().enumerate() {
            let ctx = format!("rate rule {i}");
            check_ref(&rule.target, &ctx)?;
            check_expr_indices(&rule.formula, &ctx)?;
        }
        // Verify the assignment-rule graph is acyclic.
        if !self.rules.assignments.is_empty() {
            self.rules.topo_sort()?;
        }
        Ok(())
    }

    /// The dense stoichiometry matrix `S` with `species × reactions`
    /// shape, column-major in the returned nested `Vec`.
    /// `S[i][j]` is the net change of species `i` per unit flux of
    /// reaction `j`. Used by the ODE assembler, the conservation
    /// analysis and FBA.
    pub fn stoichiometry_matrix(&self) -> Vec<Vec<f64>> {
        let ns = self.species.len();
        let nr = self.reactions.len();
        let mut s = vec![vec![0.0; nr]; ns];
        for (j, r) in self.reactions.iter().enumerate() {
            for (i, coeff) in r.net_stoichiometry() {
                s[i][j] = coeff;
            }
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A·B ⇌ C built by hand. Returns the model.
    fn ab_to_c() -> Model {
        let mut m = Model::new("ab2c");
        let a = m.add_species(Species::new("A", 10.0));
        let b = m.add_species(Species::new("B", 5.0));
        let c = m.add_species(Species::new("C", 0.0));
        m.add_reaction(Reaction {
            id: "r1".into(),
            reactants: vec![(a, 1.0), (b, 1.0)],
            products: vec![(c, 1.0)],
            rate_law: RateLaw::MassAction {
                k: 0.1,
                reactants: vec![(a, 1.0), (b, 1.0)],
            },
            reversible: false,
        });
        m
    }

    #[test]
    fn build_and_validate() {
        let m = ab_to_c();
        assert!(m.validate().is_ok());
        assert_eq!(m.species.len(), 3);
        assert_eq!(m.initial_state(), vec![10.0, 5.0, 0.0]);
    }

    #[test]
    fn validate_rejects_unphysical_kinetic_constants() {
        // A negative Michaelis constant is unphysical and would silently
        // produce wrong rates -> validation must reject it.
        let mut m = Model::new("bad-km");
        let s = m.add_species(Species::new("S", 1.0));
        let p = m.add_species(Species::new("P", 0.0));
        m.add_reaction(Reaction {
            id: "mm".into(),
            reactants: vec![(s, 1.0)],
            products: vec![(p, 1.0)],
            rate_law: RateLaw::MichaelisMenten {
                vmax: 1.0,
                km: -1.0,
                substrate: s,
            },
            reversible: false,
        });
        assert!(m.validate().is_err(), "negative Km must fail validation");

        // A non-positive Hill dissociation constant is unphysical too.
        let mut m2 = Model::new("bad-kd");
        let g = m2.add_species(Species::new("G", 1.0));
        let x = m2.add_species(Species::new("X", 0.0));
        m2.add_reaction(Reaction {
            id: "hill".into(),
            reactants: vec![],
            products: vec![(x, 1.0)],
            rate_law: RateLaw::Hill {
                vmax: 1.0,
                kd: 0.0,
                n: 2.0,
                regulator: g,
                repress: false,
            },
            reversible: false,
        });
        assert!(m2.validate().is_err(), "non-positive Hill Kd must fail validation");
    }

    #[test]
    fn stoichiometry_matrix_signs() {
        let m = ab_to_c();
        let s = m.stoichiometry_matrix();
        // species A (0): consumed -> -1; B (1): -1; C (2): +1.
        assert_eq!(s[0][0], -1.0);
        assert_eq!(s[1][0], -1.0);
        assert_eq!(s[2][0], 1.0);
    }

    #[test]
    fn validate_rejects_out_of_range_index() {
        let mut m = Model::new("bad");
        m.add_species(Species::new("A", 1.0));
        m.add_reaction(Reaction {
            id: "r".into(),
            reactants: vec![(7, 1.0)],
            products: vec![],
            rate_law: RateLaw::Constant { rate: 1.0 },
            reversible: false,
        });
        assert!(m.validate().is_err());
    }

    #[test]
    fn net_stoichiometry_merges_duplicates() {
        // 2 A -> A + B  collapses to net A: -1, B: +1.
        let r = Reaction {
            id: "r".into(),
            reactants: vec![(0, 2.0)],
            products: vec![(0, 1.0), (1, 1.0)],
            rate_law: RateLaw::Constant { rate: 1.0 },
            reversible: false,
        };
        let net = r.net_stoichiometry();
        assert_eq!(net[&0], -1.0);
        assert_eq!(net[&1], 1.0);
    }
}
