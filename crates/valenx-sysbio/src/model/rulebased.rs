//! Rule-based model expansion — feature 7 (BioNetGen-class v1).
//!
//! BioNetGen and PySB describe biochemistry with *rules* over
//! *molecules with sites and states* rather than with an explicit
//! reaction list: one rule such as
//! `"phosphorylate kinase site Y"` generates every concrete reaction
//! consistent with the pattern. The full BioNetGen language has
//! bonds, compartments, wildcards and a graph-isomorphism matcher;
//! that is a research codebase on its own.
//!
//! This module implements a **tractable v1** that captures the
//! defining BioNetGen idea — *combinatorial state expansion* — for the
//! single most common case: a set of molecules, each with a fixed list
//! of modification **sites**, every site taking one of a small set of
//! discrete **states**. The "species" of the expanded model are the
//! full state vectors; the rules are site-state transitions.
//!
//! Concretely, a [`RuleModel`] holds:
//!
//! - molecule **templates** (a name + an ordered list of sites, each
//!   site with its allowed states), and
//! - **rules**, each of which selects a molecule, requires a site to
//!   be in one state, and moves it to another — optionally gated by a
//!   second molecule acting as an enzyme.
//!
//! [`RuleModel::expand`] enumerates the reachable [`Microstate`]s by a
//! breadth-first flood from the declared seed states, instantiates one
//! [`Species`] per microstate and one mass-action [`Reaction`] per
//! `(rule, applicable microstate)` pair, and returns a flat
//! [`Model`]. That flat model then feeds every other engine in the
//! crate unchanged — exactly the BioNetGen "network generation then
//! simulate" workflow.
//!
//! ## v1 caveats
//!
//! No bonds / complexes (a microstate is one molecule, never a dimer),
//! no wildcards, no rate-law expressions beyond a per-rule constant,
//! and the reachable set is capped by [`RuleModel::max_species`] to
//! keep a combinatorial blow-up bounded. These limits are documented
//! and surfaced — [`RuleModel::expand`] returns
//! [`SysbioError::Invalid`] rather than silently truncating when the
//! cap is hit.

use std::collections::BTreeMap;

use crate::error::{Result, SysbioError};
use crate::model::{Model, RateLaw, Reaction, Species};

/// A modification site on a molecule template: a name and the discrete
/// states it may occupy (e.g. `["u", "p"]` for un/phosphorylated).
#[derive(Debug, Clone, PartialEq)]
pub struct Site {
    /// Site name (unique within its molecule).
    pub name: String,
    /// Allowed states for this site, in canonical order.
    pub states: Vec<String>,
}

impl Site {
    /// A site with the given name and state list.
    pub fn new(name: impl Into<String>, states: &[&str]) -> Self {
        Site {
            name: name.into(),
            states: states.iter().map(|s| s.to_string()).collect(),
        }
    }
}

/// A molecule template — a name plus an ordered list of sites.
#[derive(Debug, Clone, PartialEq)]
pub struct MoleculeTemplate {
    /// Molecule name.
    pub name: String,
    /// Modification sites, in canonical order.
    pub sites: Vec<Site>,
}

/// A site-state transition rule.
///
/// "On molecule `molecule`, when site `site` is in state `from`, move
/// it to state `to` at rate `rate`." If `enzyme` is set, the rule is a
/// bimolecular mass-action reaction whose rate also depends on the
/// enzyme microstate's amount.
#[derive(Debug, Clone, PartialEq)]
pub struct Rule {
    /// Rule identifier (used as a reaction-id prefix).
    pub id: String,
    /// Index of the substrate molecule in [`RuleModel::molecules`].
    pub molecule: usize,
    /// Site name on that molecule.
    pub site: String,
    /// Required state before the rule fires.
    pub from: String,
    /// State after the rule fires.
    pub to: String,
    /// Mass-action rate constant.
    pub rate: f64,
    /// Optional catalysing molecule index (`None` = unimolecular).
    pub enzyme: Option<usize>,
}

/// A concrete microstate: a molecule index plus the chosen state of
/// each of its sites (parallel to the template's `sites`).
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Microstate {
    /// Molecule template index.
    pub molecule: usize,
    /// State of each site, in the template's site order.
    pub site_states: Vec<String>,
}

/// A rule-based model awaiting network expansion.
#[derive(Debug, Clone, PartialEq)]
pub struct RuleModel {
    /// Model identifier.
    pub id: String,
    /// Molecule templates.
    pub molecules: Vec<MoleculeTemplate>,
    /// Transition rules.
    pub rules: Vec<Rule>,
    /// Seed microstates with their initial amounts.
    pub seeds: Vec<(Microstate, f64)>,
    /// Hard ceiling on the number of generated species — guards
    /// against combinatorial explosion.
    pub max_species: usize,
}

impl RuleModel {
    /// An empty rule model with a sane default species cap (10 000).
    pub fn new(id: impl Into<String>) -> Self {
        RuleModel {
            id: id.into(),
            molecules: Vec::new(),
            rules: Vec::new(),
            seeds: Vec::new(),
            max_species: 10_000,
        }
    }

    /// Append a molecule template, returning its index.
    pub fn add_molecule(&mut self, m: MoleculeTemplate) -> usize {
        self.molecules.push(m);
        self.molecules.len() - 1
    }

    /// Append a rule.
    pub fn add_rule(&mut self, r: Rule) {
        self.rules.push(r);
    }

    /// Declare a seed microstate and its initial amount.
    pub fn add_seed(&mut self, state: Microstate, amount: f64) {
        self.seeds.push((state, amount));
    }

    /// Apply a single rule to a microstate. Returns the product
    /// microstate if the rule's precondition is met, else `None`.
    fn apply(&self, rule: &Rule, state: &Microstate) -> Option<Microstate> {
        if state.molecule != rule.molecule {
            return None;
        }
        let tmpl = &self.molecules[rule.molecule];
        let site_idx = tmpl.sites.iter().position(|s| s.name == rule.site)?;
        if state.site_states[site_idx] != rule.from {
            return None;
        }
        let mut next = state.clone();
        next.site_states[site_idx] = rule.to.clone();
        Some(next)
    }

    /// Expand the rules into a flat [`Model`].
    ///
    /// Breadth-first flood from the seeds: every rule is tried against
    /// every newly discovered microstate, products are enqueued, and a
    /// mass-action reaction is recorded for each firing. Errors if the
    /// model is malformed or if the species cap is exceeded.
    pub fn expand(&self) -> Result<Model> {
        if self.molecules.is_empty() {
            return Err(SysbioError::invalid_model(
                "rule_model",
                "no molecule templates",
            ));
        }
        if self.seeds.is_empty() {
            return Err(SysbioError::invalid_model(
                "rule_model",
                "no seed microstates",
            ));
        }
        // Validate rule references up front.
        for r in &self.rules {
            let tmpl = self.molecules.get(r.molecule).ok_or_else(|| {
                SysbioError::invalid_model("rule_model", format!("rule {} bad molecule", r.id))
            })?;
            let site = tmpl
                .sites
                .iter()
                .find(|s| s.name == r.site)
                .ok_or_else(|| {
                    SysbioError::invalid_model(
                        "rule_model",
                        format!("rule {} references missing site {}", r.id, r.site),
                    )
                })?;
            for st in [&r.from, &r.to] {
                if !site.states.contains(st) {
                    return Err(SysbioError::invalid_model(
                        "rule_model",
                        format!("rule {} uses state {st} not declared on site", r.id),
                    ));
                }
            }
        }

        // BFS over reachable microstates. `index` assigns a stable
        // species index; `frontier` drives the flood.
        let mut index: BTreeMap<Microstate, usize> = BTreeMap::new();
        let mut order: Vec<Microstate> = Vec::new();
        let mut frontier: Vec<Microstate> = Vec::new();
        let mut initials: BTreeMap<Microstate, f64> = BTreeMap::new();

        for (state, amt) in &self.seeds {
            // Seed states must have the right arity.
            let tmpl = self.molecules.get(state.molecule).ok_or_else(|| {
                SysbioError::invalid_model("rule_model", "seed references missing molecule")
            })?;
            if state.site_states.len() != tmpl.sites.len() {
                return Err(SysbioError::invalid_model(
                    "rule_model",
                    "seed microstate site count mismatch",
                ));
            }
            if !index.contains_key(state) {
                index.insert(state.clone(), order.len());
                order.push(state.clone());
                frontier.push(state.clone());
            }
            *initials.entry(state.clone()).or_insert(0.0) += amt;
        }

        while let Some(state) = frontier.pop() {
            for rule in &self.rules {
                if let Some(next) = self.apply(rule, &state) {
                    if !index.contains_key(&next) {
                        if order.len() >= self.max_species {
                            return Err(SysbioError::invalid(
                                "max_species",
                                format!(
                                    "rule expansion exceeded the {} species cap",
                                    self.max_species
                                ),
                            ));
                        }
                        index.insert(next.clone(), order.len());
                        order.push(next.clone());
                        frontier.push(next);
                    }
                }
            }
        }

        // Build the flat model: one species per microstate.
        let mut model = Model::new(&self.id);
        for st in &order {
            let mol = &self.molecules[st.molecule].name;
            let label = format!("{mol}({})", st.site_states.join(","));
            let init = initials.get(st).copied().unwrap_or(0.0);
            model.add_species(Species::new(label, init));
        }
        // One reaction per (rule, applicable microstate).
        let mut rxn_counter = 0usize;
        for rule in &self.rules {
            for st in &order {
                if let Some(next) = self.apply(rule, st) {
                    let from_idx = index[st];
                    let to_idx = index[&next];
                    let mut reactants = vec![(from_idx, 1.0)];
                    let mut products = vec![(to_idx, 1.0)];
                    let mut law_reactants = vec![(from_idx, 1.0)];
                    // An enzyme catalyses: appears on both sides, and
                    // its amount enters the mass-action rate.
                    if let Some(enz_mol) = rule.enzyme {
                        // Use the first reachable microstate of the
                        // enzyme molecule as the catalytic species.
                        if let Some(enz_state) = order.iter().find(|s| s.molecule == enz_mol) {
                            let enz_idx = index[enz_state];
                            reactants.push((enz_idx, 1.0));
                            products.push((enz_idx, 1.0));
                            law_reactants.push((enz_idx, 1.0));
                        }
                    }
                    model.add_reaction(Reaction {
                        id: format!("{}_{}", rule.id, rxn_counter),
                        reactants,
                        products,
                        rate_law: RateLaw::MassAction {
                            k: rule.rate,
                            reactants: law_reactants,
                        },
                        reversible: false,
                    });
                    rxn_counter += 1;
                }
            }
        }
        model.validate()?;
        Ok(model)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A single kinase with two phosphorylation sites, each u/p.
    /// Two rules phosphorylate site_a and site_b. From the all-`u`
    /// seed the reachable set is the 2×2 grid of (a,b) states.
    fn two_site_kinase() -> RuleModel {
        let mut rm = RuleModel::new("kinase2");
        let k = rm.add_molecule(MoleculeTemplate {
            name: "K".into(),
            sites: vec![Site::new("a", &["u", "p"]), Site::new("b", &["u", "p"])],
        });
        rm.add_rule(Rule {
            id: "phos_a".into(),
            molecule: k,
            site: "a".into(),
            from: "u".into(),
            to: "p".into(),
            rate: 0.5,
            enzyme: None,
        });
        rm.add_rule(Rule {
            id: "phos_b".into(),
            molecule: k,
            site: "b".into(),
            from: "u".into(),
            to: "p".into(),
            rate: 0.3,
            enzyme: None,
        });
        rm.add_seed(
            Microstate {
                molecule: k,
                site_states: vec!["u".into(), "u".into()],
            },
            100.0,
        );
        rm
    }

    #[test]
    fn expansion_enumerates_full_state_grid() {
        let rm = two_site_kinase();
        let m = rm.expand().expect("expand");
        // 2 sites x 2 states => 4 microstates.
        assert_eq!(m.species.len(), 4);
        // From uu: 2 rules. From up & pu: 1 each. From pp: 0. => 4.
        assert_eq!(m.reactions.len(), 4);
        assert!(m.validate().is_ok());
        // Seed amount landed on the uu species.
        let uu = m
            .species
            .iter()
            .find(|s| s.id == "K(u,u)")
            .expect("uu species");
        assert_eq!(uu.initial, 100.0);
    }

    #[test]
    fn species_cap_is_enforced() {
        let mut rm = two_site_kinase();
        rm.max_species = 2;
        let err = rm.expand().unwrap_err();
        assert_eq!(err.code(), "sysbio.invalid");
    }

    #[test]
    fn missing_site_in_rule_is_rejected() {
        let mut rm = two_site_kinase();
        rm.rules[0].site = "ghost".into();
        assert!(rm.expand().is_err());
    }

    #[test]
    fn enzyme_catalysed_rule_keeps_enzyme_on_both_sides() {
        let mut rm = RuleModel::new("cat");
        let sub = rm.add_molecule(MoleculeTemplate {
            name: "S".into(),
            sites: vec![Site::new("s", &["u", "p"])],
        });
        let enz = rm.add_molecule(MoleculeTemplate {
            name: "E".into(),
            sites: vec![Site::new("x", &["on"])],
        });
        rm.add_rule(Rule {
            id: "cat".into(),
            molecule: sub,
            site: "s".into(),
            from: "u".into(),
            to: "p".into(),
            rate: 1.0,
            enzyme: Some(enz),
        });
        rm.add_seed(
            Microstate {
                molecule: sub,
                site_states: vec!["u".into()],
            },
            10.0,
        );
        rm.add_seed(
            Microstate {
                molecule: enz,
                site_states: vec!["on".into()],
            },
            1.0,
        );
        let m = rm.expand().expect("expand");
        let r = &m.reactions[0];
        // Enzyme appears once as reactant and once as product.
        assert_eq!(r.reactants.len(), 2);
        assert_eq!(r.products.len(), 2);
    }

    #[test]
    fn empty_rule_model_errors() {
        let rm = RuleModel::new("empty");
        assert!(rm.expand().is_err());
    }
}
