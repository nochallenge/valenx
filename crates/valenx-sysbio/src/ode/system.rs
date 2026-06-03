//! ODE-system assembly from a reaction network — feature 8.
//!
//! A reaction-network [`Model`] becomes a system of ordinary
//! differential equations through the fundamental systems-biology
//! identity
//!
//! ```text
//!     dy/dt = S · v(y)
//! ```
//!
//! where `y` is the species-amount vector, `S` is the
//! `species × reactions` stoichiometry matrix, and `v(y)` is the
//! reaction-rate (flux) vector obtained by evaluating every
//! [`RateLaw`](crate::model::RateLaw). Constant / boundary species
//! contribute a zero derivative regardless of `S`.
//!
//! [`OdeSystem`] precomputes the stoichiometry matrix once and stores
//! a closure-free reference to the model, so the right-hand-side
//! evaluation in an integrator's inner loop is a matrix-vector product
//! with no allocation beyond the output buffer. It also provides a
//! finite-difference [`OdeSystem::jacobian`] used by the implicit BDF
//! integrator and the steady-state solver.

use crate::model::events::{RateRule, SbmlRules, VarRef};
use crate::model::expr::Expr;
use crate::model::Model;

/// A first-order ODE system `dy/dt = f(t, y)` derived from a
/// reaction-network [`Model`].
///
/// Carries the model's stoichiometric reaction kinetics *and* (when
/// present) its SBML L3 rate-rule and assignment-rule machinery, so
/// the right-hand side passed to an integrator is the same physics
/// whether the model uses bare reactions or the full rule layer.
///
/// Parameters are exposed by reference through
/// [`OdeSystem::params_mut`] so the event driver can update parameter
/// values mid-simulation (e.g. a rate rule that drives a parameter)
/// without rebuilding the system.
#[derive(Debug, Clone)]
pub struct OdeSystem {
    /// Dense stoichiometry matrix, `species x reactions`.
    stoich: Vec<Vec<f64>>,
    /// Per-species `constant` flags (boundary species => zero
    /// derivative).
    constant: Vec<bool>,
    /// Cached rate laws, one per reaction, in reaction order.
    rate_laws: Vec<crate::model::RateLaw>,
    /// Number of species (state dimension).
    n: usize,
    /// Rate rules - additive contributions to species derivatives /
    /// parameter rate of change.
    rate_rules: Vec<RateRule>,
    /// Assignment + topology cached for the rule projector.
    rules: SbmlRules,
    /// Snapshot of model parameter values. The event driver can
    /// override this slice when it needs to.
    params: Vec<f64>,
}

impl OdeSystem {
    /// Build an ODE system from a model. The model is **not** mutated
    /// and need not outlive the system - the stoichiometry matrix and
    /// the rate laws are copied in.
    pub fn from_model(model: &Model) -> Self {
        OdeSystem {
            stoich: model.stoichiometry_matrix(),
            constant: model.species.iter().map(|s| s.constant).collect(),
            rate_laws: model.reactions.iter().map(|r| r.rate_law.clone()).collect(),
            n: model.species.len(),
            rate_rules: model.rules.rates.clone(),
            rules: model.rules.clone(),
            params: model.parameter_values(),
        }
    }

    /// Snapshot of the parameter slice this system was built with.
    pub fn params(&self) -> &[f64] {
        &self.params
    }

    /// Mutable parameter slice - the event driver uses this to apply
    /// event assignments that target parameters.
    pub fn params_mut(&mut self) -> &mut Vec<f64> {
        &mut self.params
    }

    /// Replace the cached parameter slice. Returns the previous
    /// values.
    pub fn set_params(&mut self, p: Vec<f64>) -> Vec<f64> {
        std::mem::replace(&mut self.params, p)
    }

    /// Apply the model's assignment rules to a state vector in place.
    /// Returns the (cycle-checked) rule execution order.
    pub fn project_assignments(&mut self, y: &mut [f64], t: f64) -> crate::error::Result<Vec<usize>> {
        let mut p = std::mem::take(&mut self.params);
        let res = self.rules.apply_assignments(y, &mut p, t);
        self.params = p;
        res
    }

    /// Borrow the rule set for the event driver.
    pub fn rules(&self) -> &SbmlRules {
        &self.rules
    }

    /// Whether the system has any active rate or assignment rules.
    pub fn has_rules(&self) -> bool {
        !self.rate_rules.is_empty() || !self.rules.assignments.is_empty()
    }

    /// Whether any rate rule targets the global parameter `i`. The
    /// event driver uses this to decide whether the parameter must be
    /// integrated as an extra ODE state.
    pub fn parameter_has_rate_rule(&self, i: usize) -> bool {
        self.rate_rules
            .iter()
            .any(|r| matches!(&r.target, VarRef::Parameter(p) if *p == i))
    }

    /// The set of rate-rule formulas targeting species `i` (multiple
    /// rules on the same target are summed by the driver).
    pub fn rate_rule_for_species(&self, i: usize) -> Vec<&Expr> {
        self.rate_rules
            .iter()
            .filter(|r| matches!(&r.target, VarRef::Species(s) if *s == i))
            .map(|r| &r.formula)
            .collect()
    }

    /// The set of rate-rule formulas targeting parameter `i`.
    pub fn rate_rule_for_parameter(&self, i: usize) -> Vec<&Expr> {
        self.rate_rules
            .iter()
            .filter(|r| matches!(&r.target, VarRef::Parameter(p) if *p == i))
            .map(|r| &r.formula)
            .collect()
    }

    /// State dimension (number of species).
    pub fn dim(&self) -> usize {
        self.n
    }

    /// Number of reactions.
    pub fn n_reactions(&self) -> usize {
        self.rate_laws.len()
    }

    /// The reaction-rate (flux) vector `v(y)` at state `y`.
    pub fn fluxes(&self, y: &[f64]) -> Vec<f64> {
        self.rate_laws.iter().map(|law| law.rate(y)).collect()
    }

    /// The right-hand side `f(t, y) = S * v(y) + sum(rate-rules)`.
    ///
    /// A rate rule targeting species `i` adds its formula value to
    /// the species' derivative (the COPASI convention - rate rules
    /// stack on top of any stoichiometric law). A boundary species
    /// still returns zero from its stoichiometric component, but a
    /// rate rule on a boundary species *does* take effect - that is
    /// the only way SBML L3 lets a "constant" species change.
    pub fn rhs(&self, t: f64, y: &[f64]) -> Vec<f64> {
        let v = self.fluxes(y);
        let mut dydt = vec![0.0; self.n];
        let has_rate_rules = !self.rate_rules.is_empty();
        for (i, dy) in dydt.iter_mut().enumerate() {
            let mut acc = 0.0;
            if !self.constant[i] {
                let row = &self.stoich[i];
                for (j, &s) in row.iter().enumerate() {
                    if s != 0.0 {
                        acc += s * v[j];
                    }
                }
            }
            if has_rate_rules {
                for rule in &self.rate_rules {
                    if let VarRef::Species(s) = rule.target {
                        if s == i {
                            acc += rule.formula.value(y, &self.params, t);
                        }
                    }
                }
            }
            *dy = acc;
        }
        dydt
    }

    /// A central-difference approximation of the Jacobian
    /// `∂f_i/∂y_j`, returned as an `n × n` row-major matrix.
    ///
    /// The perturbation for column `j` is scaled to the magnitude of
    /// `y_j` (relative step `sqrt(machine-eps)`), which is the
    /// standard robust choice for finite-difference Jacobians. Used by
    /// the BDF integrator's Newton solve and the steady-state solver.
    pub fn jacobian(&self, t: f64, y: &[f64]) -> Vec<Vec<f64>> {
        let n = self.n;
        let mut jac = vec![vec![0.0; n]; n];
        let base_eps = 1e-7;
        let mut yp = y.to_vec();
        for j in 0..n {
            let h = base_eps * y[j].abs().max(1.0);
            yp[j] = y[j] + h;
            let f_plus = self.rhs(t, &yp);
            yp[j] = y[j] - h;
            let f_minus = self.rhs(t, &yp);
            yp[j] = y[j];
            let inv = 1.0 / (2.0 * h);
            for (i, jrow) in jac.iter_mut().enumerate() {
                jrow[j] = (f_plus[i] - f_minus[i]) * inv;
            }
        }
        jac
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{RateLaw, Reaction, Species};

    /// Linear decay A -> 0 with rate k. dA/dt = -k A.
    fn decay(k: f64, a0: f64) -> Model {
        let mut m = Model::new("decay");
        let a = m.add_species(Species::new("A", a0));
        m.add_reaction(Reaction {
            id: "d".into(),
            reactants: vec![(a, 1.0)],
            products: vec![],
            rate_law: RateLaw::MassAction {
                k,
                reactants: vec![(a, 1.0)],
            },
            reversible: false,
        });
        m
    }

    #[test]
    fn rhs_of_linear_decay() {
        let sys = OdeSystem::from_model(&decay(0.5, 8.0));
        let dy = sys.rhs(0.0, &[8.0]);
        assert!((dy[0] - (-4.0)).abs() < 1e-12);
    }

    #[test]
    fn constant_species_has_zero_derivative() {
        let mut m = decay(0.5, 8.0);
        m.species[0].constant = true;
        let sys = OdeSystem::from_model(&m);
        assert_eq!(sys.rhs(0.0, &[8.0])[0], 0.0);
    }

    #[test]
    fn jacobian_of_linear_decay_is_minus_k() {
        let sys = OdeSystem::from_model(&decay(0.7, 5.0));
        let j = sys.jacobian(0.0, &[5.0]);
        assert!((j[0][0] - (-0.7)).abs() < 1e-5);
    }
}
