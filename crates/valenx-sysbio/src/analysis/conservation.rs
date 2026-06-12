//! Conservation analysis — feature 20.
//!
//! A *conserved moiety* is a linear combination of species amounts
//! that no reaction can change — total enzyme (free + bound), total
//! phosphate, total DNA template. Mathematically a conservation law
//! is a left null vector `g` of the stoichiometry matrix `S`:
//! `gᵀ·S = 0` implies `gᵀ·dy/dt = gᵀ·S·v = 0`, so `gᵀ·y` is constant
//! for all time.
//!
//! [`conservation_laws`] computes a basis for the left null space of
//! `S` (equivalently the right null space of `Sᵀ`) via SVD and
//! reports each law as a [`ConservationLaw`] — the coefficient vector
//! plus the conserved total evaluated at the model's initial state.
//! Knowing the moieties lets a simulator reduce the ODE dimension and
//! lets the steady-state solver understand why its Jacobian is
//! singular.

use crate::error::Result;
use crate::model::Model;
use crate::ode::linalg::null_space;

/// One conservation law of a reaction network.
#[derive(Debug, Clone, PartialEq)]
pub struct ConservationLaw {
    /// Coefficient on each species (`coefficients.len() == #species`).
    /// The combination `Σ coefficients[i]·amount[i]` is time-invariant.
    pub coefficients: Vec<f64>,
    /// The conserved total evaluated at the model's initial amounts.
    pub conserved_total: f64,
}

impl ConservationLaw {
    /// Indices of species that participate in this law (non-negligible
    /// coefficient), for a human-readable summary.
    pub fn participants(&self) -> Vec<usize> {
        self.coefficients
            .iter()
            .enumerate()
            .filter(|(_, &c)| c.abs() > 1e-9)
            .map(|(i, _)| i)
            .collect()
    }

    /// Evaluate the conserved combination at an arbitrary state.
    pub fn evaluate(&self, state: &[f64]) -> f64 {
        self.coefficients
            .iter()
            .zip(state)
            .map(|(c, x)| c * x)
            .sum()
    }
}

/// Compute the conserved moieties of `model` (feature 20).
///
/// Returns one [`ConservationLaw`] per independent conservation
/// relation. An empty result means the network conserves nothing (its
/// stoichiometry matrix has full row rank). `tol` is the relative
/// singular-value threshold passed to the null-space routine.
pub fn conservation_laws(model: &Model, tol: f64) -> Result<Vec<ConservationLaw>> {
    model.validate()?;
    let s = model.stoichiometry_matrix(); // species x reactions
    let ns = model.species.len();
    if ns == 0 {
        return Ok(Vec::new());
    }
    // Left null vectors of S = right null vectors of S^T.
    // Build S^T as reactions x species.
    let nr = model.reactions.len();
    let mut st = vec![vec![0.0; ns]; nr];
    for (i, row) in s.iter().enumerate() {
        for (j, &v) in row.iter().enumerate() {
            st[j][i] = v;
        }
    }
    // A reactionless model conserves every species independently.
    if nr == 0 {
        let init = model.initial_state();
        return Ok((0..ns)
            .map(|i| {
                let mut c = vec![0.0; ns];
                c[i] = 1.0;
                ConservationLaw {
                    coefficients: c,
                    conserved_total: init[i],
                }
            })
            .collect());
    }
    let basis = null_space(&st, tol);
    let init = model.initial_state();
    let mut laws = Vec::with_capacity(basis.len());
    for vec in basis {
        // Normalise so the largest-magnitude coefficient is +1 — gives
        // human-readable "A + B = const"-style laws.
        let max_abs = vec
            .iter()
            .cloned()
            .fold(0.0_f64, |m, v| m.max(v.abs()))
            .max(1e-300);
        let sign = vec
            .iter()
            .find(|v| v.abs() > 1e-9)
            .map(|v| v.signum())
            .unwrap_or(1.0);
        let coefficients: Vec<f64> = vec
            .iter()
            .map(|v| {
                let scaled = v / max_abs * sign;
                // Snap near-integers for readability.
                let r = scaled.round();
                if (scaled - r).abs() < 1e-6 {
                    r
                } else {
                    scaled
                }
            })
            .collect();
        let conserved_total = coefficients.iter().zip(&init).map(|(c, x)| c * x).sum();
        laws.push(ConservationLaw {
            coefficients,
            conserved_total,
        });
    }
    Ok(laws)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{RateLaw, Reaction, Species};

    /// A <-> B isomerisation: A + B is conserved.
    fn isomerise() -> Model {
        let mut m = Model::new("iso");
        let a = m.add_species(Species::new("A", 7.0));
        let b = m.add_species(Species::new("B", 3.0));
        m.add_reaction(Reaction {
            id: "fwd".into(),
            reactants: vec![(a, 1.0)],
            products: vec![(b, 1.0)],
            rate_law: RateLaw::MassAction {
                k: 1.0,
                reactants: vec![(a, 1.0)],
            },
            reversible: true,
        });
        m.add_reaction(Reaction {
            id: "rev".into(),
            reactants: vec![(b, 1.0)],
            products: vec![(a, 1.0)],
            rate_law: RateLaw::MassAction {
                k: 1.0,
                reactants: vec![(b, 1.0)],
            },
            reversible: true,
        });
        m
    }

    #[test]
    fn isomerisation_has_one_conservation_law() {
        let m = isomerise();
        let laws = conservation_laws(&m, 1e-9).unwrap();
        assert_eq!(laws.len(), 1);
        let law = &laws[0];
        // Both species participate with equal-magnitude coefficients.
        assert!((law.coefficients[0].abs() - law.coefficients[1].abs()).abs() < 1e-6);
        // Conserved total is A0 + B0 = 10 (up to the overall sign).
        assert!((law.conserved_total.abs() - 10.0).abs() < 1e-6);
    }

    #[test]
    fn conservation_law_holds_under_state_change() {
        let m = isomerise();
        let law = &conservation_laws(&m, 1e-9).unwrap()[0];
        // Moving 4 units A->B leaves the conserved total unchanged.
        let v0 = law.evaluate(&[7.0, 3.0]);
        let v1 = law.evaluate(&[3.0, 7.0]);
        assert!((v0 - v1).abs() < 1e-9);
    }

    #[test]
    fn open_system_conserves_nothing() {
        // 0 -> A -> 0 has full-rank S; no conserved moiety.
        let mut m = Model::new("open");
        let a = m.add_species(Species::new("A", 1.0));
        m.add_reaction(Reaction {
            id: "in".into(),
            reactants: vec![],
            products: vec![(a, 1.0)],
            rate_law: RateLaw::Constant { rate: 1.0 },
            reversible: false,
        });
        m.add_reaction(Reaction {
            id: "out".into(),
            reactants: vec![(a, 1.0)],
            products: vec![],
            rate_law: RateLaw::MassAction {
                k: 1.0,
                reactants: vec![(a, 1.0)],
            },
            reversible: false,
        });
        let laws = conservation_laws(&m, 1e-9).unwrap();
        assert!(laws.is_empty());
    }

    #[test]
    fn enzyme_kinetics_has_two_moieties() {
        // E + S <-> ES -> E + P. Conserved: E + ES and S + ES + P.
        let mut m = Model::new("mm");
        let e = m.add_species(Species::new("E", 1.0));
        let s = m.add_species(Species::new("S", 10.0));
        let es = m.add_species(Species::new("ES", 0.0));
        let p = m.add_species(Species::new("P", 0.0));
        // bind
        m.add_reaction(Reaction {
            id: "bind".into(),
            reactants: vec![(e, 1.0), (s, 1.0)],
            products: vec![(es, 1.0)],
            rate_law: RateLaw::MassAction {
                k: 1.0,
                reactants: vec![(e, 1.0), (s, 1.0)],
            },
            reversible: false,
        });
        // unbind
        m.add_reaction(Reaction {
            id: "unbind".into(),
            reactants: vec![(es, 1.0)],
            products: vec![(e, 1.0), (s, 1.0)],
            rate_law: RateLaw::MassAction {
                k: 0.5,
                reactants: vec![(es, 1.0)],
            },
            reversible: false,
        });
        // catalysis
        m.add_reaction(Reaction {
            id: "cat".into(),
            reactants: vec![(es, 1.0)],
            products: vec![(e, 1.0), (p, 1.0)],
            rate_law: RateLaw::MassAction {
                k: 2.0,
                reactants: vec![(es, 1.0)],
            },
            reversible: false,
        });
        let laws = conservation_laws(&m, 1e-9).unwrap();
        // 4 species, rank(S) = 2 -> 2 conservation laws.
        assert_eq!(laws.len(), 2);
    }
}
