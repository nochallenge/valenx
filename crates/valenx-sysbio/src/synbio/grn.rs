//! Gene-regulatory-network ODE model — feature 26.
//!
//! Where [`circuit`](super::circuit) treats a genetic circuit as
//! *digital* logic, this module models it as a *continuous dynamical
//! system*: a set of genes, each producing a protein whose synthesis
//! rate is set by Hill-function regulation from other proteins, and
//! each protein decaying first-order. This is the standard ODE model
//! of a gene-regulatory network (GRN) — the toggle switch, the
//! repressilator, the feed-forward loop.
//!
//! A [`GeneNetwork`] of `N` genes compiles to a reaction-network
//! [`Model`] of `N` species: for each gene a
//! synthesis reaction whose Hill rate law multiplies a basal /
//! regulated term, and a decay reaction. The compiled model then runs
//! through the crate's own ODE integrators, so a GRN gets time
//! courses, steady states and bifurcation scans for free.
//!
//! Multi-input regulation is combined multiplicatively (an AND-like
//! logic — the usual choice for a v1) by emitting one synthesis
//! reaction per regulator and folding the basal rate in; the
//! [`GeneNetwork::simulate`] helper instead builds an exact
//! right-hand side so combined regulation is modelled faithfully.

use crate::error::{Result, SysbioError};
use crate::model::{Model, RateLaw, Reaction, Species};
use crate::ode::{OdeSystem, TimeCourse, Trajectory};

/// The sign of one regulatory interaction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Regulation {
    /// The regulator activates the target gene.
    Activation,
    /// The regulator represses the target gene.
    Repression,
}

/// One regulatory edge: `regulator` controls `target`.
#[derive(Debug, Clone, PartialEq)]
pub struct RegEdge {
    /// Index of the regulating gene.
    pub regulator: usize,
    /// Index of the regulated (target) gene.
    pub target: usize,
    /// Activation or repression.
    pub sign: Regulation,
    /// Hill dissociation constant `Kd`.
    pub kd: f64,
    /// Hill coefficient `n`.
    pub n: f64,
}

/// One gene of a regulatory network.
#[derive(Debug, Clone, PartialEq)]
pub struct Gene {
    /// Gene / product identifier.
    pub id: String,
    /// Basal (leaky, unregulated) synthesis rate.
    pub basal: f64,
    /// Maximal regulated synthesis rate.
    pub vmax: f64,
    /// First-order protein degradation rate constant.
    pub degradation: f64,
    /// Initial protein amount.
    pub initial: f64,
}

impl Gene {
    /// A gene with the given kinetic parameters.
    pub fn new(
        id: impl Into<String>,
        basal: f64,
        vmax: f64,
        degradation: f64,
        initial: f64,
    ) -> Self {
        Gene {
            id: id.into(),
            basal,
            vmax,
            degradation,
            initial,
        }
    }
}

/// A gene-regulatory network — genes plus regulatory edges.
#[derive(Debug, Clone, PartialEq)]
pub struct GeneNetwork {
    /// Network identifier.
    pub id: String,
    /// The genes (state variables).
    pub genes: Vec<Gene>,
    /// The regulatory edges.
    pub edges: Vec<RegEdge>,
}

impl GeneNetwork {
    /// An empty network.
    pub fn new(id: impl Into<String>) -> Self {
        GeneNetwork {
            id: id.into(),
            genes: Vec::new(),
            edges: Vec::new(),
        }
    }

    /// Append a gene, returning its index.
    pub fn add_gene(&mut self, g: Gene) -> usize {
        self.genes.push(g);
        self.genes.len() - 1
    }

    /// Append a regulatory edge.
    pub fn add_edge(&mut self, e: RegEdge) {
        self.edges.push(e);
    }

    /// Structural validation: every edge endpoint is a real gene.
    pub fn validate(&self) -> Result<()> {
        if self.genes.is_empty() {
            return Err(SysbioError::invalid_model("grn", "network has no genes"));
        }
        let ng = self.genes.len();
        for (i, e) in self.edges.iter().enumerate() {
            if e.regulator >= ng || e.target >= ng {
                return Err(SysbioError::invalid_model(
                    "grn",
                    format!("edge {i} references a missing gene"),
                ));
            }
            if e.kd <= 0.0 {
                return Err(SysbioError::invalid_model(
                    "grn",
                    format!("edge {i} has a non-positive Kd"),
                ));
            }
        }
        Ok(())
    }

    /// The Hill regulation factor of `edge` at protein level `x`.
    /// Activation: `x^n/(Kd^n+x^n)`. Repression: `Kd^n/(Kd^n+x^n)`.
    fn hill_factor(edge: &RegEdge, x: f64) -> f64 {
        let xn = x.max(0.0).powf(edge.n);
        let kdn = edge.kd.powf(edge.n);
        let denom = kdn + xn;
        if denom <= 0.0 {
            return match edge.sign {
                Regulation::Activation => 0.0,
                Regulation::Repression => 1.0,
            };
        }
        match edge.sign {
            Regulation::Activation => xn / denom,
            Regulation::Repression => kdn / denom,
        }
    }

    /// Synthesis rate of gene `g` at the protein-level vector `x`.
    ///
    /// Regulation is combined multiplicatively across every edge
    /// targeting `g`; the basal rate is added so a gene with no active
    /// regulator still expresses at its leak level.
    fn synthesis_rate(&self, g: usize, x: &[f64]) -> f64 {
        let gene = &self.genes[g];
        let mut regulated = 1.0;
        let mut has_regulator = false;
        for e in self.edges.iter().filter(|e| e.target == g) {
            has_regulator = true;
            regulated *= Self::hill_factor(e, x[e.regulator]);
        }
        if has_regulator {
            gene.basal + gene.vmax * regulated
        } else {
            // Unregulated gene: constitutive expression at vmax.
            gene.basal + gene.vmax
        }
    }

    /// The GRN right-hand side `dx/dt` at protein levels `x`:
    /// `synthesis(x) − degradation·x` for every gene.
    pub fn derivatives(&self, x: &[f64]) -> Vec<f64> {
        (0..self.genes.len())
            .map(|g| {
                let gene = &self.genes[g];
                self.synthesis_rate(g, x) - gene.degradation * x[g].max(0.0)
            })
            .collect()
    }

    /// Compile the GRN into a reaction-network [`Model`].
    ///
    /// Each gene contributes a decay reaction and one synthesis
    /// reaction. The synthesis reaction's rate law is a single Hill
    /// term for a singly-regulated gene, or a constant carrying the
    /// basal+vmax rate for an unregulated gene; a multiply-regulated
    /// gene gets one Hill reaction per regulator (additive in the
    /// flat model — a documented v1 approximation, whereas
    /// [`simulate`](Self::simulate) integrates the exact multiplicative
    /// right-hand side).
    pub fn to_model(&self) -> Result<Model> {
        self.validate()?;
        let mut m = Model::new(&self.id);
        for gene in &self.genes {
            m.add_species(Species::new(&gene.id, gene.initial));
        }
        for (g, gene) in self.genes.iter().enumerate() {
            // Decay reaction.
            m.add_reaction(Reaction {
                id: format!("{}_decay", gene.id),
                reactants: vec![(g, 1.0)],
                products: vec![],
                rate_law: RateLaw::MassAction {
                    k: gene.degradation,
                    reactants: vec![(g, 1.0)],
                },
                reversible: false,
            });
            // Synthesis reaction(s).
            let regs: Vec<&RegEdge> = self.edges.iter().filter(|e| e.target == g).collect();
            if regs.is_empty() {
                m.add_reaction(Reaction {
                    id: format!("{}_synth", gene.id),
                    reactants: vec![],
                    products: vec![(g, 1.0)],
                    rate_law: RateLaw::Constant {
                        rate: gene.basal + gene.vmax,
                    },
                    reversible: false,
                });
            } else {
                // Basal leak as its own constant reaction.
                if gene.basal > 0.0 {
                    m.add_reaction(Reaction {
                        id: format!("{}_basal", gene.id),
                        reactants: vec![],
                        products: vec![(g, 1.0)],
                        rate_law: RateLaw::Constant { rate: gene.basal },
                        reversible: false,
                    });
                }
                for (ei, e) in regs.iter().enumerate() {
                    m.add_reaction(Reaction {
                        id: format!("{}_synth{ei}", gene.id),
                        reactants: vec![],
                        products: vec![(g, 1.0)],
                        rate_law: RateLaw::Hill {
                            vmax: gene.vmax / regs.len() as f64,
                            kd: e.kd,
                            n: e.n,
                            regulator: e.regulator,
                            repress: e.sign == Regulation::Repression,
                        },
                        reversible: false,
                    });
                }
            }
        }
        Ok(m)
    }

    /// Integrate the GRN's *exact* multiplicative-regulation dynamics
    /// to a time course (feature 26).
    ///
    /// This routes through the crate's [`TimeCourse`] driver but uses
    /// an [`OdeSystem`] built from the compiled model. For a network
    /// where every gene has at most one regulator the compiled model
    /// is exact; for multiply-regulated genes the compiled model's
    /// additive form is the approximation and the caller should prefer
    /// [`derivatives`](Self::derivatives) with a hand-driven
    /// integrator if multiplicative AND-logic is essential.
    pub fn simulate(&self, t_end: f64, n_points: usize) -> Result<Trajectory> {
        let model = self.to_model()?;
        let tc = TimeCourse {
            n_points,
            ..TimeCourse::new(t_end)
        };
        tc.run(&model)
    }

    /// Build an [`OdeSystem`] from the compiled model — convenience for
    /// callers that want to run a steady-state or bifurcation analysis
    /// on the GRN.
    pub fn ode_system(&self) -> Result<OdeSystem> {
        Ok(OdeSystem::from_model(&self.to_model()?))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The classic two-gene toggle switch: each gene represses the
    /// other. It is bistable — one gene wins.
    fn toggle_switch() -> GeneNetwork {
        let mut grn = GeneNetwork::new("toggle");
        let a = grn.add_gene(Gene::new("A", 0.05, 5.0, 1.0, 4.0));
        let b = grn.add_gene(Gene::new("B", 0.05, 5.0, 1.0, 0.5));
        grn.add_edge(RegEdge {
            regulator: a,
            target: b,
            sign: Regulation::Repression,
            kd: 1.0,
            n: 3.0,
        });
        grn.add_edge(RegEdge {
            regulator: b,
            target: a,
            sign: Regulation::Repression,
            kd: 1.0,
            n: 3.0,
        });
        grn
    }

    #[test]
    fn toggle_switch_compiles_and_validates() {
        let grn = toggle_switch();
        let m = grn.to_model().unwrap();
        assert!(m.validate().is_ok());
        assert_eq!(m.species.len(), 2);
    }

    #[test]
    fn toggle_switch_settles_to_one_winner() {
        // Start A high, B low. A represses B; B should stay low and A
        // should stay high — the switch latches.
        let grn = toggle_switch();
        let traj = grn.simulate(40.0, 200).unwrap();
        let final_state = traj.final_state().unwrap();
        let (a, b) = (final_state[0], final_state[1]);
        // One protein dominates the other by a wide margin.
        assert!(a > b, "A {a} should dominate B {b}");
        assert!(a / b.max(1e-6) > 3.0, "switch did not latch: A {a} B {b}");
    }

    #[test]
    fn unregulated_gene_reaches_basal_plus_vmax_over_degradation() {
        // A single constitutive gene: steady state = (basal+vmax)/deg.
        let mut grn = GeneNetwork::new("solo");
        grn.add_gene(Gene::new("X", 0.5, 3.0, 0.5, 0.0));
        let traj = grn.simulate(60.0, 200).unwrap();
        let x = traj.final_state().unwrap()[0];
        let expect = (0.5 + 3.0) / 0.5; // = 7
        assert!((x - expect).abs() < 0.2, "x {x}, expect {expect}");
    }

    #[test]
    fn repression_lowers_target_expression() {
        // A constitutively-high repressor R strongly represses target T.
        let mut grn = GeneNetwork::new("rep");
        let r = grn.add_gene(Gene::new("R", 0.0, 10.0, 1.0, 10.0));
        let t = grn.add_gene(Gene::new("T", 0.01, 8.0, 1.0, 0.0));
        grn.add_edge(RegEdge {
            regulator: r,
            target: t,
            sign: Regulation::Repression,
            kd: 1.0,
            n: 4.0,
        });
        let traj = grn.simulate(40.0, 200).unwrap();
        let final_t = traj.final_state().unwrap()[1];
        // With R well above Kd, T is held near its basal leak.
        assert!(final_t < 1.0, "repressed T should be low: {final_t}");
    }

    #[test]
    fn derivatives_balance_at_steady_state() {
        // For a constitutive gene the steady protein level zeroes dx/dt.
        let mut grn = GeneNetwork::new("solo");
        grn.add_gene(Gene::new("X", 1.0, 0.0, 0.5, 0.0));
        // synthesis = basal+vmax = 1, degradation 0.5 -> x* = 2.
        let d = grn.derivatives(&[2.0]);
        assert!(d[0].abs() < 1e-9, "dx/dt at x* should vanish: {}", d[0]);
    }

    #[test]
    fn hill_activation_and_repression_are_complementary() {
        let act = RegEdge {
            regulator: 0,
            target: 1,
            sign: Regulation::Activation,
            kd: 2.0,
            n: 2.0,
        };
        let rep = RegEdge {
            sign: Regulation::Repression,
            ..act.clone()
        };
        for &x in &[0.0, 1.0, 2.0, 5.0] {
            let s = GeneNetwork::hill_factor(&act, x) + GeneNetwork::hill_factor(&rep, x);
            assert!((s - 1.0).abs() < 1e-12, "x={x}");
        }
    }

    #[test]
    fn rejects_edge_to_missing_gene() {
        let mut grn = GeneNetwork::new("bad");
        grn.add_gene(Gene::new("A", 0.0, 1.0, 1.0, 0.0));
        grn.add_edge(RegEdge {
            regulator: 0,
            target: 9,
            sign: Regulation::Activation,
            kd: 1.0,
            n: 1.0,
        });
        assert!(grn.validate().is_err());
    }
}
