//! Flux-balance analysis — features 22 and 23.
//!
//! Flux-balance analysis (FBA) is the workhorse of constraint-based
//! metabolic modelling — the COBRA toolbox's central method. It does
//! *not* simulate dynamics; instead it asks: at metabolic steady
//! state, what flux distribution maximises some cellular objective
//! (classically biomass production)?
//!
//! The steady-state assumption is `S·v = 0` — every internal
//! metabolite is produced as fast as it is consumed. Together with
//! per-reaction flux bounds `lo ≤ v ≤ hi` (thermodynamics,
//! irreversibility, nutrient-uptake limits) and a linear objective
//! `cᵀ·v`, this is exactly a linear program, solved here by the
//! crate's own [`simplex`](super::simplex) solver.
//!
//! - [`FbaProblem`] — a metabolic model: a stoichiometry matrix, flux
//!   bounds and an objective vector.
//! - [`FbaProblem::optimize`] — solve the FBA LP (feature 22).
//! - [`FbaProblem::flux_variability`] — flux variability analysis:
//!   the min and max each reaction can carry at (near-)optimal growth
//!   (feature 23).
//! - [`FbaProblem::parsimonious`] — parsimonious FBA (pFBA): among all
//!   flux distributions achieving the optimal objective, the one of
//!   minimum total absolute flux — the "least enzyme" solution
//!   (feature 23).

use crate::error::{Result, SysbioError};
use crate::fba::simplex::{ConstraintSense, LinearProgram, LpStatus};
use crate::model::Model;

/// A constraint-based metabolic model ready for FBA.
#[derive(Debug, Clone, PartialEq)]
pub struct FbaProblem {
    /// Metabolite (row) names.
    pub metabolites: Vec<String>,
    /// Reaction (column) names.
    pub reactions: Vec<String>,
    /// Stoichiometry matrix `S`, row-major (`metabolites × reactions`).
    pub stoich: Vec<Vec<f64>>,
    /// Per-reaction lower flux bounds.
    pub lb: Vec<f64>,
    /// Per-reaction upper flux bounds.
    pub ub: Vec<f64>,
    /// Objective coefficients (usually a 1 on the biomass reaction).
    pub objective: Vec<f64>,
}

/// The result of an FBA optimisation.
#[derive(Debug, Clone, PartialEq)]
pub struct FbaSolution {
    /// `true` if an optimum was found.
    pub feasible: bool,
    /// Optimal objective value.
    pub objective_value: f64,
    /// Optimal flux through each reaction.
    pub fluxes: Vec<f64>,
}

/// One row of a flux-variability table.
#[derive(Debug, Clone, PartialEq)]
pub struct FluxRange {
    /// Reaction name.
    pub reaction: String,
    /// Minimum flux at the (near-)optimal objective.
    pub min: f64,
    /// Maximum flux at the (near-)optimal objective.
    pub max: f64,
}

impl FbaProblem {
    /// Build an FBA problem from a reaction-network [`Model`].
    ///
    /// The model's stoichiometry matrix becomes `S`; every reaction is
    /// given the bounds `[-1000, 1000]` (reversible) or `[0, 1000]`
    /// (irreversible) — the COBRA default flux scale. The objective is
    /// left empty for the caller to set with
    /// [`FbaProblem::set_objective`].
    pub fn from_model(model: &Model) -> Result<Self> {
        model.validate()?;
        let stoich = model.stoichiometry_matrix();
        let nr = model.reactions.len();
        let lb = model
            .reactions
            .iter()
            .map(|r| if r.reversible { -1000.0 } else { 0.0 })
            .collect();
        Ok(FbaProblem {
            metabolites: model.species.iter().map(|s| s.id.clone()).collect(),
            reactions: model.reactions.iter().map(|r| r.id.clone()).collect(),
            stoich,
            lb,
            ub: vec![1000.0; nr],
            objective: vec![0.0; nr],
        })
    }

    /// Set the objective to maximise flux through one named reaction.
    pub fn set_objective(&mut self, reaction: &str) -> Result<()> {
        let idx = self
            .reactions
            .iter()
            .position(|r| r == reaction)
            .ok_or_else(|| SysbioError::invalid("reaction", "unknown objective reaction"))?;
        self.objective = vec![0.0; self.reactions.len()];
        self.objective[idx] = 1.0;
        Ok(())
    }

    /// Constrain reaction `reaction`'s flux to `[lo, hi]`.
    pub fn set_bounds(&mut self, reaction: &str, lo: f64, hi: f64) -> Result<()> {
        let idx = self
            .reactions
            .iter()
            .position(|r| r == reaction)
            .ok_or_else(|| SysbioError::invalid("reaction", "unknown reaction"))?;
        if lo > hi {
            return Err(SysbioError::invalid("bounds", "lower bound exceeds upper"));
        }
        self.lb[idx] = lo;
        self.ub[idx] = hi;
        Ok(())
    }

    /// Structural check: matrix shape consistent with the bound /
    /// objective vectors.
    fn check_shape(&self) -> Result<()> {
        let nr = self.reactions.len();
        if self.lb.len() != nr || self.ub.len() != nr || self.objective.len() != nr {
            return Err(SysbioError::invalid_model(
                "fba",
                "bound / objective length disagrees with reaction count",
            ));
        }
        if self.stoich.len() != self.metabolites.len()
            || self.stoich.iter().any(|r| r.len() != nr)
        {
            return Err(SysbioError::invalid_model(
                "fba",
                "stoichiometry matrix shape mismatch",
            ));
        }
        Ok(())
    }

    /// The base LP: `S·v = 0`, `lb ≤ v ≤ ub`.
    fn base_lp(&self, objective: &[f64]) -> LinearProgram {
        LinearProgram {
            c: objective.to_vec(),
            a: self.stoich.clone(),
            b: vec![0.0; self.metabolites.len()],
            sense: vec![ConstraintSense::Eq; self.metabolites.len()],
            lo: self.lb.clone(),
            hi: self.ub.clone(),
        }
    }

    /// Solve the FBA LP — maximise the objective at steady state
    /// (feature 22).
    pub fn optimize(&self) -> Result<FbaSolution> {
        self.check_shape()?;
        let lp = self.base_lp(&self.objective);
        let sol = lp.solve();
        match sol.status {
            LpStatus::Optimal => Ok(FbaSolution {
                feasible: true,
                objective_value: sol.objective,
                fluxes: sol.x,
            }),
            LpStatus::Infeasible => Ok(FbaSolution {
                feasible: false,
                objective_value: 0.0,
                fluxes: Vec::new(),
            }),
            LpStatus::Unbounded => Err(SysbioError::not_converged(
                "simplex",
                "FBA objective is unbounded — check the flux bounds",
            )),
        }
    }

    /// Flux variability analysis (feature 23).
    ///
    /// First solves the FBA LP for the optimal objective `z*`, then —
    /// for every reaction — solves two further LPs (minimise and
    /// maximise that reaction's flux) under the added constraint that
    /// the objective stays at least `gamma · z*`. `gamma = 1.0` pins
    /// growth to the exact optimum; `gamma = 0.9` explores the flux
    /// space within 90 % of optimal growth.
    pub fn flux_variability(&self, gamma: f64) -> Result<Vec<FluxRange>> {
        self.check_shape()?;
        if !(0.0..=1.0).contains(&gamma) {
            return Err(SysbioError::invalid("gamma", "fraction must be in [0, 1]"));
        }
        let opt = self.optimize()?;
        if !opt.feasible {
            return Err(SysbioError::not_converged(
                "simplex",
                "FVA base problem is infeasible",
            ));
        }
        let nr = self.reactions.len();
        let z_star = opt.objective_value;

        // Augment the base LP with the objective-floor row.
        let mut a = self.stoich.clone();
        let mut b = vec![0.0; self.metabolites.len()];
        let mut sense = vec![ConstraintSense::Eq; self.metabolites.len()];
        a.push(self.objective.clone());
        // c·v >= gamma z*  (handle a negative z* by flipping sense).
        let floor = gamma * z_star;
        if z_star >= 0.0 {
            b.push(floor);
            sense.push(ConstraintSense::Ge);
        } else {
            b.push(floor);
            sense.push(ConstraintSense::Le);
        }

        let mut ranges = Vec::with_capacity(nr);
        for j in 0..nr {
            let mut c = vec![0.0; nr];
            c[j] = 1.0;
            // Maximise reaction j.
            let lp_max = LinearProgram {
                c: c.clone(),
                a: a.clone(),
                b: b.clone(),
                sense: sense.clone(),
                lo: self.lb.clone(),
                hi: self.ub.clone(),
            };
            let max = match lp_max.solve().status {
                LpStatus::Optimal => lp_max.solve().objective,
                _ => self.ub[j],
            };
            // Minimise reaction j (maximise -j).
            let neg_c: Vec<f64> = c.iter().map(|v| -v).collect();
            let lp_min = LinearProgram {
                c: neg_c,
                a: a.clone(),
                b: b.clone(),
                sense: sense.clone(),
                lo: self.lb.clone(),
                hi: self.ub.clone(),
            };
            let min = match lp_min.solve().status {
                LpStatus::Optimal => -lp_min.solve().objective,
                _ => self.lb[j],
            };
            ranges.push(FluxRange {
                reaction: self.reactions[j].clone(),
                min,
                max,
            });
        }
        Ok(ranges)
    }

    /// Parsimonious FBA — pFBA (feature 23).
    ///
    /// Among all flux distributions that achieve the optimal objective
    /// `z*`, pFBA returns the one minimising the **total absolute
    /// flux** `Σ|v_j|` — biologically, the solution that uses the
    /// least total enzyme. The L1 minimisation is linearised in the
    /// standard way: each reaction's flux `v_j` is split into a
    /// non-negative forward part `v⁺` and reverse part `v⁻` with
    /// `v_j = v⁺ − v⁻`, and `Σ(v⁺ + v⁻)` is minimised subject to the
    /// objective staying at `z*`.
    pub fn parsimonious(&self) -> Result<FbaSolution> {
        self.check_shape()?;
        let opt = self.optimize()?;
        if !opt.feasible {
            return Ok(FbaSolution {
                feasible: false,
                objective_value: 0.0,
                fluxes: Vec::new(),
            });
        }
        let nr = self.reactions.len();
        let nm = self.metabolites.len();
        let z_star = opt.objective_value;

        // Variables: [v_1..v_nr, p_1..p_nr, m_1..m_nr]
        //   v_j is the flux (kept with its own bounds),
        //   p_j >= 0, m_j >= 0 are the split magnitudes.
        let total_vars = 3 * nr;

        // Objective: minimise sum(p + m) -> maximise -(sum p + m).
        let mut c = vec![0.0; total_vars];
        for j in 0..nr {
            c[nr + j] = -1.0; // p
            c[2 * nr + j] = -1.0; // m
        }

        let mut a: Vec<Vec<f64>> = Vec::new();
        let mut b: Vec<f64> = Vec::new();
        let mut sense: Vec<ConstraintSense> = Vec::new();

        // Steady state on v: S·v = 0.
        for (i, srow) in self.stoich.iter().enumerate() {
            let mut row = vec![0.0; total_vars];
            row[..nr].copy_from_slice(&srow[..nr]);
            a.push(row);
            b.push(0.0);
            sense.push(ConstraintSense::Eq);
            let _ = i;
        }
        // Link rows: v_j - p_j + m_j = 0.
        for j in 0..nr {
            let mut row = vec![0.0; total_vars];
            row[j] = 1.0;
            row[nr + j] = -1.0;
            row[2 * nr + j] = 1.0;
            a.push(row);
            b.push(0.0);
            sense.push(ConstraintSense::Eq);
        }
        // Objective floor: c_obj · v = z*  (pin growth exactly).
        let mut obj_row = vec![0.0; total_vars];
        obj_row[..nr].copy_from_slice(&self.objective[..nr]);
        a.push(obj_row);
        b.push(z_star);
        sense.push(ConstraintSense::Eq);
        let _ = nm;

        // Bounds: v keeps lb/ub; p, m in [0, large].
        let big = 1e6;
        let mut lo = vec![0.0; total_vars];
        let mut hi = vec![big; total_vars];
        lo[..nr].copy_from_slice(&self.lb[..nr]);
        hi[..nr].copy_from_slice(&self.ub[..nr]);

        let lp = LinearProgram {
            c,
            a,
            b,
            sense,
            lo,
            hi,
        };
        let sol = lp.solve();
        match sol.status {
            LpStatus::Optimal => Ok(FbaSolution {
                feasible: true,
                objective_value: z_star,
                fluxes: sol.x[..nr].to_vec(),
            }),
            _ => Ok(FbaSolution {
                feasible: false,
                objective_value: 0.0,
                fluxes: Vec::new(),
            }),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A tiny linear pathway:  Aex --R1--> A --R2--> B --R3--> Bex.
    /// One internal metabolite each for A and B; R1 imports, R3
    /// exports. With R1 capped at 10, all fluxes equal 10 at optimum
    /// and biomass-as-R3 is maximised at 10.
    fn linear_pathway() -> FbaProblem {
        // metabolites: A, B.  reactions: R1 (->A), R2 (A->B), R3 (B->).
        FbaProblem {
            metabolites: vec!["A".into(), "B".into()],
            reactions: vec!["R1".into(), "R2".into(), "R3".into()],
            stoich: vec![
                // A:  +1 from R1, -1 from R2,  0 from R3
                vec![1.0, -1.0, 0.0],
                // B:   0 from R1, +1 from R2, -1 from R3
                vec![0.0, 1.0, -1.0],
            ],
            lb: vec![0.0, 0.0, 0.0],
            ub: vec![10.0, 1000.0, 1000.0],
            objective: vec![0.0, 0.0, 1.0], // maximise export R3
        }
    }

    #[test]
    fn fba_maximises_through_bottleneck() {
        let p = linear_pathway();
        let sol = p.optimize().unwrap();
        assert!(sol.feasible);
        assert!((sol.objective_value - 10.0).abs() < 1e-6);
        // Steady state forces every flux to the import cap.
        for f in &sol.fluxes {
            assert!((f - 10.0).abs() < 1e-6, "flux {f}");
        }
    }

    #[test]
    fn fba_from_model_builds_consistent_problem() {
        use crate::model::{RateLaw, Reaction, Species};
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
        let mut p = FbaProblem::from_model(&m).unwrap();
        p.set_bounds("imp", 0.0, 5.0).unwrap();
        p.set_objective("exp").unwrap();
        let sol = p.optimize().unwrap();
        assert!((sol.objective_value - 5.0).abs() < 1e-6);
    }

    #[test]
    fn fva_pins_bottleneck_but_frees_nothing_extra() {
        let p = linear_pathway();
        let ranges = p.flux_variability(1.0).unwrap();
        assert_eq!(ranges.len(), 3);
        // At gamma = 1 every reaction is pinned to 10 (the only flux
        // distribution achieving optimal export).
        for r in &ranges {
            assert!((r.min - 10.0).abs() < 1e-5, "{} min {}", r.reaction, r.min);
            assert!((r.max - 10.0).abs() < 1e-5, "{} max {}", r.reaction, r.max);
        }
    }

    #[test]
    fn fva_relaxed_growth_widens_ranges() {
        let p = linear_pathway();
        // At 50 % of optimal growth the import flux may drop below 10.
        let ranges = p.flux_variability(0.5).unwrap();
        let r1 = ranges.iter().find(|r| r.reaction == "R1").unwrap();
        assert!(r1.min < 10.0 - 1e-6, "R1 min should relax: {}", r1.min);
        assert!((r1.max - 10.0).abs() < 1e-5);
    }

    #[test]
    fn parsimonious_fba_matches_optimum_and_is_minimal() {
        let p = linear_pathway();
        let p_sol = p.parsimonious().unwrap();
        assert!(p_sol.feasible);
        assert!((p_sol.objective_value - 10.0).abs() < 1e-6);
        // The pathway is linear so the parsimonious solution is the
        // same all-10 distribution.
        for f in &p_sol.fluxes {
            assert!((f - 10.0).abs() < 1e-5, "flux {f}");
        }
    }

    #[test]
    fn parsimonious_prunes_a_wasteful_cycle() {
        // A pathway with a redundant internal cycle: pFBA should
        // carry zero flux through the cycle while FBA may not.
        // metabolites A, B. reactions: imp(->A,<=10), conv(A->B),
        // exp(B->), cyc_f(A->B), cyc_r(B->A).
        let p = FbaProblem {
            metabolites: vec!["A".into(), "B".into()],
            reactions: vec![
                "imp".into(),
                "conv".into(),
                "exp".into(),
                "cyc_f".into(),
                "cyc_r".into(),
            ],
            stoich: vec![
                // A
                vec![1.0, -1.0, 0.0, -1.0, 1.0],
                // B
                vec![0.0, 1.0, -1.0, 1.0, -1.0],
            ],
            lb: vec![0.0; 5],
            ub: vec![10.0, 1000.0, 1000.0, 1000.0, 1000.0],
            objective: vec![0.0, 0.0, 1.0, 0.0, 0.0],
        };
        let par = p.parsimonious().unwrap();
        assert!(par.feasible);
        assert!((par.objective_value - 10.0).abs() < 1e-6);
        // The futile cycle carries (near) zero net enzyme in pFBA.
        assert!(par.fluxes[3] < 1e-4, "cyc_f flux {}", par.fluxes[3]);
        assert!(par.fluxes[4] < 1e-4, "cyc_r flux {}", par.fluxes[4]);
    }

    #[test]
    fn rejects_unknown_objective_reaction() {
        let mut p = linear_pathway();
        assert!(p.set_objective("ghost").is_err());
    }
}
