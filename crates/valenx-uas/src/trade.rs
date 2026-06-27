//! Trade studies: sweep design parameters, evaluate performance, and extract
//! the **Pareto (non-dominated) front** over chosen objectives.
//!
//! A trade study is deliberately generic over the design and the objectives:
//! a [`DesignPoint`] is a named bag of design parameters, the caller supplies
//! an *evaluator* `Fn(&DesignPoint) -> Result<Vec<f64>, E>` that turns each
//! point into a vector of objective values, and [`TradeStudy::run`] returns
//! every evaluated point plus the Pareto front.
//!
//! ## Pareto dominance
//!
//! Each objective is tagged [`Objective::Maximize`] or [`Objective::Minimize`].
//! Design `a` **dominates** `b` when `a` is no worse than `b` on *every*
//! objective and strictly better on *at least one*. The **Pareto front** is the
//! set of points that no other point dominates — the designs where you cannot
//! improve one objective without sacrificing another. A 2-design study where one
//! design is worse on every objective drops the dominated design and keeps only
//! the winner (pinned in the tests).
//!
//! Non-finite objective values are treated as *dominated by everything*: a
//! design that fails to evaluate to a finite objective never appears on the
//! front (it cannot dominate, and any finite design dominates it). This keeps a
//! single bad evaluation from poisoning the front.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::error::UasError;

/// The optimisation sense of one objective.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Objective {
    /// Larger is better (e.g. endurance, range, payload).
    Maximize,
    /// Smaller is better (e.g. mass, cost, hover power).
    Minimize,
}

/// A single design point: a label plus its named scalar design parameters.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DesignPoint {
    /// A human-readable label for this design (e.g. `"R=0.15,m=1.5"`).
    pub label: String,
    /// Named design parameters (e.g. `rotor_radius_m`, `mass_kg`).
    pub params: BTreeMap<String, f64>,
}

impl DesignPoint {
    /// Build a design point from a label and `(name, value)` parameters.
    pub fn new(label: impl Into<String>, params: impl IntoIterator<Item = (String, f64)>) -> Self {
        Self {
            label: label.into(),
            params: params.into_iter().collect(),
        }
    }

    /// Look up a design parameter by name.
    pub fn get(&self, name: &str) -> Option<f64> {
        self.params.get(name).copied()
    }
}

/// One evaluated design: the design point and its objective values, in the same
/// order as the study's [`Objective`] list.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EvaluatedDesign {
    /// The design that was evaluated.
    pub design: DesignPoint,
    /// Objective values, aligned with the study's objective senses.
    pub objectives: Vec<f64>,
}

impl EvaluatedDesign {
    /// True if every objective value is finite.
    fn all_finite(&self) -> bool {
        self.objectives.iter().all(|v| v.is_finite())
    }
}

/// The result of a trade study: all evaluated designs and the indices /
/// designs on the Pareto front.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ParetoFront {
    /// Every design that was evaluated (sweep order preserved).
    pub all: Vec<EvaluatedDesign>,
    /// Indices into [`all`](Self::all) of the non-dominated designs.
    pub front_indices: Vec<usize>,
}

impl ParetoFront {
    /// The non-dominated designs themselves (a convenience over
    /// [`front_indices`](Self::front_indices)).
    pub fn front(&self) -> Vec<&EvaluatedDesign> {
        self.front_indices.iter().map(|&i| &self.all[i]).collect()
    }

    /// Number of designs on the Pareto front.
    pub fn front_len(&self) -> usize {
        self.front_indices.len()
    }
}

/// A trade study over a set of design points and a set of objectives.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TradeStudy {
    /// The objective senses (length must match each evaluation's objective
    /// vector).
    pub objectives: Vec<Objective>,
    /// The design points to sweep.
    pub designs: Vec<DesignPoint>,
}

impl TradeStudy {
    /// Build a trade study from its objectives and design points.
    pub fn new(objectives: Vec<Objective>, designs: Vec<DesignPoint>) -> Self {
        Self {
            objectives,
            designs,
        }
    }

    /// Evaluate every design with `evaluator`, then compute the Pareto front.
    ///
    /// `evaluator` maps a [`DesignPoint`] to its objective values, in the same
    /// order as [`objectives`](Self::objectives). Its error type `E` is
    /// converted into [`UasError`] (so a vehicle-performance evaluator that
    /// already returns [`UasError`] composes directly).
    ///
    /// # Errors
    ///
    /// - [`UasError::EmptyStudy`] if there are no designs or no objectives.
    /// - Any error returned by `evaluator` for a design.
    /// - [`UasError::NotFinite`] is *not* raised for non-finite objectives —
    ///   those are kept in `all` but excluded from the front (see the module
    ///   docs); only an evaluator-side failure aborts the run.
    pub fn run<E>(
        &self,
        mut evaluator: impl FnMut(&DesignPoint) -> Result<Vec<f64>, E>,
    ) -> Result<ParetoFront, UasError>
    where
        UasError: From<E>,
    {
        if self.designs.is_empty() {
            return Err(UasError::EmptyStudy {
                what: "design points",
            });
        }
        if self.objectives.is_empty() {
            return Err(UasError::EmptyStudy { what: "objectives" });
        }

        let mut all = Vec::with_capacity(self.designs.len());
        for design in &self.designs {
            let objectives = evaluator(design)?;
            // A mismatched objective count is a programming error in the
            // evaluator; treat it as a non-finite (always-dominated) row by
            // padding/truncation would hide the bug, so reject loudly.
            if objectives.len() != self.objectives.len() {
                return Err(UasError::NotFinite {
                    quantity: "objective-vector length mismatch",
                    value: objectives.len() as f64,
                });
            }
            all.push(EvaluatedDesign {
                design: design.clone(),
                objectives,
            });
        }

        let front_indices = self.pareto_indices(&all);
        Ok(ParetoFront { all, front_indices })
    }

    /// Compute the indices of the non-dominated designs.
    fn pareto_indices(&self, all: &[EvaluatedDesign]) -> Vec<usize> {
        let mut front = Vec::new();
        for (i, cand) in all.iter().enumerate() {
            // A design with any non-finite objective can never be on the front.
            if !cand.all_finite() {
                continue;
            }
            let dominated = all.iter().enumerate().any(|(j, other)| {
                j != i && other.all_finite() && self.dominates(&other.objectives, &cand.objectives)
            });
            if !dominated {
                front.push(i);
            }
        }
        front
    }

    /// Does objective vector `a` dominate `b`? (`a` no worse on all, strictly
    /// better on at least one, per the per-objective sense). Both are assumed
    /// finite and the same length as `self.objectives`.
    fn dominates(&self, a: &[f64], b: &[f64]) -> bool {
        let mut strictly_better_somewhere = false;
        for ((sense, &av), &bv) in self.objectives.iter().zip(a.iter()).zip(b.iter()) {
            let (better, worse) = match sense {
                Objective::Maximize => (av > bv, av < bv),
                Objective::Minimize => (av < bv, av > bv),
            };
            if worse {
                return false; // a is worse on this objective -> cannot dominate
            }
            if better {
                strictly_better_somewhere = true;
            }
        }
        strictly_better_somewhere
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vehicle::{Battery, MultirotorUas, SEA_LEVEL_AIR_DENSITY};

    fn dp(label: &str, pairs: &[(&str, f64)]) -> DesignPoint {
        DesignPoint::new(label, pairs.iter().map(|(k, v)| ((*k).to_string(), *v)))
    }

    // ---- BENCHMARK PIN: a 2-design sweep drops the dominated design --------
    #[test]
    fn two_design_sweep_drops_the_dominated_design() {
        // Objectives: maximize endurance, maximize payload.
        // Design A: endurance 30, payload 2.  Design B: 20, 1 (worse on both).
        // A dominates B, so the front is {A} only.
        let study = TradeStudy::new(
            vec![Objective::Maximize, Objective::Maximize],
            vec![dp("A", &[]), dp("B", &[])],
        );
        let table = [("A", vec![30.0, 2.0]), ("B", vec![20.0, 1.0])];
        let front = study
            .run::<UasError>(|d| Ok(table.iter().find(|(l, _)| *l == d.label).unwrap().1.clone()))
            .unwrap();
        assert_eq!(front.front_len(), 1, "exactly one non-dominated design");
        assert_eq!(front.front()[0].design.label, "A");
        assert_eq!(front.all.len(), 2, "both designs are still recorded");
    }

    #[test]
    fn genuine_tradeoff_keeps_both_designs() {
        // A: high endurance, low payload. B: low endurance, high payload.
        // Neither dominates -> both on the front.
        let study = TradeStudy::new(
            vec![Objective::Maximize, Objective::Maximize],
            vec![dp("A", &[]), dp("B", &[])],
        );
        let table = [("A", vec![40.0, 1.0]), ("B", vec![20.0, 3.0])];
        let front = study
            .run::<UasError>(|d| Ok(table.iter().find(|(l, _)| *l == d.label).unwrap().1.clone()))
            .unwrap();
        assert_eq!(front.front_len(), 2);
    }

    #[test]
    fn minimize_objective_sense_is_respected() {
        // Maximize range, MINIMIZE mass.
        // A: range 100, mass 2. B: range 100, mass 5 (heavier, same range).
        // A dominates B (same range, lighter). Front = {A}.
        let study = TradeStudy::new(
            vec![Objective::Maximize, Objective::Minimize],
            vec![dp("A", &[]), dp("B", &[])],
        );
        let table = [("A", vec![100.0, 2.0]), ("B", vec![100.0, 5.0])];
        let front = study
            .run::<UasError>(|d| Ok(table.iter().find(|(l, _)| *l == d.label).unwrap().1.clone()))
            .unwrap();
        assert_eq!(front.front_len(), 1);
        assert_eq!(front.front()[0].design.label, "A");
    }

    #[test]
    fn non_finite_objective_is_excluded_from_front() {
        // A finite design and a design that evaluates to NaN: the NaN design is
        // never on the front, and does not knock the finite one off.
        let study = TradeStudy::new(
            vec![Objective::Maximize],
            vec![dp("good", &[]), dp("bad", &[])],
        );
        let front = study
            .run::<UasError>(|d| {
                Ok(if d.label == "good" {
                    vec![10.0]
                } else {
                    vec![f64::NAN]
                })
            })
            .unwrap();
        assert_eq!(front.front_len(), 1);
        assert_eq!(front.front()[0].design.label, "good");
    }

    #[test]
    fn empty_study_is_rejected() {
        let no_designs = TradeStudy::new(vec![Objective::Maximize], vec![]);
        assert!(no_designs.run::<UasError>(|_| Ok(vec![1.0])).is_err());
        let no_objectives = TradeStudy::new(vec![], vec![dp("A", &[])]);
        assert!(no_objectives.run::<UasError>(|_| Ok(vec![])).is_err());
    }

    #[test]
    fn objective_length_mismatch_is_rejected() {
        let study = TradeStudy::new(
            vec![Objective::Maximize, Objective::Minimize],
            vec![dp("A", &[])],
        );
        // Evaluator returns one value but two objectives are declared.
        assert!(study.run::<UasError>(|_| Ok(vec![1.0])).is_err());
    }

    #[test]
    fn end_to_end_rotor_radius_sweep_real_vehicles() {
        // Sweep multirotor rotor radius; objectives: maximize hover endurance,
        // minimize disk loading. Drives the real composed performance model
        // and confirms the front is non-empty and a strict subset is possible.
        let radii = [0.10, 0.15, 0.20, 0.25];
        let designs: Vec<DesignPoint> = radii
            .iter()
            .map(|r| dp(&format!("R={r}"), &[("rotor_radius_m", *r)]))
            .collect();
        let study = TradeStudy::new(vec![Objective::Maximize, Objective::Minimize], designs);
        let front = study
            .run::<UasError>(|d| {
                let r = d.get("rotor_radius_m").unwrap();
                let uas = MultirotorUas::new(
                    4,
                    r,
                    1.5,
                    0.7,
                    SEA_LEVEL_AIR_DENSITY,
                    Battery::new(100.0, 0.8).unwrap(),
                    0.3,
                    0.75,
                )?;
                Ok(vec![uas.hover_endurance_s()?, uas.rotor.disk_loading()])
            })
            .unwrap();
        assert_eq!(front.all.len(), 4);
        assert!(front.front_len() >= 1);
        // Bigger rotors -> more endurance AND lower disk loading, so larger
        // radius dominates smaller on BOTH objectives: the front collapses to
        // just the largest radius.
        assert_eq!(front.front_len(), 1);
        assert_eq!(front.front()[0].design.label, "R=0.25");
    }
}
