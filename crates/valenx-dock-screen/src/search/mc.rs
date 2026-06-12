//! Feature 10 — Monte-Carlo / simulated-annealing pose search.
//!
//! AutoDock's original (pre-LGA) search was simulated annealing: a
//! single walker perturbs its pose, scores the candidate, and accepts
//! it by the Metropolis criterion
//!
//! ```text
//! accept  if  ΔE < 0
//! accept  with probability exp(-ΔE / kT)  otherwise
//! ```
//!
//! while the temperature `T` is slowly lowered. High temperature early
//! on lets the walker escape local minima; low temperature late on
//! settles it into a basin.
//!
//! This module implements both the constant-temperature Metropolis
//! Monte-Carlo case and the cooled (annealed) case via the
//! [`McSchedule`] enum. The pose perturbation reuses
//! [`valenx_dock::search::mc::perturb`] — the same translation /
//! rotation / torsion kick the dock crate's own search applies.

use nalgebra::Vector3;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use valenx_dock::pose::Pose;
use valenx_dock::search::mc::perturb;

use crate::error::{DockScreenError, Result};
use crate::search::ga::coordinate_descent;
use crate::search::objective::PoseObjective;

/// The temperature schedule for a Monte-Carlo run.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum McSchedule {
    /// Constant-temperature Metropolis Monte Carlo — `kT` is fixed for
    /// the whole run.
    ConstantTemperature {
        /// The fixed `kT` value (energy units).
        kt: f64,
    },
    /// Geometric simulated annealing — `kT` starts at `kt_start` and is
    /// multiplied by `cooling` (< 1) each step until `kt_end`.
    SimulatedAnnealing {
        /// Starting `kT`.
        kt_start: f64,
        /// Final `kT`.
        kt_end: f64,
        /// Per-step geometric cooling factor (`0 < cooling < 1`).
        cooling: f64,
    },
}

impl McSchedule {
    /// A reasonable default annealing schedule.
    pub fn annealing() -> Self {
        McSchedule::SimulatedAnnealing {
            kt_start: 5.0,
            kt_end: 0.1,
            cooling: 0.95,
        }
    }

    /// The `kT` for step `i` of `total` steps.
    fn kt_at(&self, i: usize) -> f64 {
        match *self {
            McSchedule::ConstantTemperature { kt } => kt,
            McSchedule::SimulatedAnnealing {
                kt_start,
                kt_end,
                cooling,
            } => (kt_start * cooling.powi(i as i32)).max(kt_end),
        }
    }

    fn validate(&self) -> Result<()> {
        match *self {
            McSchedule::ConstantTemperature { kt } => {
                if !kt.is_finite() || kt <= 0.0 {
                    return Err(DockScreenError::invalid("kt", "kT must be positive"));
                }
            }
            McSchedule::SimulatedAnnealing {
                kt_start,
                kt_end,
                cooling,
            } => {
                if !kt_start.is_finite() || kt_start <= 0.0 {
                    return Err(DockScreenError::invalid("kt_start", "must be positive"));
                }
                if !kt_end.is_finite() || kt_end <= 0.0 || kt_end > kt_start {
                    return Err(DockScreenError::invalid(
                        "kt_end",
                        "must be positive and ≤ kt_start",
                    ));
                }
                if !(0.0..1.0).contains(&cooling) || cooling <= 0.0 {
                    return Err(DockScreenError::invalid(
                        "cooling",
                        "cooling factor must be in (0, 1)",
                    ));
                }
            }
        }
        Ok(())
    }
}

/// Tunable parameters for a Monte-Carlo / simulated-annealing search.
#[derive(Clone, Copy, Debug)]
pub struct McParams {
    /// Number of Metropolis steps.
    pub steps: usize,
    /// Maximum translation perturbation per step (Å).
    pub translation_step: f64,
    /// Maximum rotation perturbation per step (radians).
    pub rotation_step: f64,
    /// The temperature schedule.
    pub schedule: McSchedule,
    /// If `true`, run a local-search refinement on the best pose at
    /// the end of the walk.
    pub polish: bool,
}

impl Default for McParams {
    fn default() -> Self {
        McParams {
            steps: 250,
            translation_step: 2.0,
            rotation_step: std::f64::consts::FRAC_PI_6,
            schedule: McSchedule::annealing(),
            polish: true,
        }
    }
}

impl McParams {
    /// A small, fast configuration for previews and unit tests.
    pub fn fast() -> Self {
        McParams {
            steps: 60,
            translation_step: 2.0,
            rotation_step: std::f64::consts::FRAC_PI_6,
            schedule: McSchedule::ConstantTemperature { kt: 1.0 },
            polish: true,
        }
    }
}

/// The outcome of a Monte-Carlo / annealing run.
#[derive(Clone, Debug)]
pub struct McResult {
    /// The best pose visited.
    pub best_pose: Pose,
    /// Score of [`McResult::best_pose`].
    pub best_score: f64,
    /// Fraction `[0,1]` of proposed moves that were accepted.
    pub acceptance_ratio: f64,
}

/// Run a Monte-Carlo / simulated-annealing pose search from `start`.
///
/// Returns the best pose visited (refined by a local search if
/// `params.polish` is set), plus the move-acceptance ratio.
pub fn monte_carlo_search(
    objective: &PoseObjective,
    start: &Pose,
    params: &McParams,
    seed: u64,
) -> Result<McResult> {
    params.schedule.validate()?;
    if params.steps == 0 {
        return Err(DockScreenError::invalid("steps", "step count must be ≥ 1"));
    }
    let mut rng = StdRng::seed_from_u64(seed);
    let mut current = start.clone();
    let mut current_score = objective.score(&current);
    let mut best = current.clone();
    let mut best_score = current_score;
    let mut accepted = 0usize;

    for i in 0..params.steps {
        let kt = params.schedule.kt_at(i);
        let candidate = perturb(
            &current,
            &mut rng,
            params.translation_step,
            params.rotation_step,
        );
        let candidate_score = objective.score(&candidate);
        let delta = candidate_score - current_score;
        let accept = delta < 0.0 || rng.gen::<f64>() < (-delta / kt).exp();
        if accept {
            current = candidate;
            current_score = candidate_score;
            accepted += 1;
            if current_score < best_score {
                best = current.clone();
                best_score = current_score;
            }
        }
    }

    if params.polish {
        let (refined, refined_score) = coordinate_descent(objective, &best, 30);
        if refined_score < best_score {
            best = refined;
            best_score = refined_score;
        }
    }

    Ok(McResult {
        best_pose: best,
        best_score,
        acceptance_ratio: accepted as f64 / params.steps as f64,
    })
}

/// Convenience: run a Monte-Carlo search from a random start inside
/// the box `center` ± `size/2`.
pub fn monte_carlo_from_box(
    objective: &PoseObjective,
    center: Vector3<f64>,
    size: Vector3<f64>,
    params: &McParams,
    seed: u64,
) -> Result<McResult> {
    let mut rng = StdRng::seed_from_u64(seed);
    let half = size / 2.0;
    let mut start = Pose::identity(objective.n_torsions());
    start.translation = Vector3::new(
        center.x + rng.gen_range(-half.x..half.x),
        center.y + rng.gen_range(-half.y..half.y),
        center.z + rng.gen_range(-half.z..half.z),
    );
    monte_carlo_search(objective, &start, params, seed.wrapping_add(1))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prep::gridbox::GridBox;
    use crate::score::gridmap::{AffinityMapSet, MapKind};
    use valenx_dock::atom_type::Ad4AtomType;
    use valenx_dock::ligand::Ligand;
    use valenx_dock::receptor::{Receptor, ReceptorAtom};

    fn setup() -> (Ligand, AffinityMapSet) {
        let receptor = Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::zeros(),
                ad4_type: Ad4AtomType::C,
                partial_charge: 0.0,
            }],
        };
        let pdbqt = "ROOT
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C
ENDROOT
TORSDOF 0
";
        let lig = Ligand::from_pdbqt(pdbqt).unwrap();
        let grid = GridBox::with_spacing([0.0; 3], [16.0; 3], 0.375).unwrap();
        let maps =
            AffinityMapSet::precompute(&receptor, &[Ad4AtomType::C], &grid, MapKind::Vina).unwrap();
        (lig, maps)
    }

    #[test]
    fn schedule_cools_down_over_steps() {
        let s = McSchedule::SimulatedAnnealing {
            kt_start: 10.0,
            kt_end: 0.5,
            cooling: 0.9,
        };
        assert!(s.kt_at(0) > s.kt_at(20));
        // It never drops below kt_end.
        assert!(s.kt_at(1000) >= 0.5 - 1e-12);
    }

    #[test]
    fn constant_schedule_is_flat() {
        let s = McSchedule::ConstantTemperature { kt: 2.0 };
        assert_eq!(s.kt_at(0), 2.0);
        assert_eq!(s.kt_at(500), 2.0);
    }

    #[test]
    fn schedule_validation_catches_bad_values() {
        assert!(McSchedule::ConstantTemperature { kt: 0.0 }
            .validate()
            .is_err());
        assert!(McSchedule::SimulatedAnnealing {
            kt_start: 1.0,
            kt_end: 2.0, // end > start
            cooling: 0.9,
        }
        .validate()
        .is_err());
        assert!(McSchedule::SimulatedAnnealing {
            kt_start: 5.0,
            kt_end: 0.1,
            cooling: 1.5, // cooling >= 1
        }
        .validate()
        .is_err());
    }

    #[test]
    fn mc_rejects_zero_steps() {
        let (lig, maps) = setup();
        let obj = PoseObjective::uncharged(&lig, &maps);
        let params = McParams {
            steps: 0,
            ..McParams::fast()
        };
        assert!(monte_carlo_search(&obj, &Pose::identity(0), &params, 1).is_err());
    }

    #[test]
    fn mc_finds_an_attractive_pose_from_a_far_start() {
        let (lig, maps) = setup();
        let obj = PoseObjective::uncharged(&lig, &maps);
        let mut start = Pose::identity(0);
        start.translation = Vector3::new(6.5, 0.0, 0.0);
        let result = monte_carlo_search(&obj, &start, &McParams::fast(), 42).unwrap();
        assert!(
            result.best_score < -0.005,
            "MC best = {}",
            result.best_score
        );
        assert!((0.0..=1.0).contains(&result.acceptance_ratio));
    }

    #[test]
    fn annealing_run_finds_a_basin() {
        let (lig, maps) = setup();
        let obj = PoseObjective::uncharged(&lig, &maps);
        let params = McParams {
            steps: 120,
            schedule: McSchedule::annealing(),
            ..McParams::default()
        };
        let result = monte_carlo_from_box(
            &obj,
            Vector3::zeros(),
            Vector3::new(10.0, 10.0, 10.0),
            &params,
            3,
        )
        .unwrap();
        assert!(
            result.best_score < 0.0,
            "annealing best = {}",
            result.best_score
        );
    }

    #[test]
    fn mc_is_deterministic_for_a_fixed_seed() {
        let (lig, maps) = setup();
        let obj = PoseObjective::uncharged(&lig, &maps);
        let mut start = Pose::identity(0);
        start.translation = Vector3::new(5.0, 0.0, 0.0);
        let a = monte_carlo_search(&obj, &start, &McParams::fast(), 77).unwrap();
        let b = monte_carlo_search(&obj, &start, &McParams::fast(), 77).unwrap();
        assert!((a.best_score - b.best_score).abs() < 1e-12);
    }
}
