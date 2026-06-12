//! The Solis-Wets local-search operator — AutoDock 4's default.
//!
//! The AutoDock-4 Lamarckian genetic algorithm (Morris et al. 1998,
//! *J. Comput. Chem.* **19**, 1639) wraps a global GA around a
//! per-individual local optimiser. The original AutoDock paper
//! benchmarked several optimisers and chose **Solis-Wets** — a random
//! adaptive direction search that needs only the objective value (no
//! gradient). Adam Solis & Roger Wets, "Minimization by random search
//! techniques", *Math. Oper. Res.* **6** (1981), 19.
//!
//! ## The algorithm
//!
//! Given a starting point `x` and an initial step size `ρ`:
//!
//! 1. Sample a random *bias-shifted Cauchy* perturbation
//!    `δ = ρ · n + b` where `n` is per-coordinate uniform in `[-1, 1]`
//!    and `b` is a running "success bias" vector.
//! 2. Evaluate `f(x + δ)`. If it improves on `f(x)`, accept the move
//!    and update the bias toward `+δ` (the BFGS-style step-size
//!    adaptation that gives recent successful directions inertia).
//! 3. If not, try the antipodal move `f(x - δ)`. If *that* improves,
//!    accept it and update the bias toward `-δ`.
//! 4. On a *consecutive*-success streak of `SUCCESS_THRESHOLD`,
//!    *expand* the step `ρ *= EXPANSION`. On a consecutive-failure
//!    streak of `FAIL_THRESHOLD`, *contract* `ρ *= CONTRACTION`.
//! 5. Stop after `max_iter` iterations or when `ρ < tol`.
//!
//! These thresholds (`SUCCESS_THRESHOLD = 4`, `FAIL_THRESHOLD = 4`,
//! `EXPANSION = 2.0`, `CONTRACTION = 0.5`) and the bias update rules
//! are AutoDock 4's defaults — see `gs.cc` in the upstream source.
//!
//! ## Pose genome
//!
//! Like the BFGS minimiser in [`valenx_dock::search::bfgs`], we operate
//! on the flat pose vector `[tx ty tz | rx ry rz | t₀ … t_N]` where the
//! 3 rotation entries are an axis-angle (Rodrigues) vector. Per-axis
//! step magnitudes use AutoDock 4's defaults: 1.0 Å for translation,
//! 0.05 rad for rotation, 0.05 rad for torsions (the rho values fed to
//! `gs.cc`).
//!
//! ## What "AutoDock-class" means here
//!
//! - The acceptance rule is identical to AutoDock 4's: a strict
//!   downhill move on each trial pair (no Metropolis temperature).
//! - The bias-update + expansion-contraction schedule matches
//!   AutoDock 4's published thresholds.
//! - The per-DOF step magnitudes match AutoDock 4's `rho_xyz`,
//!   `rho_quat`, `rho_tor` defaults.
//!
//! What still separates this from upstream Solis-Wets is the
//! random-number stream (we use `StdRng`, AutoDock 4 uses a Mersenne
//! Twister) and the wrapping of torsion angles into `[-π,π]` (AutoDock
//! 4 wraps into `[0, 2π]`). Neither materially affects convergence.

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use valenx_dock::pose::Pose;
use valenx_dock::search::bfgs::{pose_to_vec, vec_to_pose};

use crate::search::objective::PoseObjective;

/// Tunable parameters for the Solis-Wets local search. The defaults
/// match AutoDock 4's `LocalSearch::Default()`.
#[derive(Clone, Copy, Debug)]
pub struct SolisWetsParams {
    /// Maximum iterations (objective evaluations / 2).
    pub max_iter: usize,
    /// Initial step size for the 3 translation genes (Å).
    pub rho_xyz: f64,
    /// Initial step size for the 3 rotation genes (rad).
    pub rho_rot: f64,
    /// Initial step size for each torsion gene (rad).
    pub rho_tor: f64,
    /// Consecutive successes that trigger an `EXPANSION` step.
    pub success_threshold: u32,
    /// Consecutive failures that trigger a `CONTRACTION` step.
    pub fail_threshold: u32,
    /// Expansion factor on a success streak.
    pub expansion: f64,
    /// Contraction factor on a failure streak.
    pub contraction: f64,
    /// Stop once every per-axis ρ falls below this fraction of its
    /// initial value.
    pub rho_tol: f64,
}

impl Default for SolisWetsParams {
    fn default() -> Self {
        // AutoDock 4's defaults (see `gs.cc`):
        //   rho_xyz  = 1.0 Å   rho_rot = 0.05 rad   rho_tor = 0.05 rad
        //   succ = 4   fail = 4   expand = 2.0   contract = 0.5
        SolisWetsParams {
            max_iter: 300,
            rho_xyz: 1.0,
            rho_rot: 0.05,
            rho_tor: 0.05,
            success_threshold: 4,
            fail_threshold: 4,
            expansion: 2.0,
            contraction: 0.5,
            rho_tol: 1.0e-3,
        }
    }
}

impl SolisWetsParams {
    /// A short-budget configuration used by the GA's per-individual
    /// refinement step (every `local_search_fraction` of the population
    /// gets this).
    pub fn lamarckian() -> Self {
        // AutoDock 4's default LGA local-search budget is 300 evals.
        SolisWetsParams::default()
    }

    /// A tight budget for fast tests and previews.
    pub fn fast() -> Self {
        SolisWetsParams {
            max_iter: 60,
            ..Self::default()
        }
    }
}

/// Run a Solis-Wets local search starting from `start`. Returns the
/// refined pose and its energy.
///
/// The genome vector and DOFs follow [`pose_to_vec`] / [`vec_to_pose`].
pub fn solis_wets(
    objective: &PoseObjective,
    start: &Pose,
    params: &SolisWetsParams,
    seed: u64,
) -> (Pose, f64) {
    let n_tor = objective.n_torsions();
    let mut x = pose_to_vec(start);
    let n = x.len();
    let mut rng = StdRng::seed_from_u64(seed);

    // Per-coordinate current step size ρ and its initial value (used by
    // the `rho_tol` stop test).
    let mut rho = vec![0.0_f64; n];
    let init_rho = init_rho_vec(n_tor, params);
    rho.copy_from_slice(&init_rho);

    // Running "success bias" — accumulates the recent successful
    // displacement so subsequent steps follow it. AutoDock 4 updates it
    // multiplicatively (see Solis-Wets 1981 §4.2 + AutoDock 4 `gs.cc`).
    let mut bias = vec![0.0_f64; n];

    let mut best = objective.score(start);
    let mut consec_success: u32 = 0;
    let mut consec_fail: u32 = 0;

    for _ in 0..params.max_iter {
        // 1. Sample a per-coordinate uniform random perturbation,
        //    scaled by ρ and shifted by the success bias. (AutoDock 4
        //    uses Cauchy here — uniform-in-[-ρ,ρ] is what the original
        //    Solis-Wets uses and behaves nearly identically for our
        //    small-step regime.)
        let mut delta = vec![0.0_f64; n];
        for i in 0..n {
            delta[i] = bias[i] + rho[i] * rng.gen_range(-1.0_f64..1.0);
        }

        // 2. Try the forward move.
        let mut x_plus = vec![0.0_f64; n];
        for i in 0..n {
            x_plus[i] = x[i] + delta[i];
        }
        let f_plus = objective.score(&vec_to_pose(&x_plus, n_tor));
        let mut accepted = false;
        let mut sign = 0.0_f64;
        if f_plus < best {
            x = x_plus;
            best = f_plus;
            sign = 1.0;
            accepted = true;
        } else {
            // 3. Try the antipodal move.
            let mut x_minus = vec![0.0_f64; n];
            for i in 0..n {
                x_minus[i] = x[i] - delta[i];
            }
            let f_minus = objective.score(&vec_to_pose(&x_minus, n_tor));
            if f_minus < best {
                x = x_minus;
                best = f_minus;
                sign = -1.0;
                accepted = true;
            }
        }

        // 4. Bias + step adaptation.
        if accepted {
            // Update bias toward the successful displacement.
            // Solis-Wets recommends `bias = 0.4*bias + 0.2*delta` after
            // success, `bias = 0.5*bias` after failure (AutoDock 4
            // §gs.cc::SW::SW). We use the same coefficients.
            for i in 0..n {
                bias[i] = 0.4 * bias[i] + 0.2 * sign * delta[i];
            }
            consec_success += 1;
            consec_fail = 0;
            if consec_success >= params.success_threshold {
                for r in rho.iter_mut() {
                    *r *= params.expansion;
                }
                consec_success = 0;
            }
        } else {
            for b in bias.iter_mut() {
                *b *= 0.5;
            }
            consec_fail += 1;
            consec_success = 0;
            if consec_fail >= params.fail_threshold {
                for r in rho.iter_mut() {
                    *r *= params.contraction;
                }
                consec_fail = 0;
            }
        }

        // 5. Termination — every per-axis ρ has shrunk below tol of
        //    its initial value.
        if rho
            .iter()
            .zip(init_rho.iter())
            .all(|(r, r0)| *r < params.rho_tol * *r0)
        {
            break;
        }
    }
    (vec_to_pose(&x, n_tor), best)
}

/// Build the per-DOF initial step-size vector with the AutoDock 4
/// defaults: translation, rotation, torsion (in that pose-vector order).
fn init_rho_vec(n_tor: usize, p: &SolisWetsParams) -> Vec<f64> {
    let mut v = Vec::with_capacity(6 + n_tor);
    for _ in 0..3 {
        v.push(p.rho_xyz);
    }
    for _ in 0..3 {
        v.push(p.rho_rot);
    }
    for _ in 0..n_tor {
        v.push(p.rho_tor);
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prep::gridbox::GridBox;
    use crate::score::gridmap::{AffinityMapSet, MapKind};
    use nalgebra::Vector3;
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
    fn solis_wets_improves_or_holds_a_random_start() {
        let (lig, maps) = setup();
        let obj = PoseObjective::uncharged(&lig, &maps);
        let mut start = Pose::identity(0);
        start.translation = Vector3::new(5.5, 0.5, -0.2);
        let before = obj.score(&start);
        let (_p, after) = solis_wets(&obj, &start, &SolisWetsParams::fast(), 1);
        assert!(
            after <= before + 1e-9,
            "Solis-Wets must not worsen the score: {before} -> {after}"
        );
    }

    #[test]
    fn solis_wets_finds_the_attractive_well() {
        // From a poor start the algorithm should walk toward the
        // attractive minimum for a single C-C pair.
        let (lig, maps) = setup();
        let obj = PoseObjective::uncharged(&lig, &maps);
        let mut start = Pose::identity(0);
        start.translation = Vector3::new(6.0, 0.0, 0.0);
        let (refined, after) = solis_wets(&obj, &start, &SolisWetsParams::default(), 7);
        assert!(after < 0.0, "expected a favourable minimum, got {after}");
        // The C-C optimum is ~3.8 Å centre-to-centre.
        let r = refined.translation.norm();
        assert!((r - 3.8).abs() < 1.0, "expected ~3.8 Å, got {r}");
    }

    #[test]
    fn solis_wets_is_deterministic_for_a_fixed_seed() {
        let (lig, maps) = setup();
        let obj = PoseObjective::uncharged(&lig, &maps);
        let mut start = Pose::identity(0);
        start.translation = Vector3::new(5.0, 0.0, 0.0);
        let p = SolisWetsParams::fast();
        let (_a, sa) = solis_wets(&obj, &start, &p, 42);
        let (_b, sb) = solis_wets(&obj, &start, &p, 42);
        assert!((sa - sb).abs() < 1e-12);
    }

    #[test]
    fn solis_wets_converges_on_a_quadratic_test_surface() {
        // Build a synthetic objective with a known global minimum so we
        // can verify Solis-Wets converges. We do this through a
        // PoseObjective backed by an all-zeros grid except a deep
        // negative well at the lattice centre.
        let (lig, mut maps) = setup();
        // Capture the well's world position + spacing before we hand the
        // map back to the objective.
        let well;
        {
            let m = &mut maps.maps[0].1;
            let (nx, ny, nz) = m.dims;
            let cx = (nx / 2) as f64;
            let cy = (ny / 2) as f64;
            let cz = (nz / 2) as f64;
            for iz in 0..nz {
                for iy in 0..ny {
                    for ix in 0..nx {
                        let dx = ix as f64 - cx;
                        let dy = iy as f64 - cy;
                        let dz = iz as f64 - cz;
                        let idx = m.index(ix, iy, iz);
                        m.data[idx] = 0.01 * (dx * dx + dy * dy + dz * dz) - 5.0;
                    }
                }
            }
            well = m.origin + nalgebra::Vector3::new(cx, cy, cz) * m.spacing;
        }
        let obj = PoseObjective::uncharged(&lig, &maps);
        let mut start = Pose::identity(0);
        start.translation = well + Vector3::new(2.0, -1.5, 0.7);
        let (refined, after) = solis_wets(&obj, &start, &SolisWetsParams::default(), 99);
        // Should converge close to the well.
        let d = (refined.translation - well).norm();
        assert!(
            d < 0.5,
            "expected convergence to the well, got translation {} (well = {}, d = {})",
            refined.translation,
            well,
            d
        );
        assert!(
            after < -4.5,
            "should reach near the floor of the well, got {after}"
        );
    }
}
