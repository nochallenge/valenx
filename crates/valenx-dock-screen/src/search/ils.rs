//! Feature 11 — Vina-class iterated local search.
//!
//! AutoDock Vina's search is *iterated local search* (ILS), also
//! called basin hopping:
//!
//! 1. start from a pose, run a local optimisation (Vina uses BFGS);
//! 2. perturb the optimised pose;
//! 3. run the local optimisation again;
//! 4. accept the new local minimum by the Metropolis criterion;
//! 5. repeat.
//!
//! The local-optimisation step is what turns a coarse random walk into
//! an efficient search — every accepted point is already at the bottom
//! of a basin.
//!
//! For the Vina-class scoring function this module **reuses
//! [`valenx_dock`]'s ILS verbatim** — `valenx_dock`'s
//! [`valenx_dock::search::mc::iterated_local_search`]
//! and its BFGS minimiser are exactly the Vina-class search, so the
//! Vina path here just builds the dock crate's grid bundle and calls
//! it. For an arbitrary [`crate::search::objective::PoseObjective`]
//! (e.g. an AutoDock4 grid) ILS falls back to a self-contained
//! perturb + coordinate-descent + Metropolis loop.

use nalgebra::Vector3;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use valenx_dock::grid::GridBundle;
use valenx_dock::ligand::Ligand;
use valenx_dock::pose::Pose;
use valenx_dock::receptor::Receptor;
use valenx_dock::search::mc::{iterated_local_search, perturb};

use crate::error::{DockScreenError, Result};
use crate::prep::gridbox::GridBox;
use crate::search::ga::coordinate_descent;
use crate::search::objective::PoseObjective;

/// Tunable parameters for an iterated-local-search run.
#[derive(Clone, Copy, Debug)]
pub struct IlsParams {
    /// Number of outer ILS iterations (perturb + local-optimise).
    pub iterations: usize,
    /// Metropolis temperature `kT` for accepting a new local minimum.
    pub kt: f64,
    /// Maximum translation perturbation between basins (Å).
    pub translation_step: f64,
    /// Maximum rotation perturbation between basins (radians).
    pub rotation_step: f64,
}

impl Default for IlsParams {
    fn default() -> Self {
        IlsParams {
            iterations: 50,
            kt: 1.0,
            translation_step: 2.0,
            rotation_step: std::f64::consts::FRAC_PI_6,
        }
    }
}

impl IlsParams {
    /// A small, fast configuration for previews and unit tests.
    pub fn fast() -> Self {
        IlsParams {
            iterations: 12,
            kt: 1.0,
            translation_step: 2.0,
            rotation_step: std::f64::consts::FRAC_PI_6,
        }
    }
}

/// The outcome of an iterated-local-search run.
#[derive(Clone, Debug)]
pub struct IlsResult {
    /// The best pose found.
    pub best_pose: Pose,
    /// Score of [`IlsResult::best_pose`].
    pub best_score: f64,
}

/// Vina-class ILS using `valenx_dock`'s own BFGS-backed iterated local
/// search. This is the exact AutoDock Vina search — it builds the dock
/// crate's grid bundle for the receptor and ligand and runs its
/// `iterated_local_search`.
///
/// Use this when the scoring function is Vina-class; for other
/// objectives use [`iterated_local_search_objective`].
pub fn iterated_local_search_vina(
    receptor: &Receptor,
    ligand: &Ligand,
    grid: &GridBox,
    start: &Pose,
    params: &IlsParams,
    seed: u64,
) -> Result<IlsResult> {
    if params.iterations == 0 {
        return Err(DockScreenError::invalid(
            "iterations",
            "ILS iteration count must be ≥ 1",
        ));
    }
    // Build the dock crate's own grid bundle so the search is
    // bit-identical to a native valenx-dock run.
    let bundle = GridBundle::build(
        receptor,
        ligand,
        grid.origin(),
        grid.spacing,
        grid.dims(),
    );
    let (best_pose, best_score) =
        iterated_local_search(ligand, start, &bundle, params.iterations, seed);
    Ok(IlsResult {
        best_pose,
        best_score,
    })
}

/// Generic ILS over an arbitrary [`PoseObjective`] — used when the
/// scoring function is *not* Vina-class (for instance an AutoDock4
/// affinity-map objective). Perturb, run a coordinate-descent local
/// optimisation, accept by Metropolis.
pub fn iterated_local_search_objective(
    objective: &PoseObjective,
    start: &Pose,
    params: &IlsParams,
    seed: u64,
) -> Result<IlsResult> {
    if params.iterations == 0 {
        return Err(DockScreenError::invalid(
            "iterations",
            "ILS iteration count must be ≥ 1",
        ));
    }
    if !params.kt.is_finite() || params.kt <= 0.0 {
        return Err(DockScreenError::invalid("kt", "kT must be positive"));
    }
    let mut rng = StdRng::seed_from_u64(seed);
    let (mut current, mut current_score) = coordinate_descent(objective, start, 30);
    let mut best = current.clone();
    let mut best_score = current_score;

    for _ in 0..params.iterations {
        let candidate = perturb(
            &current,
            &mut rng,
            params.translation_step,
            params.rotation_step,
        );
        let (refined, refined_score) = coordinate_descent(objective, &candidate, 30);
        let delta = refined_score - current_score;
        let accept = delta < 0.0 || rng.gen::<f64>() < (-delta / params.kt).exp();
        if accept {
            current = refined;
            current_score = refined_score;
            if current_score < best_score {
                best = current.clone();
                best_score = current_score;
            }
        }
    }
    Ok(IlsResult {
        best_pose: best,
        best_score,
    })
}

/// Convenience: Vina-class ILS from a random start inside the box.
pub fn ils_vina_from_box(
    receptor: &Receptor,
    ligand: &Ligand,
    grid: &GridBox,
    params: &IlsParams,
    seed: u64,
) -> Result<IlsResult> {
    let mut rng = StdRng::seed_from_u64(seed);
    let half = grid.size / 2.0;
    let mut start = Pose::identity(ligand.n_torsions());
    start.translation = Vector3::new(
        grid.center.x + rng.gen_range(-half.x..half.x),
        grid.center.y + rng.gen_range(-half.y..half.y),
        grid.center.z + rng.gen_range(-half.z..half.z),
    );
    iterated_local_search_vina(receptor, ligand, grid, &start, params, seed.wrapping_add(1))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::score::gridmap::{AffinityMapSet, MapKind};
    use valenx_dock::atom_type::Ad4AtomType;
    use valenx_dock::receptor::ReceptorAtom;

    fn carbon_receptor() -> Receptor {
        Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::zeros(),
                ad4_type: Ad4AtomType::C,
                partial_charge: 0.0,
            }],
        }
    }

    fn one_carbon_ligand() -> Ligand {
        let pdbqt = "ROOT
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C
ENDROOT
TORSDOF 0
";
        Ligand::from_pdbqt(pdbqt).unwrap()
    }

    #[test]
    fn vina_ils_rejects_zero_iterations() {
        let r = carbon_receptor();
        let lig = one_carbon_ligand();
        let grid = GridBox::with_spacing([0.0; 3], [12.0; 3], 0.5).unwrap();
        let params = IlsParams {
            iterations: 0,
            ..IlsParams::fast()
        };
        assert!(
            iterated_local_search_vina(&r, &lig, &grid, &Pose::identity(0), &params, 1).is_err()
        );
    }

    #[test]
    fn vina_ils_finds_attractive_well() {
        let r = carbon_receptor();
        let lig = one_carbon_ligand();
        let grid = GridBox::with_spacing([0.0; 3], [12.0; 3], 0.5).unwrap();
        let mut start = Pose::identity(0);
        start.translation = Vector3::new(5.0, 0.0, 0.0);
        let result =
            iterated_local_search_vina(&r, &lig, &grid, &start, &IlsParams::fast(), 42).unwrap();
        assert!(
            result.best_score < -0.005,
            "Vina ILS best = {}",
            result.best_score
        );
    }

    #[test]
    fn objective_ils_improves_a_far_start() {
        let r = carbon_receptor();
        let lig = one_carbon_ligand();
        let grid = GridBox::with_spacing([0.0; 3], [14.0; 3], 0.375).unwrap();
        let maps =
            AffinityMapSet::precompute(&r, &[Ad4AtomType::C], &grid, MapKind::Vina).unwrap();
        let obj = PoseObjective::uncharged(&lig, &maps);
        let mut start = Pose::identity(0);
        start.translation = Vector3::new(6.0, 0.0, 0.0);
        let before = obj.score(&start);
        let result =
            iterated_local_search_objective(&obj, &start, &IlsParams::fast(), 7).unwrap();
        assert!(result.best_score <= before, "ILS must not worsen the score");
    }

    #[test]
    fn objective_ils_rejects_bad_kt() {
        let r = carbon_receptor();
        let lig = one_carbon_ligand();
        let grid = GridBox::with_spacing([0.0; 3], [14.0; 3], 0.5).unwrap();
        let maps =
            AffinityMapSet::precompute(&r, &[Ad4AtomType::C], &grid, MapKind::Vina).unwrap();
        let obj = PoseObjective::uncharged(&lig, &maps);
        let params = IlsParams {
            kt: 0.0,
            ..IlsParams::fast()
        };
        assert!(iterated_local_search_objective(&obj, &Pose::identity(0), &params, 1).is_err());
    }

    #[test]
    fn ils_from_box_runs() {
        let r = carbon_receptor();
        let lig = one_carbon_ligand();
        let grid = GridBox::with_spacing([0.0; 3], [12.0; 3], 0.5).unwrap();
        let result = ils_vina_from_box(&r, &lig, &grid, &IlsParams::fast(), 5).unwrap();
        assert!(result.best_score.is_finite());
    }
}
