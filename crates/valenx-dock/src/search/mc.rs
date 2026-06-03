//! Monte Carlo / iterated local search over poses.

use nalgebra::{UnitQuaternion, Vector3};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use crate::grid::GridBundle;
use crate::ligand::Ligand;
use crate::pose::Pose;
use crate::search::bfgs::minimize_bfgs;

/// Apply a random perturbation: small translation, small rotation,
/// each torsion +/- up to 30 deg.
pub fn perturb(pose: &Pose, rng: &mut StdRng, tr_step: f64, rot_step: f64) -> Pose {
    let trans = Vector3::new(
        rng.gen_range(-tr_step..tr_step),
        rng.gen_range(-tr_step..tr_step),
        rng.gen_range(-tr_step..tr_step),
    );
    let axis = nalgebra::Unit::new_normalize(Vector3::new(
        rng.gen_range(-1.0..1.0),
        rng.gen_range(-1.0..1.0),
        rng.gen_range(-1.0..1.0),
    ));
    let angle = rng.gen_range(-rot_step..rot_step);
    let drot = UnitQuaternion::from_axis_angle(&axis, angle);
    let mut torsions = pose.torsions.clone();
    for t in torsions.iter_mut() {
        *t += rng.gen_range(-0.5..0.5);
    }
    Pose {
        translation: pose.translation + trans,
        orientation: drot * pose.orientation,
        torsions,
    }
}

/// Run iterated local search for `n_iter` outer iterations, starting
/// from `start`. Returns the best pose ever seen.
pub fn iterated_local_search(
    ligand: &Ligand,
    start: &Pose,
    grids: &GridBundle,
    n_iter: usize,
    seed: u64,
) -> (Pose, f64) {
    let mut rng = StdRng::seed_from_u64(seed);
    let (mut current, mut current_score) = minimize_bfgs(ligand, start, grids, 50, 1e-3);
    let mut best = current.clone();
    let mut best_score = current_score;
    let kt = 1.0;
    for _ in 0..n_iter {
        let candidate = perturb(&current, &mut rng, 2.0, std::f64::consts::PI / 6.0);
        let (refined, score) = minimize_bfgs(ligand, &candidate, grids, 50, 1e-3);
        let dscore = score - current_score;
        let accept = dscore < 0.0 || rng.gen::<f64>() < (-dscore / kt).exp();
        if accept {
            current = refined;
            current_score = score;
            if current_score < best_score {
                best = current.clone();
                best_score = current_score;
            }
        }
    }
    (best, best_score)
}

use rayon::prelude::*;

/// Run `n_chains` ILS chains in parallel from random starts within
/// the search box (`center` +/- `size/2`). Returns all `(Pose, score)`
/// sorted ascending by score.
pub fn exhaustiveness_search(
    ligand: &Ligand,
    grids: &GridBundle,
    center: Vector3<f64>,
    size: Vector3<f64>,
    n_chains: usize,
    n_iter_per_chain: usize,
    seed: u64,
) -> Vec<(Pose, f64)> {
    let half = size / 2.0;
    let mut results: Vec<(Pose, f64)> = (0..n_chains)
        .into_par_iter()
        .map(|i| {
            // splitmix64 finalizer. The previous
            // `seed.wrapping_add(i).wrapping_mul(2654435761)` mixer
            // had a degenerate seed=0,i=0 → chain_seed=0 case and
            // poorly decorrelated nearby (i, i+1) seeds. splitmix64
            // is the canonical fix and is what `StdRng`'s own
            // `next_u64` jump table uses internally.
            let mut x = seed.wrapping_add((i as u64).wrapping_mul(0x9E3779B97F4A7C15));
            x = (x ^ (x >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
            let chain_seed = x ^ (x >> 31);
            let mut rng = StdRng::seed_from_u64(chain_seed);
            let mut start = Pose::identity(ligand.n_torsions());
            start.translation = Vector3::new(
                center.x + rng.gen_range(-half.x..half.x),
                center.y + rng.gen_range(-half.y..half.y),
                center.z + rng.gen_range(-half.z..half.z),
            );
            // Random initial orientation: pick an axis on the unit
            // sphere + uniform angle.
            let axis = nalgebra::Unit::new_normalize(Vector3::new(
                rng.gen_range(-1.0..1.0),
                rng.gen_range(-1.0..1.0),
                rng.gen_range(-1.0..1.0),
            ));
            let angle = rng.gen_range(-std::f64::consts::PI..std::f64::consts::PI);
            start.orientation = UnitQuaternion::from_axis_angle(&axis, angle);
            // Random initial torsions.
            for t in start.torsions.iter_mut() {
                *t = rng.gen_range(-std::f64::consts::PI..std::f64::consts::PI);
            }
            iterated_local_search(ligand, &start, grids, n_iter_per_chain, chain_seed)
        })
        .collect();
    results.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom_type::Ad4AtomType;
    use crate::pose::Pose;
    use crate::receptor::{Receptor, ReceptorAtom};

    #[test]
    fn ils_finds_attractive_well_from_far_start() {
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
        let grids = GridBundle::build(
            &receptor,
            &lig,
            Vector3::new(-5.0, -5.0, -5.0),
            0.5,
            (21, 21, 21),
        );
        let mut start = Pose::identity(0);
        start.translation = Vector3::new(4.5, 0.0, 0.0);
        let (_best, best_score) = iterated_local_search(&lig, &start, &grids, 20, 42);
        // Any reasonable pose in this attractive well should score
        // well below zero given a single C-C pair.
        assert!(best_score < -0.005, "best_score = {best_score}");
    }
    #[test]
    fn parallel_exhaustiveness_returns_best_of_all_chains() {
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
        let grids = GridBundle::build(
            &receptor,
            &lig,
            Vector3::new(-5.0, -5.0, -5.0),
            0.5,
            (21, 21, 21),
        );
        let center = Vector3::zeros();
        let size = Vector3::new(8.0, 8.0, 8.0);
        let results = exhaustiveness_search(&lig, &grids, center, size, 4, 10, 42);
        assert_eq!(results.len(), 4);
        // Sorted ascending by score — the best (lowest) is first.
        for w in results.windows(2) {
            assert!(w[0].1 <= w[1].1, "results not sorted: {results:?}");
        }
    }
}
