//! Features 12 & 13 — rigid and flexible docking drivers.
//!
//! The drivers tie the preparation, scoring and search layers
//! together: they take a parsed receptor + ligand and a search box,
//! build affinity maps, run a chosen search, and return ranked poses.
//!
//! - [`rigid_dock`] (feature 12) — the standard case: a rigid
//!   receptor, a flexible ligand. Builds the affinity-map set once and
//!   runs the search.
//! - [`flexible_dock`] (feature 13) — flexible-receptor docking. A set
//!   of receptor side chains (chosen by [`crate::prep::flex`]) is
//!   allowed to move. v1 models this by *splitting the flexible
//!   side-chain atoms out of the rigid receptor* — they no longer
//!   contribute to the static affinity grid — and re-scoring the
//!   moved side chains explicitly against the ligand for each of a
//!   small set of side-chain conformations, keeping the best. This
//!   captures the dominant induced-fit effect (a clashing side chain
//!   rotating away) without a full receptor torsion tree.
//!
//! Both drivers are generic over [`SearchAlgorithm`].

use nalgebra::Vector3;

use valenx_dock::atom_type::Ad4AtomType;
use valenx_dock::ligand::Ligand;
use valenx_dock::pose::Pose;
use valenx_dock::receptor::{Receptor, ReceptorAtom};

use crate::error::{DockScreenError, Result};
use crate::prep::flex::FlexSelection;
use crate::prep::gridbox::GridBox;
use crate::score::gridmap::{AffinityMapSet, MapKind};
use crate::score::vina::score_complex as vina_score_complex;
use crate::search::flex_pose::{induced_fit_dock, ChiRotation, InducedFitResult};
use crate::search::ga::{GaParams, LamarckianGa};
use crate::search::ils::{iterated_local_search_objective, IlsParams};
use crate::search::mc::{monte_carlo_from_box, McParams};
use crate::search::objective::PoseObjective;

/// Which search algorithm a docking driver runs.
#[derive(Clone, Copy, Debug)]
pub enum SearchAlgorithm {
    /// The Lamarckian genetic algorithm ([`crate::search::ga`]).
    Genetic(GaParams),
    /// Monte-Carlo / simulated annealing ([`crate::search::mc`]).
    MonteCarlo(McParams),
    /// Iterated local search ([`crate::search::ils`]).
    IteratedLocal(IlsParams),
}

impl SearchAlgorithm {
    /// A small, fast genetic-algorithm configuration — the default for
    /// previews and tests.
    pub fn fast() -> Self {
        SearchAlgorithm::Genetic(GaParams::fast())
    }

    /// A short label naming the algorithm.
    pub fn label(&self) -> &'static str {
        match self {
            SearchAlgorithm::Genetic(_) => "lamarckian-ga",
            SearchAlgorithm::MonteCarlo(_) => "monte-carlo",
            SearchAlgorithm::IteratedLocal(_) => "iterated-local-search",
        }
    }
}

/// A scored pose returned by a docking driver.
#[derive(Clone, Debug)]
pub struct ScoredPose {
    /// The ligand pose.
    pub pose: Pose,
    /// The pose's score in kcal/mol (lower is better).
    pub score: f64,
}

/// The result of a single docking driver run: the ranked poses plus
/// the affinity-map set that produced them (kept so callers can
/// rescore or cluster without rebuilding the grids).
#[derive(Clone, Debug)]
pub struct DockRun {
    /// Poses sorted best (lowest score) first.
    pub poses: Vec<ScoredPose>,
    /// The receptor affinity maps used.
    pub maps: AffinityMapSet,
}

impl DockRun {
    /// The best (lowest-scoring) pose, if any.
    pub fn best(&self) -> Option<&ScoredPose> {
        self.poses.first()
    }
}

/// Feature 12 — dock a flexible ligand against a rigid receptor.
///
/// `n_runs` independent searches are launched from different seeds;
/// every search's best pose is collected and the pool is sorted by
/// score. The caller typically clusters the result (see
/// [`crate::screen::cluster`]).
///
/// Returns [`DockScreenError`] if the inputs are degenerate or a
/// search parameter is out of range.
pub fn rigid_dock(
    receptor: &Receptor,
    ligand: &Ligand,
    grid: &GridBox,
    algorithm: SearchAlgorithm,
    n_runs: usize,
    seed: u64,
) -> Result<DockRun> {
    if receptor.atoms.is_empty() {
        return Err(DockScreenError::invalid_receptor("receptor has no atoms"));
    }
    if ligand.atoms.is_empty() {
        return Err(DockScreenError::invalid_ligand("ligand has no atoms"));
    }
    if n_runs == 0 {
        return Err(DockScreenError::invalid("n_runs", "run count must be ≥ 1"));
    }
    let ligand_types: Vec<Ad4AtomType> = ligand.atoms.iter().map(|a| a.ad4_type).collect();
    let maps = AffinityMapSet::precompute(receptor, &ligand_types, grid, MapKind::Vina)?;
    let charges: Vec<f64> = ligand.atoms.iter().map(|a| a.partial_charge).collect();
    let objective = PoseObjective::new(ligand, &maps, charges);

    let mut poses: Vec<ScoredPose> = Vec::with_capacity(n_runs);
    for run in 0..n_runs {
        let run_seed = mix_seed(seed, run as u64);
        let (pose, score) = run_one(&objective, receptor, ligand, grid, &algorithm, run_seed)?;
        poses.push(ScoredPose { pose, score });
    }
    poses.sort_by(|a, b| {
        a.score
            .partial_cmp(&b.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Ok(DockRun { poses, maps })
}

/// The flexible-receptor side of a [`flexible_dock`] call: which
/// side-chain residues flex, which receptor atoms belong to them, and
/// the alternate side-chain conformations to try.
#[derive(Clone, Debug, Default)]
pub struct FlexibleSpec {
    /// The selected flexible side chains (from
    /// [`crate::prep::flex::select_flexible`]).
    pub selection: FlexSelection,
    /// Receptor-atom indices belonging to the flexible side chains.
    pub atom_indices: Vec<usize>,
    /// Alternate side-chain conformations to try. Each is the world
    /// positions of the `atom_indices` atoms, in the same order.
    pub conformations: Vec<Vec<Vector3<f64>>>,
}

impl FlexibleSpec {
    /// Build a flexible-receptor spec.
    pub fn new(
        selection: FlexSelection,
        atom_indices: Vec<usize>,
        conformations: Vec<Vec<Vector3<f64>>>,
    ) -> Self {
        FlexibleSpec {
            selection,
            atom_indices,
            conformations,
        }
    }
}

/// Feature 13 — dock a flexible ligand against a receptor with
/// flexible side chains.
///
/// The `flex` spec (built from [`crate::prep::flex::select_flexible`])
/// names which receptor side chains flex, the receptor-atom indices
/// belonging to them, and the alternate side-chain conformations to
/// try.
///
/// v1 strategy: the flexible side-chain atoms are removed from the
/// rigid receptor before the affinity grid is built, so the grid
/// reflects only the fixed part of the receptor. After the ligand
/// search converges, each candidate side-chain conformation is scored
/// explicitly against the best ligand pose with the Vina-class
/// function; the conformation that minimises the total (grid +
/// side-chain) score is chosen. This captures the dominant induced-fit
/// move — a clashing side chain swinging out of the way — without a
/// full receptor torsion tree.
pub fn flexible_dock(
    receptor: &Receptor,
    ligand: &Ligand,
    grid: &GridBox,
    flex: &FlexibleSpec,
    algorithm: SearchAlgorithm,
    n_runs: usize,
    seed: u64,
) -> Result<FlexibleDockRun> {
    let flexible = &flex.selection;
    let flex_atom_indices = flex.atom_indices.as_slice();
    let sidechain_conformations = flex.conformations.as_slice();
    if flexible.is_empty() {
        // Nothing flexible — fall straight back to rigid docking.
        let run = rigid_dock(receptor, ligand, grid, algorithm, n_runs, seed)?;
        return Ok(FlexibleDockRun {
            poses: run.poses,
            chosen_conformation: 0,
            flexible_residues: 0,
        });
    }
    // Validate the flexible-atom indices.
    for &idx in flex_atom_indices {
        if idx >= receptor.atoms.len() {
            return Err(DockScreenError::invalid_receptor(format!(
                "flexible-sidechain atom index {idx} out of range ({} receptor atoms)",
                receptor.atoms.len()
            )));
        }
    }
    // 1. Build a rigid receptor with the flexible side-chain atoms
    //    removed — the static grid sees only the fixed core.
    let flex_set: std::collections::BTreeSet<usize> =
        flex_atom_indices.iter().copied().collect();
    let rigid_core = Receptor {
        atoms: receptor
            .atoms
            .iter()
            .enumerate()
            .filter(|(i, _)| !flex_set.contains(i))
            .map(|(_, a)| a.clone())
            .collect(),
    };
    if rigid_core.atoms.is_empty() {
        return Err(DockScreenError::invalid_receptor(
            "flexible selection removed every receptor atom",
        ));
    }
    // 2. Dock the ligand against the rigid core.
    let core_run = rigid_dock(&rigid_core, ligand, grid, algorithm, n_runs, seed)?;
    // 3. For the best ligand pose, pick the side-chain conformation
    //    that minimises the explicit ligand / side-chain energy.
    let flex_atom_types: Vec<Ad4AtomType> = flex_atom_indices
        .iter()
        .map(|&i| receptor.atoms[i].ad4_type)
        .collect();
    let mut best_conf = 0usize;
    let mut adjusted_poses = core_run.poses.clone();
    if let Some(top) = core_run.poses.first() {
        let ligand_world = ligand.apply_pose(&top.pose);
        let ligand_for_sidechain: Vec<(Vector3<f64>, Ad4AtomType)> = ligand_world
            .iter()
            .zip(ligand.atoms.iter())
            .map(|(p, a)| (*p, a.ad4_type))
            .collect();
        // Conformation 0 is the input side-chain geometry.
        let mut confs: Vec<Vec<Vector3<f64>>> = Vec::new();
        confs.push(
            flex_atom_indices
                .iter()
                .map(|&i| receptor.atoms[i].position)
                .collect(),
        );
        confs.extend(sidechain_conformations.iter().cloned());

        let mut best_total = f64::INFINITY;
        for (ci, conf) in confs.iter().enumerate() {
            if conf.len() != flex_atom_indices.len() {
                // Skip a malformed conformation rather than failing
                // the whole run.
                continue;
            }
            // Score the moved side-chain atoms as a "mini receptor"
            // against the ligand.
            let sidechain_receptor = Receptor {
                atoms: conf
                    .iter()
                    .zip(flex_atom_types.iter())
                    .map(|(p, t)| ReceptorAtom {
                        position: *p,
                        ad4_type: *t,
                        partial_charge: 0.0,
                    })
                    .collect(),
            };
            let sidechain_energy =
                vina_score_complex(&sidechain_receptor, &ligand_for_sidechain, 0)
                    .intermolecular();
            let total = top.score + sidechain_energy;
            if total < best_total {
                best_total = total;
                best_conf = ci;
            }
        }
        // Re-rank: every pose's score gains the chosen side-chain term
        // for the *top* pose's geometry (a v1 approximation — see the
        // module note; a full treatment would re-evaluate per pose).
        if best_total.is_finite() {
            let delta = best_total - top.score;
            for p in adjusted_poses.iter_mut() {
                p.score += delta;
            }
        }
    }
    adjusted_poses.sort_by(|a, b| {
        a.score
            .partial_cmp(&b.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Ok(FlexibleDockRun {
        poses: adjusted_poses,
        chosen_conformation: best_conf,
        flexible_residues: flexible.len(),
    })
}

/// The result of a flexible-receptor docking run.
#[derive(Clone, Debug)]
pub struct FlexibleDockRun {
    /// Poses sorted best-first, scores adjusted for the chosen
    /// side-chain conformation.
    pub poses: Vec<ScoredPose>,
    /// Index into the conformation list of the side-chain
    /// conformation that was selected (0 = the input geometry).
    pub chosen_conformation: usize,
    /// Number of receptor side chains treated as flexible.
    pub flexible_residues: usize,
}

/// True induced-fit docking driver: dock a flexible ligand against a
/// receptor with flexible side-chain χ angles co-optimised *during* the
/// search.
///
/// Unlike [`flexible_dock`] (the post-search rescoring approach), this
/// runs the GA / local-search jointly over `(ligand_pose ∪ χ_angles)`,
/// so the receptor truly relaxes around the ligand during the search.
/// Optionally rescores the final complex with the MM-GBSA-class
/// generalized-Born model.
///
/// `chi_rotations` lists which side-chain χ angles the search treats as
/// free variables. Use [`crate::search::flex_pose::chi_from_axis_atoms`]
/// to build them from receptor atom indices.
pub fn induced_fit_dock_driver(
    receptor: &Receptor,
    ligand: &Ligand,
    grid: &GridBox,
    chi_rotations: Vec<ChiRotation>,
    n_restarts: usize,
    seed: u64,
    rescore_with_mmgbsa: bool,
) -> Result<InducedFitResult> {
    induced_fit_dock(
        receptor,
        ligand,
        grid,
        chi_rotations,
        n_restarts,
        seed,
        rescore_with_mmgbsa,
    )
}

/// Run one search with the chosen algorithm and return `(pose, score)`.
fn run_one(
    objective: &PoseObjective,
    receptor: &Receptor,
    ligand: &Ligand,
    grid: &GridBox,
    algorithm: &SearchAlgorithm,
    seed: u64,
) -> Result<(Pose, f64)> {
    match algorithm {
        SearchAlgorithm::Genetic(p) => {
            let res = LamarckianGa::new(*p).run(objective, grid.center, grid.size, seed)?;
            Ok((res.best_pose, res.best_score))
        }
        SearchAlgorithm::MonteCarlo(p) => {
            let res = monte_carlo_from_box(objective, grid.center, grid.size, p, seed)?;
            Ok((res.best_pose, res.best_score))
        }
        SearchAlgorithm::IteratedLocal(p) => {
            // ILS from a random box start, using the generic objective
            // path (works for any scoring function / grid).
            let _ = (receptor, ligand); // signature symmetry with other arms
            let start = random_start(objective, grid, seed);
            let res = iterated_local_search_objective(objective, &start, p, seed)?;
            Ok((res.best_pose, res.best_score))
        }
    }
}

/// A random starting pose inside the box.
fn random_start(objective: &PoseObjective, grid: &GridBox, seed: u64) -> Pose {
    use rand::{Rng, SeedableRng};
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    let half = grid.size / 2.0;
    let mut p = Pose::identity(objective.n_torsions());
    p.translation = Vector3::new(
        grid.center.x + rng.gen_range(-half.x..half.x),
        grid.center.y + rng.gen_range(-half.y..half.y),
        grid.center.z + rng.gen_range(-half.z..half.z),
    );
    p
}

/// splitmix64-style seed mixer so per-run seeds are well decorrelated.
fn mix_seed(seed: u64, run: u64) -> u64 {
    let mut x = seed.wrapping_add(run.wrapping_mul(0x9E37_79B9_7F4A_7C15));
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prep::flex::FlexibleSidechain;

    fn carbon_receptor() -> Receptor {
        Receptor {
            atoms: vec![
                ReceptorAtom {
                    position: Vector3::zeros(),
                    ad4_type: Ad4AtomType::C,
                    partial_charge: 0.0,
                },
                ReceptorAtom {
                    position: Vector3::new(2.0, 0.0, 0.0),
                    ad4_type: Ad4AtomType::C,
                    partial_charge: 0.0,
                },
            ],
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
    fn rigid_dock_rejects_degenerate_inputs() {
        let r = carbon_receptor();
        let lig = one_carbon_ligand();
        let grid = GridBox::cubic([1.0, 0.0, 0.0], 10.0).unwrap();
        // Empty receptor.
        assert!(rigid_dock(&Receptor::default(), &lig, &grid, SearchAlgorithm::fast(), 1, 1).is_err());
        // Zero runs.
        assert!(rigid_dock(&r, &lig, &grid, SearchAlgorithm::fast(), 0, 1).is_err());
    }

    #[test]
    fn rigid_dock_returns_sorted_poses() {
        let r = carbon_receptor();
        let lig = one_carbon_ligand();
        let grid = GridBox::with_spacing([1.0, 0.0, 0.0], [12.0; 3], 0.5).unwrap();
        let run = rigid_dock(&r, &lig, &grid, SearchAlgorithm::fast(), 4, 42).unwrap();
        assert_eq!(run.poses.len(), 4);
        for w in run.poses.windows(2) {
            assert!(w[0].score <= w[1].score, "poses not sorted");
        }
        assert!(run.best().is_some());
    }

    #[test]
    fn rigid_dock_works_with_each_algorithm() {
        let r = carbon_receptor();
        let lig = one_carbon_ligand();
        let grid = GridBox::with_spacing([1.0, 0.0, 0.0], [12.0; 3], 0.5).unwrap();
        for algo in [
            SearchAlgorithm::Genetic(GaParams::fast()),
            SearchAlgorithm::MonteCarlo(McParams::fast()),
            SearchAlgorithm::IteratedLocal(IlsParams::fast()),
        ] {
            let run = rigid_dock(&r, &lig, &grid, algo, 2, 7).unwrap();
            assert!(!run.poses.is_empty(), "{} produced no poses", algo.label());
            assert!(run.best().unwrap().score.is_finite());
        }
    }

    fn ser_selection() -> FlexSelection {
        FlexSelection {
            residues: vec![FlexibleSidechain {
                chain_id: "A".into(),
                residue_name: "SER".into(),
                seq_num: 1,
                chi_count: 1,
                min_distance: 2.0,
            }],
            cutoff: 5.0,
        }
    }

    #[test]
    fn flexible_dock_with_empty_selection_falls_back_to_rigid() {
        let r = carbon_receptor();
        let lig = one_carbon_ligand();
        let grid = GridBox::with_spacing([1.0, 0.0, 0.0], [12.0; 3], 0.5).unwrap();
        let spec = FlexibleSpec::default();
        let run = flexible_dock(&r, &lig, &grid, &spec, SearchAlgorithm::fast(), 2, 1).unwrap();
        assert!(!run.poses.is_empty());
        assert_eq!(run.flexible_residues, 0);
    }

    #[test]
    fn flexible_dock_picks_a_sidechain_conformation() {
        let r = carbon_receptor();
        let lig = one_carbon_ligand();
        let grid = GridBox::with_spacing([1.0, 0.0, 0.0], [12.0; 3], 0.5).unwrap();
        // Treat receptor atom 1 as a flexible side-chain atom, with one
        // alternate conformation that moves it far away.
        let spec = FlexibleSpec::new(
            ser_selection(),
            vec![1],
            vec![vec![Vector3::new(50.0, 50.0, 50.0)]],
        );
        let run = flexible_dock(&r, &lig, &grid, &spec, SearchAlgorithm::fast(), 2, 3).unwrap();
        assert_eq!(run.flexible_residues, 1);
        assert!(!run.poses.is_empty());
        // chosen_conformation is 0 (input) or 1 (the alternate).
        assert!(run.chosen_conformation <= 1);
    }

    #[test]
    fn flexible_dock_rejects_out_of_range_flex_atom() {
        let r = carbon_receptor();
        let lig = one_carbon_ligand();
        let grid = GridBox::with_spacing([1.0, 0.0, 0.0], [12.0; 3], 0.5).unwrap();
        // Atom index 99 does not exist.
        let spec = FlexibleSpec::new(ser_selection(), vec![99], vec![]);
        assert!(flexible_dock(&r, &lig, &grid, &spec, SearchAlgorithm::fast(), 1, 1).is_err());
    }
}
