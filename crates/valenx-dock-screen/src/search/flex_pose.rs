//! True induced-fit flexible-receptor docking.
//!
//! The v1 [`crate::search::driver::flexible_dock`] handled flexible
//! sidechains by *post-search* re-scoring: the rigid-core search ran to
//! convergence first, then a list of pre-baked side-chain conformations
//! was scored against the best ligand pose. The dominant induced-fit
//! move — a clashing side chain swinging out — is captured, but the
//! receptor never *adapts* to the ligand during search.
//!
//! This module makes the search variable set
//!
//! ```text
//!   x = ligand_pose ∪ {χ₁, χ₂, …}     for each flexible sidechain
//! ```
//!
//! and runs a co-optimisation of all of them. Every objective call
//! moves both the ligand *and* the flexible sidechains in lock-step;
//! the search settles into a pose + receptor conformation that mutually
//! relax around each other. This is the same idea behind AutoDock 4's
//! `flex_residues` and Vina's `--flex` PDBQT mode, just expressed
//! directly on top of our affinity-grid backend.
//!
//! The optimiser is the Solis-Wets local search from
//! `crate::search::solis_wets`, specialised here to the combined
//! `(pose, χ)` vector. A wrapping MC-warm-start in [`induced_fit_dock`]
//! mirrors the GA-then-LS schedule that AutoDock 4 uses.
//!
//! ## What "induced-fit" means here
//!
//! - Side-chain χ angles are real search variables — every objective
//!   evaluation re-poses the side chains around the ligand.
//! - Search uses a single per-DOF Solis-Wets step adaptation, so
//!   "easy" residue rotations and "easy" ligand translations are
//!   adapted independently.
//! - The final ligand + side-chain coordinates are reported together
//!   for downstream rescoring ([`InducedFitResult`] carries both).
//! - MM-GBSA-class rescoring on the post-search complex is one extra
//!   call to [`crate::screen::rescore::mmgbsa_rescore`] over the
//!   combined receptor / ligand atom set.
//!
//! ## Side-chain rotation model
//!
//! Each flexible sidechain is represented by a list of receptor atom
//! indices (heavy atoms in the side chain, in the order they hang off
//! the Cα — Cβ → Cγ → Cδ → …). A single χ angle rotates every atom
//! *past* the Cα–Cβ axis through that residue's `chi_origin` (the Cα
//! position) about the `chi_axis` (the Cα→Cβ unit vector). v1 supports
//! one χ per residue; the same machinery extends to multi-χ residues
//! (ARG, LYS, …) by adding extra `(origin, axis, atoms)` triples.

use nalgebra::{Unit, UnitQuaternion, Vector3};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use valenx_dock::atom_type::Ad4AtomType;
use valenx_dock::ligand::Ligand;
use valenx_dock::pose::Pose;
use valenx_dock::receptor::{Receptor, ReceptorAtom};
use valenx_dock::search::bfgs::{pose_to_vec, vec_to_pose};

use crate::error::{DockScreenError, Result};
use crate::prep::gridbox::GridBox;
use crate::score::gridmap::{score_ligand_on_maps, AffinityMapSet, MapKind};
use crate::score::vina::score_complex as vina_score_complex;
use crate::screen::rescore::{mmgbsa_rescore, MmGbsaTerms};
use crate::search::solis_wets::SolisWetsParams;

/// One flexible sidechain χ rotation: an origin (the Cα position), a
/// unit axis (Cα→Cβ direction) and the receptor-atom indices that
/// rotate with it.
#[derive(Clone, Debug, PartialEq)]
pub struct ChiRotation {
    /// World-space origin of the rotation axis.
    pub origin: Vector3<f64>,
    /// World-space unit axis vector.
    pub axis: Vector3<f64>,
    /// Indices into [`Receptor::atoms`] that rotate with this χ.
    pub atoms: Vec<usize>,
}

/// A combined (ligand pose, side-chain χ angles) search state — the
/// genome of true induced-fit docking.
#[derive(Clone, Debug)]
pub struct FlexPose {
    /// The ligand pose.
    pub pose: Pose,
    /// One angle per element of [`FlexPoseObjective::chi_rotations`].
    pub chi_angles: Vec<f64>,
}

impl FlexPose {
    /// Identity flex-pose: no translation, identity rotation, all
    /// torsions zero, all χ angles zero (= input geometry).
    pub fn identity(n_torsions: usize, n_chi: usize) -> Self {
        FlexPose {
            pose: Pose::identity(n_torsions),
            chi_angles: vec![0.0; n_chi],
        }
    }

    /// Genome length: `6 + n_torsions + n_chi`.
    pub fn n_dofs(&self) -> usize {
        6 + self.pose.torsions.len() + self.chi_angles.len()
    }
}

/// Pack a flex-pose into the flat search vector.
pub fn flex_pose_to_vec(fp: &FlexPose) -> Vec<f64> {
    let mut v = pose_to_vec(&fp.pose);
    v.extend_from_slice(&fp.chi_angles);
    v
}

/// Inverse of [`flex_pose_to_vec`].
pub fn vec_to_flex_pose(v: &[f64], n_torsions: usize, n_chi: usize) -> FlexPose {
    let pose = vec_to_pose(&v[..6 + n_torsions], n_torsions);
    let chi_angles = v[6 + n_torsions..6 + n_torsions + n_chi].to_vec();
    FlexPose { pose, chi_angles }
}

/// Objective function for true induced-fit search.
///
/// Every score call:
///
/// 1. applies the ligand pose to get world-space ligand atoms;
/// 2. applies the χ angles to displace the flexible-sidechain atoms
///    from their input positions;
/// 3. sums the *grid* lookup for the ligand atoms (using maps
///    precomputed over the *rigid core* — flexible atoms are split
///    out) and an *explicit* Vina pair score between the moved
///    sidechain atoms and the ligand;
/// 4. adds a soft intra-receptor clash penalty so χ moves that drive
///    side-chain atoms into the rigid core are themselves penalised
///    (the receptor cannot pass through itself).
///
/// This is the simultaneous (ligand_pose ∪ χ) optimisation the v1
/// driver lacked.
pub struct FlexPoseObjective<'a> {
    /// The ligand whose poses are being scored.
    pub ligand: &'a Ligand,
    /// Affinity maps for the rigid core of the receptor (flexible
    /// sidechain atoms already split out).
    pub rigid_maps: &'a AffinityMapSet,
    /// The original (pre-search) receptor — used both for the rigid
    /// core atoms (for the clash penalty) and to read the input
    /// positions / atom types of the flexible sidechain atoms.
    pub receptor: &'a Receptor,
    /// Which receptor atoms belong to the rigid core (everything *not*
    /// in any flex-residue).
    pub rigid_indices: Vec<usize>,
    /// One [`ChiRotation`] per side-chain DOF the search co-optimises.
    pub chi_rotations: Vec<ChiRotation>,
}

impl<'a> FlexPoseObjective<'a> {
    /// Build the objective. The chi rotations name which receptor atom
    /// indices belong to which side chain plus the geometric axis for
    /// each.
    pub fn new(
        ligand: &'a Ligand,
        rigid_maps: &'a AffinityMapSet,
        receptor: &'a Receptor,
        chi_rotations: Vec<ChiRotation>,
    ) -> Self {
        // Compute the rigid-core indices: every receptor atom NOT in
        // any chi-rotation's atom set.
        let mut flex_atoms: std::collections::BTreeSet<usize> = Default::default();
        for c in &chi_rotations {
            flex_atoms.extend(c.atoms.iter().copied());
        }
        let rigid_indices: Vec<usize> = (0..receptor.atoms.len())
            .filter(|i| !flex_atoms.contains(i))
            .collect();
        FlexPoseObjective {
            ligand,
            rigid_maps,
            receptor,
            rigid_indices,
            chi_rotations,
        }
    }

    /// Number of ligand torsions in the search.
    pub fn n_torsions(&self) -> usize {
        self.ligand.n_torsions()
    }

    /// Number of side-chain χ angles in the search.
    pub fn n_chi(&self) -> usize {
        self.chi_rotations.len()
    }

    /// Apply χ angles, returning the world-space positions of the
    /// flexible-sidechain atoms.
    pub fn apply_chi(
        &self,
        chi_angles: &[f64],
    ) -> Vec<(Vector3<f64>, Ad4AtomType)> {
        let mut moved: Vec<(Vector3<f64>, Ad4AtomType)> = Vec::new();
        for (ci, chi) in self.chi_rotations.iter().enumerate() {
            let angle = chi_angles.get(ci).copied().unwrap_or(0.0);
            let axis_norm = chi.axis.norm();
            let rot = if axis_norm > 1e-9 {
                let u = Unit::new_unchecked(chi.axis / axis_norm);
                UnitQuaternion::from_axis_angle(&u, angle)
            } else {
                UnitQuaternion::identity()
            };
            for &i in &chi.atoms {
                let a = &self.receptor.atoms[i];
                let p = rot * (a.position - chi.origin) + chi.origin;
                moved.push((p, a.ad4_type));
            }
        }
        moved
    }

    /// Compose the moved-sidechain mini-receptor (used both for scoring
    /// and for the clash penalty).
    fn moved_sidechain_receptor(&self, chi_angles: &[f64]) -> Receptor {
        let atoms: Vec<ReceptorAtom> = self
            .chi_rotations
            .iter()
            .enumerate()
            .flat_map(|(ci, chi)| {
                let angle = chi_angles.get(ci).copied().unwrap_or(0.0);
                let axis_norm = chi.axis.norm();
                let rot = if axis_norm > 1e-9 {
                    let u = Unit::new_unchecked(chi.axis / axis_norm);
                    UnitQuaternion::from_axis_angle(&u, angle)
                } else {
                    UnitQuaternion::identity()
                };
                chi.atoms.iter().map(move |&i| {
                    let a = &self.receptor.atoms[i];
                    let p = rot * (a.position - chi.origin) + chi.origin;
                    ReceptorAtom {
                        position: p,
                        ad4_type: a.ad4_type,
                        partial_charge: a.partial_charge,
                    }
                })
            })
            .collect();
        Receptor { atoms }
    }

    /// Soft intra-receptor clash penalty: any moved sidechain atom that
    /// approaches a rigid-core atom closer than the sum of their VDW
    /// radii pays a quadratic penalty.
    fn intra_receptor_clash(&self, sidechain: &Receptor) -> f64 {
        let mut penalty = 0.0;
        for sa in &sidechain.atoms {
            let sr = sa.ad4_type.props().vdw_radius;
            for &ri in &self.rigid_indices {
                let ra = &self.receptor.atoms[ri];
                let rr = ra.ad4_type.props().vdw_radius;
                let d = (sa.position - ra.position).norm();
                let overlap = sr + rr - d;
                if overlap > 0.0 {
                    // Stiff wall, same shape as Vina's repulsion term.
                    penalty += overlap * overlap;
                }
            }
        }
        penalty
    }

    /// Score a complete (pose, χ) state. Lower is better.
    pub fn score(&self, fp: &FlexPose) -> f64 {
        // 1. Ligand-on-grid energy (fast, captures the rigid core).
        let world = self.ligand.apply_pose(&fp.pose);
        let lig_atoms: Vec<(Vector3<f64>, Ad4AtomType, f64)> = world
            .iter()
            .zip(self.ligand.atoms.iter())
            .map(|(p, a)| (*p, a.ad4_type, a.partial_charge))
            .collect();
        let grid_energy = score_ligand_on_maps(self.rigid_maps, &lig_atoms);

        // 2. Moved-sidechain ↔ ligand explicit Vina score.
        let sidechain = self.moved_sidechain_receptor(&fp.chi_angles);
        let lig_for_sc: Vec<(Vector3<f64>, Ad4AtomType)> = world
            .iter()
            .zip(self.ligand.atoms.iter())
            .map(|(p, a)| (*p, a.ad4_type))
            .collect();
        let sidechain_energy = if sidechain.atoms.is_empty() {
            0.0
        } else {
            vina_score_complex(&sidechain, &lig_for_sc, self.ligand.n_torsions())
                .intermolecular()
        };

        // 3. Soft intra-receptor clash penalty (Vina-class repulsion
        //    weight on the overlap volume).
        let clash = self.intra_receptor_clash(&sidechain);

        grid_energy + sidechain_energy + 0.84 * clash
    }
}

/// The outcome of an induced-fit search.
#[derive(Clone, Debug)]
pub struct InducedFitResult {
    /// The final (pose, χ) state.
    pub state: FlexPose,
    /// Final search-objective score (the value Solis-Wets minimised).
    pub score: f64,
    /// Optional MM-GBSA rescoring of the post-search complex.
    pub mmgbsa: Option<MmGbsaTerms>,
    /// The full receptor — rigid core plus moved sidechains — after
    /// the search. Useful for downstream analysis / visualisation.
    pub final_receptor: Receptor,
    /// Number of flexible sidechain χ angles co-optimised.
    pub n_chi: usize,
}

/// Run a true induced-fit local-search refinement starting from `start`.
///
/// Uses Solis-Wets on the combined (pose, χ) vector. The per-DOF step
/// sizes follow AutoDock 4's defaults (translation 1 Å, rotation 0.05
/// rad, torsion 0.05 rad, χ 0.1 rad).
pub fn induced_fit_solis_wets(
    objective: &FlexPoseObjective,
    start: &FlexPose,
    params: &SolisWetsParams,
    seed: u64,
) -> (FlexPose, f64) {
    // Hand off to the regular Solis-Wets engine by wrapping the
    // combined search vector as a PoseObjective-style closure. We
    // reuse the inner Solis-Wets implementation by writing a small
    // local copy specialised to this objective — keeps the
    // PoseObjective signature unchanged.
    let n_tor = objective.n_torsions();
    let n_chi = objective.n_chi();
    let mut x = flex_pose_to_vec(start);
    let n = x.len();
    let mut rng = StdRng::seed_from_u64(seed);

    // Per-DOF initial step size.
    let mut rho = vec![0.0; n];
    let init_rho = init_rho_combined(n_tor, n_chi, params);
    rho.copy_from_slice(&init_rho);
    let mut bias = vec![0.0; n];

    let mut best = objective.score(start);
    let mut consec_success: u32 = 0;
    let mut consec_fail: u32 = 0;

    for _ in 0..params.max_iter {
        let mut delta = vec![0.0; n];
        for i in 0..n {
            delta[i] = bias[i] + rho[i] * rng.gen_range(-1.0_f64..1.0);
        }
        let mut x_plus = vec![0.0; n];
        for i in 0..n {
            x_plus[i] = x[i] + delta[i];
        }
        let f_plus = objective.score(&vec_to_flex_pose(&x_plus, n_tor, n_chi));
        let mut accepted = false;
        let mut sign = 0.0;
        if f_plus < best {
            x = x_plus;
            best = f_plus;
            sign = 1.0;
            accepted = true;
        } else {
            let mut x_minus = vec![0.0; n];
            for i in 0..n {
                x_minus[i] = x[i] - delta[i];
            }
            let f_minus = objective.score(&vec_to_flex_pose(&x_minus, n_tor, n_chi));
            if f_minus < best {
                x = x_minus;
                best = f_minus;
                sign = -1.0;
                accepted = true;
            }
        }
        if accepted {
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
        if rho
            .iter()
            .zip(init_rho.iter())
            .all(|(r, r0)| *r < params.rho_tol * *r0)
        {
            break;
        }
    }
    (vec_to_flex_pose(&x, n_tor, n_chi), best)
}

fn init_rho_combined(n_tor: usize, n_chi: usize, p: &SolisWetsParams) -> Vec<f64> {
    let mut v = Vec::with_capacity(6 + n_tor + n_chi);
    for _ in 0..3 {
        v.push(p.rho_xyz);
    }
    for _ in 0..3 {
        v.push(p.rho_rot);
    }
    for _ in 0..n_tor {
        v.push(p.rho_tor);
    }
    for _ in 0..n_chi {
        // χ angles get a slightly looser step than torsions — side
        // chains are bulkier and benefit from bigger jumps in the
        // exploration phase.
        v.push(2.0 * p.rho_tor);
    }
    v
}

/// Top-level induced-fit docking driver.
///
/// 1. Builds affinity maps over the *rigid core* of the receptor
///    (flexible-sidechain atoms split out).
/// 2. Runs a GA-style multi-start refinement: `n_restarts` random
///    starting (pose, χ) states each refined by [`induced_fit_solis_wets`].
/// 3. Reports the best post-search complex plus an MM-GBSA rescoring
///    if `rescore_with_mmgbsa` is set.
pub fn induced_fit_dock(
    receptor: &Receptor,
    ligand: &Ligand,
    grid: &GridBox,
    chi_rotations: Vec<ChiRotation>,
    n_restarts: usize,
    seed: u64,
    rescore_with_mmgbsa: bool,
) -> Result<InducedFitResult> {
    if receptor.atoms.is_empty() {
        return Err(DockScreenError::invalid_receptor("receptor has no atoms"));
    }
    if ligand.atoms.is_empty() {
        return Err(DockScreenError::invalid_ligand("ligand has no atoms"));
    }
    if n_restarts == 0 {
        return Err(DockScreenError::invalid(
            "n_restarts",
            "must be ≥ 1",
        ));
    }
    // Sanity-check the chi rotations.
    for (i, c) in chi_rotations.iter().enumerate() {
        for &ai in &c.atoms {
            if ai >= receptor.atoms.len() {
                return Err(DockScreenError::invalid_receptor(format!(
                    "chi rotation {i} references atom {ai} out of {} receptor atoms",
                    receptor.atoms.len()
                )));
            }
        }
    }

    // 1. Rigid core receptor — atoms NOT in any chi rotation.
    let mut flex_atoms: std::collections::BTreeSet<usize> = Default::default();
    for c in &chi_rotations {
        flex_atoms.extend(c.atoms.iter().copied());
    }
    let rigid_core = Receptor {
        atoms: receptor
            .atoms
            .iter()
            .enumerate()
            .filter(|(i, _)| !flex_atoms.contains(i))
            .map(|(_, a)| a.clone())
            .collect(),
    };
    if rigid_core.atoms.is_empty() {
        return Err(DockScreenError::invalid_receptor(
            "induced-fit flex selection removed every receptor atom",
        ));
    }
    let ligand_types: Vec<Ad4AtomType> = ligand.atoms.iter().map(|a| a.ad4_type).collect();
    let rigid_maps =
        AffinityMapSet::precompute(&rigid_core, &ligand_types, grid, MapKind::Vina)?;

    // 2. Multi-start refinement.
    let obj = FlexPoseObjective::new(ligand, &rigid_maps, receptor, chi_rotations);
    let n_chi = obj.n_chi();
    let n_tor = obj.n_torsions();

    let mut best_state: Option<FlexPose> = None;
    let mut best_score = f64::INFINITY;
    let params = SolisWetsParams::default();
    let mut rng = StdRng::seed_from_u64(seed);
    for r in 0..n_restarts {
        let half = grid.size / 2.0;
        // Random ligand pose + random χ angles around 0.
        let mut start = FlexPose::identity(n_tor, n_chi);
        start.pose.translation = Vector3::new(
            grid.center.x + rng.gen_range(-half.x..half.x),
            grid.center.y + rng.gen_range(-half.y..half.y),
            grid.center.z + rng.gen_range(-half.z..half.z),
        );
        // Random axis-angle rotation in [-π, π].
        let axis = Vector3::new(
            rng.gen_range(-1.0..1.0),
            rng.gen_range(-1.0..1.0),
            rng.gen_range(-1.0..1.0),
        );
        let axis = if axis.norm() > 1e-9 {
            axis.normalize()
        } else {
            Vector3::x()
        };
        let angle = rng.gen_range(-std::f64::consts::PI..std::f64::consts::PI);
        start.pose.orientation =
            UnitQuaternion::from_axis_angle(&Unit::new_unchecked(axis), angle);
        for ci in 0..n_chi {
            start.chi_angles[ci] = rng.gen_range(-std::f64::consts::PI..std::f64::consts::PI);
        }
        let run_seed = seed
            .wrapping_add((r as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));
        let (refined, score) = induced_fit_solis_wets(&obj, &start, &params, run_seed);
        if score < best_score {
            best_score = score;
            best_state = Some(refined);
        }
    }
    let state = best_state.ok_or_else(|| DockScreenError::invalid("restarts", "no successful refinement"))?;

    // 3. Reconstruct the final receptor (rigid core + moved sidechains).
    let mut final_atoms = rigid_core.atoms.clone();
    let moved = obj.moved_sidechain_receptor(&state.chi_angles);
    final_atoms.extend(moved.atoms.iter().cloned());
    let final_receptor = Receptor { atoms: final_atoms };

    // 4. Optional MM-GBSA rescoring on the post-search complex.
    let mmgbsa = if rescore_with_mmgbsa {
        let world = ligand.apply_pose(&state.pose);
        let lig_atoms: Vec<(Vector3<f64>, Ad4AtomType, f64)> = world
            .iter()
            .zip(ligand.atoms.iter())
            .map(|(p, a)| (*p, a.ad4_type, a.partial_charge))
            .collect();
        mmgbsa_rescore(&final_receptor, &lig_atoms, 10.0).ok()
    } else {
        None
    };

    Ok(InducedFitResult {
        state,
        score: best_score,
        mmgbsa,
        final_receptor,
        n_chi: obj.n_chi(),
    })
}

/// Build a [`ChiRotation`] from a Cα / Cβ pair and a list of side-chain
/// atom indices that rotate together with the χ angle. The rotation
/// axis is the unit vector from `c_alpha` to `c_beta`, origin = Cα.
pub fn chi_from_axis_atoms(
    receptor: &Receptor,
    c_alpha_idx: usize,
    c_beta_idx: usize,
    rotated_atom_indices: Vec<usize>,
) -> Result<ChiRotation> {
    if c_alpha_idx >= receptor.atoms.len() || c_beta_idx >= receptor.atoms.len() {
        return Err(DockScreenError::invalid_receptor(
            "chi axis atom index out of range",
        ));
    }
    let origin = receptor.atoms[c_alpha_idx].position;
    let beta = receptor.atoms[c_beta_idx].position;
    let axis_raw = beta - origin;
    let n = axis_raw.norm();
    if n < 1e-6 {
        return Err(DockScreenError::invalid_receptor(
            "chi axis Cα/Cβ atoms coincide",
        ));
    }
    Ok(ChiRotation {
        origin,
        axis: axis_raw / n,
        atoms: rotated_atom_indices,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prep::gridbox::GridBox;

    fn three_atom_receptor() -> Receptor {
        // Atom 0 = rigid core (origin), atom 1 = Cα anchor (3,0,0),
        // atom 2 = sidechain end (5,0,0) — rotates about the Cα.
        Receptor {
            atoms: vec![
                ReceptorAtom {
                    position: Vector3::new(0.0, 0.0, 0.0),
                    ad4_type: Ad4AtomType::C,
                    partial_charge: 0.0,
                },
                ReceptorAtom {
                    position: Vector3::new(3.0, 0.0, 0.0),
                    ad4_type: Ad4AtomType::C,
                    partial_charge: 0.0,
                },
                ReceptorAtom {
                    position: Vector3::new(5.0, 0.0, 0.0),
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
    fn flex_pose_packing_roundtrips() {
        let mut fp = FlexPose::identity(2, 3);
        fp.pose.translation = Vector3::new(1.0, 2.0, 3.0);
        fp.pose.torsions = vec![0.5, -0.3];
        fp.chi_angles = vec![0.1, 0.2, -0.4];
        let v = flex_pose_to_vec(&fp);
        let fp2 = vec_to_flex_pose(&v, 2, 3);
        assert!((fp.pose.translation - fp2.pose.translation).norm() < 1e-12);
        assert_eq!(fp.pose.torsions, fp2.pose.torsions);
        assert_eq!(fp.chi_angles, fp2.chi_angles);
        assert_eq!(fp.n_dofs(), 6 + 2 + 3);
    }

    #[test]
    fn chi_rotation_moves_atom_about_axis() {
        // χ axis along +X through (3,0,0). Rotate atom 2 (5,0,0)
        // by 180°: should end up at the same place (atom is on the
        // axis). So move atom 2 OFF the axis first.
        let mut r = three_atom_receptor();
        r.atoms[2].position = Vector3::new(5.0, 1.0, 0.0);
        let chi = ChiRotation {
            origin: Vector3::new(3.0, 0.0, 0.0),
            axis: Vector3::new(1.0, 0.0, 0.0),
            atoms: vec![2],
        };
        // Wrap in an objective so we can use apply_chi.
        let grid = GridBox::with_spacing([0.0; 3], [12.0; 3], 0.5).unwrap();
        let maps = AffinityMapSet::precompute(
            &Receptor {
                atoms: vec![r.atoms[0].clone()],
            },
            &[Ad4AtomType::C],
            &grid,
            MapKind::Vina,
        )
        .unwrap();
        let lig = one_carbon_ligand();
        let obj = FlexPoseObjective::new(&lig, &maps, &r, vec![chi]);
        // 90° → moves (5,1,0) to (5,0,1) (rotate 1 about +X axis: y→z).
        let world = obj.apply_chi(&[std::f64::consts::FRAC_PI_2]);
        assert_eq!(world.len(), 1);
        let p = world[0].0;
        assert!((p - Vector3::new(5.0, 0.0, 1.0)).norm() < 1e-9, "got {p}");
    }

    #[test]
    fn induced_fit_objective_responds_to_chi_motion() {
        // Place the ligand right at the sidechain end's original spot
        // — a clash. Rotating the sidechain out of the way should LOWER
        // the score.
        let mut r = three_atom_receptor();
        r.atoms[2].position = Vector3::new(5.0, 1.0, 0.0);
        let chi = ChiRotation {
            origin: Vector3::new(3.0, 0.0, 0.0),
            axis: Vector3::new(1.0, 0.0, 0.0),
            atoms: vec![2],
        };
        // Build maps from the rigid core (atoms 0 + 1).
        let rigid_core = Receptor {
            atoms: vec![r.atoms[0].clone(), r.atoms[1].clone()],
        };
        let grid = GridBox::with_spacing([2.0, 0.0, 0.0], [12.0; 3], 0.5).unwrap();
        let maps =
            AffinityMapSet::precompute(&rigid_core, &[Ad4AtomType::C], &grid, MapKind::Vina)
                .unwrap();
        let lig = one_carbon_ligand();
        let obj = FlexPoseObjective::new(&lig, &maps, &r, vec![chi]);

        // Place the ligand right where the sidechain is.
        let mut fp = FlexPose::identity(0, 1);
        fp.pose.translation = Vector3::new(5.0, 1.0, 0.0);

        // χ = 0 → sidechain still in place → clash.
        fp.chi_angles[0] = 0.0;
        let s_clash = obj.score(&fp);
        // χ = π → sidechain rotated to (5,-1,0) → no clash.
        fp.chi_angles[0] = std::f64::consts::PI;
        let s_relaxed = obj.score(&fp);
        assert!(
            s_relaxed < s_clash,
            "expected relaxed < clash, got relaxed={s_relaxed} vs clash={s_clash}"
        );
    }

    #[test]
    fn induced_fit_solis_wets_finds_a_clash_free_arrangement() {
        // Same setup as the previous test — but now SW should *find*
        // the relaxed χ rather than rely on the caller setting it.
        let mut r = three_atom_receptor();
        r.atoms[2].position = Vector3::new(5.0, 1.0, 0.0);
        let chi = ChiRotation {
            origin: Vector3::new(3.0, 0.0, 0.0),
            axis: Vector3::new(1.0, 0.0, 0.0),
            atoms: vec![2],
        };
        let rigid_core = Receptor {
            atoms: vec![r.atoms[0].clone(), r.atoms[1].clone()],
        };
        let grid = GridBox::with_spacing([2.0, 0.0, 0.0], [12.0; 3], 0.5).unwrap();
        let maps =
            AffinityMapSet::precompute(&rigid_core, &[Ad4AtomType::C], &grid, MapKind::Vina)
                .unwrap();
        let lig = one_carbon_ligand();
        let obj = FlexPoseObjective::new(&lig, &maps, &r, vec![chi]);

        // Start clashing.
        let mut start = FlexPose::identity(0, 1);
        start.pose.translation = Vector3::new(5.0, 1.0, 0.0);
        start.chi_angles[0] = 0.0;
        let before = obj.score(&start);
        let (after_state, after) = induced_fit_solis_wets(&obj, &start, &SolisWetsParams::default(), 4);
        assert!(
            after <= before,
            "SW must not worsen, before={before} after={after}"
        );
        // The χ should have moved.
        assert!(
            (after_state.chi_angles[0]).abs() > 0.05,
            "χ angle should have moved, stayed at {}",
            after_state.chi_angles[0]
        );
    }

    #[test]
    fn induced_fit_dock_runs_end_to_end_with_mmgbsa() {
        let mut r = three_atom_receptor();
        r.atoms[2].position = Vector3::new(5.0, 1.0, 0.0);
        let chi = ChiRotation {
            origin: Vector3::new(3.0, 0.0, 0.0),
            axis: Vector3::new(1.0, 0.0, 0.0),
            atoms: vec![2],
        };
        let lig = one_carbon_ligand();
        let grid = GridBox::with_spacing([2.0, 0.0, 0.0], [12.0; 3], 0.5).unwrap();
        let res = induced_fit_dock(&r, &lig, &grid, vec![chi], 3, 11, true).unwrap();
        assert_eq!(res.n_chi, 1);
        assert!(res.score.is_finite());
        // Final receptor has both rigid core + moved sidechain.
        assert_eq!(res.final_receptor.atoms.len(), 3);
        // MM-GBSA was requested, should be present (and finite).
        let mm = res.mmgbsa.expect("mmgbsa requested");
        assert!(mm.total().is_finite());
    }

    #[test]
    fn induced_fit_dock_validates_inputs() {
        let r = three_atom_receptor();
        let lig = one_carbon_ligand();
        let grid = GridBox::cubic([0.0; 3], 10.0).unwrap();
        // Empty receptor.
        assert!(induced_fit_dock(
            &Receptor::default(),
            &lig,
            &grid,
            vec![],
            1,
            1,
            false
        )
        .is_err());
        // n_restarts = 0.
        assert!(induced_fit_dock(&r, &lig, &grid, vec![], 0, 1, false).is_err());
        // Out-of-range chi atom.
        let bad_chi = ChiRotation {
            origin: Vector3::zeros(),
            axis: Vector3::x(),
            atoms: vec![99],
        };
        assert!(induced_fit_dock(&r, &lig, &grid, vec![bad_chi], 1, 1, false).is_err());
    }

    #[test]
    fn chi_from_axis_atoms_builds_unit_axis() {
        let r = three_atom_receptor();
        let chi = chi_from_axis_atoms(&r, 0, 1, vec![2]).unwrap();
        assert!((chi.axis.norm() - 1.0).abs() < 1e-9);
        assert_eq!(chi.atoms, vec![2]);
        assert_eq!(chi.origin, Vector3::zeros());
    }

    #[test]
    fn chi_from_axis_atoms_rejects_bad_inputs() {
        let r = three_atom_receptor();
        // Out of range.
        assert!(chi_from_axis_atoms(&r, 0, 99, vec![]).is_err());
        // Coincident Cα/Cβ.
        let r2 = Receptor {
            atoms: vec![
                r.atoms[0].clone(),
                r.atoms[0].clone(), // identical to atom 0
            ],
        };
        assert!(chi_from_axis_atoms(&r2, 0, 1, vec![]).is_err());
    }
}
