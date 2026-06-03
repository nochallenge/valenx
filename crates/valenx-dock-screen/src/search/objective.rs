//! The pose objective function shared by every search algorithm.
//!
//! Each search strategy ([`crate::search::ga`], [`crate::search::mc`],
//! [`crate::search::ils`]) needs the same thing: a way to take a
//! candidate [`valenx_dock::pose::Pose`] and return its energy.
//! That energy must be cheap — the search inner loop calls it
//! thousands of times.
//!
//! [`PoseObjective`] bundles a ligand and a precomputed
//! [`AffinityMapSet`] so the energy of a pose is: apply the pose to
//! get world-space atom positions, then sum a trilinear grid lookup
//! per atom. No receptor-pair loop, no closures captured across thread
//! boundaries — the search algorithms hold a `&PoseObjective` and call
//! [`PoseObjective::score`].

use valenx_dock::atom_type::Ad4AtomType;
use valenx_dock::ligand::Ligand;
use valenx_dock::pose::Pose;

use crate::score::gridmap::{score_ligand_on_maps, AffinityMapSet};

/// A grid-backed pose energy function: a ligand plus the receptor
/// affinity maps its atoms are scored against.
pub struct PoseObjective<'a> {
    /// The ligand whose poses are being scored.
    pub ligand: &'a Ligand,
    /// Precomputed per-atom-type affinity maps for the receptor.
    pub maps: &'a AffinityMapSet,
    /// Per-atom partial charges, parallel to `ligand.atoms`. Used only
    /// by the AutoDock4 electrostatic map; an all-zeros vector is fine
    /// for the Vina-class function.
    pub charges: Vec<f64>,
}

impl<'a> PoseObjective<'a> {
    /// Build an objective. `charges` must have one entry per ligand
    /// atom; if it does not, it is padded / truncated to fit (a short
    /// vector is treated as all-zeros tail — the Vina function ignores
    /// charges anyway).
    pub fn new(ligand: &'a Ligand, maps: &'a AffinityMapSet, charges: Vec<f64>) -> Self {
        let mut q = charges;
        q.resize(ligand.atoms.len(), 0.0);
        PoseObjective {
            ligand,
            maps,
            charges: q,
        }
    }

    /// Build an objective with all-zero charges — the common case for
    /// the Vina-class scoring function.
    pub fn uncharged(ligand: &'a Ligand, maps: &'a AffinityMapSet) -> Self {
        let q = vec![0.0; ligand.atoms.len()];
        PoseObjective {
            ligand,
            maps,
            charges: q,
        }
    }

    /// Energy of `pose` in scoring-function units (kcal/mol). Lower is
    /// better. This is the value every search algorithm minimises.
    pub fn score(&self, pose: &Pose) -> f64 {
        let world = self.ligand.apply_pose(pose);
        let atoms: Vec<(nalgebra::Vector3<f64>, Ad4AtomType, f64)> = world
            .iter()
            .zip(self.ligand.atoms.iter())
            .enumerate()
            .map(|(i, (p, a))| (*p, a.ad4_type, self.charges[i]))
            .collect();
        score_ligand_on_maps(self.maps, &atoms)
    }

    /// Number of torsional degrees of freedom of the ligand.
    pub fn n_torsions(&self) -> usize {
        self.ligand.n_torsions()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prep::gridbox::GridBox;
    use crate::score::gridmap::MapKind;
    use nalgebra::Vector3;
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
        let grid = GridBox::with_spacing([0.0; 3], [12.0; 3], 0.375).unwrap();
        let maps =
            AffinityMapSet::precompute(&receptor, &[Ad4AtomType::C], &grid, MapKind::Vina).unwrap();
        (lig, maps)
    }

    #[test]
    fn closer_pose_scores_lower_than_far_pose() {
        let (lig, maps) = setup();
        let obj = PoseObjective::uncharged(&lig, &maps);
        let mut far = Pose::identity(0);
        far.translation = Vector3::new(6.0, 0.0, 0.0);
        let mut near = Pose::identity(0);
        near.translation = Vector3::new(3.8, 0.0, 0.0);
        assert!(
            obj.score(&near) < obj.score(&far),
            "near pose should score better"
        );
    }

    #[test]
    fn charges_resize_to_atom_count() {
        let (lig, maps) = setup();
        // An over-long charge vector is truncated; a short one padded.
        let obj = PoseObjective::new(&lig, &maps, vec![0.1, 0.2, 0.3]);
        assert_eq!(obj.charges.len(), lig.atoms.len());
        let obj2 = PoseObjective::new(&lig, &maps, vec![]);
        assert_eq!(obj2.charges.len(), lig.atoms.len());
    }

    #[test]
    fn n_torsions_passes_through() {
        let (lig, maps) = setup();
        let obj = PoseObjective::uncharged(&lig, &maps);
        assert_eq!(obj.n_torsions(), 0);
    }
}
