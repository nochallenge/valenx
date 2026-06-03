//! Pose evaluation: transform ligand, look up in receptor grids, sum.

use crate::grid::GridBundle;
use crate::ligand::Ligand;
use crate::pose::Pose;

/// Inter-molecular score for `pose`, using precomputed grids.
pub fn inter_score(ligand: &Ligand, pose: &Pose, grids: &GridBundle) -> f64 {
    let world = ligand.apply_pose(pose);
    let mut sum = 0.0;
    for (atom, p) in ligand.atoms.iter().zip(world.iter()) {
        let Some(g) = grids.grids.get(&atom.ad4_type) else {
            // Should never happen if bundle was built from this ligand.
            continue;
        };
        sum += g.sample(*p);
    }
    sum
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom_type::Ad4AtomType;
    use crate::grid::GridBundle;
    use crate::ligand::Ligand;
    use crate::pose::Pose;
    use crate::receptor::{Receptor, ReceptorAtom};
    use nalgebra::Vector3;

    #[test]
    fn moving_attractive_pair_toward_each_other_lowers_score() {
        // Receptor: one C at origin. Ligand: one C, varied translation.
        let receptor = Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::zeros(),
                ad4_type: Ad4AtomType::C,
                partial_charge: 0.0,
            }],
        };
        let pdbqt = "\
ROOT
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
        // "close" sits near the Vina attractive minimum (surface
        // distance ≈ 0 for two C atoms with 1.9 Å VDW each, so the
        // centroid-to-centroid sweet spot is ≈ 3.8 Å). "far" sits in
        // the attractive tail. Anything closer than 3.8 Å walks into
        // the repulsion wall and the score climbs sharply.
        let mut far = Pose::identity(0);
        far.translation = Vector3::new(6.0, 0.0, 0.0);
        let mut close = Pose::identity(0);
        close.translation = Vector3::new(4.0, 0.0, 0.0);
        let s_far = inter_score(&lig, &far, &grids);
        let s_close = inter_score(&lig, &close, &grids);
        assert!(s_close < s_far, "close={s_close} should be < far={s_far}");
    }
}
