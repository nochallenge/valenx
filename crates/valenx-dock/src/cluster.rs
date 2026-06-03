//! Pose clustering by heavy-atom RMSD.

use crate::ligand::Ligand;
use crate::pose::Pose;

/// Heavy-atom RMSD between two poses of the same ligand. Hydrogens
/// (AD4 type H or HD) are excluded â same convention as reference Vina.
pub fn rmsd(ligand: &Ligand, a: &Pose, b: &Pose) -> f64 {
    use crate::atom_type::Ad4AtomType;
    let pa = ligand.apply_pose(a);
    let pb = ligand.apply_pose(b);
    let mut sum_sq = 0.0;
    let mut n = 0;
    for (i, atom) in ligand.atoms.iter().enumerate() {
        if matches!(atom.ad4_type, Ad4AtomType::H | Ad4AtomType::HD) {
            continue;
        }
        sum_sq += (pa[i] - pb[i]).norm_squared();
        n += 1;
    }
    if n == 0 {
        return 0.0;
    }
    (sum_sq / n as f64).sqrt()
}

/// Group poses so that within each cluster every pair has RMSD <= `tol`.
/// Returns clusters sorted by their best (lowest) score; within each
/// cluster the lowest-score pose is the representative.
///
/// Simple greedy single-linkage: take the lowest-energy unclustered
/// pose as the next centroid, swallow everything within `tol` of it.
pub fn cluster_poses(ligand: &Ligand, poses: &[(Pose, f64)], tol: f64) -> Vec<(Pose, f64)> {
    let mut sorted: Vec<&(Pose, f64)> = poses.iter().collect();
    sorted.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    let mut taken = vec![false; sorted.len()];
    let mut out = Vec::new();
    for i in 0..sorted.len() {
        if taken[i] {
            continue;
        }
        out.push(sorted[i].clone());
        taken[i] = true;
        for j in (i + 1)..sorted.len() {
            if taken[j] {
                continue;
            }
            if rmsd(ligand, &sorted[i].0, &sorted[j].0) <= tol {
                taken[j] = true;
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    fn one_atom_ligand() -> Ligand {
        let pdbqt = "ROOT
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C 
ENDROOT
TORSDOF 0
";
        Ligand::from_pdbqt(pdbqt).unwrap()
    }

    #[test]
    fn rmsd_identity_is_zero() {
        let lig = one_atom_ligand();
        let p = Pose::identity(0);
        assert_eq!(rmsd(&lig, &p, &p), 0.0);
    }

    #[test]
    fn rmsd_pure_translation_equals_distance() {
        let lig = one_atom_ligand();
        let mut a = Pose::identity(0);
        let mut b = Pose::identity(0);
        a.translation = Vector3::zeros();
        b.translation = Vector3::new(3.0, 4.0, 0.0);
        let d = rmsd(&lig, &a, &b);
        assert!((d - 5.0).abs() < 1e-9);
    }

    #[test]
    fn clustering_merges_nearby_poses() {
        let lig = one_atom_ligand();
        let mut p1 = Pose::identity(0);
        p1.translation = Vector3::new(0.0, 0.0, 0.0);
        let mut p2 = Pose::identity(0);
        p2.translation = Vector3::new(0.5, 0.0, 0.0);
        let mut p3 = Pose::identity(0);
        p3.translation = Vector3::new(5.0, 0.0, 0.0);
        let poses = vec![
            (p1, -2.0),
            (p2, -1.5), // within 1.0 of p1 -> swallowed
            (p3, -1.0), // far -> kept
        ];
        let clusters = cluster_poses(&lig, &poses, 1.0);
        assert_eq!(clusters.len(), 2);
        // Best score first.
        assert!((clusters[0].1 - -2.0).abs() < 1e-9);
        assert!((clusters[1].1 - -1.0).abs() < 1e-9);
    }
}
