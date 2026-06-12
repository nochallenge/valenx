//! Feature 15 — pose clustering (RMSD-based) and ranking.
//!
//! A docking search returns dozens of poses, but many are tiny
//! variations on the same binding mode. Clustering collapses them:
//! poses within an RMSD cutoff of each other form a cluster, and the
//! lowest-energy member is the cluster representative. The number of
//! members is a useful confidence signal — a deep, well-populated
//! cluster is more trustworthy than a lone low-energy pose.
//!
//! The RMSD itself is heavy-atom RMSD from
//! [`valenx_dock::cluster::rmsd`] (hydrogens excluded, the AutoDock
//! convention). This module adds the cluster *grouping* and a richer
//! [`PoseCluster`] type carrying the member list and population.
//!
//! The algorithm is greedy single-linkage, the same one AutoDock uses
//! by default: take the lowest-energy unclustered pose as a new
//! cluster seed, swallow every pose within `rmsd_cutoff` of it,
//! repeat.

use valenx_dock::cluster::rmsd;
use valenx_dock::ligand::Ligand;
use valenx_dock::pose::Pose;

use crate::error::{DockScreenError, Result};
use crate::search::driver::ScoredPose;

/// One cluster of docked poses sharing a binding mode.
#[derive(Clone, Debug)]
pub struct PoseCluster {
    /// The lowest-energy pose in the cluster — its representative.
    pub representative: Pose,
    /// Score of the representative (the cluster's best score).
    pub best_score: f64,
    /// Every pose in the cluster, including the representative, sorted
    /// best-first.
    pub members: Vec<ScoredPose>,
}

impl PoseCluster {
    /// Number of poses in the cluster.
    pub fn population(&self) -> usize {
        self.members.len()
    }

    /// Mean score across the cluster's members.
    pub fn mean_score(&self) -> f64 {
        if self.members.is_empty() {
            return 0.0;
        }
        self.members.iter().map(|m| m.score).sum::<f64>() / self.members.len() as f64
    }
}

/// Cluster a set of scored poses by heavy-atom RMSD.
///
/// Returns clusters sorted by representative score (best first). Each
/// cluster's lowest-energy pose is its representative; every pair
/// inside a cluster need not be within the cutoff (single-linkage),
/// but every member is within `rmsd_cutoff` of the seed.
///
/// Returns [`DockScreenError::Invalid`] if `rmsd_cutoff` is
/// non-positive.
pub fn cluster_poses(
    ligand: &Ligand,
    poses: &[ScoredPose],
    rmsd_cutoff: f64,
) -> Result<Vec<PoseCluster>> {
    if !rmsd_cutoff.is_finite() || rmsd_cutoff <= 0.0 {
        return Err(DockScreenError::invalid(
            "rmsd_cutoff",
            format!("clustering RMSD cutoff must be positive, got {rmsd_cutoff}"),
        ));
    }
    // Sort poses ascending by score so the lowest-energy pose seeds
    // each new cluster.
    let mut order: Vec<usize> = (0..poses.len()).collect();
    order.sort_by(|&a, &b| {
        poses[a]
            .score
            .partial_cmp(&poses[b].score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut taken = vec![false; poses.len()];
    let mut clusters: Vec<PoseCluster> = Vec::new();
    for &seed in &order {
        if taken[seed] {
            continue;
        }
        taken[seed] = true;
        let mut members = vec![poses[seed].clone()];
        for &other in &order {
            if taken[other] {
                continue;
            }
            if rmsd(ligand, &poses[seed].pose, &poses[other].pose) <= rmsd_cutoff {
                taken[other] = true;
                members.push(poses[other].clone());
            }
        }
        members.sort_by(|a, b| {
            a.score
                .partial_cmp(&b.score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        clusters.push(PoseCluster {
            representative: members[0].pose.clone(),
            best_score: members[0].score,
            members,
        });
    }
    // Clusters are already discovered best-first (seeds came from the
    // score-sorted order), but sort defensively.
    clusters.sort_by(|a, b| {
        a.best_score
            .partial_cmp(&b.best_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Ok(clusters)
}

/// Return the cluster representatives only — one best pose per binding
/// mode, ranked. The common "give me the top N distinct poses" call.
pub fn ranked_representatives(
    ligand: &Ligand,
    poses: &[ScoredPose],
    rmsd_cutoff: f64,
) -> Result<Vec<ScoredPose>> {
    let clusters = cluster_poses(ligand, poses, rmsd_cutoff)?;
    Ok(clusters
        .into_iter()
        .map(|c| ScoredPose {
            pose: c.representative,
            score: c.best_score,
        })
        .collect())
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

    fn pose_at(x: f64, score: f64) -> ScoredPose {
        let mut p = Pose::identity(0);
        p.translation = Vector3::new(x, 0.0, 0.0);
        ScoredPose { pose: p, score }
    }

    #[test]
    fn rejects_bad_cutoff() {
        let lig = one_atom_ligand();
        assert!(cluster_poses(&lig, &[], 0.0).is_err());
        assert!(cluster_poses(&lig, &[], -1.0).is_err());
    }

    #[test]
    fn nearby_poses_merge_into_one_cluster() {
        let lig = one_atom_ligand();
        let poses = vec![
            pose_at(0.0, -8.0),
            pose_at(0.5, -7.5),  // 0.5 Å from the first — same cluster
            pose_at(10.0, -6.0), // far — its own cluster
        ];
        let clusters = cluster_poses(&lig, &poses, 1.0).unwrap();
        assert_eq!(clusters.len(), 2);
        // The lowest-energy cluster comes first; it has two members.
        assert_eq!(clusters[0].population(), 2);
        assert!((clusters[0].best_score - -8.0).abs() < 1e-9);
        assert_eq!(clusters[1].population(), 1);
    }

    #[test]
    fn representative_is_lowest_energy_member() {
        let lig = one_atom_ligand();
        // Two near poses; the SECOND has the better score.
        let poses = vec![pose_at(0.0, -5.0), pose_at(0.3, -9.0)];
        let clusters = cluster_poses(&lig, &poses, 1.0).unwrap();
        assert_eq!(clusters.len(), 1);
        assert!((clusters[0].best_score - -9.0).abs() < 1e-9);
        // The representative pose is the -9.0 one (at x = 0.3).
        assert!((clusters[0].representative.translation.x - 0.3).abs() < 1e-9);
    }

    #[test]
    fn ranked_representatives_returns_one_per_cluster() {
        let lig = one_atom_ligand();
        let poses = vec![pose_at(0.0, -8.0), pose_at(0.4, -7.0), pose_at(20.0, -9.0)];
        let reps = ranked_representatives(&lig, &poses, 1.0).unwrap();
        assert_eq!(reps.len(), 2);
        // Best-first ranking: the -9.0 pose leads.
        assert!((reps[0].score - -9.0).abs() < 1e-9);
    }

    #[test]
    fn mean_score_averages_members() {
        let lig = one_atom_ligand();
        let poses = vec![pose_at(0.0, -8.0), pose_at(0.2, -6.0)];
        let clusters = cluster_poses(&lig, &poses, 1.0).unwrap();
        assert!((clusters[0].mean_score() - -7.0).abs() < 1e-9);
    }

    #[test]
    fn empty_input_yields_no_clusters() {
        let lig = one_atom_ligand();
        let clusters = cluster_poses(&lig, &[], 2.0).unwrap();
        assert!(clusters.is_empty());
    }
}
