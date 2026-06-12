//! **Feature 12 — decoy generation, clustering and model selection.**
//!
//! A single ab-initio trajectory is noisy: the lowest-energy model
//! from one run may not be the best fold. The classical remedy
//! (Rosetta's `cluster`) is to generate **many decoys** — independent
//! trajectories from different random seeds — then **cluster** them by
//! structural similarity. The largest cluster is the most
//! reproducibly-sampled fold; its centre is the predicted model. This
//! works because the native fold has the widest energy funnel and so
//! is sampled most often.
//!
//! This module:
//!
//! - generates decoys by running [`crate::abinitio::fragment_assembly`]
//!   over a range of seeds;
//! - clusters them by Cα-RMSD with a leader / threshold-clustering
//!   algorithm;
//! - selects a model — by largest cluster (consensus) or by lowest
//!   energy.

use serde::{Deserialize, Serialize};

use crate::abinitio::assemble::{fragment_assembly, AssemblyOptions, AssemblyResult};
use crate::abinitio::fragments::FragmentLibrary;
use crate::error::{Result, StructPredictError};
use crate::model::ProteinModel;
use crate::refine::superpose::ca_rmsd_superposed;

/// A cluster of structurally-similar decoys.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DecoyCluster {
    /// Indices (into the decoy list) of this cluster's members.
    pub members: Vec<usize>,
    /// Index of the cluster centre — the decoy closest on average to
    /// the others in the cluster.
    pub center: usize,
    /// The cluster centre's knowledge-score total.
    pub center_score: f64,
}

impl DecoyCluster {
    /// Cluster size (member count).
    pub fn size(&self) -> usize {
        self.members.len()
    }
}

/// Generates a set of decoys for a sequence.
///
/// Runs [`fragment_assembly`] `count` times, with seeds derived from
/// `base_options.seed` so the whole set is reproducible. Returns the
/// per-run [`AssemblyResult`]s.
///
/// # Errors
/// [`StructPredictError::Invalid`] for `count == 0`; propagates
/// [`fragment_assembly`] errors.
pub fn generate_decoys(
    sequence: &str,
    library: &FragmentLibrary,
    base_options: AssemblyOptions,
    count: usize,
) -> Result<Vec<AssemblyResult>> {
    if count == 0 {
        return Err(StructPredictError::invalid("count", "must be at least 1"));
    }
    let mut decoys = Vec::with_capacity(count);
    for k in 0..count {
        let mut opts = base_options;
        // Distinct, deterministic seed per decoy.
        opts.seed = base_options
            .seed
            .wrapping_add((k as u64).wrapping_mul(0x9E3779B97F4A7C15));
        decoys.push(fragment_assembly(sequence, library, opts)?);
    }
    Ok(decoys)
}

/// Clusters a set of decoy models by Cα-RMSD.
///
/// Uses **leader clustering**: walk the decoys; each decoy joins the
/// first existing cluster whose leader is within `rmsd_threshold`
/// ångström, or founds a new cluster. Each cluster's centre is then
/// re-chosen as the member with the smallest mean RMSD to the others.
/// Clusters are returned largest-first.
///
/// `scores` must be one knowledge-score total per model (same order
/// and length as `models`).
///
/// # Errors
/// [`StructPredictError::Invalid`] for an empty model list, a
/// `scores` length mismatch, or a non-positive threshold.
pub fn cluster_decoys(
    models: &[ProteinModel],
    scores: &[f64],
    rmsd_threshold: f64,
) -> Result<Vec<DecoyCluster>> {
    if models.is_empty() {
        return Err(StructPredictError::invalid(
            "models",
            "no decoys to cluster",
        ));
    }
    if scores.len() != models.len() {
        return Err(StructPredictError::invalid(
            "scores",
            "length disagrees with the model list",
        ));
    }
    if !(rmsd_threshold.is_finite() && rmsd_threshold > 0.0) {
        return Err(StructPredictError::invalid(
            "rmsd_threshold",
            "must be finite and positive",
        ));
    }

    // Leader clustering.
    let mut leaders: Vec<usize> = Vec::new();
    let mut clusters: Vec<Vec<usize>> = Vec::new();
    for (i, model) in models.iter().enumerate() {
        let mut joined = false;
        for (c, &leader) in leaders.iter().enumerate() {
            let rmsd = ca_rmsd_superposed(model, &models[leader]).unwrap_or(f64::INFINITY);
            if rmsd <= rmsd_threshold {
                clusters[c].push(i);
                joined = true;
                break;
            }
        }
        if !joined {
            leaders.push(i);
            clusters.push(vec![i]);
        }
    }

    // Re-choose each cluster centre: the member with the smallest mean
    // RMSD to the other members.
    let mut out = Vec::with_capacity(clusters.len());
    for members in clusters {
        let center = if members.len() == 1 {
            members[0]
        } else {
            let mut best = members[0];
            let mut best_mean = f64::INFINITY;
            for &a in &members {
                let mut sum = 0.0;
                for &b in &members {
                    if a != b {
                        sum += ca_rmsd_superposed(&models[a], &models[b]).unwrap_or(f64::INFINITY);
                    }
                }
                let mean = sum / (members.len() - 1) as f64;
                if mean < best_mean {
                    best_mean = mean;
                    best = a;
                }
            }
            best
        };
        out.push(DecoyCluster {
            center,
            center_score: scores[center],
            members,
        });
    }
    // Largest cluster first.
    out.sort_by_key(|c| std::cmp::Reverse(c.size()));
    Ok(out)
}

/// How a final model is chosen from clustered decoys.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum SelectionPolicy {
    /// Pick the centre of the largest cluster — the most
    /// reproducibly-sampled fold (Rosetta's default).
    LargestCluster,
    /// Pick the single lowest-energy decoy, regardless of clustering.
    LowestEnergy,
}

/// Selects a final model from a clustered decoy set.
///
/// Returns the chosen decoy's index. With
/// [`SelectionPolicy::LargestCluster`] this is the largest cluster's
/// centre; with [`SelectionPolicy::LowestEnergy`] the global
/// minimum-score decoy.
///
/// # Errors
/// [`StructPredictError::Invalid`] for an empty cluster list or empty
/// scores.
pub fn select_model(
    clusters: &[DecoyCluster],
    scores: &[f64],
    policy: SelectionPolicy,
) -> Result<usize> {
    if clusters.is_empty() {
        return Err(StructPredictError::invalid("clusters", "empty"));
    }
    match policy {
        SelectionPolicy::LargestCluster => {
            // Clusters are already largest-first.
            Ok(clusters[0].center)
        }
        SelectionPolicy::LowestEnergy => {
            if scores.is_empty() {
                return Err(StructPredictError::invalid("scores", "empty"));
            }
            let (idx, _) = scores
                .iter()
                .enumerate()
                .min_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
                .expect("non-empty");
            Ok(idx)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Point3;

    /// A linear Cα model offset by `dx` from the origin.
    fn shifted_model(n: usize, dx: f64) -> ProteinModel {
        let seq = "A".repeat(n);
        let mut m = ProteinModel::from_sequence(&seq).expect("model");
        for (i, r) in m.residues.iter_mut().enumerate() {
            r.ca = Some(Point3::new(i as f64 * 3.8 + dx, 0.0, 0.0));
        }
        m
    }

    /// A bent Cα model — structurally distinct from the linear one.
    fn bent_model(n: usize) -> ProteinModel {
        let seq = "A".repeat(n);
        let mut m = ProteinModel::from_sequence(&seq).expect("model");
        for (i, r) in m.residues.iter_mut().enumerate() {
            let x = i as f64 * 3.8;
            r.ca = Some(Point3::new(x, (x * 0.3).sin() * 10.0, 0.0));
        }
        m
    }

    #[test]
    fn translation_does_not_split_a_cluster() {
        // Three translated copies of the same shape — RMSD-superposed
        // they are identical, so they form one cluster.
        let models = vec![
            shifted_model(8, 0.0),
            shifted_model(8, 5.0),
            shifted_model(8, -3.0),
        ];
        let scores = vec![1.0, 2.0, 3.0];
        let clusters = cluster_decoys(&models, &scores, 1.0).expect("cluster");
        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters[0].size(), 3);
    }

    #[test]
    fn distinct_shapes_form_separate_clusters() {
        let models = vec![shifted_model(8, 0.0), bent_model(8), shifted_model(8, 2.0)];
        let scores = vec![1.0, 2.0, 1.5];
        let clusters = cluster_decoys(&models, &scores, 1.0).expect("cluster");
        // Two linear + one bent → two clusters; the linear one is bigger.
        assert_eq!(clusters.len(), 2);
        assert_eq!(clusters[0].size(), 2);
    }

    #[test]
    fn selection_policies_differ() {
        let models = vec![shifted_model(8, 0.0), bent_model(8), shifted_model(8, 2.0)];
        let scores = vec![5.0, -9.0, 4.0]; // the bent decoy has the best score
        let clusters = cluster_decoys(&models, &scores, 1.0).expect("cluster");
        let consensus =
            select_model(&clusters, &scores, SelectionPolicy::LargestCluster).expect("consensus");
        let lowest =
            select_model(&clusters, &scores, SelectionPolicy::LowestEnergy).expect("lowest");
        // Largest cluster is the linear pair; lowest energy is the bent one.
        assert_ne!(consensus, lowest);
        assert_eq!(lowest, 1);
    }

    #[test]
    fn empty_inputs_rejected() {
        assert!(cluster_decoys(&[], &[], 1.0).is_err());
        assert!(select_model(&[], &[], SelectionPolicy::LowestEnergy).is_err());
    }
}
