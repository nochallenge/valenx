//! **Feature 16 — model quality assessment.**
//!
//! A predicted model needs an honest, *reference-free* quality score
//! — a number you can compute without knowing the true structure,
//! the way MolProbity scores a crystal structure. The classical
//! quality signals are:
//!
//! - **Clash score** — the count of severe inter-atomic overlaps per
//!   1000 atoms. A good structure has almost none; clashes mean a
//!   physically impossible model.
//! - **Ramachandran statistics** — the fraction of residues with
//!   backbone φ/ψ in allowed regions, and the outlier fraction. Real
//!   structures have > 98 % of residues in allowed regions.
//! - **Packing** — globular proteins have a compact, well-packed
//!   core. A loosely-packed or hollow model (too large a radius of
//!   gyration for its length, too few contacts) is suspect.
//!
//! This module computes all three from a [`crate::model::ProteinModel`]
//! and blends them into a single [`QualityReport`] with an overall
//! `0-100` score. It is a genuine classical model-quality estimator —
//! the same signals MolProbity / Verify3D use — though, like every
//! reference-free score, it can be fooled and is not a substitute for
//! experimental validation.

use serde::{Deserialize, Serialize};

use crate::error::{Result, StructPredictError};
use crate::model::ProteinModel;
use crate::refine::ramachandran::{is_allowed, model_phi_psi};

/// A classical reference-free model-quality report.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct QualityReport {
    /// Severe clashes per 1000 atoms (Cα-level). Lower is better; a
    /// good model is near 0.
    pub clash_score: f64,
    /// Fraction of (non-terminal) residues in allowed Ramachandran
    /// regions, `[0, 1]`. Higher is better.
    pub ramachandran_favored: f64,
    /// Fraction of (non-terminal) residues that are Ramachandran
    /// outliers, `[0, 1]`. Lower is better.
    pub ramachandran_outliers: f64,
    /// Packing score in `[0, 1]` — how close the model's radius of
    /// gyration and contact density are to a well-folded globular
    /// protein. Higher is better.
    pub packing: f64,
    /// Overall quality in `[0, 100]` — a blend of the above.
    pub overall: f64,
}

impl QualityReport {
    /// A coarse verdict string for display.
    pub fn verdict(&self) -> &'static str {
        if self.overall >= 75.0 {
            "good"
        } else if self.overall >= 50.0 {
            "acceptable"
        } else {
            "poor"
        }
    }
}

/// Assesses the quality of a model.
///
/// Computes the clash score, Ramachandran statistics and packing
/// score from the model's Cα trace (and full backbone, where present)
/// and blends them into an overall `0-100` score.
///
/// # Errors
/// [`StructPredictError::Invalid`] if the model has fewer than 3
/// residues with Cα atoms.
pub fn assess_quality(model: &ProteinModel) -> Result<QualityReport> {
    let trace = model.ca_trace();
    let n = trace.len();
    if n < 3 {
        return Err(StructPredictError::invalid(
            "model",
            "need at least 3 residues with Cα atoms to assess",
        ));
    }

    // --- Clash score ---------------------------------------------------
    // Count Cα pairs (non-adjacent) closer than a severe-overlap
    // radius. Adjacent Cα atoms are legitimately ~3.8 Å apart.
    const SEVERE: f64 = 3.0;
    let mut clashes = 0usize;
    for i in 0..n {
        for j in (i + 2)..n {
            if (trace[i] - trace[j]).norm() < SEVERE {
                clashes += 1;
            }
        }
    }
    let clash_score = clashes as f64 / n as f64 * 1000.0;

    // --- Ramachandran --------------------------------------------------
    let (favored, outliers) = if model.is_complete() {
        let pp = model_phi_psi(model);
        let mut fav = 0usize;
        let mut out = 0usize;
        let mut counted = 0usize;
        for (i, &(phi, psi)) in pp.iter().enumerate() {
            if i == 0 || i + 1 == pp.len() {
                continue; // terminal residues have partial φ/ψ
            }
            counted += 1;
            if is_allowed(phi, psi) {
                fav += 1;
            } else {
                out += 1;
            }
        }
        if counted > 0 {
            (fav as f64 / counted as f64, out as f64 / counted as f64)
        } else {
            (1.0, 0.0)
        }
    } else {
        // No full backbone — Ramachandran cannot be assessed; report
        // a neutral 1.0 favored so it does not unfairly sink the
        // score of a Cα-only model.
        (1.0, 0.0)
    };

    // --- Packing -------------------------------------------------------
    // Radius of gyration vs the expected Rg of a folded globular
    // protein, plus the mean contact count.
    let centroid = {
        let mut c = nalgebra::Vector3::zeros();
        for p in &trace {
            c += p.coords;
        }
        c / n as f64
    };
    let rg = {
        let mut s = 0.0;
        for p in &trace {
            s += (p.coords - centroid).norm_squared();
        }
        (s / n as f64).sqrt()
    };
    let expected_rg = 2.2 * (n as f64).powf(0.38);
    // Rg score: 1.0 when rg ≈ expected, decaying away from it.
    let rg_ratio = rg / expected_rg;
    let rg_score = (1.0 - (rg_ratio - 1.0).abs()).clamp(0.0, 1.0);
    // Contact density: count Cα–Cα contacts within 8 Å.
    let mut contacts = 0usize;
    for i in 0..n {
        for j in (i + 2)..n {
            if (trace[i] - trace[j]).norm() < 8.0 {
                contacts += 1;
            }
        }
    }
    // A well-folded protein has ~3-6 such contacts per residue.
    let contact_density = contacts as f64 / n as f64;
    let contact_score = (contact_density / 4.0).clamp(0.0, 1.0);
    let packing = 0.5 * rg_score + 0.5 * contact_score;

    // --- Overall -------------------------------------------------------
    // A clash-score → [0,1] term: 1 at zero clashes, decaying.
    let clash_term = (1.0 - clash_score / 60.0).clamp(0.0, 1.0);
    let overall = 100.0
        * (0.35 * clash_term + 0.35 * favored + 0.30 * packing - 0.0 * outliers).clamp(0.0, 1.0);

    Ok(QualityReport {
        clash_score,
        ramachandran_favored: favored,
        ramachandran_outliers: outliers,
        packing,
        overall,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::build_backbone_from_torsions;
    use nalgebra::Point3;

    #[test]
    fn clean_helix_scores_well() {
        let mut m = ProteinModel::from_sequence("AAAAAAAAAAAAAAAAAAAA").expect("model");
        build_backbone_from_torsions(&mut m, &vec![(-63.0, -42.0); 20]).expect("build");
        let q = assess_quality(&m).expect("assess");
        assert!(q.clash_score < 5.0, "clash score {}", q.clash_score);
        assert!(
            q.ramachandran_favored > 0.9,
            "favored {}",
            q.ramachandran_favored
        );
    }

    #[test]
    fn clashing_model_scores_poorly() {
        let mut m = ProteinModel::from_sequence("AAAAAAAA").expect("model");
        // Pile every residue near the origin → huge clash score.
        for (i, r) in m.residues.iter_mut().enumerate() {
            r.ca = Some(Point3::new((i as f64) * 0.4, 0.0, 0.0));
        }
        let q = assess_quality(&m).expect("assess");
        assert!(q.clash_score > 50.0, "clash score {}", q.clash_score);
        assert!(q.overall < 50.0, "overall {}", q.overall);
        assert_eq!(q.verdict(), "poor");
    }

    #[test]
    fn ramachandran_outliers_are_counted() {
        let mut m = ProteinModel::from_sequence("AAAAAAAAAA").expect("model");
        let mut torsions = vec![(-63.0, -42.0); 10];
        torsions[3] = (20.0, 100.0); // disallowed
        torsions[6] = (40.0, 110.0); // disallowed
        build_backbone_from_torsions(&mut m, &torsions).expect("build");
        let q = assess_quality(&m).expect("assess");
        assert!(q.ramachandran_outliers > 0.0, "outliers detected");
    }

    #[test]
    fn too_short_rejected() {
        let m = ProteinModel::from_sequence("A").expect("model");
        assert!(assess_quality(&m).is_err());
    }
}
