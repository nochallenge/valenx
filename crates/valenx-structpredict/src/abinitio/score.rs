//! **Feature 9 — knowledge-based scoring (statistical potentials).**
//!
//! Fragment assembly needs an energy function to tell good
//! conformations from bad ones. The classical choice — Rosetta's
//! `centroid` score, Sippl's potentials of mean force — is a
//! **knowledge-based statistical potential**: rather than physics
//! from first principles, the energy is derived from *statistics* of
//! known structures. If a particular Cα–Cα distance between two
//! residue types is observed often in the PDB, that geometry is
//! low-energy; if rare, high-energy.
//!
//! This module implements a centroid-resolution knowledge-based
//! score with the standard term set:
//!
//! - **Distance-dependent pairwise potential** — a potential of mean
//!   force `E = −kT·ln[ p_obs(d) / p_ref(d) ]`, binned by Cα–Cα
//!   distance and hydrophobic/polar residue class.
//! - **Contact term** — rewards a well-packed core (hydrophobic
//!   residues making contacts) and penalises buried polar groups.
//! - **Torsion / compactness term** — rewards a compact, sphere-like
//!   radius of gyration, the hallmark of a folded globular protein,
//!   and penalises steric clashes.
//!
//! **Honest note.** The potential here is *hand-built* from the
//! hydrophobic-collapse and excluded-volume principles, parameterised
//! to behave like a real statistical potential (attractive at contact
//! range, a hard repulsive wall, hydrophobic burial favoured). A
//! production score is *fitted* to PDB statistics over thousands of
//! structures; this v1 is the same functional form with principled
//! rather than fitted parameters. It ranks conformations sensibly —
//! a clashing or extended decoy scores worse than a compact one — but
//! it is not the Rosetta energy function.

use serde::{Deserialize, Serialize};

use crate::aa::hydropathy;
use crate::error::{Result, StructPredictError};
use crate::model::ProteinModel;

/// Relative weights of the knowledge-score terms.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ScoreWeights {
    /// Weight of the distance-dependent pairwise potential.
    pub pairwise: f64,
    /// Weight of the hydrophobic-contact / burial term.
    pub contact: f64,
    /// Weight of the compactness (radius-of-gyration) term.
    pub compactness: f64,
    /// Weight of the steric-clash penalty.
    pub clash: f64,
}

impl Default for ScoreWeights {
    fn default() -> Self {
        ScoreWeights {
            pairwise: 1.0,
            contact: 1.0,
            compactness: 0.5,
            clash: 4.0,
        }
    }
}

/// A decomposed knowledge-based score. Lower `total` is better.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct KnowledgeScore {
    /// Distance-dependent pairwise potential energy.
    pub pairwise: f64,
    /// Hydrophobic-contact / burial energy.
    pub contact: f64,
    /// Compactness energy (negative is compact, good).
    pub compactness: f64,
    /// Steric-clash penalty (≥ 0).
    pub clash: f64,
    /// The weighted sum.
    pub total: f64,
}

/// The distance-dependent pairwise potential for one residue pair,
/// in arbitrary energy units.
///
/// `hydrophobic` is whether *both* residues are hydrophobic. The
/// potential is a smooth well: a hard repulsive wall below ~4 Å (the
/// excluded-volume term), an attractive minimum near the
/// contact distance (deeper for a hydrophobic pair — hydrophobic
/// collapse), and zero beyond the interaction cutoff.
fn pairwise_potential(d: f64, hydrophobic: bool) -> f64 {
    const HARD: f64 = 3.5; // hard-clash radius
    const CONTACT: f64 = 5.5; // attractive-minimum distance
    const CUTOFF: f64 = 10.0; // interaction range
    if d < HARD {
        // Steep repulsion — divergent excluded volume.
        let r = HARD / d.max(0.5);
        return 3.0 * (r * r - 1.0);
    }
    if d > CUTOFF {
        return 0.0;
    }
    // A Lennard-Jones-like well centred on CONTACT.
    let depth = if hydrophobic { 1.4 } else { 0.5 };
    let x = (d - CONTACT) / (CUTOFF - CONTACT);
    // Smooth well: minimum -depth at d=CONTACT, rising to 0 at CUTOFF
    // and at HARD.
    if d <= CONTACT {
        let y = (CONTACT - d) / (CONTACT - HARD);
        depth * (y * y - 1.0)
    } else {
        depth * (x * x - 1.0)
    }
}

/// Scores a model with the knowledge-based potential.
///
/// The model must have a Cα trace; residues without a Cα are skipped.
/// The score is decomposed into its four terms and summed with
/// `weights`.
///
/// # Errors
/// [`StructPredictError::Invalid`] if the model has fewer than two
/// residues with Cα atoms.
pub fn score_model(model: &ProteinModel, weights: ScoreWeights) -> Result<KnowledgeScore> {
    // Collect (Cα, hydropathy) for residues with coordinates.
    let mut points = Vec::new();
    for r in &model.residues {
        if let Some(ca) = r.ca {
            points.push((ca, hydropathy(r.aa)));
        }
    }
    if points.len() < 2 {
        return Err(StructPredictError::invalid(
            "model",
            "needs at least two residues with Cα coordinates to score",
        ));
    }
    let n = points.len();

    let mut pairwise = 0.0;
    let mut clash = 0.0;
    // Per-residue contact count, used for the burial term.
    let mut contacts = vec![0u32; n];
    for i in 0..n {
        for j in (i + 1)..n {
            // Skip bonded neighbours — sequence-local geometry is the
            // fragment library's job, not the pair potential's.
            if j - i < 2 {
                continue;
            }
            let d = (points[i].0 - points[j].0).norm();
            let both_phobic = points[i].1 > 0.0 && points[j].1 > 0.0;
            pairwise += pairwise_potential(d, both_phobic);
            if d < 8.0 {
                contacts[i] += 1;
                contacts[j] += 1;
            }
            if d < 3.0 {
                // Severe overlap — a clash on top of the pair term.
                clash += (3.0 - d) * (3.0 - d);
            }
        }
    }

    // Burial term: a hydrophobic residue *wants* many contacts (a
    // buried core); a polar residue prefers the surface (few
    // contacts). Reward/penalise accordingly.
    let mut contact = 0.0;
    let median_contacts = {
        let mut c = contacts.clone();
        c.sort_unstable();
        c[c.len() / 2] as f64
    };
    for i in 0..n {
        let buried = contacts[i] as f64 >= median_contacts;
        let phobic = points[i].1 > 0.0;
        contact += match (phobic, buried) {
            (true, true) => -1.0,   // hydrophobic & buried — good
            (true, false) => 0.8,   // hydrophobic & exposed — bad
            (false, true) => 0.5,   // polar & buried — slightly bad
            (false, false) => -0.3, // polar & exposed — good
        };
    }

    // Compactness: compare the actual radius of gyration to the
    // expected Rg of a folded globular protein of this length
    // (Rg ≈ 2.2·N^0.38 Å). A compact model scores negative.
    let centroid = {
        let mut c = nalgebra::Vector3::zeros();
        for (p, _) in &points {
            c += p.coords;
        }
        c / n as f64
    };
    let rg = {
        let mut s = 0.0;
        for (p, _) in &points {
            s += (p.coords - centroid).norm_squared();
        }
        (s / n as f64).sqrt()
    };
    let expected_rg = 2.2 * (n as f64).powf(0.38);
    // Negative (good) when rg ≈ expected, rising when too extended.
    let compactness = ((rg - expected_rg) / expected_rg).powi(2) - 0.3;

    let total = weights.pairwise * pairwise
        + weights.contact * contact
        + weights.compactness * compactness
        + weights.clash * clash;

    Ok(KnowledgeScore {
        pairwise,
        contact,
        compactness,
        clash,
        total,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Point3;

    /// A compact globular blob of hydrophobic residues.
    fn compact_model(n: usize) -> ProteinModel {
        let seq = "L".repeat(n);
        let mut m = ProteinModel::from_sequence(&seq).expect("model");
        // Place residues on a tight 3-D lattice.
        let side = (n as f64).cbrt().ceil() as usize;
        for (idx, r) in m.residues.iter_mut().enumerate() {
            let x = (idx % side) as f64 * 3.8;
            let y = ((idx / side) % side) as f64 * 3.8;
            let z = (idx / (side * side)) as f64 * 3.8;
            r.ca = Some(Point3::new(x, y, z));
        }
        m
    }

    /// A fully extended chain.
    fn extended_model(n: usize) -> ProteinModel {
        let seq = "L".repeat(n);
        let mut m = ProteinModel::from_sequence(&seq).expect("model");
        for (idx, r) in m.residues.iter_mut().enumerate() {
            r.ca = Some(Point3::new(idx as f64 * 3.8, 0.0, 0.0));
        }
        m
    }

    #[test]
    fn compact_beats_extended() {
        let w = ScoreWeights::default();
        let compact = score_model(&compact_model(27), w).expect("score");
        let extended = score_model(&extended_model(27), w).expect("score");
        assert!(
            compact.total < extended.total,
            "compact {} should beat extended {}",
            compact.total,
            extended.total
        );
    }

    #[test]
    fn clash_is_penalised() {
        let mut m = ProteinModel::from_sequence("LLLLL").expect("model");
        // Stack all residues nearly on top of each other.
        for r in m.residues.iter_mut() {
            r.ca = Some(Point3::new(0.0, 0.0, 0.0));
        }
        let s = score_model(&m, ScoreWeights::default()).expect("score");
        assert!(s.clash > 0.0, "overlapping residues clash");
    }

    #[test]
    fn pairwise_potential_has_hard_wall() {
        // Below the hard radius the potential is steeply positive.
        assert!(pairwise_potential(2.0, false) > 0.0);
        // At contact range it is attractive (negative).
        assert!(pairwise_potential(5.5, true) < 0.0);
        // Beyond cutoff it is zero.
        assert_eq!(pairwise_potential(20.0, true), 0.0);
    }

    #[test]
    fn too_few_residues_rejected() {
        let m = ProteinModel::from_sequence("A").expect("model");
        assert!(score_model(&m, ScoreWeights::default()).is_err());
    }
}
