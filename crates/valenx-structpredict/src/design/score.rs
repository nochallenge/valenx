//! **Feature 20 — sequence-design scoring.**
//!
//! The design search needs an energy function that says how well a
//! *candidate sequence* fits a *fixed backbone*. The classical
//! design score combines several physics / statistics terms; this
//! module implements the standard set at centroid resolution:
//!
//! - **Burial / hydropathy term** — the dominant signal in protein
//!   design. A buried position (one whose Cα makes many contacts)
//!   should carry a hydrophobic residue; an exposed position a polar
//!   one. This is the hydrophobic-core / polar-surface principle that
//!   makes a globular protein fold.
//! - **Packing term** — the chosen sidechain must fit the available
//!   cavity. A big sidechain in a tight pocket clashes; a tiny one in
//!   a large pocket leaves a void. Both are penalised.
//! - **Secondary-structure-propensity term** — a residue in a helix
//!   should be a good helix former (Ala, Glu, Leu …); a residue in a
//!   strand a good strand former. Uses the same Chou-Fasman
//!   propensities as the SS predictor.
//! - **Reference / amino-acid-composition term** — a small per-amino
//!   acid reference energy that keeps the designed sequence's
//!   composition realistic (without it the search piles up whatever
//!   single residue scores best).
//!
//! Lower energy is a better-designed sequence. Like the ab-initio
//! score, the parameters are principled rather than PDB-fitted — a
//! real, honest classical design potential, not the Rosetta energy
//! function.

use serde::{Deserialize, Serialize};

use crate::aa::{hydropathy, sidechain_volume};
use crate::abinitio::ss::SecondaryStructure;
use crate::error::{Result, StructPredictError};
use crate::model::ProteinModel;

/// Relative weights of the design-score terms.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DesignScoreWeights {
    /// Weight of the burial / hydropathy term.
    pub burial: f64,
    /// Weight of the cavity-packing term.
    pub packing: f64,
    /// Weight of the secondary-structure-propensity term.
    pub ss_propensity: f64,
    /// Weight of the amino-acid reference term.
    pub reference: f64,
}

impl Default for DesignScoreWeights {
    fn default() -> Self {
        DesignScoreWeights {
            burial: 1.0,
            packing: 0.6,
            ss_propensity: 0.8,
            reference: 0.3,
        }
    }
}

/// A decomposed design score. Lower `total` is a better sequence.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DesignScore {
    /// Burial / hydropathy energy.
    pub burial: f64,
    /// Cavity-packing energy.
    pub packing: f64,
    /// Secondary-structure-propensity energy.
    pub ss_propensity: f64,
    /// Amino-acid reference energy.
    pub reference: f64,
    /// The weighted sum.
    pub total: f64,
}

/// The Chou-Fasman helix / strand propensity of an amino acid,
/// reused from the SS predictor's statistics.
fn ss_propensity(aa: char, state: SecondaryStructure) -> f64 {
    // Re-derive the two propensities (helix, strand) inline.
    let (ph, pe): (f64, f64) = match aa {
        'E' => (1.51, 0.37),
        'M' => (1.45, 1.05),
        'A' => (1.42, 0.83),
        'L' => (1.21, 1.30),
        'K' => (1.16, 0.74),
        'F' => (1.13, 1.38),
        'Q' => (1.11, 1.10),
        'W' => (1.08, 1.37),
        'I' => (1.08, 1.60),
        'V' => (1.06, 1.70),
        'D' => (1.01, 0.54),
        'H' => (1.00, 0.87),
        'R' => (0.98, 0.93),
        'T' => (0.83, 1.19),
        'S' => (0.77, 0.75),
        'C' => (0.70, 1.19),
        'Y' => (0.69, 1.47),
        'N' => (0.67, 0.89),
        'P' => (0.57, 0.55),
        'G' => (0.57, 0.75),
        _ => (1.0, 1.0),
    };
    match state {
        // Energy = −ln(propensity): a good former is low-energy.
        SecondaryStructure::Helix => -(ph.max(1e-3)).ln(),
        SecondaryStructure::Strand => -(pe.max(1e-3)).ln(),
        SecondaryStructure::Coil => 0.0,
    }
}

/// A small per-amino-acid reference energy. Positive for the rarest
/// residues (Trp, Cys) so the design search does not over-use them;
/// near zero for the common ones.
fn reference_energy(aa: char) -> f64 {
    match aa {
        'W' => 1.2,
        'C' => 0.9,
        'M' => 0.6,
        'H' => 0.5,
        'F' | 'Y' => 0.3,
        'A' | 'L' | 'G' | 'V' | 'E' | 'S' | 'K' => -0.1,
        _ => 0.0,
    }
}

/// Per-residue burial — the number of Cα–Cα contacts within 10 Å.
/// A buried residue has many contacts; an exposed one few.
pub fn residue_burial(model: &ProteinModel) -> Vec<usize> {
    let trace: Vec<_> = model
        .residues
        .iter()
        .map(|r| r.ca)
        .collect::<Vec<_>>();
    let n = trace.len();
    let mut burial = vec![0usize; n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            if let (Some(ci), Some(cj)) = (trace[i], trace[j]) {
                if (ci - cj).norm() < 10.0 {
                    burial[i] += 1;
                }
            }
        }
    }
    burial
}

/// Scores a candidate sequence on a fixed backbone.
///
/// `model` supplies the backbone (its Cα trace and residue count);
/// `sequence` is the candidate amino-acid sequence (same length);
/// `ss` is the per-residue secondary structure (same length, or
/// empty to skip the SS term).
///
/// # Errors
/// [`StructPredictError::Invalid`] for a length mismatch between the
/// model, the sequence and (if non-empty) `ss`, or a model with
/// fewer than 2 Cα atoms.
pub fn design_score(
    model: &ProteinModel,
    sequence: &str,
    ss: &[SecondaryStructure],
    weights: DesignScoreWeights,
) -> Result<DesignScore> {
    let chars: Vec<char> = sequence.trim().to_ascii_uppercase().chars().collect();
    if chars.len() != model.residues.len() {
        return Err(StructPredictError::invalid(
            "sequence",
            format!(
                "length {} disagrees with the backbone's {} residues",
                chars.len(),
                model.residues.len()
            ),
        ));
    }
    if !ss.is_empty() && ss.len() != chars.len() {
        return Err(StructPredictError::invalid(
            "ss",
            "length disagrees with the sequence",
        ));
    }
    let burial_counts = residue_burial(model);
    let with_ca = burial_counts.len();
    if model.ca_trace().len() < 2 {
        return Err(StructPredictError::invalid(
            "model",
            "need at least 2 Cα atoms to score a design",
        ));
    }

    // Median burial splits buried from exposed.
    let median = {
        let mut c = burial_counts.clone();
        c.sort_unstable();
        c[c.len() / 2]
    };

    let mut burial = 0.0;
    let mut packing = 0.0;
    let mut ss_prop = 0.0;
    let mut reference = 0.0;
    for i in 0..with_ca {
        let aa = chars[i];
        let buried = burial_counts[i] >= median;
        let phobic = hydropathy(aa) > 0.0;
        // Burial term: hydrophobic-buried & polar-exposed are good.
        burial += match (phobic, buried) {
            (true, true) => -1.0,
            (true, false) => 1.0,
            (false, true) => 0.6,
            (false, false) => -0.4,
        };
        // Packing: a buried position wants a residue whose volume
        // matches the local cavity. Approximate the cavity size from
        // the contact count (more contacts → tighter), and penalise a
        // mismatch.
        let cavity = (12 - burial_counts[i].min(12)) as f64 * 12.0; // Å³-ish
        let vol = sidechain_volume(aa);
        let mismatch = (vol - cavity).abs() / 100.0;
        packing += mismatch * mismatch;
        // SS propensity.
        if let Some(&state) = ss.get(i) {
            ss_prop += ss_propensity(aa, state);
        }
        reference += reference_energy(aa);
    }

    let total = weights.burial * burial
        + weights.packing * packing
        + weights.ss_propensity * ss_prop
        + weights.reference * reference;

    Ok(DesignScore {
        burial,
        packing,
        ss_propensity: ss_prop,
        reference,
        total,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Point3;

    /// A compact globular blob — interior residues are buried,
    /// surface residues exposed.
    fn blob(n: usize) -> ProteinModel {
        let seq = "A".repeat(n);
        let mut m = ProteinModel::from_sequence(&seq).expect("model");
        let side = (n as f64).cbrt().ceil() as usize;
        for (idx, r) in m.residues.iter_mut().enumerate() {
            let x = (idx % side) as f64 * 3.8;
            let y = ((idx / side) % side) as f64 * 3.8;
            let z = (idx / (side * side)) as f64 * 3.8;
            r.ca = Some(Point3::new(x, y, z));
        }
        m
    }

    #[test]
    fn hydrophobic_core_beats_hydrophilic_core() {
        let m = blob(27);
        let n = m.residues.len();
        // All-hydrophobic vs all-hydrophilic sequences.
        let hydrophobic = "L".repeat(n);
        let hydrophilic = "K".repeat(n);
        let w = DesignScoreWeights::default();
        let phobic = design_score(&m, &hydrophobic, &[], w).expect("score");
        let philic = design_score(&m, &hydrophilic, &[], w).expect("score");
        // A blob is mostly buried → hydrophobic burial is favoured.
        assert!(
            phobic.burial < philic.burial,
            "burial: phobic {} vs philic {}",
            phobic.burial,
            philic.burial
        );
    }

    #[test]
    fn helix_propensity_rewards_helix_formers() {
        let m = blob(9);
        let n = m.residues.len();
        let ss = vec![SecondaryStructure::Helix; n];
        let w = DesignScoreWeights::default();
        // Glu/Ala are strong helix formers; Gly/Pro are not.
        let good = "E".repeat(n);
        let bad = "P".repeat(n);
        let sg = design_score(&m, &good, &ss, w).expect("score");
        let sb = design_score(&m, &bad, &ss, w).expect("score");
        assert!(
            sg.ss_propensity < sb.ss_propensity,
            "ss: good {} vs bad {}",
            sg.ss_propensity,
            sb.ss_propensity
        );
    }

    #[test]
    fn length_mismatch_rejected() {
        let m = blob(9);
        assert!(design_score(&m, "ACDEF", &[], DesignScoreWeights::default()).is_err());
    }
}
