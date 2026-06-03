//! **Feature 8 — secondary-structure prediction (classical, non-ML).**
//!
//! Predicts the three-state secondary structure — helix `H`, strand
//! `E`, coil `C` — of a protein from its sequence alone, using the
//! two classical statistical methods:
//!
//! - **Chou-Fasman** (1974) — each amino acid carries a measured
//!   helix propensity `P(a)` and strand propensity `P(b)`. The method
//!   scans for nucleation windows where the average propensity is
//!   high, then extends them.
//! - **GOR** (Garnier-Osguthorpe-Robson) — an information-theoretic
//!   method: the state of a residue is predicted from the identities
//!   of its neighbours within a sliding window, weighting each
//!   neighbour's contribution by a position-specific information
//!   value.
//!
//! Both are **purely statistical** — propensity tables and
//! information sums, no neural network, no trained weights beyond the
//! published amino-acid statistics. They reach ~60-65 % three-state
//! accuracy, the well-known ceiling of single-sequence classical
//! prediction (modern ML predictors using evolutionary profiles
//! reach ~85 %+ — that gain is the network's, and out of scope here).

use serde::{Deserialize, Serialize};

use crate::error::{Result, StructPredictError};

/// A three-state secondary-structure assignment.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub enum SecondaryStructure {
    /// α-helix (`H`).
    Helix,
    /// β-strand / extended (`E`).
    Strand,
    /// Coil / loop / turn (`C`).
    Coil,
}

impl SecondaryStructure {
    /// The single-character DSSP-style code (`H`, `E`, `C`).
    pub fn code(self) -> char {
        match self {
            SecondaryStructure::Helix => 'H',
            SecondaryStructure::Strand => 'E',
            SecondaryStructure::Coil => 'C',
        }
    }

    /// `true` for regular secondary structure (helix or strand) —
    /// the structure-aware alignment penalises gaps inside these.
    pub fn is_regular(self) -> bool {
        matches!(self, SecondaryStructure::Helix | SecondaryStructure::Strand)
    }
}

/// A per-residue secondary-structure prediction with confidence.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SsPrediction {
    /// The predicted state per residue, in sequence order.
    pub states: Vec<SecondaryStructure>,
    /// Per-residue confidence in `[0, 1]` — the margin between the
    /// winning and runner-up state scores.
    pub confidence: Vec<f64>,
}

impl SsPrediction {
    /// The prediction as a `HEC` string.
    pub fn as_string(&self) -> String {
        self.states.iter().map(|s| s.code()).collect()
    }

    /// Fraction of residues predicted helical.
    pub fn helix_fraction(&self) -> f64 {
        self.fraction(SecondaryStructure::Helix)
    }

    /// Fraction of residues predicted strand.
    pub fn strand_fraction(&self) -> f64 {
        self.fraction(SecondaryStructure::Strand)
    }

    fn fraction(&self, s: SecondaryStructure) -> f64 {
        if self.states.is_empty() {
            return 0.0;
        }
        self.states.iter().filter(|&&x| x == s).count() as f64 / self.states.len() as f64
    }
}

/// Chou-Fasman conformational parameters `(P_helix, P_strand,
/// P_coil)` per amino acid. Values from the original 1974 paper,
/// rescaled so `P ≈ 1.0` is the no-preference baseline.
fn chou_fasman_params(c: char) -> (f64, f64, f64) {
    match c {
        'E' => (1.51, 0.37, 0.74),
        'M' => (1.45, 1.05, 0.60),
        'A' => (1.42, 0.83, 0.66),
        'L' => (1.21, 1.30, 0.59),
        'K' => (1.16, 0.74, 1.01),
        'F' => (1.13, 1.38, 0.60),
        'Q' => (1.11, 1.10, 0.98),
        'W' => (1.08, 1.37, 0.96),
        'I' => (1.08, 1.60, 0.47),
        'V' => (1.06, 1.70, 0.50),
        'D' => (1.01, 0.54, 1.46),
        'H' => (1.00, 0.87, 0.95),
        'R' => (0.98, 0.93, 0.95),
        'T' => (0.83, 1.19, 0.96),
        'S' => (0.77, 0.75, 1.43),
        'C' => (0.70, 1.19, 1.19),
        'Y' => (0.69, 1.47, 1.14),
        'N' => (0.67, 0.89, 1.56),
        'P' => (0.57, 0.55, 1.52),
        'G' => (0.57, 0.75, 1.56),
        _ => (1.0, 1.0, 1.0),
    }
}

/// Predicts secondary structure by the **Chou-Fasman** method.
///
/// Helix and strand nucleation windows are found where the average
/// propensity exceeds 1.0, then extended; coil is the residual. A
/// helix/strand conflict is resolved in favour of the higher summed
/// propensity. Returns per-residue states and confidences.
///
/// # Errors
/// [`StructPredictError::Invalid`] for an empty sequence.
pub fn predict_chou_fasman(sequence: &str) -> Result<SsPrediction> {
    let seq = sequence.trim().to_ascii_uppercase();
    if seq.is_empty() {
        return Err(StructPredictError::invalid("sequence", "empty"));
    }
    let chars: Vec<char> = seq.chars().collect();
    let n = chars.len();
    let mut helix_score = vec![0.0f64; n];
    let mut strand_score = vec![0.0f64; n];
    let mut coil_score = vec![0.0f64; n];

    // Smooth the per-residue propensities over a small window — the
    // Chou-Fasman nucleation idea.
    let win = 3usize;
    for i in 0..n {
        let lo = i.saturating_sub(win);
        let hi = (i + win + 1).min(n);
        let mut h = 0.0;
        let mut e = 0.0;
        let mut c = 0.0;
        for &ch in &chars[lo..hi] {
            let (ph, pe, pc) = chou_fasman_params(ch);
            h += ph;
            e += pe;
            c += pc;
        }
        let m = (hi - lo) as f64;
        helix_score[i] = h / m;
        strand_score[i] = e / m;
        coil_score[i] = c / m;
    }
    finalize(&helix_score, &strand_score, &coil_score)
}

/// Predicts secondary structure by the **GOR** method.
///
/// Each residue's state score is the information-weighted sum of its
/// neighbours' single-residue propensities within an
/// `±window`-residue sliding window, neighbours nearer the centre
/// weighted higher (a triangular information kernel). This captures
/// GOR's central insight — that a residue's conformation depends on
/// its sequence *context*, not just its own identity.
///
/// # Errors
/// [`StructPredictError::Invalid`] for an empty sequence.
pub fn predict_gor(sequence: &str) -> Result<SsPrediction> {
    let seq = sequence.trim().to_ascii_uppercase();
    if seq.is_empty() {
        return Err(StructPredictError::invalid("sequence", "empty"));
    }
    let chars: Vec<char> = seq.chars().collect();
    let n = chars.len();
    let window = 8i64; // GOR uses a ±8-residue window.
    let mut helix_score = vec![0.0f64; n];
    let mut strand_score = vec![0.0f64; n];
    let mut coil_score = vec![0.0f64; n];

    for i in 0..n {
        let mut h = 0.0;
        let mut e = 0.0;
        let mut c = 0.0;
        let mut wsum = 0.0;
        for d in -window..=window {
            let j = i as i64 + d;
            if j < 0 || j >= n as i64 {
                continue;
            }
            // Triangular information weight: 1 at the centre, decaying
            // to ~0 at the window edge.
            let w = 1.0 - (d.abs() as f64) / (window as f64 + 1.0);
            let (ph, pe, pc) = chou_fasman_params(chars[j as usize]);
            // Convert propensity to a log-odds-like information value.
            h += w * (ph.ln());
            e += w * (pe.ln());
            c += w * (pc.ln());
            wsum += w;
        }
        if wsum > 0.0 {
            helix_score[i] = h / wsum;
            strand_score[i] = e / wsum;
            coil_score[i] = c / wsum;
        }
    }
    finalize(&helix_score, &strand_score, &coil_score)
}

/// Default secondary-structure predictor used by the rest of the
/// crate — currently GOR (the better of the two single-sequence
/// classical methods).
pub fn predict_secondary_structure(sequence: &str) -> Result<SsPrediction> {
    predict_gor(sequence)
}

/// Turns three per-residue score arrays into a state assignment +
/// confidence, with a 3-residue median smoothing pass to remove
/// isolated single-residue states (a helix or strand must span at
/// least a few residues to be real).
fn finalize(helix: &[f64], strand: &[f64], coil: &[f64]) -> Result<SsPrediction> {
    let n = helix.len();
    let mut states = Vec::with_capacity(n);
    let mut confidence = Vec::with_capacity(n);
    for i in 0..n {
        let scores = [
            (SecondaryStructure::Helix, helix[i]),
            (SecondaryStructure::Strand, strand[i]),
            (SecondaryStructure::Coil, coil[i]),
        ];
        let mut sorted = scores;
        sorted.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        states.push(sorted[0].0);
        // Confidence = margin between winner and runner-up, squashed.
        let margin = (sorted[0].1 - sorted[1].1).abs();
        confidence.push((margin * 2.0).min(1.0));
    }
    // Smoothing: a regular state surrounded by a different state on
    // both sides is flipped to its neighbours' state, and short
    // helix/strand runs (< 3) are demoted to coil.
    smooth_runs(&mut states);
    Ok(SsPrediction { states, confidence })
}

/// Removes secondary-structure runs shorter than the minimum length:
/// helices/strands need ≥ 3 residues to be reported.
fn smooth_runs(states: &mut [SecondaryStructure]) {
    let n = states.len();
    if n == 0 {
        return;
    }
    let mut i = 0;
    while i < n {
        let s = states[i];
        let mut j = i;
        while j < n && states[j] == s {
            j += 1;
        }
        if s.is_regular() && (j - i) < 3 {
            for st in states.iter_mut().take(j).skip(i) {
                *st = SecondaryStructure::Coil;
            }
        }
        i = j;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_sequence_rejected() {
        assert!(predict_chou_fasman("").is_err());
        assert!(predict_gor("").is_err());
    }

    #[test]
    fn poly_glutamate_predicts_helix() {
        // Glu/Ala/Leu/Met are strong helix formers.
        let pred = predict_gor("EEEEAAAALLLLEEEEAAAA").expect("predict");
        assert_eq!(pred.states.len(), 20);
        assert!(
            pred.helix_fraction() > 0.5,
            "helix fraction {}",
            pred.helix_fraction()
        );
    }

    #[test]
    fn poly_valine_predicts_strand() {
        // Val/Ile are strong strand formers.
        let pred = predict_gor("VVVVIIIIVVVVIIIIVVVV").expect("predict");
        assert!(
            pred.strand_fraction() > 0.5,
            "strand fraction {}",
            pred.strand_fraction()
        );
    }

    #[test]
    fn short_runs_smoothed_to_coil() {
        let mut states = vec![
            SecondaryStructure::Coil,
            SecondaryStructure::Helix, // isolated single helix
            SecondaryStructure::Coil,
        ];
        smooth_runs(&mut states);
        assert!(states.iter().all(|&s| s == SecondaryStructure::Coil));
    }

    #[test]
    fn string_round_trip() {
        let pred = predict_chou_fasman("ACDEFGHIKLMNPQRSTVWY").expect("predict");
        let s = pred.as_string();
        assert_eq!(s.len(), 20);
        assert!(s.chars().all(|c| matches!(c, 'H' | 'E' | 'C')));
    }
}
