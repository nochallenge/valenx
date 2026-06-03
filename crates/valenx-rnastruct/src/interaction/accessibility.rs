//! Accessibility — unpaired-region probability profiles.
//!
//! A region of an RNA can only be bound by another molecule (a
//! miRNA, an antisense oligo, a protein motif) if it is *accessible* —
//! not already locked inside intramolecular structure. Accessibility
//! is quantified from the McCaskill ensemble
//! ([`crate::ensemble::partition`]):
//!
//! - the **per-base unpaired probability** `p_u(i)`;
//! - the probability that a whole **window** `[i, i+w)` is unpaired,
//!   and the corresponding **opening free energy** `ΔG_open =
//!   −RT·ln P(window unpaired)` — the energetic price of clearing
//!   structure off that window.
//!
//! ## v1 scope
//!
//! The exact joint probability that *every* base in a window is
//! simultaneously unpaired requires a constrained partition function
//! per window. This v1 uses the standard tractable estimate — the
//! product of per-base unpaired probabilities — which is exact when
//! the window bases are structurally independent and a close upper
//! estimate otherwise. The approximation is documented here and on
//! [`AccessibilityProfile::window_unpaired_probability`].

use crate::ensemble::partition::{partition_function, PartitionFunction};
use crate::error::{Result, RnaStructError};
use crate::fold::energy::GAS_CONSTANT;
use crate::rna::RnaSeq;

/// A per-base accessibility profile.
#[derive(Clone, Debug)]
pub struct AccessibilityProfile {
    /// `unpaired[i]` = probability base `i` is unpaired, in `[0, 1]`.
    unpaired: Vec<f64>,
    /// The folding temperature, kelvin.
    temperature_k: f64,
}

impl AccessibilityProfile {
    /// Sequence length.
    pub fn len(&self) -> usize {
        self.unpaired.len()
    }

    /// `true` if the profile is empty.
    pub fn is_empty(&self) -> bool {
        self.unpaired.is_empty()
    }

    /// The unpaired probability of base `i` (0.0 out of range).
    pub fn unpaired(&self, i: usize) -> f64 {
        self.unpaired.get(i).copied().unwrap_or(0.0)
    }

    /// The whole per-base unpaired-probability vector.
    pub fn as_slice(&self) -> &[f64] {
        &self.unpaired
    }

    /// Estimated probability that the window `[start, start+len)` is
    /// entirely unpaired.
    ///
    /// Computed as the product of the per-base unpaired probabilities
    /// — see the module note on the independence approximation.
    /// Returns 0.0 for an out-of-range window.
    pub fn window_unpaired_probability(&self, start: usize, len: usize) -> f64 {
        if len == 0 || start + len > self.unpaired.len() {
            return 0.0;
        }
        self.unpaired[start..start + len].iter().product()
    }

    /// The opening free energy of the window `[start, start+len)`:
    /// `ΔG_open = −RT·ln P(window unpaired)`, kcal/mol.
    ///
    /// A large positive value means the window is buried (expensive to
    /// expose); near zero means it is already accessible. Returns
    /// `None` for an out-of-range window or a window whose probability
    /// underflows to zero.
    pub fn opening_energy(&self, start: usize, len: usize) -> Option<f64> {
        let p = self.window_unpaired_probability(start, len);
        if p <= 0.0 {
            return None;
        }
        Some(-GAS_CONSTANT * self.temperature_k * p.ln())
    }

    /// The most accessible window of width `w` — the `(start, ΔG_open)`
    /// minimising the opening energy. Returns `None` if `w` exceeds the
    /// sequence length.
    pub fn most_accessible_window(&self, w: usize) -> Option<(usize, f64)> {
        let n = self.unpaired.len();
        if w == 0 || w > n {
            return None;
        }
        let mut best: Option<(usize, f64)> = None;
        for start in 0..=(n - w) {
            if let Some(dg) = self.opening_energy(start, w) {
                if best.map(|(_, b)| dg < b).unwrap_or(true) {
                    best = Some((start, dg));
                }
            }
        }
        best
    }
}

/// Computes the accessibility profile of `seq` at 37 °C.
pub fn accessibility(seq: &RnaSeq) -> Result<AccessibilityProfile> {
    let pf = partition_function(seq)?;
    Ok(from_partition_function(&pf))
}

/// Builds an [`AccessibilityProfile`] from an existing partition
/// function (avoids recomputing the ensemble when the caller already
/// has one).
pub fn from_partition_function(pf: &PartitionFunction) -> AccessibilityProfile {
    let n = pf.len();
    let unpaired: Vec<f64> = (0..n).map(|i| pf.unpaired_probability(i)).collect();
    AccessibilityProfile {
        unpaired,
        temperature_k: pf.temperature_k(),
    }
}

/// Convenience: the per-base unpaired-probability vector of `seq`.
///
/// # Errors
/// Propagates partition-function errors.
pub fn unpaired_profile(seq: &RnaSeq) -> Result<Vec<f64>> {
    Ok(accessibility(seq)?.unpaired)
}

/// Returns the indices whose unpaired probability is at least
/// `threshold` — the "accessible" positions.
///
/// # Errors
/// [`RnaStructError::Invalid`] if `threshold` is outside `[0, 1]`.
pub fn accessible_positions(seq: &RnaSeq, threshold: f64) -> Result<Vec<usize>> {
    if !(0.0..=1.0).contains(&threshold) {
        return Err(RnaStructError::invalid(
            "threshold",
            "must be a probability in [0, 1]",
        ));
    }
    let prof = accessibility(seq)?;
    Ok((0..prof.len())
        .filter(|&i| prof.unpaired(i) >= threshold)
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn profile_length_matches_sequence() {
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let p = accessibility(&seq).unwrap();
        assert_eq!(p.len(), seq.len());
    }

    #[test]
    fn loop_region_is_more_accessible_than_stem() {
        // In a hairpin the loop bases should be more often unpaired
        // than the stem bases.
        let seq = RnaSeq::parse("GGGGGGGAAAACCCCCCC").unwrap();
        let p = accessibility(&seq).unwrap();
        let loop_acc = p.unpaired(9); // middle of the AAAA loop
        let stem_acc = p.unpaired(1); // inside the G stem
        assert!(
            loop_acc >= stem_acc,
            "loop ({loop_acc}) should be at least as accessible as stem ({stem_acc})"
        );
    }

    #[test]
    fn unpairable_rna_is_fully_accessible() {
        let seq = RnaSeq::parse("AAAAAAAAAAAA").unwrap();
        let p = accessibility(&seq).unwrap();
        for i in 0..seq.len() {
            assert!((p.unpaired(i) - 1.0).abs() < 1e-6);
        }
        // a fully-unpaired window has ~zero opening energy
        let dg = p.opening_energy(0, 5).unwrap();
        assert!(dg.abs() < 1e-3, "opening energy should be ~0, got {dg}");
    }

    #[test]
    fn opening_energy_is_nonnegative() {
        let seq = RnaSeq::parse("GGGGGGGAAAACCCCCCC").unwrap();
        let p = accessibility(&seq).unwrap();
        for start in 0..(seq.len() - 4) {
            if let Some(dg) = p.opening_energy(start, 4) {
                assert!(dg >= -1e-6, "opening energy {dg} must be >= 0");
            }
        }
    }

    #[test]
    fn most_accessible_window_found() {
        let seq = RnaSeq::parse("GGGGGGGAAAACCCCCCC").unwrap();
        let p = accessibility(&seq).unwrap();
        let (start, dg) = p.most_accessible_window(4).unwrap();
        assert!(start + 4 <= seq.len());
        assert!(dg.is_finite());
        // The loop window should beat a stem window. The stem at
        // positions 0..4 may be so deeply buried that its window
        // unpaired probability rounds to zero, in which case
        // `opening_energy` returns `None` — that is an *infinitely*
        // expensive window to open, so the loop trivially wins.
        match p.opening_energy(0, 4) {
            Some(stem_dg) => assert!(dg <= stem_dg + 1e-6),
            None => { /* stem is fully buried: ΔG_open = +∞ */ }
        }
    }

    #[test]
    fn accessible_positions_rejects_bad_threshold() {
        let seq = RnaSeq::parse("GGGGCCCC").unwrap();
        assert!(accessible_positions(&seq, 1.5).is_err());
        assert!(accessible_positions(&seq, -0.1).is_err());
        assert!(accessible_positions(&seq, 0.5).is_ok());
    }
}
