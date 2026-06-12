//! SHAPE-reactivity-constrained folding (Deigan model).
//!
//! SHAPE (Selective 2′-Hydroxyl Acylation analysed by Primer
//! Extension) chemistry probes which nucleotides are flexible:
//! flexible bases — typically unpaired — react strongly, structured
//! bases react weakly. Deigan *et al.* (2009) showed a per-nucleotide
//! pseudo-energy derived from the reactivity, added to the folding
//! free energy, dramatically improves prediction accuracy.
//!
//! This module is a thin façade: it converts a reactivity vector into
//! a [`FoldConstraints`] soft term via [`FoldConstraints::from_shape`]
//! and runs the Zuker MFE folder. The interesting code is the Deigan
//! transform in [`crate::fold::constraint`].

use crate::error::{Result, RnaStructError};
use crate::fold::constraint::FoldConstraints;
use crate::fold::zuker::{mfe_constrained, MfeResult};
use crate::rna::RnaSeq;

/// Parameters for SHAPE-directed folding.
#[derive(Copy, Clone, Debug)]
pub struct ShapeParams {
    /// Deigan slope `m` (kcal/mol).
    pub slope: f64,
    /// Deigan intercept `b` (kcal/mol).
    pub intercept: f64,
}

impl Default for ShapeParams {
    fn default() -> Self {
        ShapeParams {
            slope: crate::fold::constraint::DEIGAN_SLOPE,
            intercept: crate::fold::constraint::DEIGAN_INTERCEPT,
        }
    }
}

/// Folds `seq` with SHAPE reactivities using the default Deigan
/// parameters.
///
/// `reactivities[i]` is the normalised SHAPE reactivity of base `i`;
/// a negative value marks missing data (and contributes no
/// pseudo-energy). The vector length must equal the sequence length.
///
/// # Errors
/// [`RnaStructError::Invalid`] on a length mismatch.
pub fn fold_with_shape(seq: &RnaSeq, reactivities: &[f64]) -> Result<MfeResult> {
    fold_with_shape_params(seq, reactivities, ShapeParams::default())
}

/// [`fold_with_shape`] with explicit [`ShapeParams`].
///
/// # Errors
/// [`RnaStructError::Invalid`] on a length mismatch.
pub fn fold_with_shape_params(
    seq: &RnaSeq,
    reactivities: &[f64],
    params: ShapeParams,
) -> Result<MfeResult> {
    if reactivities.len() != seq.len() {
        return Err(RnaStructError::invalid(
            "reactivities",
            format!(
                "{} SHAPE values for a length-{} sequence",
                reactivities.len(),
                seq.len()
            ),
        ));
    }
    let cons =
        FoldConstraints::from_shape_with(seq.len(), reactivities, params.slope, params.intercept)?;
    mfe_constrained(seq, &cons)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fold::zuker::mfe;

    #[test]
    fn shape_length_must_match() {
        let seq = RnaSeq::parse("GGGGAAAACCCC").unwrap();
        assert!(fold_with_shape(&seq, &[1.0, 2.0]).is_err());
    }

    #[test]
    fn flat_low_reactivity_resembles_plain_fold() {
        // All-low reactivity adds a near-uniform constant; the MFE
        // structure should still be a sensible fold.
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let react = vec![0.0; seq.len()];
        let shaped = fold_with_shape(&seq, &react).unwrap();
        assert!(shaped.structure.is_nested());
        // plain fold still works
        let plain = mfe(&seq).unwrap();
        assert!(plain.structure.is_nested());
    }

    #[test]
    fn strong_reactivity_discourages_pairing() {
        // Make the whole stem highly reactive -> SHAPE should push the
        // fold towards fewer pairs than the plain MFE.
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let plain = mfe(&seq).unwrap();
        let mut react = vec![0.0; seq.len()];
        for r in react.iter_mut() {
            *r = 4.0; // very reactive everywhere
        }
        let shaped = fold_with_shape(&seq, &react).unwrap();
        assert!(
            shaped.structure.n_pairs() <= plain.structure.n_pairs(),
            "high SHAPE reactivity should not increase pairing"
        );
    }

    #[test]
    fn missing_data_is_tolerated() {
        let seq = RnaSeq::parse("GGGGAAAACCCC").unwrap();
        let react = vec![-1.0; seq.len()]; // all missing
        let r = fold_with_shape(&seq, &react).unwrap();
        assert!(r.structure.is_nested());
    }
}
