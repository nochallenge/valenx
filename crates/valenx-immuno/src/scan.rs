//! Sliding-window epitope scanning and summary scores.

use serde::{Deserialize, Serialize};

use crate::aa::aa_index;
use crate::error::ImmunoError;
use crate::matrix::Pssm;

/// One scored window produced by [`scan`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EpitopeHit {
    /// Zero-based start position of the window within the protein.
    pub start: usize,
    /// The window's residues.
    pub peptide: String,
    /// The window's PSSM score.
    pub score: f64,
}

/// Score every length-`pssm.length()` window of `protein`, in sequence order.
///
/// Validates that the protein is non-empty, contains only standard amino acids
/// (failing fast with the offending position), and is at least one window long.
pub fn scan(pssm: &Pssm, protein: &str) -> Result<Vec<EpitopeHit>, ImmunoError> {
    let window = pssm.length();
    let bytes = protein.as_bytes();
    if bytes.is_empty() {
        return Err(ImmunoError::Empty { what: "protein" });
    }
    for (pos, &b) in bytes.iter().enumerate() {
        if aa_index(b).is_none() {
            return Err(ImmunoError::InvalidResidue {
                residue: char::from(b),
                pos,
            });
        }
    }
    if bytes.len() < window {
        return Err(ImmunoError::ProteinTooShort {
            protein: bytes.len(),
            window,
        });
    }
    let mut hits = Vec::with_capacity(bytes.len() - window + 1);
    for start in 0..=(bytes.len() - window) {
        let peptide = &protein[start..start + window];
        let score = pssm.score(peptide)?;
        hits.push(EpitopeHit {
            start,
            peptide: peptide.to_string(),
            score,
        });
    }
    Ok(hits)
}

/// All windows scoring at or above `threshold`, in sequence order.
pub fn scan_threshold(
    pssm: &Pssm,
    protein: &str,
    threshold: f64,
) -> Result<Vec<EpitopeHit>, ImmunoError> {
    if !threshold.is_finite() {
        return Err(ImmunoError::NonFiniteThreshold(threshold));
    }
    let mut hits = scan(pssm, protein)?;
    hits.retain(|h| h.score >= threshold);
    Ok(hits)
}

/// The `n` highest-scoring hits, score-descending, ties broken by earliest
/// start. Consumes `hits`; pass the output of [`scan`].
pub fn top_n(mut hits: Vec<EpitopeHit>, n: usize) -> Vec<EpitopeHit> {
    hits.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.start.cmp(&b.start))
    });
    hits.truncate(n);
    hits
}

/// The fraction of windows scoring at or above `threshold`, in `[0, 1]`.
///
/// A coarse, single-number immunogenicity flag: the predicted-epitope density.
pub fn epitope_density(pssm: &Pssm, protein: &str, threshold: f64) -> Result<f64, ImmunoError> {
    if !threshold.is_finite() {
        return Err(ImmunoError::NonFiniteThreshold(threshold));
    }
    let hits = scan(pssm, protein)?;
    // `scan` guarantees at least one window, so the denominator is never zero.
    let above = hits.iter().filter(|h| h.score >= threshold).count();
    Ok(above as f64 / hits.len() as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aa::N_AA;

    /// Matrix of length 3 that rewards `K` at the middle position by 5.0.
    fn mid_k() -> Pssm {
        let mut rows = vec![[0.0; N_AA]; 3];
        rows[1][aa_index(b'K').unwrap()] = 5.0;
        Pssm::new("mid_k", rows).unwrap()
    }

    #[test]
    fn scan_yields_one_window_per_offset() {
        let hits = scan(&mid_k(), "AAKAA").unwrap();
        assert_eq!(hits.len(), 5 - 3 + 1);
        assert_eq!(hits[0].start, 0);
        assert_eq!(hits[1].start, 1);
        assert_eq!(hits[2].start, 2);
        assert_eq!(hits[0].peptide, "AAK");
        assert_eq!(hits[1].peptide, "AKA");
        assert_eq!(hits[2].peptide, "KAA");
        // 'K' sits at the rewarded middle position only in window 1.
        assert!((hits[1].score - 5.0).abs() < 1e-12);
        assert!((hits[0].score - 0.0).abs() < 1e-12);
    }

    #[test]
    fn scan_rejects_short_protein() {
        let err = scan(&mid_k(), "AA").unwrap_err();
        assert_eq!(err.code(), "protein_too_short");
    }

    #[test]
    fn scan_rejects_empty_and_invalid() {
        assert_eq!(scan(&mid_k(), "").unwrap_err().code(), "empty");
        let err = scan(&mid_k(), "AABAA").unwrap_err();
        assert_eq!(err.code(), "invalid_residue"); // 'B' is not standard
    }

    #[test]
    fn threshold_filters_and_validates() {
        let hits = scan_threshold(&mid_k(), "AAKAA", 1.0).unwrap();
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].peptide, "AKA");
        assert_eq!(
            scan_threshold(&mid_k(), "AAKAA", f64::NAN)
                .unwrap_err()
                .code(),
            "non_finite_threshold"
        );
    }

    #[test]
    fn top_n_sorts_descending_then_by_start() {
        let hits = scan(&mid_k(), "AKAKA").unwrap(); // two K-in-middle windows
        let top = top_n(hits, 2);
        assert!(top[0].score >= top[1].score);
        // both top windows score 5.0; tie broken by earliest start
        assert!(top[0].start < top[1].start);
    }

    #[test]
    fn density_is_bounded_fraction() {
        // 3 windows; one above threshold 1.0 -> 1/3.
        let d = epitope_density(&mid_k(), "AAKAA", 1.0).unwrap();
        assert!((d - 1.0 / 3.0).abs() < 1e-12);
        // threshold below every score -> all windows count -> 1.0
        let all = epitope_density(&mid_k(), "AAKAA", -1.0).unwrap();
        assert!((all - 1.0).abs() < 1e-12);
        // threshold above every score -> 0.0
        let none = epitope_density(&mid_k(), "AAKAA", 100.0).unwrap();
        assert!(none.abs() < 1e-12);
    }
}
