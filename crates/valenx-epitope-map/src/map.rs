//! Sliding-window propensity profile and linear-epitope region calls.

use crate::error::EpitopeError;
use crate::scale::{aa_index, PropensityScale};

fn validate(seq: &str) -> Result<(), EpitopeError> {
    if seq.is_empty() {
        return Err(EpitopeError::EmptySequence);
    }
    for (pos, b) in seq.bytes().enumerate() {
        if aa_index(b).is_none() {
            return Err(EpitopeError::InvalidResidue {
                residue: char::from(b),
                pos,
            });
        }
    }
    Ok(())
}

/// Mean propensity of each sliding window of length `window`. The result has
/// `len - window + 1` entries.
pub fn propensity_profile(
    seq: &str,
    scale: &PropensityScale,
    window: usize,
) -> Result<Vec<f64>, EpitopeError> {
    validate(seq)?;
    let bytes = seq.as_bytes();
    if window == 0 || window > bytes.len() {
        return Err(EpitopeError::BadWindow {
            window,
            len: bytes.len(),
        });
    }
    let v: Vec<f64> = bytes.iter().map(|&b| scale.value(b).unwrap()).collect();
    let w = window as f64;
    Ok((0..=bytes.len() - window)
        .map(|i| v[i..i + window].iter().sum::<f64>() / w)
        .collect())
}

/// Contiguous residue spans `(start, end)` (end exclusive) whose windowed
/// propensity is at or above `threshold` — the candidate linear epitopes.
pub fn linear_epitope_regions(
    seq: &str,
    scale: &PropensityScale,
    window: usize,
    threshold: f64,
) -> Result<Vec<(usize, usize)>, EpitopeError> {
    if !threshold.is_finite() {
        return Err(EpitopeError::NonFinite { what: "threshold" });
    }
    let profile = propensity_profile(seq, scale, window)?;
    let mut regions = Vec::new();
    let mut start: Option<usize> = None;
    for (i, &m) in profile.iter().enumerate() {
        if m >= threshold {
            start.get_or_insert(i);
        } else if let Some(s) = start.take() {
            regions.push((s, i + window - 1));
        }
    }
    if let Some(s) = start.take() {
        regions.push((s, profile.len() - 1 + window));
    }
    Ok(regions)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scale::hydrophilicity_kd;

    #[test]
    fn profile_length_and_mean() {
        let s = hydrophilicity_kd();
        let p = propensity_profile("DEKRDE", &s, 2).unwrap();
        assert_eq!(p.len(), 6 - 2 + 1);
        // D and E both have hydrophilicity 3.5 -> window mean 3.5
        assert!((p[0] - 3.5).abs() < 1e-9);
    }

    #[test]
    fn hydrophilic_flagged_hydrophobic_not() {
        let s = hydrophilicity_kd();
        let hot = linear_epitope_regions("DEKRDEKR", &s, 3, 0.0).unwrap();
        assert_eq!(hot.len(), 1);
        assert_eq!(hot[0], (0, 8)); // whole hydrophilic stretch
        let cold = linear_epitope_regions("IILLVVFF", &s, 3, 0.0).unwrap();
        assert!(cold.is_empty());
    }

    #[test]
    fn threshold_filters() {
        let s = hydrophilicity_kd();
        // a very high threshold flags nothing
        assert!(linear_epitope_regions("DEKRDEKR", &s, 3, 100.0)
            .unwrap()
            .is_empty());
    }

    #[test]
    fn rejects_bad_input() {
        let s = hydrophilicity_kd();
        assert_eq!(
            propensity_profile("", &s, 3).unwrap_err().code(),
            "empty_sequence"
        );
        assert_eq!(
            propensity_profile("DEX", &s, 2).unwrap_err().code(),
            "invalid_residue"
        );
        assert_eq!(
            propensity_profile("DE", &s, 5).unwrap_err().code(),
            "bad_window"
        );
        assert_eq!(
            linear_epitope_regions("DE", &s, 2, f64::NAN)
                .unwrap_err()
                .code(),
            "non_finite"
        );
    }
}
