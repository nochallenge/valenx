//! Hydrophobicity: GRAVY, windowed profile, and aggregation-prone regions.

use crate::aa::hydropathy;
use crate::error::DevelopabilityError;

fn validate(seq: &str) -> Result<(), DevelopabilityError> {
    if seq.is_empty() {
        return Err(DevelopabilityError::EmptySequence);
    }
    if let Some((pos, residue)) = crate::aa::first_invalid(seq) {
        return Err(DevelopabilityError::InvalidResidue { residue, pos });
    }
    Ok(())
}

/// Grand average of hydropathy (GRAVY): the mean Kyte-Doolittle value over the
/// sequence. Positive is overall hydrophobic.
pub fn gravy(seq: &str) -> Result<f64, DevelopabilityError> {
    validate(seq)?;
    let sum: f64 = seq.bytes().map(|b| hydropathy(b).unwrap()).sum();
    Ok(sum / seq.len() as f64)
}

/// Mean hydropathy of each sliding window of length `window`. The result has
/// `len - window + 1` entries.
pub fn hydropathy_profile(seq: &str, window: usize) -> Result<Vec<f64>, DevelopabilityError> {
    validate(seq)?;
    let bytes = seq.as_bytes();
    if window == 0 || window > bytes.len() {
        return Err(DevelopabilityError::BadWindow {
            window,
            len: bytes.len(),
        });
    }
    let h: Vec<f64> = bytes.iter().map(|&b| hydropathy(b).unwrap()).collect();
    let w = window as f64;
    Ok((0..=bytes.len() - window)
        .map(|i| h[i..i + window].iter().sum::<f64>() / w)
        .collect())
}

/// Aggregation-prone regions: maximal runs of `window`-windows whose mean
/// hydropathy is at or above `threshold`, returned as `(start, end)` (end
/// exclusive) residue spans. An illustrative hydrophobicity proxy, not a
/// validated aggregation predictor.
pub fn aggregation_prone_regions(
    seq: &str,
    window: usize,
    threshold: f64,
) -> Result<Vec<(usize, usize)>, DevelopabilityError> {
    if !threshold.is_finite() {
        return Err(DevelopabilityError::NonFinite { what: "threshold" });
    }
    let profile = hydropathy_profile(seq, window)?;
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

    #[test]
    fn gravy_is_mean_hydropathy() {
        // I = 4.5, V = 4.2 -> mean 4.35
        assert!((gravy("IV").unwrap() - 4.35).abs() < 1e-9);
        // D = E = -3.5 -> mean -3.5
        assert!((gravy("DE").unwrap() + 3.5).abs() < 1e-9);
    }

    #[test]
    fn hydrophobic_positive_charged_negative() {
        assert!(gravy("IIVVLL").unwrap() > 0.0);
        assert!(gravy("DDEEKK").unwrap() < 0.0);
    }

    #[test]
    fn profile_has_right_length() {
        let p = hydropathy_profile("IIIIVV", 3).unwrap();
        assert_eq!(p.len(), 6 - 3 + 1);
    }

    #[test]
    fn apr_flags_hydrophobic_not_charged() {
        // hydrophobic stretch -> one region spanning it
        let hot = aggregation_prone_regions("IIIIVVVV", 3, 1.5).unwrap();
        assert_eq!(hot.len(), 1);
        assert_eq!(hot[0], (0, 8));
        // hydrophilic stretch -> none
        let cold = aggregation_prone_regions("DDDDEEEE", 3, 1.5).unwrap();
        assert!(cold.is_empty());
    }

    #[test]
    fn rejects_bad_input() {
        assert_eq!(gravy("").unwrap_err().code(), "empty_sequence");
        assert_eq!(gravy("IXV").unwrap_err().code(), "invalid_residue");
        assert_eq!(
            hydropathy_profile("AA", 5).unwrap_err().code(),
            "bad_window"
        );
    }
}
