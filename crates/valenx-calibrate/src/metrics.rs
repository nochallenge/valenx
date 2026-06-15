//! Calibration metrics: Brier score, expected calibration error, and the
//! reliability diagram they summarise.

use serde::{Deserialize, Serialize};

use crate::error::CalibrateError;

fn check_probs_labels(probs: &[f64], labels: &[u8]) -> Result<(), CalibrateError> {
    if probs.is_empty() {
        return Err(CalibrateError::Empty { what: "probs" });
    }
    if probs.len() != labels.len() {
        return Err(CalibrateError::LengthMismatch {
            a: "probs",
            a_len: probs.len(),
            b: "labels",
            b_len: labels.len(),
        });
    }
    for &p in probs {
        if !p.is_finite() {
            return Err(CalibrateError::NonFinite {
                what: "probability",
            });
        }
        if !(0.0..=1.0).contains(&p) {
            return Err(CalibrateError::ProbOutOfRange { value: p });
        }
    }
    for &y in labels {
        if y > 1 {
            return Err(CalibrateError::LabelNotBinary { value: y });
        }
    }
    Ok(())
}

/// The Brier score: mean squared error between predicted probabilities and
/// `0/1` outcomes. Lower is better; `0` is perfect, `0.25` is the score of a
/// constant `0.5` predictor on balanced data.
pub fn brier_score(probs: &[f64], labels: &[u8]) -> Result<f64, CalibrateError> {
    check_probs_labels(probs, labels)?;
    let sum: f64 = probs
        .iter()
        .zip(labels)
        .map(|(&p, &y)| {
            let d = p - f64::from(y);
            d * d
        })
        .sum();
    Ok(sum / probs.len() as f64)
}

/// One equal-width confidence bin of a [`reliability_diagram`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ReliabilityBin {
    /// Inclusive lower edge of the bin's confidence range.
    pub lower: f64,
    /// Exclusive upper edge (inclusive for the final bin).
    pub upper: f64,
    /// Number of predictions that fell in this bin.
    pub count: usize,
    /// Mean predicted confidence of the bin's predictions (`NaN`-free: `0` when
    /// the bin is empty).
    pub mean_confidence: f64,
    /// Observed accuracy (mean label) of the bin's predictions (`0` when
    /// empty).
    pub accuracy: f64,
}

/// Partition `[0, 1]` into `n_bins` equal-width bins and summarise the
/// predictions' mean confidence vs observed accuracy in each. The basis of the
/// reliability diagram and of [`expected_calibration_error`].
pub fn reliability_diagram(
    probs: &[f64],
    labels: &[u8],
    n_bins: usize,
) -> Result<Vec<ReliabilityBin>, CalibrateError> {
    check_probs_labels(probs, labels)?;
    if n_bins == 0 {
        return Err(CalibrateError::ZeroBins);
    }
    let nb = n_bins as f64;
    let mut conf_sum = vec![0.0_f64; n_bins];
    let mut acc_sum = vec![0.0_f64; n_bins];
    let mut count = vec![0usize; n_bins];
    for (&p, &y) in probs.iter().zip(labels) {
        // floor(p * n_bins), clamped so p == 1.0 lands in the last bin.
        let mut idx = (p * nb).floor() as usize;
        if idx >= n_bins {
            idx = n_bins - 1;
        }
        conf_sum[idx] += p;
        acc_sum[idx] += f64::from(y);
        count[idx] += 1;
    }
    Ok((0..n_bins)
        .map(|i| {
            let c = count[i];
            let (mean_confidence, accuracy) = if c == 0 {
                (0.0, 0.0)
            } else {
                (conf_sum[i] / c as f64, acc_sum[i] / c as f64)
            };
            ReliabilityBin {
                lower: i as f64 / nb,
                upper: (i + 1) as f64 / nb,
                count: c,
                mean_confidence,
                accuracy,
            }
        })
        .collect())
}

/// Expected Calibration Error: the sample-weighted mean absolute gap between
/// confidence and accuracy across `n_bins` equal-width bins. `0` is perfectly
/// calibrated.
pub fn expected_calibration_error(
    probs: &[f64],
    labels: &[u8],
    n_bins: usize,
) -> Result<f64, CalibrateError> {
    let bins = reliability_diagram(probs, labels, n_bins)?;
    let n = probs.len() as f64;
    Ok(bins
        .iter()
        .filter(|b| b.count > 0)
        .map(|b| (b.count as f64 / n) * (b.accuracy - b.mean_confidence).abs())
        .sum())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn brier_perfect_and_chance() {
        assert!(brier_score(&[1.0, 0.0, 1.0], &[1, 0, 1]).unwrap().abs() < 1e-12);
        // constant 0.5 on balanced labels -> 0.25
        let b = brier_score(&[0.5, 0.5, 0.5, 0.5], &[1, 0, 1, 0]).unwrap();
        assert!((b - 0.25).abs() < 1e-12);
    }

    #[test]
    fn ece_zero_when_perfectly_calibrated() {
        // confidence exactly equals accuracy in every occupied bin.
        let probs = [1.0, 1.0, 0.0, 0.0];
        let labels = [1u8, 1, 0, 0];
        assert!(
            expected_calibration_error(&probs, &labels, 10)
                .unwrap()
                .abs()
                < 1e-12
        );
    }

    #[test]
    fn ece_positive_when_overconfident() {
        // all predicted 0.8, but only half are correct -> gap 0.3.
        let probs = [0.8, 0.8, 0.8, 0.8];
        let labels = [1u8, 1, 0, 0];
        let ece = expected_calibration_error(&probs, &labels, 10).unwrap();
        assert!((ece - 0.3).abs() < 1e-12);
    }

    #[test]
    fn reliability_bins_partition_counts() {
        let bins = reliability_diagram(&[0.05, 0.15, 0.95], &[0, 0, 1], 10).unwrap();
        assert_eq!(bins.len(), 10);
        assert_eq!(bins[0].count, 1); // 0.05 -> bin 0
        assert_eq!(bins[1].count, 1); // 0.15 -> bin 1
        assert_eq!(bins[9].count, 1); // 0.95 -> bin 9
        let total: usize = bins.iter().map(|b| b.count).sum();
        assert_eq!(total, 3);
    }

    #[test]
    fn prob_one_lands_in_last_bin() {
        let bins = reliability_diagram(&[1.0], &[1], 5).unwrap();
        assert_eq!(bins[4].count, 1);
    }

    #[test]
    fn rejects_bad_input() {
        assert_eq!(
            brier_score(&[1.5], &[1]).unwrap_err().code(),
            "prob_out_of_range"
        );
        assert_eq!(
            brier_score(&[0.5, 0.5], &[1]).unwrap_err().code(),
            "length_mismatch"
        );
        assert_eq!(
            reliability_diagram(&[0.5], &[1], 0).unwrap_err().code(),
            "zero_bins"
        );
    }
}
