//! Isotonic (monotonic) calibration via pool-adjacent-violators (PAV).

use serde::{Deserialize, Serialize};

use crate::error::CalibrateError;

/// Pool-adjacent-violators: the least-squares **non-decreasing** fit to
/// `values` with positive `weights`. The returned vector is the same length as
/// the input and is monotonically non-decreasing.
pub fn pav(values: &[f64], weights: &[f64]) -> Result<Vec<f64>, CalibrateError> {
    if values.is_empty() {
        return Err(CalibrateError::Empty { what: "values" });
    }
    if values.len() != weights.len() {
        return Err(CalibrateError::LengthMismatch {
            a: "values",
            a_len: values.len(),
            b: "weights",
            b_len: weights.len(),
        });
    }
    for (&v, &w) in values.iter().zip(weights) {
        if !v.is_finite() {
            return Err(CalibrateError::NonFinite { what: "value" });
        }
        if !w.is_finite() || w <= 0.0 {
            return Err(CalibrateError::NonFinite { what: "weight" });
        }
    }

    // Each block: (total weight, weighted sum of values, element count).
    let mut blocks: Vec<(f64, f64, usize)> = Vec::with_capacity(values.len());
    for (&v, &w) in values.iter().zip(weights) {
        let mut cur = (w, w * v, 1usize);
        while let Some(&last) = blocks.last() {
            let last_mean = last.1 / last.0;
            let cur_mean = cur.1 / cur.0;
            if last_mean > cur_mean {
                blocks.pop();
                cur = (last.0 + cur.0, last.1 + cur.1, last.2 + cur.2);
            } else {
                break;
            }
        }
        blocks.push(cur);
    }

    let mut out = Vec::with_capacity(values.len());
    for &(wsum, wysum, count) in &blocks {
        let mean = wysum / wsum;
        for _ in 0..count {
            out.push(mean);
        }
    }
    Ok(out)
}

/// A fitted isotonic calibrator: a monotone step function from score to
/// probability, learned from held-out `(score, label)` pairs.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct IsotonicCalibrator {
    /// Sorted scores (the step breakpoints).
    x: Vec<f64>,
    /// The PAV-fitted, non-decreasing calibrated probability at each breakpoint.
    y: Vec<f64>,
}

impl IsotonicCalibrator {
    /// Fit by sorting on score and running [`pav`] over the binary labels.
    pub fn fit(scores: &[f64], labels: &[u8]) -> Result<Self, CalibrateError> {
        crate::error::check_scores_labels(scores, labels)?;
        let mut idx: Vec<usize> = (0..scores.len()).collect();
        idx.sort_by(|&i, &j| {
            scores[i]
                .partial_cmp(&scores[j])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let x: Vec<f64> = idx.iter().map(|&i| scores[i]).collect();
        let sorted_labels: Vec<f64> = idx.iter().map(|&i| f64::from(labels[i])).collect();
        let weights = vec![1.0; sorted_labels.len()];
        let y = pav(&sorted_labels, &weights)?;
        Ok(Self { x, y })
    }

    /// Predict a calibrated probability for `score` by the monotone step
    /// function: the fitted value of the largest breakpoint not exceeding
    /// `score`, clamped to the endpoints outside the fitted range.
    pub fn calibrate(&self, score: f64) -> f64 {
        if score < self.x[0] {
            return self.y[0];
        }
        // number of breakpoints <= score
        let k = self.x.partition_point(|&b| b <= score);
        self.y[k - 1]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pav_pools_violators() {
        assert_eq!(
            pav(&[1.0, 2.0, 4.0, 3.0, 5.0], &[1.0; 5]).unwrap(),
            vec![1.0, 2.0, 3.5, 3.5, 5.0]
        );
    }

    #[test]
    fn pav_leaves_monotone_untouched() {
        assert_eq!(
            pav(&[1.0, 2.0, 3.0], &[1.0; 3]).unwrap(),
            vec![1.0, 2.0, 3.0]
        );
    }

    #[test]
    fn pav_pools_all_decreasing_to_mean() {
        assert_eq!(
            pav(&[3.0, 2.0, 1.0], &[1.0; 3]).unwrap(),
            vec![2.0, 2.0, 2.0]
        );
    }

    #[test]
    fn pav_output_is_non_decreasing() {
        let out = pav(&[0.0, 1.0, 0.0, 1.0, 1.0], &[1.0; 5]).unwrap();
        assert_eq!(out, vec![0.0, 0.5, 0.5, 1.0, 1.0]);
        assert!(out.windows(2).all(|w| w[0] <= w[1] + 1e-12));
    }

    #[test]
    fn isotonic_calibrator_is_monotone() {
        let cal = IsotonicCalibrator::fit(&[1.0, 2.0, 3.0, 4.0], &[0, 0, 1, 1]).unwrap();
        assert!((cal.calibrate(0.5) - 0.0).abs() < 1e-12); // below range -> first
        assert!((cal.calibrate(4.5) - 1.0).abs() < 1e-12); // above range -> last
                                                           // monotone across the range
        let a = cal.calibrate(1.5);
        let b = cal.calibrate(3.5);
        assert!(a <= b);
    }

    #[test]
    fn pav_rejects_bad_input() {
        assert_eq!(pav(&[], &[]).unwrap_err().code(), "empty");
        assert_eq!(
            pav(&[1.0], &[1.0, 2.0]).unwrap_err().code(),
            "length_mismatch"
        );
        assert_eq!(pav(&[1.0], &[0.0]).unwrap_err().code(), "non_finite"); // weight <= 0
    }
}
