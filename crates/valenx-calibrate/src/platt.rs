//! Platt scaling: logistic recalibration `sigmoid(a * score + b)`.

use serde::{Deserialize, Serialize};

use crate::error::CalibrateError;
use crate::temperature::sigmoid;

/// Negative log-likelihood of a Platt fit — used by the tests to confirm a fit
/// beats the uninformative baseline. (The fit itself works in probability
/// space, so this is test-only.)
#[cfg(test)]
const EPS: f64 = 1e-12;

#[cfg(test)]
fn nll(scores: &[f64], labels: &[u8], a: f64, b: f64) -> f64 {
    scores
        .iter()
        .zip(labels)
        .map(|(&s, &y)| {
            let p = sigmoid(a * s + b).clamp(EPS, 1.0 - EPS);
            -(f64::from(y) * p.ln() + (1.0 - f64::from(y)) * (1.0 - p).ln())
        })
        .sum()
}

/// A fitted Platt (logistic) calibrator. `calibrate(s) = sigmoid(a*s + b)`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PlattScaler {
    a: f64,
    b: f64,
}

impl PlattScaler {
    /// Wrap known parameters (both must be finite).
    pub fn new(a: f64, b: f64) -> Result<Self, CalibrateError> {
        if !a.is_finite() || !b.is_finite() {
            return Err(CalibrateError::NonFinite { what: "parameter" });
        }
        Ok(Self { a, b })
    }

    /// Fit `a, b` by gradient descent on the negative log-likelihood of the
    /// held-out `(score, label)` pairs. Requires both classes to be present.
    pub fn fit(scores: &[f64], labels: &[u8]) -> Result<Self, CalibrateError> {
        crate::error::check_scores_labels(scores, labels)?;
        let pos = labels.iter().filter(|&&y| y == 1).count();
        if pos == 0 || pos == labels.len() {
            return Err(CalibrateError::SingleClass);
        }
        let n = scores.len() as f64;
        let (mut a, mut b) = (0.0_f64, 0.0_f64);
        let lr = 0.5;
        for _ in 0..5000 {
            let mut grad_a = 0.0;
            let mut grad_b = 0.0;
            for (&s, &y) in scores.iter().zip(labels) {
                let p = sigmoid(a * s + b);
                let r = p - f64::from(y);
                grad_a += r * s;
                grad_b += r;
            }
            a -= lr * grad_a / n;
            b -= lr * grad_b / n;
        }
        Self::new(a, b)
    }

    /// The fitted slope `a`.
    pub fn slope(&self) -> f64 {
        self.a
    }

    /// The fitted intercept `b`.
    pub fn intercept(&self) -> f64 {
        self.b
    }

    /// Map a raw score to a calibrated probability.
    pub fn calibrate(&self, score: f64) -> f64 {
        sigmoid(self.a * score + self.b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_parameters_recover_sigmoid() {
        let p = PlattScaler::new(1.0, 0.0).unwrap();
        assert!((p.calibrate(2.0) - sigmoid(2.0)).abs() < 1e-12);
    }

    #[test]
    fn fit_separable_data_has_positive_slope_and_lowers_nll() {
        let scores = [-2.0, -1.0, 1.0, 2.0];
        let labels = [0u8, 0, 1, 1];
        let p = PlattScaler::fit(&scores, &labels).unwrap();
        // higher score -> higher probability
        assert!(p.slope() > 0.0);
        assert!(p.calibrate(2.0) > p.calibrate(-2.0));
        // fit beats the uninformative a=b=0 baseline (constant 0.5)
        let baseline = nll(&scores, &labels, 0.0, 0.0);
        let fitted = nll(&scores, &labels, p.slope(), p.intercept());
        assert!(fitted < baseline);
    }

    #[test]
    fn fit_rejects_single_class() {
        assert_eq!(
            PlattScaler::fit(&[1.0, 2.0, 3.0], &[0, 0, 0])
                .unwrap_err()
                .code(),
            "single_class"
        );
    }
}
