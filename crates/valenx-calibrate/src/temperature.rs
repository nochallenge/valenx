//! Temperature scaling: a single-parameter recalibration of logits
//! (Guo *et al.*, 2017).

use serde::{Deserialize, Serialize};

use crate::error::CalibrateError;

/// The logistic sigmoid `1 / (1 + e^-z)`.
pub fn sigmoid(z: f64) -> f64 {
    1.0 / (1.0 + (-z).exp())
}

const EPS: f64 = 1e-12;

fn nll(logits: &[f64], labels: &[u8], t: f64) -> f64 {
    logits
        .iter()
        .zip(labels)
        .map(|(&z, &y)| {
            let p = sigmoid(z / t).clamp(EPS, 1.0 - EPS);
            -(f64::from(y) * p.ln() + (1.0 - f64::from(y)) * (1.0 - p).ln())
        })
        .sum()
}

/// A fitted temperature. Divide a logit by [`TemperatureScaler::temperature`]
/// before the sigmoid: `T > 1` softens overconfident predictions toward `0.5`,
/// `T < 1` sharpens them.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct TemperatureScaler {
    t: f64,
}

impl TemperatureScaler {
    /// Wrap a known temperature (`> 0`).
    pub fn new(t: f64) -> Result<Self, CalibrateError> {
        if !t.is_finite() || t <= 0.0 {
            return Err(CalibrateError::NonPositiveTemperature { t });
        }
        Ok(Self { t })
    }

    /// Fit the temperature that minimises negative log-likelihood on held-out
    /// `(logit, label)` pairs, by golden-section search over `[0.05, 100]`.
    ///
    /// Requires both classes to be present (a single-class set has no
    /// finite-temperature optimum).
    pub fn fit(logits: &[f64], labels: &[u8]) -> Result<Self, CalibrateError> {
        crate::error::check_scores_labels(logits, labels)?;
        let pos = labels.iter().filter(|&&y| y == 1).count();
        if pos == 0 || pos == labels.len() {
            return Err(CalibrateError::SingleClass);
        }
        // Golden-section minimisation of a well-behaved 1-D objective.
        let gr = (5.0_f64.sqrt() - 1.0) / 2.0;
        let (mut a, mut b) = (0.05_f64, 100.0_f64);
        let mut c = b - gr * (b - a);
        let mut d = a + gr * (b - a);
        for _ in 0..200 {
            if nll(logits, labels, c) < nll(logits, labels, d) {
                b = d;
            } else {
                a = c;
            }
            c = b - gr * (b - a);
            d = a + gr * (b - a);
        }
        Self::new((a + b) / 2.0)
    }

    /// The fitted (or supplied) temperature.
    pub fn temperature(&self) -> f64 {
        self.t
    }

    /// Apply the temperature to a logit, returning a calibrated probability.
    pub fn calibrate(&self, logit: f64) -> f64 {
        sigmoid(logit / self.t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sigmoid_known_values() {
        assert!((sigmoid(0.0) - 0.5).abs() < 1e-12);
        assert!(sigmoid(50.0) > 0.999);
        assert!(sigmoid(-50.0) < 0.001);
    }

    #[test]
    fn temperature_one_is_identity() {
        let s = TemperatureScaler::new(1.0).unwrap();
        assert!((s.calibrate(2.0) - sigmoid(2.0)).abs() < 1e-12);
    }

    #[test]
    fn rejects_non_positive_temperature() {
        assert_eq!(
            TemperatureScaler::new(0.0).unwrap_err().code(),
            "non_positive_temperature"
        );
        assert_eq!(
            TemperatureScaler::new(-1.0).unwrap_err().code(),
            "non_positive_temperature"
        );
    }

    #[test]
    fn fit_softens_overconfident_logits() {
        // Large-magnitude logits but only chance-level accuracy -> the model is
        // overconfident, so the best temperature is > 1 (push toward 0.5).
        let logits = [6.0, 6.0, 6.0, 6.0, -6.0, -6.0, -6.0, -6.0];
        let labels = [1u8, 1, 0, 0, 0, 0, 1, 1]; // half wrong
        let s = TemperatureScaler::fit(&logits, &labels).unwrap();
        assert!(s.temperature() > 1.0, "T = {}", s.temperature());
        // And the fit does not increase NLL versus T = 1.
        assert!(nll(&logits, &labels, s.temperature()) <= nll(&logits, &labels, 1.0) + 1e-9);
    }

    #[test]
    fn fit_rejects_single_class() {
        assert_eq!(
            TemperatureScaler::fit(&[1.0, 2.0], &[1, 1])
                .unwrap_err()
                .code(),
            "single_class"
        );
    }
}
