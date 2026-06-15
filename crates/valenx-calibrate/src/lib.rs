//! # valenx-calibrate
//!
//! Turn raw predictive scores into **calibrated confidences** — and measure how
//! calibrated they are. This is the "never an automatic green light" layer of a
//! biologic-design pipeline: a model that says 0.9 should be right 90% of the
//! time, and this crate is what makes (and checks) that promise.
//!
//! ## What
//!
//! - [`metrics::brier_score`] — mean-squared error of probabilistic predictions
//!   (a proper scoring rule).
//! - [`metrics::expected_calibration_error`] / [`metrics::reliability_diagram`]
//!   — bin predictions and compare confidence against observed accuracy.
//! - [`platt::PlattScaler`] — logistic (Platt) calibration `sigmoid(a*s + b)`.
//! - [`temperature::TemperatureScaler`] — single-parameter temperature scaling
//!   of logits (Guo *et al.* 2017).
//! - [`isotonic::pav`] / [`isotonic::IsotonicCalibrator`] — non-parametric
//!   monotonic calibration by pool-adjacent-violators.
//! - [`conformal::split_conformal_quantile`] / [`conformal::conformal_interval`]
//!   — distribution-free prediction intervals with a finite-sample coverage
//!   guarantee (split conformal).
//!
//! ## Model
//!
//! Calibration learns a monotone map from a model's raw score to a probability
//! (or, for conformal, a coverage-calibrated interval) using a held-out
//! calibration set. Platt fits a two-parameter sigmoid; temperature scaling
//! fits one parameter on logits; isotonic fits an arbitrary monotone step
//! function. The metrics quantify the residual gap between stated confidence
//! and empirical accuracy. Split conformal converts held-out nonconformity
//! scores into an interval that covers a new point with probability at least
//! `1 - alpha`, assuming exchangeability.
//!
//! ## Honest scope
//!
//! Research/educational grade. These are the standard, citable calibration
//! methods, implemented transparently and checked against analytic ground
//! truth. But calibration is only as good as the calibration set: it assumes
//! the held-out data is representative and exchangeable with deployment, and it
//! does **not** create information the underlying model lacks. Nothing here is a
//! validated clinical or regulatory confidence statement — a calibrated 0.9 is
//! a statement about a dataset, not a guarantee about a new biological reality.
//!
//! ## Example
//!
//! ```
//! use valenx_calibrate::metrics::brier_score;
//! use valenx_calibrate::isotonic::pav;
//!
//! // A perfectly confident, perfectly correct predictor scores Brier 0.
//! let probs = [1.0, 0.0, 1.0];
//! let labels = [1u8, 0, 1];
//! assert!(brier_score(&probs, &labels).unwrap().abs() < 1e-12);
//!
//! // Pool-adjacent-violators makes a sequence non-decreasing.
//! let fitted = pav(&[1.0, 2.0, 4.0, 3.0, 5.0], &[1.0; 5]).unwrap();
//! assert_eq!(fitted, vec![1.0, 2.0, 3.5, 3.5, 5.0]);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod conformal;
pub mod error;
pub mod isotonic;
pub mod metrics;
pub mod platt;
pub mod temperature;

pub use conformal::{conformal_interval, split_conformal_quantile};
pub use error::CalibrateError;
pub use isotonic::{pav, IsotonicCalibrator};
pub use metrics::{brier_score, expected_calibration_error, reliability_diagram, ReliabilityBin};
pub use platt::PlattScaler;
pub use temperature::{sigmoid, TemperatureScaler};
