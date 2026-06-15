//! Error taxonomy for calibration.

use thiserror::Error;

/// Errors raised while fitting or evaluating a calibrator.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum CalibrateError {
    /// A required input slice was empty.
    #[error("empty {what}")]
    Empty {
        /// What was empty (e.g. `"scores"`, `"residuals"`).
        what: &'static str,
    },

    /// Two paired slices had different lengths.
    #[error("length mismatch: {a} has {a_len} elements but {b} has {b_len}")]
    LengthMismatch {
        /// Name of the first slice.
        a: &'static str,
        /// Length of the first slice.
        a_len: usize,
        /// Name of the second slice.
        b: &'static str,
        /// Length of the second slice.
        b_len: usize,
    },

    /// A value that must be finite was `NaN` or infinite.
    #[error("non-finite {what}")]
    NonFinite {
        /// What was non-finite.
        what: &'static str,
    },

    /// A probability fell outside `[0, 1]`.
    #[error("probability {value} is outside [0, 1]")]
    ProbOutOfRange {
        /// The offending probability.
        value: f64,
    },

    /// A label was not `0` or `1`.
    #[error("label {value} is not binary (expected 0 or 1)")]
    LabelNotBinary {
        /// The offending label.
        value: u8,
    },

    /// A miscoverage level `alpha` was not in the open interval `(0, 1)`.
    #[error("alpha {alpha} is not in the open interval (0, 1)")]
    AlphaOutOfRange {
        /// The offending alpha.
        alpha: f64,
    },

    /// A temperature was not strictly positive.
    #[error("temperature {t} must be > 0")]
    NonPositiveTemperature {
        /// The offending temperature.
        t: f64,
    },

    /// A bin count was zero.
    #[error("number of bins must be >= 1")]
    ZeroBins,

    /// The calibration set contained only one class, so a parametric fit is
    /// undefined.
    #[error("calibration labels are all the same class; cannot fit")]
    SingleClass,
}

impl CalibrateError {
    /// A short, stable machine-readable code for this error.
    pub fn code(&self) -> &'static str {
        match self {
            CalibrateError::Empty { .. } => "empty",
            CalibrateError::LengthMismatch { .. } => "length_mismatch",
            CalibrateError::NonFinite { .. } => "non_finite",
            CalibrateError::ProbOutOfRange { .. } => "prob_out_of_range",
            CalibrateError::LabelNotBinary { .. } => "label_not_binary",
            CalibrateError::AlphaOutOfRange { .. } => "alpha_out_of_range",
            CalibrateError::NonPositiveTemperature { .. } => "non_positive_temperature",
            CalibrateError::ZeroBins => "zero_bins",
            CalibrateError::SingleClass => "single_class",
        }
    }
}

/// Validate paired score/label slices: non-empty, equal length, finite scores,
/// binary labels. Shared by the parametric calibrators.
pub(crate) fn check_scores_labels(scores: &[f64], labels: &[u8]) -> Result<(), CalibrateError> {
    if scores.is_empty() {
        return Err(CalibrateError::Empty { what: "scores" });
    }
    if scores.len() != labels.len() {
        return Err(CalibrateError::LengthMismatch {
            a: "scores",
            a_len: scores.len(),
            b: "labels",
            b_len: labels.len(),
        });
    }
    for &s in scores {
        if !s.is_finite() {
            return Err(CalibrateError::NonFinite { what: "score" });
        }
    }
    for &y in labels {
        if y > 1 {
            return Err(CalibrateError::LabelNotBinary { value: y });
        }
    }
    Ok(())
}
