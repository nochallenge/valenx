//! Error taxonomy for scoring.

use thiserror::Error;

/// Errors raised by the energy primitives, the endpoint estimator, or the
/// comparable-score aggregator.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum ScoreError {
    /// A value that must be finite was `NaN` or infinite.
    #[error("non-finite {what}")]
    NonFinite {
        /// What was non-finite.
        what: &'static str,
    },

    /// A quantity that must be strictly positive (distance, sigma, epsilon,
    /// Born radius) was zero or negative.
    #[error("{what} must be > 0, got {value}")]
    NonPositive {
        /// What was non-positive.
        what: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A dielectric constant was below 1 (vacuum).
    #[error("dielectric {value} must be >= 1")]
    DielectricTooSmall {
        /// The offending dielectric.
        value: f64,
    },

    /// A confidence value (ipTM / pLDDT) was outside `[0, 1]`.
    #[error("{what} {value} is outside [0, 1]")]
    ConfidenceOutOfRange {
        /// Which confidence channel.
        what: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A comparable score was requested but no evidence channel was present.
    #[error("no scoring components present")]
    NoComponents,
}

impl ScoreError {
    /// A short, stable machine-readable code for this error.
    pub fn code(&self) -> &'static str {
        match self {
            ScoreError::NonFinite { .. } => "non_finite",
            ScoreError::NonPositive { .. } => "non_positive",
            ScoreError::DielectricTooSmall { .. } => "dielectric_too_small",
            ScoreError::ConfidenceOutOfRange { .. } => "confidence_out_of_range",
            ScoreError::NoComponents => "no_components",
        }
    }
}

pub(crate) fn require_positive(what: &'static str, value: f64) -> Result<(), ScoreError> {
    if !value.is_finite() {
        return Err(ScoreError::NonFinite { what });
    }
    if value <= 0.0 {
        return Err(ScoreError::NonPositive { what, value });
    }
    Ok(())
}

pub(crate) fn require_finite(what: &'static str, value: f64) -> Result<(), ScoreError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(ScoreError::NonFinite { what })
    }
}
