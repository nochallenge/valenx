//! Error taxonomy for binder scoring.

use thiserror::Error;

/// Errors raised while computing a binder-quality score.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum BinderError {
    /// No channel was present (or all had zero weight).
    #[error("no scoring components present")]
    NoComponents,

    /// A value that must be finite was `NaN` or infinite.
    #[error("non-finite {what}")]
    NonFinite {
        /// What was non-finite.
        what: &'static str,
    },

    /// The developability score was outside `[0, 1]`.
    #[error("developability {value} is outside [0, 1]")]
    DevelopabilityOutOfRange {
        /// The offending value.
        value: f64,
    },

    /// The safety severity was outside `0..=4`.
    #[error("safety severity {value} is outside 0..=4")]
    SeverityOutOfRange {
        /// The offending severity.
        value: u8,
    },

    /// A weight was negative or non-finite.
    #[error("weight {value} for {what} must be finite and >= 0")]
    BadWeight {
        /// Which weight.
        what: &'static str,
        /// The offending value.
        value: f64,
    },
}

impl BinderError {
    /// A short, stable machine-readable code for this error.
    pub fn code(&self) -> &'static str {
        match self {
            BinderError::NoComponents => "no_components",
            BinderError::NonFinite { .. } => "non_finite",
            BinderError::DevelopabilityOutOfRange { .. } => "developability_out_of_range",
            BinderError::SeverityOutOfRange { .. } => "severity_out_of_range",
            BinderError::BadWeight { .. } => "bad_weight",
        }
    }
}
