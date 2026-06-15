//! Bolted-joint error taxonomy.
//!
//! Every fallible constructor in this crate funnels through
//! [`BoltError`]. The variants are deliberately specific so a caller can
//! tell a bad nut factor from a bad stiffness or a negative load apart,
//! and each carries the offending value for diagnostics.

use thiserror::Error;

/// Errors raised while building or evaluating a bolted joint.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum BoltError {
    /// A length, diameter or area that must be strictly positive was
    /// given as zero or negative.
    #[error("`{name}` must be > 0, got {value}")]
    NonPositive {
        /// Name of the offending parameter.
        name: &'static str,
        /// The value that was supplied.
        value: f64,
    },

    /// A value that must be finite (not NaN / not infinite) was not.
    #[error("`{name}` must be finite, got {value}")]
    NotFinite {
        /// Name of the offending parameter.
        name: &'static str,
        /// The value that was supplied.
        value: f64,
    },

    /// The nut factor `K` (a dimensionless friction coefficient) was
    /// outside the physically meaningful open interval `(0, 1)`.
    /// Typical values are ~0.10 (lubricated) to ~0.30 (dry, plated).
    #[error("nut factor K must be in (0, 1), got {value}")]
    NutFactorRange {
        /// The value that was supplied.
        value: f64,
    },

    /// A stiffness ratio / load-share fraction that must lie in the open
    /// interval `(0, 1)` was out of range. The bolt and the clamped
    /// members are springs in parallel, so neither can carry all nor
    /// none of an external load.
    #[error("stiffness ratio C must be in (0, 1), got {value}")]
    StiffnessRatioRange {
        /// The value that was supplied.
        value: f64,
    },

    /// An external service load was negative. The convention here is a
    /// tensile (joint-opening) external load is non-negative; a negative
    /// value is rejected so sign errors surface early.
    #[error("external load `{name}` must be >= 0, got {value}")]
    NegativeLoad {
        /// Name of the offending parameter.
        name: &'static str,
        /// The value that was supplied.
        value: f64,
    },
}

/// Coarse error category, useful for telemetry / UI grouping.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The error stems from an out-of-domain user input.
    Input,
    /// The error stems from a non-finite / NaN value leaking in.
    Numeric,
}

impl BoltError {
    /// Build a [`BoltError::NonPositive`] after also rejecting non-finite
    /// input, so callers get the most specific error for a bad value.
    pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64, BoltError> {
        if !value.is_finite() {
            return Err(BoltError::NotFinite { name, value });
        }
        if value <= 0.0 {
            return Err(BoltError::NonPositive { name, value });
        }
        Ok(value)
    }

    /// Reject a negative or non-finite external load.
    pub(crate) fn require_non_negative(name: &'static str, value: f64) -> Result<f64, BoltError> {
        if !value.is_finite() {
            return Err(BoltError::NotFinite { name, value });
        }
        if value < 0.0 {
            return Err(BoltError::NegativeLoad { name, value });
        }
        Ok(value)
    }

    /// Stable kebab-cased identifier for logs / dashboards.
    pub fn code(&self) -> &'static str {
        match self {
            BoltError::NonPositive { .. } => "bolt.non-positive",
            BoltError::NotFinite { .. } => "bolt.not-finite",
            BoltError::NutFactorRange { .. } => "bolt.nut-factor-range",
            BoltError::StiffnessRatioRange { .. } => "bolt.stiffness-ratio-range",
            BoltError::NegativeLoad { .. } => "bolt.negative-load",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            BoltError::NotFinite { .. } => ErrorCategory::Numeric,
            _ => ErrorCategory::Input,
        }
    }
}
