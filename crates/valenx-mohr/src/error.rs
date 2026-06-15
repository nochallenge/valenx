//! Error taxonomy for 2D stress transformation.

use thiserror::Error;

/// Errors raised when constructing or transforming a plane-stress state.
#[derive(Debug, Error)]
pub enum MohrError {
    /// A supplied stress component was not a finite number (`NaN` or
    /// infinite). The closed-form transformation is only meaningful for
    /// finite inputs, so non-finite components are rejected at the
    /// constructor boundary.
    #[error("non-finite stress component `{name}`: {value}")]
    NonFinite {
        /// Name of the offending component (`sx`, `sy`, or `txy`).
        name: &'static str,
        /// The non-finite value that was supplied.
        value: f64,
    },

    /// A supplied plane angle was not a finite number (`NaN` or
    /// infinite).
    #[error("non-finite angle `{name}`: {value} rad")]
    NonFiniteAngle {
        /// Name of the offending angle parameter.
        name: &'static str,
        /// The non-finite value that was supplied.
        value: f64,
    },
}

/// Coarse category for a [`MohrError`], mirroring the sibling crates'
/// error-classification convention.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The error originates from caller-supplied input values.
    Input,
}

impl MohrError {
    /// Stable, kebab-cased identifier suitable for logging or telemetry.
    pub fn code(&self) -> &'static str {
        match self {
            MohrError::NonFinite { .. } => "mohr.non_finite",
            MohrError::NonFiniteAngle { .. } => "mohr.non_finite_angle",
        }
    }

    /// Coarse category for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            MohrError::NonFinite { .. } | MohrError::NonFiniteAngle { .. } => ErrorCategory::Input,
        }
    }
}
