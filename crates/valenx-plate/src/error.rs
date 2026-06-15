//! Error taxonomy for the plate-bending workbench.
//!
//! Every fallible constructor and analysis function in this crate returns
//! [`Result<_, PlateError>`]. The variants distinguish a bad scalar
//! parameter (the common case — a non-positive thickness, a Poisson ratio
//! outside the physically admissible open interval, etc.) from a thin-plate
//! modelling-assumption violation (an aspect ratio so small the plate is no
//! longer "thin" and Kirchhoff-Love theory no longer applies).

use thiserror::Error;

/// Errors raised while building plate inputs or evaluating plate models.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum PlateError {
    /// A scalar input was outside its admissible range.
    ///
    /// `name` is the offending parameter (e.g. `"thickness"`,
    /// `"poisson_ratio"`); `reason` explains the violated constraint.
    #[error("invalid parameter `{name}`: {reason} (got {value})")]
    InvalidParameter {
        /// Name of the offending parameter.
        name: &'static str,
        /// Human-readable description of the violated constraint.
        reason: String,
        /// The supplied (rejected) value.
        value: f64,
    },

    /// The thin-plate (Kirchhoff-Love) assumption is violated.
    ///
    /// Raised when the radius-to-thickness ratio `a / t` is below the
    /// small threshold beyond which transverse-shear and
    /// thickness-stretch effects this crate ignores stop being
    /// negligible, so the closed-form results would be misleading.
    #[error(
        "thin-plate assumption violated: radius/thickness ratio {ratio} \
         is below the minimum {min} for Kirchhoff-Love theory"
    )]
    NotThin {
        /// The supplied radius-to-thickness ratio `a / t`.
        ratio: f64,
        /// The minimum admissible ratio.
        min: f64,
    },
}

/// Coarse category for telemetry / triage, mirroring sibling crates.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The caller supplied an out-of-range input value.
    Input,
    /// A modelling assumption of the chosen theory was violated.
    Algorithm,
}

impl PlateError {
    /// Construct an [`InvalidParameter`](PlateError::InvalidParameter) error.
    ///
    /// This is the single internal helper the validated constructors use so
    /// the `name` / `reason` / `value` triple is assembled consistently.
    pub fn invalid(name: &'static str, reason: impl Into<String>, value: f64) -> Self {
        PlateError::InvalidParameter {
            name,
            reason: reason.into(),
            value,
        }
    }

    /// Stable, kebab-cased identifier suitable for logs and error tables.
    pub fn code(&self) -> &'static str {
        match self {
            PlateError::InvalidParameter { .. } => "plate.invalid_parameter",
            PlateError::NotThin { .. } => "plate.not_thin",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            PlateError::InvalidParameter { .. } => ErrorCategory::Input,
            PlateError::NotThin { .. } => ErrorCategory::Algorithm,
        }
    }
}

/// Validate that `value` is finite and strictly positive.
///
/// Returns `value` unchanged on success, or an
/// [`InvalidParameter`](PlateError::InvalidParameter) error naming `name`.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64, PlateError> {
    if !value.is_finite() {
        return Err(PlateError::invalid(name, "must be a finite number", value));
    }
    if value <= 0.0 {
        return Err(PlateError::invalid(
            name,
            "must be strictly positive",
            value,
        ));
    }
    Ok(value)
}
