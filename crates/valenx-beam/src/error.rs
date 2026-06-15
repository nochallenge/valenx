//! Beam-bending error taxonomy.
//!
//! A single [`BeamError`] enum covers the two ways a calculation can be
//! rejected: a structurally invalid input (a non-positive or non-finite
//! geometric / material / load quantity) and a degenerate section whose
//! second moment of area works out to zero.

use thiserror::Error;

/// Errors raised by beam, section and stress calculations.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum BeamError {
    /// A scalar input was non-positive or non-finite where a strictly
    /// positive, finite value is required (length, load, modulus,
    /// dimension, …).
    #[error("bad parameter `{name}`: expected a finite value > 0, got {value}")]
    BadParameter {
        /// Name of the offending parameter.
        name: &'static str,
        /// The value that was supplied.
        value: f64,
    },

    /// The supplied geometry has a zero (or numerically vanishing)
    /// second moment of area, so stress `M*c/I` is undefined.
    #[error("degenerate section: second moment of area is zero ({reason})")]
    DegenerateSection {
        /// Human-readable reason.
        reason: &'static str,
    },
}

/// Coarse category for an error, useful for UI grouping and metrics.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// Caused by a user-supplied input value.
    Input,
    /// Caused by a degenerate / unsupported geometric configuration.
    Geometry,
}

impl BeamError {
    /// Construct a [`BeamError::BadParameter`] for `name`/`value`.
    ///
    /// Intended for use inside validated constructors so the call site
    /// reads as a guard clause.
    pub fn bad_parameter(name: &'static str, value: f64) -> Self {
        BeamError::BadParameter { name, value }
    }

    /// Validate that `value` is finite and strictly positive, returning
    /// it unchanged on success or a [`BeamError::BadParameter`] tagged
    /// with `name` on failure.
    ///
    /// ```
    /// use valenx_beam::error::BeamError;
    /// assert_eq!(BeamError::require_positive("len", 2.0).unwrap(), 2.0);
    /// assert!(BeamError::require_positive("len", 0.0).is_err());
    /// assert!(BeamError::require_positive("len", -1.0).is_err());
    /// assert!(BeamError::require_positive("len", f64::NAN).is_err());
    /// ```
    pub fn require_positive(name: &'static str, value: f64) -> Result<f64, BeamError> {
        if value.is_finite() && value > 0.0 {
            Ok(value)
        } else {
            Err(BeamError::bad_parameter(name, value))
        }
    }

    /// Stable, kebab-cased identifier for the error variant.
    pub fn code(&self) -> &'static str {
        match self {
            BeamError::BadParameter { .. } => "beam.bad-parameter",
            BeamError::DegenerateSection { .. } => "beam.degenerate-section",
        }
    }

    /// Coarse [`ErrorCategory`] for the error variant.
    pub fn category(&self) -> ErrorCategory {
        match self {
            BeamError::BadParameter { .. } => ErrorCategory::Input,
            BeamError::DegenerateSection { .. } => ErrorCategory::Geometry,
        }
    }
}
