//! Error taxonomy for the thick-walled-cylinder solver.
//!
//! All public entry points validate their arguments up front and return
//! [`ThickCylinderError`] rather than panicking or silently producing
//! NaN/Inf. The constructors here centralise the domain checks so every
//! module rejects the same bad input the same way.

use thiserror::Error;

/// Errors raised when constructing or evaluating a thick-walled cylinder.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum ThickCylinderError {
    /// A supplied value was not finite (NaN or +/-infinity).
    #[error("parameter `{name}` must be finite, got {value}")]
    NotFinite {
        /// Name of the offending parameter.
        name: &'static str,
        /// The non-finite value that was supplied.
        value: f64,
    },

    /// A radius (or wall thickness) was zero or negative where physics
    /// demands a strictly positive length.
    #[error("parameter `{name}` must be > 0, got {value}")]
    NonPositive {
        /// Name of the offending parameter.
        name: &'static str,
        /// The non-positive value that was supplied.
        value: f64,
    },

    /// The outer radius was not strictly greater than the inner radius,
    /// so the wall has zero or negative thickness.
    #[error("outer radius ({outer}) must be > inner radius ({inner})")]
    BadRadii {
        /// Inner radius `a`.
        inner: f64,
        /// Outer radius `b`.
        outer: f64,
    },

    /// An evaluation radius `r` fell outside the closed wall interval
    /// `[a, b]`, where the Lame solution is not defined.
    #[error("evaluation radius {r} is outside the wall [{inner}, {outer}]")]
    OutsideWall {
        /// Inner radius `a`.
        inner: f64,
        /// Outer radius `b`.
        outer: f64,
        /// The requested evaluation radius.
        r: f64,
    },
}

impl ThickCylinderError {
    /// Stable, kebab-cased identifier for logging / matching in tests.
    ///
    /// The string is part of the public contract and will not change for
    /// a given variant.
    pub fn code(&self) -> &'static str {
        match self {
            ThickCylinderError::NotFinite { .. } => "thickcylinder.not-finite",
            ThickCylinderError::NonPositive { .. } => "thickcylinder.non-positive",
            ThickCylinderError::BadRadii { .. } => "thickcylinder.bad-radii",
            ThickCylinderError::OutsideWall { .. } => "thickcylinder.outside-wall",
        }
    }
}

/// Reject a non-finite scalar.
///
/// Returns the value unchanged when it is finite, otherwise
/// [`ThickCylinderError::NotFinite`].
pub fn finite(name: &'static str, value: f64) -> Result<f64, ThickCylinderError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ThickCylinderError::NotFinite { name, value })
    }
}

/// Require a finite, strictly positive scalar (a length / pressure that
/// physics forbids from being zero or negative).
pub fn positive(name: &'static str, value: f64) -> Result<f64, ThickCylinderError> {
    let value = finite(name, value)?;
    if value > 0.0 {
        Ok(value)
    } else {
        Err(ThickCylinderError::NonPositive { name, value })
    }
}

/// Require a finite scalar that may be zero or positive but not negative
/// (e.g. an external pressure, which may legitimately be zero).
pub fn non_negative(name: &'static str, value: f64) -> Result<f64, ThickCylinderError> {
    let value = finite(name, value)?;
    if value >= 0.0 {
        Ok(value)
    } else {
        Err(ThickCylinderError::NonPositive { name, value })
    }
}
