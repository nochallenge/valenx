//! Curved-beam error taxonomy.
//!
//! Every public constructor in this crate routes its inputs through the
//! [`CurvedBeamError`] validators below, so a successfully built value is
//! guaranteed finite and physically admissible (positive radii, positive
//! dimensions, inner radius strictly less than outer radius). Downstream
//! formula code can therefore assume clean inputs.

use thiserror::Error;

/// Errors raised while building or evaluating a curved-beam model.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum CurvedBeamError {
    /// A numeric input was NaN or infinite.
    #[error("parameter `{name}` must be finite, got {value}")]
    NotFinite {
        /// Parameter name.
        name: &'static str,
        /// Offending value.
        value: f64,
    },

    /// A quantity that physics requires to be strictly positive was not.
    #[error("parameter `{name}` must be > 0, got {value}")]
    NotPositive {
        /// Parameter name.
        name: &'static str,
        /// Offending value.
        value: f64,
    },

    /// Inner radius was not strictly less than the outer radius.
    #[error("inner radius r_i = {r_i} must be < outer radius r_o = {r_o}")]
    BadRadii {
        /// Inner-fibre radius.
        r_i: f64,
        /// Outer-fibre radius.
        r_o: f64,
    },

    /// A solid-circular section does not fit: the centroidal radius must
    /// exceed the cross-section radius (`R > c`) so the bore stays clear
    /// of the centre of curvature.
    #[error("centroidal radius R = {r_bar} must be > section radius c = {c}")]
    SectionTooLarge {
        /// Centroidal radius.
        r_bar: f64,
        /// Section (circle) radius.
        c: f64,
    },
}

/// Coarse error category for callers that triage by class rather than code.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// Non-finite numeric input (NaN / infinity).
    Domain,
    /// Finite but physically inadmissible input.
    Input,
}

impl CurvedBeamError {
    /// Stable kebab-cased identifier, suitable for logs and tests.
    pub fn code(&self) -> &'static str {
        match self {
            CurvedBeamError::NotFinite { .. } => "curvedbeam.not_finite",
            CurvedBeamError::NotPositive { .. } => "curvedbeam.not_positive",
            CurvedBeamError::BadRadii { .. } => "curvedbeam.bad_radii",
            CurvedBeamError::SectionTooLarge { .. } => "curvedbeam.section_too_large",
        }
    }

    /// Coarse category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            CurvedBeamError::NotFinite { .. } => ErrorCategory::Domain,
            CurvedBeamError::NotPositive { .. }
            | CurvedBeamError::BadRadii { .. }
            | CurvedBeamError::SectionTooLarge { .. } => ErrorCategory::Input,
        }
    }
}

/// Reject NaN / infinite values, returning the finite value on success.
///
/// # Errors
///
/// Returns [`CurvedBeamError::NotFinite`] when `value` is not finite.
pub fn require_finite(name: &'static str, value: f64) -> Result<f64, CurvedBeamError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(CurvedBeamError::NotFinite { name, value })
    }
}

/// Reject non-finite and non-positive values, returning the value on success.
///
/// # Errors
///
/// Returns [`CurvedBeamError::NotFinite`] when `value` is not finite, or
/// [`CurvedBeamError::NotPositive`] when `value <= 0`.
pub fn require_positive(name: &'static str, value: f64) -> Result<f64, CurvedBeamError> {
    let v = require_finite(name, value)?;
    if v > 0.0 {
        Ok(v)
    } else {
        Err(CurvedBeamError::NotPositive { name, value: v })
    }
}
