//! Error taxonomy for `valenx-fluid-statics`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, FluidStaticsError>`]. The variants are intentionally
//! coarse — a hydrostatics caller usually only cares about three things:
//!
//! 1. Did the caller pass a physically impossible argument — a negative
//!    density, a non-positive `g`, an out-of-range geometry
//!    ([`FluidStaticsError::Invalid`])?
//! 2. Do two inputs disagree on a constraint — a vertical extent that
//!    pokes above the free surface, a depth shallower than the top edge
//!    ([`FluidStaticsError::Geometry`])?
//! 3. Would the requested quantity divide by zero — a centre of pressure
//!    asked for a plate whose centroid sits exactly on the free surface
//!    ([`FluidStaticsError::Singular`])?
//!
//! Use [`FluidStaticsError::code`] for stable log / telemetry tagging
//! and [`FluidStaticsError::category`] to bucket failures without
//! matching every variant. The pattern mirrors `valenx-popgen`'s
//! `PopgenError` and `valenx-astro`'s error type.

use thiserror::Error;

/// Errors produced by `valenx-fluid-statics`.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum FluidStaticsError {
    /// Caller passed an argument the model cannot accept: a negative
    /// density, a non-positive gravitational acceleration, a negative
    /// depth or area, a non-finite value. A property of a single
    /// argument, independent of any other input.
    #[error("invalid `{what}`: {reason}")]
    Invalid {
        /// Logical parameter name (e.g. `"density"`, `"gravity"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// Two inputs are individually valid but jointly inconsistent — a
    /// plate that extends above the free surface, a submerged depth
    /// shallower than half the plate height, a manometer column whose
    /// legs disagree. A property of the *combination* of arguments.
    #[error("inconsistent geometry for {context}: {reason}")]
    Geometry {
        /// Short context label (e.g. `"submerged plate"`).
        context: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// The requested quantity is mathematically undefined for the given
    /// inputs — typically a division by zero, such as a centre of
    /// pressure for a plate whose centroidal depth is zero (its centroid
    /// lies on the free surface, so the resultant force is zero and its
    /// line of action is indeterminate).
    #[error("undefined result for {context}: {reason}")]
    Singular {
        /// Short context label (e.g. `"centre of pressure"`).
        context: &'static str,
        /// Human-readable reason.
        reason: String,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this rather
/// than on every error variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A single argument is out of its physical domain.
    Input,
    /// Otherwise-valid arguments are jointly inconsistent.
    Geometry,
    /// The requested quantity is mathematically undefined.
    Singular,
}

impl FluidStaticsError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"fluid_statics.<sub_id>"`. Codes never change
    /// across minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            FluidStaticsError::Invalid { .. } => "fluid_statics.invalid",
            FluidStaticsError::Geometry { .. } => "fluid_statics.geometry",
            FluidStaticsError::Singular { .. } => "fluid_statics.singular",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            FluidStaticsError::Invalid { .. } => "input",
            FluidStaticsError::Geometry { .. } => "geometry",
            FluidStaticsError::Singular { .. } => "singular",
        }
    }

    /// Typed category enum (for callers that want to `match` instead of
    /// comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            FluidStaticsError::Invalid { .. } => ErrorCategory::Input,
            FluidStaticsError::Geometry { .. } => ErrorCategory::Geometry,
            FluidStaticsError::Singular { .. } => ErrorCategory::Singular,
        }
    }

    /// Convenience constructor for [`FluidStaticsError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        FluidStaticsError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`FluidStaticsError::Geometry`].
    pub fn geometry(context: &'static str, reason: impl Into<String>) -> Self {
        FluidStaticsError::Geometry {
            context,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`FluidStaticsError::Singular`].
    pub fn singular(context: &'static str, reason: impl Into<String>) -> Self {
        FluidStaticsError::Singular {
            context,
            reason: reason.into(),
        }
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, FluidStaticsError>;

// --- Shared validation helpers ---------------------------------------

/// Validate that a value is finite (neither `NaN` nor infinite),
/// returning [`FluidStaticsError::Invalid`] otherwise.
pub(crate) fn require_finite(what: &'static str, value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(FluidStaticsError::invalid(
            what,
            format!("must be finite, got {value}"),
        ))
    }
}

/// Validate that a value is finite and strictly positive.
pub(crate) fn require_positive(what: &'static str, value: f64) -> Result<f64> {
    let value = require_finite(what, value)?;
    if value > 0.0 {
        Ok(value)
    } else {
        Err(FluidStaticsError::invalid(
            what,
            format!("must be strictly positive, got {value}"),
        ))
    }
}

/// Validate that a value is finite and non-negative.
pub(crate) fn require_non_negative(what: &'static str, value: f64) -> Result<f64> {
    let value = require_finite(what, value)?;
    if value >= 0.0 {
        Ok(value)
    } else {
        Err(FluidStaticsError::invalid(
            what,
            format!("must be non-negative, got {value}"),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = FluidStaticsError::invalid("density", "must be positive");
        assert_eq!(err.code(), "fluid_statics.invalid");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = FluidStaticsError::geometry("submerged plate", "above surface");
        assert_eq!(err.code(), "fluid_statics.geometry");
        assert_eq!(err.category(), "geometry");
        assert_eq!(err.category_enum(), ErrorCategory::Geometry);

        let err = FluidStaticsError::singular("centre of pressure", "zero depth");
        assert_eq!(err.code(), "fluid_statics.singular");
        assert_eq!(err.category(), "singular");
        assert_eq!(err.category_enum(), ErrorCategory::Singular);
    }

    #[test]
    fn display_is_informative() {
        let msg = FluidStaticsError::invalid("gravity", "must be positive").to_string();
        assert!(msg.contains("gravity"), "got: {msg}");
        assert!(msg.contains("must be positive"), "got: {msg}");

        let msg = FluidStaticsError::geometry("plate", "extends above surface").to_string();
        assert!(msg.contains("plate"), "got: {msg}");
        assert!(msg.contains("extends above surface"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> =
            Box::new(FluidStaticsError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }

    #[test]
    fn require_positive_rejects_zero_and_negative() {
        assert!(require_positive("g", 9.81).is_ok());
        assert!(require_positive("g", 0.0).is_err());
        assert!(require_positive("g", -1.0).is_err());
        assert!(require_positive("g", f64::NAN).is_err());
        assert!(require_positive("g", f64::INFINITY).is_err());
    }

    #[test]
    fn require_non_negative_accepts_zero_rejects_negative() {
        assert!(require_non_negative("h", 0.0).is_ok());
        assert!(require_non_negative("h", 2.5).is_ok());
        assert!(require_non_negative("h", -0.1).is_err());
        assert!(require_non_negative("h", f64::NAN).is_err());
    }
}
