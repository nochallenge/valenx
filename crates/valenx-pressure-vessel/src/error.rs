//! Error taxonomy for `valenx-pressure-vessel`.
//!
//! Every fallible constructor in this crate returns
//! [`Result<_, VesselError>`]. The pressure-vessel formulae are total
//! functions of their inputs *once the geometry is physically valid*, so
//! the only ways a call can fail are:
//!
//! 1. A dimension or pressure is non-finite, zero, or negative when the
//!    model requires it to be strictly positive
//!    ([`VesselError::BadParameter`]).
//! 2. The geometry is inconsistent — an outer radius that is not strictly
//!    larger than the inner radius, or an evaluation radius outside the
//!    `[r_inner, r_outer]` wall ([`VesselError::Geometry`]).
//!
//! Use [`VesselError::code`] for stable log / telemetry tagging and
//! [`VesselError::category`] to bucket failures without matching every
//! variant. The pattern mirrors `valenx-springs`' `SpringsError`.

use thiserror::Error;

/// Errors raised by pressure-vessel construction and evaluation.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum VesselError {
    /// A scalar parameter is outside its valid domain — non-finite
    /// (`NaN` / infinity), or not strictly positive where the model
    /// demands it (a wall thickness, radius, or pressure of zero or
    /// less). A property of the individual *argument*.
    #[error("bad parameter `{name}`: {reason} (got {value})")]
    BadParameter {
        /// Logical parameter name (e.g. `"thickness"`, `"r_outer"`).
        name: &'static str,
        /// Human-readable reason the value is rejected.
        reason: &'static str,
        /// The offending value, echoed for diagnostics.
        value: f64,
    },

    /// The geometry is internally inconsistent: the outer radius is not
    /// strictly greater than the inner radius, or an evaluation radius
    /// falls outside the closed wall interval `[r_inner, r_outer]`. A
    /// property of how several arguments *relate*, not of one in
    /// isolation.
    #[error("inconsistent geometry: {reason}")]
    Geometry {
        /// Human-readable reason.
        reason: String,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this rather
/// than on the full variant set.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A single supplied argument is out of range.
    Input,
    /// Several arguments are mutually inconsistent (geometry).
    Geometry,
}

impl VesselError {
    /// Stable kebab/dotted error code for log / telemetry tagging.
    /// Format: `"pressure_vessel.<sub_id>"`; codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            VesselError::BadParameter { .. } => "pressure_vessel.bad_parameter",
            VesselError::Geometry { .. } => "pressure_vessel.geometry",
        }
    }

    /// Coarse [`ErrorCategory`] for routing without matching every
    /// variant.
    pub fn category(&self) -> ErrorCategory {
        match self {
            VesselError::BadParameter { .. } => ErrorCategory::Input,
            VesselError::Geometry { .. } => ErrorCategory::Geometry,
        }
    }

    /// Construct a [`VesselError::Geometry`] from any string-like reason.
    pub fn geometry(reason: impl Into<String>) -> Self {
        VesselError::Geometry {
            reason: reason.into(),
        }
    }

    /// Validate that `value` is finite and strictly positive, returning
    /// it on success or a [`VesselError::BadParameter`] otherwise.
    ///
    /// This is the single gate every public constructor routes radii,
    /// thicknesses and pressures through, so the rejection message is
    /// identical everywhere.
    pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64> {
        if !value.is_finite() {
            return Err(VesselError::BadParameter {
                name,
                reason: "must be a finite number",
                value,
            });
        }
        if value <= 0.0 {
            return Err(VesselError::BadParameter {
                name,
                reason: "must be strictly positive",
                value,
            });
        }
        Ok(value)
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, VesselError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = VesselError::require_positive("thickness", -1.0).unwrap_err();
        assert_eq!(err.code(), "pressure_vessel.bad_parameter");
        assert_eq!(err.category(), ErrorCategory::Input);

        let err = VesselError::geometry("r_outer <= r_inner");
        assert_eq!(err.code(), "pressure_vessel.geometry");
        assert_eq!(err.category(), ErrorCategory::Geometry);
    }

    #[test]
    fn require_positive_accepts_and_rejects() {
        assert!(VesselError::require_positive("p", 2.5).is_ok());

        // Zero is rejected.
        let err = VesselError::require_positive("p", 0.0).unwrap_err();
        assert!(matches!(err, VesselError::BadParameter { .. }));

        // NaN and infinity are rejected with the finiteness message.
        let err = VesselError::require_positive("p", f64::NAN).unwrap_err();
        match err {
            VesselError::BadParameter { reason, .. } => {
                assert!(reason.contains("finite"), "got reason: {reason}");
            }
            other => panic!("expected BadParameter, got {other:?}"),
        }
        assert!(VesselError::require_positive("p", f64::INFINITY).is_err());
    }

    #[test]
    fn display_is_informative() {
        let msg = VesselError::require_positive("r_outer", -3.0)
            .unwrap_err()
            .to_string();
        assert!(msg.contains("r_outer"), "got: {msg}");
        assert!(msg.contains("-3"), "got: {msg}");

        let msg = VesselError::geometry("radius 9 outside wall").to_string();
        assert!(msg.contains("radius 9 outside wall"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(VesselError::geometry("x"));
        assert!(err.to_string().contains('x'));
    }
}
