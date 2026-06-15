//! Hydraulics error taxonomy.
//!
//! A single [`HydraulicsError`] enum covers the failure modes of the
//! whole crate. Construct it through the validated helpers
//! ([`HydraulicsError::non_positive`], [`HydraulicsError::negative`],
//! and [`HydraulicsError::geometry`]) rather than building the
//! variants by hand, so the messages stay uniform.

use thiserror::Error;

/// Errors raised while validating hydraulic inputs.
#[derive(Debug, Error)]
pub enum HydraulicsError {
    /// A quantity that must be strictly positive was zero or negative.
    #[error("parameter `{name}` must be > 0, got {value}")]
    NonPositive {
        /// Offending parameter name.
        name: &'static str,
        /// Offending value.
        value: f64,
    },

    /// A quantity that must be non-negative was negative.
    #[error("parameter `{name}` must be >= 0, got {value}")]
    Negative {
        /// Offending parameter name.
        name: &'static str,
        /// Offending value.
        value: f64,
    },

    /// The geometry is physically impossible (e.g. a rod at least as
    /// thick as the bore leaves no annulus area on the rod side).
    #[error("invalid geometry: {0}")]
    Geometry(String),
}

/// Coarse bucket a [`HydraulicsError`] falls into, for callers that
/// want to branch on the kind of failure without matching every
/// variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// A scalar input was out of its allowed range.
    Input,
    /// The combination of inputs is geometrically inconsistent.
    Geometry,
}

impl HydraulicsError {
    /// Build a [`HydraulicsError::NonPositive`] for `name`/`value`.
    ///
    /// Use this when a parameter (a pressure, area, diameter, flow
    /// coefficient, ...) must be strictly greater than zero.
    pub fn non_positive(name: &'static str, value: f64) -> Self {
        HydraulicsError::NonPositive { name, value }
    }

    /// Build a [`HydraulicsError::Negative`] for `name`/`value`.
    ///
    /// Use this when a parameter may be zero but never negative (for
    /// instance a pressure drop across a fully closed valve).
    pub fn negative(name: &'static str, value: f64) -> Self {
        HydraulicsError::Negative { name, value }
    }

    /// Build a [`HydraulicsError::Geometry`] from a free-form message.
    pub fn geometry(reason: impl Into<String>) -> Self {
        HydraulicsError::Geometry(reason.into())
    }

    /// Validate that `value` is finite and strictly positive,
    /// returning it unchanged or a [`HydraulicsError::NonPositive`].
    ///
    /// `NaN` and infinities are rejected as non-positive so downstream
    /// arithmetic never propagates them.
    pub fn require_positive(name: &'static str, value: f64) -> Result<f64, Self> {
        if value.is_finite() && value > 0.0 {
            Ok(value)
        } else {
            Err(Self::non_positive(name, value))
        }
    }

    /// Validate that `value` is finite and non-negative, returning it
    /// unchanged or a [`HydraulicsError::Negative`].
    pub fn require_non_negative(name: &'static str, value: f64) -> Result<f64, Self> {
        if value.is_finite() && value >= 0.0 {
            Ok(value)
        } else {
            Err(Self::negative(name, value))
        }
    }

    /// Stable kebab-cased identifier, handy for logs and tests.
    pub fn code(&self) -> &'static str {
        match self {
            HydraulicsError::NonPositive { .. } => "hydraulics.non_positive",
            HydraulicsError::Negative { .. } => "hydraulics.negative",
            HydraulicsError::Geometry(_) => "hydraulics.geometry",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            HydraulicsError::NonPositive { .. } | HydraulicsError::Negative { .. } => {
                ErrorCategory::Input
            }
            HydraulicsError::Geometry(_) => ErrorCategory::Geometry,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn require_positive_accepts_positive() {
        assert_eq!(HydraulicsError::require_positive("p", 3.0).unwrap(), 3.0);
    }

    #[test]
    fn require_positive_rejects_zero_negative_and_nonfinite() {
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let err = HydraulicsError::require_positive("p", bad).unwrap_err();
            assert_eq!(err.code(), "hydraulics.non_positive");
            assert_eq!(err.category(), ErrorCategory::Input);
        }
    }

    #[test]
    fn require_non_negative_accepts_zero_but_not_negative() {
        assert_eq!(
            HydraulicsError::require_non_negative("dp", 0.0).unwrap(),
            0.0
        );
        let err = HydraulicsError::require_non_negative("dp", -0.5).unwrap_err();
        assert_eq!(err.code(), "hydraulics.negative");
        assert_eq!(err.category(), ErrorCategory::Input);
    }

    #[test]
    fn geometry_category_and_code() {
        let err = HydraulicsError::geometry("rod thicker than bore");
        assert_eq!(err.code(), "hydraulics.geometry");
        assert_eq!(err.category(), ErrorCategory::Geometry);
    }

    #[test]
    fn display_includes_name_and_value() {
        let msg = HydraulicsError::non_positive("bore_diameter", -2.0).to_string();
        assert!(msg.contains("bore_diameter"));
        assert!(msg.contains("-2"));
    }
}
