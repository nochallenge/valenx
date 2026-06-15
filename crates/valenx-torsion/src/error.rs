//! Error taxonomy for circular-shaft torsion calculations.
//!
//! Every fallible constructor in this crate returns a [`TorsionError`].
//! The variants distinguish *which* input was rejected and *why*, so a
//! caller (CLI, GUI, or another crate) can surface an actionable message
//! and a stable [`TorsionError::code`] for logging or tests.

use thiserror::Error;

/// Errors raised when validating torsion inputs.
///
/// All physical inputs to the closed-form models must be strictly
/// positive and finite. The geometric variants additionally guard the
/// hollow-shaft annulus (outer diameter strictly greater than the bore)
/// and the radius-within-section invariant used by the shear-stress
/// query.
#[derive(Debug, Error)]
pub enum TorsionError {
    /// A scalar parameter was not strictly positive (or was NaN/∞).
    ///
    /// Diameters, length, shear modulus, applied torque magnitude and
    /// angular speed are all required to be `> 0` and finite. A value of
    /// exactly `0.0` is rejected because it produces a degenerate section
    /// (`J = 0`) or a meaningless query.
    #[error("parameter `{name}` must be finite and strictly positive, got {value}")]
    NonPositive {
        /// Name of the offending parameter (stable, kebab-free identifier).
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A hollow shaft's bore was not strictly smaller than its outer diameter.
    ///
    /// A hollow round section requires `inner_diameter < outer_diameter`;
    /// an equal or inverted pair has no material and yields a non-positive
    /// polar second moment of area.
    #[error(
        "hollow shaft requires inner diameter < outer diameter, got inner = {inner}, outer = {outer}"
    )]
    InvertedAnnulus {
        /// Inner (bore) diameter that was supplied.
        inner: f64,
        /// Outer diameter that was supplied.
        outer: f64,
    },

    /// A radius supplied to a stress query lay outside the cross-section.
    ///
    /// The shear stress varies linearly with radius, so a query radius
    /// must satisfy `inner_radius <= r <= outer_radius`. For a solid
    /// shaft the inner radius is zero.
    #[error("query radius {radius} is outside the section (allowed {min_radius}..={max_radius})")]
    RadiusOutOfRange {
        /// The radius that was queried.
        radius: f64,
        /// Smallest admissible radius (bore radius; zero for solid shafts).
        min_radius: f64,
        /// Largest admissible radius (outer radius).
        max_radius: f64,
    },
}

/// Coarse classification of a [`TorsionError`], useful for routing
/// (e.g. colouring an input field versus reporting an internal bug).
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The user supplied an out-of-domain value (non-positive, NaN, ∞).
    Input,
    /// The geometry is self-inconsistent (e.g. inverted annulus, radius
    /// outside the section).
    Geometry,
}

impl TorsionError {
    /// Stable, dotted identifier for this error.
    ///
    /// Intended for structured logs and assertions; the string is part of
    /// the crate's contract and will not change for a given variant.
    pub fn code(&self) -> &'static str {
        match self {
            TorsionError::NonPositive { .. } => "torsion.non_positive",
            TorsionError::InvertedAnnulus { .. } => "torsion.inverted_annulus",
            TorsionError::RadiusOutOfRange { .. } => "torsion.radius_out_of_range",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            TorsionError::NonPositive { .. } => ErrorCategory::Input,
            TorsionError::InvertedAnnulus { .. } | TorsionError::RadiusOutOfRange { .. } => {
                ErrorCategory::Geometry
            }
        }
    }
}

/// Validate that `value` is finite and strictly positive.
///
/// Returns the value unchanged on success, or [`TorsionError::NonPositive`]
/// naming `name` on failure. This is the single choke-point used by every
/// constructor so the positivity rule is enforced identically everywhere.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64, TorsionError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(TorsionError::NonPositive { name, value })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn require_positive_accepts_positive_finite() {
        assert_eq!(require_positive("d", 3.5).unwrap(), 3.5);
    }

    #[test]
    fn require_positive_rejects_zero_negative_and_nonfinite() {
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let err = require_positive("d", bad).unwrap_err();
            assert!(matches!(err, TorsionError::NonPositive { name: "d", .. }));
        }
    }

    #[test]
    fn codes_and_categories_are_stable() {
        let np = TorsionError::NonPositive {
            name: "d",
            value: 0.0,
        };
        assert_eq!(np.code(), "torsion.non_positive");
        assert_eq!(np.category(), ErrorCategory::Input);

        let inv = TorsionError::InvertedAnnulus {
            inner: 5.0,
            outer: 5.0,
        };
        assert_eq!(inv.code(), "torsion.inverted_annulus");
        assert_eq!(inv.category(), ErrorCategory::Geometry);

        let oor = TorsionError::RadiusOutOfRange {
            radius: 9.0,
            min_radius: 0.0,
            max_radius: 5.0,
        };
        assert_eq!(oor.code(), "torsion.radius_out_of_range");
        assert_eq!(oor.category(), ErrorCategory::Geometry);
    }

    #[test]
    fn display_mentions_parameter_name() {
        let err = TorsionError::NonPositive {
            name: "shear_modulus",
            value: -2.0,
        };
        let text = format!("{err}");
        assert!(text.contains("shear_modulus"), "got: {text}");
    }
}
