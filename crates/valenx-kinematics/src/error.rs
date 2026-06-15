//! Kinematics error taxonomy.
//!
//! A single [`KinematicsError`] enum covers the two failure modes the
//! planar models can hit: an out-of-range numeric input (a non-positive
//! link length, a negative rise duration) and a geometric assembly
//! failure (a crank angle at which the four-bar loop cannot close
//! because the dyad triangle is non-realisable). Construction of the
//! validated input types funnels through this enum so a caller never
//! gets a silently wrong answer from a degenerate mechanism.

use thiserror::Error;

/// Errors raised by the planar-kinematics models.
#[derive(Debug, Error)]
pub enum KinematicsError {
    /// A scalar input fell outside its valid domain (e.g. a link
    /// length that is not strictly positive, or a rise angle/lift that
    /// is negative).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Name of the offending parameter.
        name: &'static str,
        /// Human-readable explanation of why it was rejected.
        reason: String,
    },

    /// The four-bar linkage cannot be assembled at the requested crank
    /// angle: the coupler/rocker dyad would have to form a triangle
    /// whose closing side is longer than the sum (or shorter than the
    /// absolute difference) of the other two, so no real solution
    /// exists. The diagonal length is reported for diagnosis.
    #[error(
        "linkage cannot close at crank angle {crank_rad} rad: \
         diagonal {diagonal} outside reachable range \
         [{reach_min}, {reach_max}]"
    )]
    CannotClose {
        /// Crank angle (radians) at which closure failed.
        crank_rad: f64,
        /// Length of the ground-to-coupler-pin diagonal.
        diagonal: f64,
        /// Minimum reachable diagonal (`|coupler - rocker|`).
        reach_min: f64,
        /// Maximum reachable diagonal (`coupler + rocker`).
        reach_max: f64,
    },
}

/// Coarse category for routing/telemetry, mirroring the sibling
/// workbench crates.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// Caller-supplied input was invalid.
    Input,
    /// The algorithm hit a domain limit (non-realisable geometry).
    Algorithm,
}

impl KinematicsError {
    /// Stable kebab-cased identifier, handy for logs and tests.
    pub fn code(&self) -> &'static str {
        match self {
            KinematicsError::BadParameter { .. } => "kinematics.bad_parameter",
            KinematicsError::CannotClose { .. } => "kinematics.cannot_close",
        }
    }

    /// Coarse category for the error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            KinematicsError::BadParameter { .. } => ErrorCategory::Input,
            KinematicsError::CannotClose { .. } => ErrorCategory::Algorithm,
        }
    }
}

/// Internal helper: reject a value that must be strictly positive.
///
/// The `is_finite` guard runs first so a `NaN`/`±∞` input is rejected
/// before the magnitude check (a bare `value <= 0.0` would let `NaN`
/// through, since every `NaN` comparison is `false`).
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<(), KinematicsError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(KinematicsError::BadParameter {
            name,
            reason: format!("must be a finite value > 0, got {value}"),
        });
    }
    Ok(())
}

/// Internal helper: reject a value that must be non-negative.
///
/// Same `is_finite`-first ordering as [`require_positive`] so a
/// non-finite input cannot slip past the `< 0.0` magnitude check.
pub(crate) fn require_non_negative(name: &'static str, value: f64) -> Result<(), KinematicsError> {
    if !value.is_finite() || value < 0.0 {
        return Err(KinematicsError::BadParameter {
            name,
            reason: format!("must be a finite value >= 0, got {value}"),
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_and_categories_are_stable() {
        let bad = KinematicsError::BadParameter {
            name: "crank",
            reason: "x".into(),
        };
        assert_eq!(bad.code(), "kinematics.bad_parameter");
        assert_eq!(bad.category(), ErrorCategory::Input);

        let close = KinematicsError::CannotClose {
            crank_rad: 0.0,
            diagonal: 10.0,
            reach_min: 1.0,
            reach_max: 5.0,
        };
        assert_eq!(close.code(), "kinematics.cannot_close");
        assert_eq!(close.category(), ErrorCategory::Algorithm);
    }

    #[test]
    fn require_positive_rejects_zero_negative_and_nonfinite() {
        assert!(require_positive("l", 1.0).is_ok());
        assert!(require_positive("l", 0.0).is_err());
        assert!(require_positive("l", -2.0).is_err());
        assert!(require_positive("l", f64::NAN).is_err());
        assert!(require_positive("l", f64::INFINITY).is_err());
    }

    #[test]
    fn require_non_negative_allows_zero_rejects_negative() {
        assert!(require_non_negative("lift", 0.0).is_ok());
        assert!(require_non_negative("lift", 3.5).is_ok());
        assert!(require_non_negative("lift", -0.1).is_err());
        assert!(require_non_negative("lift", f64::NAN).is_err());
    }
}
