//! Power-screw error taxonomy.
//!
//! Every fallible constructor in this crate returns
//! [`ScrewThreadError`]. Construction is validated up front so that
//! the analytic routines in [`crate::screw`] can assume a physically
//! meaningful geometry (positive diameters, a non-negative friction
//! coefficient, at least one thread start, and a denominator that does
//! not pass through zero in the raising-torque expression).

use thiserror::Error;

/// Errors raised while building or analysing a power screw.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum ScrewThreadError {
    /// A scalar parameter fell outside its physically admissible range
    /// (for example a non-positive diameter or pitch, or a negative
    /// friction coefficient).
    #[error("invalid parameter `{name}`: {reason} (got {value})")]
    InvalidParameter {
        /// Name of the offending parameter.
        name: &'static str,
        /// Why the value is rejected.
        reason: &'static str,
        /// The rejected value, echoed back for diagnostics.
        value: f64,
    },

    /// The raising-torque denominator `pi*dm - mu*l` is zero (or
    /// numerically indistinguishable from zero), which makes the
    /// raise-torque expression singular. Physically this means the
    /// lead is so large relative to the mean diameter, scaled by the
    /// friction coefficient, that the screw model breaks down.
    #[error(
        "singular raise-torque denominator: pi*dm - mu*l = {denominator} \
         (dm = {mean_diameter}, lead = {lead}, mu = {friction})"
    )]
    SingularDenominator {
        /// Computed value of `pi*dm - mu*l`.
        denominator: f64,
        /// Mean (pitch-line) diameter that produced it.
        mean_diameter: f64,
        /// Lead that produced it.
        lead: f64,
        /// Coefficient of friction that produced it.
        friction: f64,
    },
}

/// Coarse classification of a [`ScrewThreadError`], handy for routing
/// the error to a UI channel (user-input banner vs. internal log).
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The caller supplied an out-of-range input value.
    Input,
    /// The geometry is admissible per-parameter but drives an analytic
    /// expression into a singularity.
    Domain,
}

impl ScrewThreadError {
    /// Stable, kebab-cased identifier suitable for logs and tests.
    ///
    /// The string is part of the crate's public contract and will not
    /// change for a given variant across patch releases.
    pub fn code(&self) -> &'static str {
        match self {
            ScrewThreadError::InvalidParameter { .. } => "screwthread.invalid_parameter",
            ScrewThreadError::SingularDenominator { .. } => "screwthread.singular_denominator",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            ScrewThreadError::InvalidParameter { .. } => ErrorCategory::Input,
            ScrewThreadError::SingularDenominator { .. } => ErrorCategory::Domain,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_are_stable_and_distinct() {
        let a = ScrewThreadError::InvalidParameter {
            name: "pitch",
            reason: "must be positive",
            value: -1.0,
        };
        let b = ScrewThreadError::SingularDenominator {
            denominator: 0.0,
            mean_diameter: 1.0,
            lead: 1.0,
            friction: 1.0,
        };
        assert_eq!(a.code(), "screwthread.invalid_parameter");
        assert_eq!(b.code(), "screwthread.singular_denominator");
        assert_ne!(a.code(), b.code());
    }

    #[test]
    fn categories_route_correctly() {
        let a = ScrewThreadError::InvalidParameter {
            name: "pitch",
            reason: "must be positive",
            value: -1.0,
        };
        let b = ScrewThreadError::SingularDenominator {
            denominator: 0.0,
            mean_diameter: 1.0,
            lead: 1.0,
            friction: 1.0,
        };
        assert_eq!(a.category(), ErrorCategory::Input);
        assert_eq!(b.category(), ErrorCategory::Domain);
    }

    #[test]
    fn display_includes_value() {
        let a = ScrewThreadError::InvalidParameter {
            name: "mean_diameter_mm",
            reason: "must be positive",
            value: -3.0,
        };
        let text = a.to_string();
        assert!(text.contains("mean_diameter_mm"));
        assert!(text.contains("-3"));
    }
}
