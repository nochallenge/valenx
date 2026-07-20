//! Error taxonomy for `valenx-radiation`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, RadiationError>`]. The variants are intentionally
//! coarse — a radiation-heat-transfer caller usually only cares about
//! three failure modes:
//!
//! 1. A value is not finite — a `NaN` or an infinity slipped in
//!    ([`RadiationError::NotFinite`]).
//! 2. A value is outside its physical domain — a negative absolute
//!    temperature, a non-positive area, an emissivity outside `[0, 1]`,
//!    a view factor outside `[0, 1]` ([`RadiationError::OutOfDomain`]).
//! 3. A modelling invariant is violated — e.g. reciprocity would force
//!    a non-physical complementary area ([`RadiationError::Inconsistent`]).
//!
//! Use [`RadiationError::code`] for stable log / telemetry tagging.
//! The pattern mirrors `valenx-popgen`'s `PopgenError` and the other
//! sim crates' error enums.

use thiserror::Error;

/// Errors produced by `valenx-radiation`.
///
/// All constructors are *validating*: they reject non-finite inputs and
/// values outside their physical domain, so a successfully constructed
/// value (or a function that took validated inputs) is always in range.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum RadiationError {
    /// A supplied value was `NaN` or `±∞`. `what` names the logical
    /// parameter (e.g. `"temperature_k"`).
    #[error("value `{what}` is not finite ({value})")]
    NotFinite {
        /// Logical parameter name.
        what: &'static str,
        /// The offending value, formatted for the message.
        value: f64,
    },

    /// A finite value fell outside its allowed physical domain — a
    /// negative absolute temperature, a non-positive area, an emissivity
    /// or view factor outside `[0, 1]`, etc.
    #[error("value `{what}` = {value} is out of domain: {reason}")]
    OutOfDomain {
        /// Logical parameter name.
        what: &'static str,
        /// The offending value.
        value: f64,
        /// Human-readable statement of the allowed range.
        reason: &'static str,
    },

    /// A relationship between several already-valid inputs is
    /// physically inconsistent — e.g. a requested reciprocity that would
    /// imply a non-positive complementary area.
    #[error("inconsistent radiation model: {reason}")]
    Inconsistent {
        /// Human-readable reason.
        reason: String,
    },
}

/// Coarse error category for routing / display, stable across versions.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A value was non-finite or out of its physical domain — i.e. the
    /// caller passed bad input.
    Input,
    /// The model's own invariants were violated.
    Model,
}

impl RadiationError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"radiation.<sub_id>"`. Codes never change
    /// across minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            RadiationError::NotFinite { .. } => "radiation.not_finite",
            RadiationError::OutOfDomain { .. } => "radiation.out_of_domain",
            RadiationError::Inconsistent { .. } => "radiation.inconsistent",
        }
    }

    /// Coarse [`ErrorCategory`] for callers that want to bucket failures
    /// without matching every variant.
    pub fn category(&self) -> ErrorCategory {
        match self {
            RadiationError::NotFinite { .. } | RadiationError::OutOfDomain { .. } => {
                ErrorCategory::Input
            }
            RadiationError::Inconsistent { .. } => ErrorCategory::Model,
        }
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, RadiationError>;

/// Reject a non-finite value, returning it unchanged when finite.
///
/// Used internally by the domain-checking helpers; exposed so callers
/// can validate raw inputs with the same rule the crate uses.
pub fn finite(what: &'static str, value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(RadiationError::NotFinite { what, value })
    }
}

/// Validate an absolute temperature in kelvin: finite and `>= 0`.
///
/// Absolute zero is permitted (it makes the emissive power exactly
/// zero); negative absolute temperatures are rejected.
pub fn check_temperature(what: &'static str, t_k: f64) -> Result<f64> {
    let t = finite(what, t_k)?;
    if t < 0.0 {
        Err(RadiationError::OutOfDomain {
            what,
            value: t,
            reason: "absolute temperature must be >= 0 K",
        })
    } else {
        Ok(t)
    }
}

/// Validate an area: finite and strictly positive.
pub fn check_area(what: &'static str, area: f64) -> Result<f64> {
    let a = finite(what, area)?;
    if a <= 0.0 {
        Err(RadiationError::OutOfDomain {
            what,
            value: a,
            reason: "area must be > 0",
        })
    } else {
        Ok(a)
    }
}

/// Validate a dimensionless quantity that must lie in the closed unit
/// interval `[0, 1]` — emissivity, absorptivity, or a view factor.
pub fn check_unit_interval(what: &'static str, x: f64) -> Result<f64> {
    let v = finite(what, x)?;
    if !(0.0..=1.0).contains(&v) {
        Err(RadiationError::OutOfDomain {
            what,
            value: v,
            reason: "must lie in [0, 1]",
        })
    } else {
        Ok(v)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let e = RadiationError::NotFinite {
            what: "temperature_k",
            value: f64::NAN,
        };
        assert_eq!(e.code(), "radiation.not_finite");
        assert_eq!(e.category(), ErrorCategory::Input);

        let e = RadiationError::OutOfDomain {
            what: "emissivity",
            value: 1.5,
            reason: "must lie in [0, 1]",
        };
        assert_eq!(e.code(), "radiation.out_of_domain");
        assert_eq!(e.category(), ErrorCategory::Input);

        let e = RadiationError::Inconsistent {
            reason: "bad".into(),
        };
        assert_eq!(e.code(), "radiation.inconsistent");
        assert_eq!(e.category(), ErrorCategory::Model);
    }

    #[test]
    fn finite_rejects_nan_and_inf() {
        assert!(finite("x", f64::NAN).is_err());
        assert!(finite("x", f64::INFINITY).is_err());
        assert!(finite("x", -f64::INFINITY).is_err());
        assert_eq!(finite("x", 3.0).unwrap(), 3.0);
    }

    #[test]
    fn temperature_domain() {
        assert_eq!(check_temperature("t", 0.0).unwrap(), 0.0);
        assert_eq!(check_temperature("t", 300.0).unwrap(), 300.0);
        let e = check_temperature("t", -1.0).unwrap_err();
        assert_eq!(e.code(), "radiation.out_of_domain");
        assert!(check_temperature("t", f64::NAN).is_err());
    }

    #[test]
    fn area_domain() {
        assert_eq!(check_area("a", 2.0).unwrap(), 2.0);
        assert!(check_area("a", 0.0).is_err());
        assert!(check_area("a", -1.0).is_err());
    }

    #[test]
    fn unit_interval_domain() {
        assert_eq!(check_unit_interval("e", 0.0).unwrap(), 0.0);
        assert_eq!(check_unit_interval("e", 1.0).unwrap(), 1.0);
        assert_eq!(check_unit_interval("e", 0.5).unwrap(), 0.5);
        assert!(check_unit_interval("e", -0.01).is_err());
        assert!(check_unit_interval("e", 1.01).is_err());
    }

    #[test]
    fn error_is_std_error() {
        let err: Box<dyn std::error::Error> = Box::new(RadiationError::Inconsistent {
            reason: "boom".into(),
        });
        assert!(err.to_string().contains("boom"));
    }
}
