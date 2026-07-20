//! Error taxonomy for `valenx-buoyancy`.
//!
//! Every fallible public function in this crate returns
//! [`Result<T, BuoyancyError>`]. Two concerns dominate hydrostatic
//! input validation, so the enum keeps to two variants:
//!
//! 1. A value is not a finite number (`NaN`, `+inf`, `-inf`) — a
//!    [`BuoyancyError::NotFinite`]. Non-finite inputs poison every
//!    downstream comparison, so they are rejected at the door.
//! 2. A value is finite but outside the physical domain — a density,
//!    volume, gravity or dimension that must be strictly positive but
//!    is zero or negative ([`BuoyancyError::OutOfDomain`]).
//!
//! Construct errors through the validated helpers
//! [`BuoyancyError::finite`] and [`BuoyancyError::positive`] rather than
//! building variants by hand: they perform the check and return
//! `Ok(value)` on success, so call sites read as
//! `let rho = BuoyancyError::positive("rho_fluid", rho_fluid)?;`.
//!
//! Use [`BuoyancyError::code`] for stable log / telemetry tagging — the
//! returned string never changes across minor versions.

use thiserror::Error;

/// Errors produced by `valenx-buoyancy`.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum BuoyancyError {
    /// A supplied quantity was not a finite number (`NaN` or infinity).
    /// `name` is the logical parameter (e.g. `"rho_fluid"`); `value` is
    /// the offending float, surfaced verbatim for debugging.
    #[error("parameter `{name}` must be finite, got {value}")]
    NotFinite {
        /// Logical parameter name.
        name: &'static str,
        /// The non-finite value that was supplied.
        value: f64,
    },

    /// A quantity was finite but outside its physical domain — a
    /// density, volume, acceleration or length that the model requires
    /// to be strictly positive but which was zero or negative.
    #[error("parameter `{name}` is out of domain: {reason} (got {value})")]
    OutOfDomain {
        /// Logical parameter name.
        name: &'static str,
        /// The offending value.
        value: f64,
        /// Human-readable statement of the domain that was violated.
        reason: &'static str,
    },
}

/// Coarse category for routing / display. Stable across crate versions:
/// match on this rather than on the individual error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A value was `NaN` or infinite.
    NotFinite,
    /// A finite value violated a physical domain constraint.
    Domain,
}

impl BuoyancyError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"buoyancy.<sub_id>"`. Codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            BuoyancyError::NotFinite { .. } => "buoyancy.not_finite",
            BuoyancyError::OutOfDomain { .. } => "buoyancy.out_of_domain",
        }
    }

    /// Coarse [`ErrorCategory`] for callers that want to `match` on a
    /// small stable enum instead of every variant.
    pub fn category(&self) -> ErrorCategory {
        match self {
            BuoyancyError::NotFinite { .. } => ErrorCategory::NotFinite,
            BuoyancyError::OutOfDomain { .. } => ErrorCategory::Domain,
        }
    }

    /// Validate that `value` is finite, returning it unchanged on
    /// success. Rejects `NaN`, `+inf` and `-inf` with
    /// [`BuoyancyError::NotFinite`].
    ///
    /// ```
    /// use valenx_buoyancy::BuoyancyError;
    /// assert_eq!(BuoyancyError::finite("g", 9.81).unwrap(), 9.81);
    /// assert!(BuoyancyError::finite("g", f64::NAN).is_err());
    /// ```
    pub fn finite(name: &'static str, value: f64) -> std::result::Result<f64, BuoyancyError> {
        if value.is_finite() {
            Ok(value)
        } else {
            Err(BuoyancyError::NotFinite { name, value })
        }
    }

    /// Validate that `value` is finite and strictly positive (`> 0`),
    /// returning it unchanged on success. A finite but non-positive
    /// value yields [`BuoyancyError::OutOfDomain`]; a non-finite value
    /// yields [`BuoyancyError::NotFinite`] (the finiteness check runs
    /// first).
    ///
    /// ```
    /// use valenx_buoyancy::BuoyancyError;
    /// assert_eq!(BuoyancyError::positive("rho", 1000.0).unwrap(), 1000.0);
    /// assert!(BuoyancyError::positive("rho", 0.0).is_err());
    /// assert!(BuoyancyError::positive("rho", -1.0).is_err());
    /// ```
    pub fn positive(name: &'static str, value: f64) -> std::result::Result<f64, BuoyancyError> {
        let value = BuoyancyError::finite(name, value)?;
        if value > 0.0 {
            Ok(value)
        } else {
            Err(BuoyancyError::OutOfDomain {
                name,
                value,
                reason: "must be strictly positive",
            })
        }
    }

    /// Validate that `value` is finite and non-negative (`>= 0`),
    /// returning it unchanged on success. Used for quantities that may
    /// legitimately be zero (e.g. a vertical separation `BG = 0`).
    ///
    /// ```
    /// use valenx_buoyancy::BuoyancyError;
    /// assert_eq!(BuoyancyError::non_negative("bg", 0.0).unwrap(), 0.0);
    /// assert!(BuoyancyError::non_negative("bg", -0.5).is_err());
    /// ```
    pub fn non_negative(name: &'static str, value: f64) -> std::result::Result<f64, BuoyancyError> {
        let value = BuoyancyError::finite(name, value)?;
        if value >= 0.0 {
            Ok(value)
        } else {
            Err(BuoyancyError::OutOfDomain {
                name,
                value,
                reason: "must be non-negative",
            })
        }
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, BuoyancyError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finite_accepts_normal_and_rejects_nonfinite() {
        assert_eq!(BuoyancyError::finite("x", 3.5).unwrap(), 3.5);
        assert_eq!(BuoyancyError::finite("x", -2.0).unwrap(), -2.0);
        assert_eq!(BuoyancyError::finite("x", 0.0).unwrap(), 0.0);
        assert!(BuoyancyError::finite("x", f64::NAN).is_err());
        assert!(BuoyancyError::finite("x", f64::INFINITY).is_err());
        assert!(BuoyancyError::finite("x", f64::NEG_INFINITY).is_err());
    }

    #[test]
    fn positive_rejects_zero_negative_and_nonfinite() {
        assert_eq!(BuoyancyError::positive("x", 1.0).unwrap(), 1.0);
        assert!(BuoyancyError::positive("x", 0.0).is_err());
        assert!(BuoyancyError::positive("x", -1.0).is_err());
        // Non-finite is caught as NotFinite, not OutOfDomain.
        assert_eq!(
            BuoyancyError::positive("x", f64::NAN)
                .unwrap_err()
                .category(),
            ErrorCategory::NotFinite
        );
        assert_eq!(
            BuoyancyError::positive("x", -1.0).unwrap_err().category(),
            ErrorCategory::Domain
        );
    }

    #[test]
    fn non_negative_allows_zero_but_not_negative() {
        assert_eq!(BuoyancyError::non_negative("x", 0.0).unwrap(), 0.0);
        assert_eq!(BuoyancyError::non_negative("x", 2.0).unwrap(), 2.0);
        assert!(BuoyancyError::non_negative("x", -0.001).is_err());
        assert!(BuoyancyError::non_negative("x", f64::NEG_INFINITY).is_err());
    }

    #[test]
    fn code_and_category_are_stable() {
        let nf = BuoyancyError::finite("g", f64::NAN).unwrap_err();
        assert_eq!(nf.code(), "buoyancy.not_finite");
        assert_eq!(nf.category(), ErrorCategory::NotFinite);

        let dom = BuoyancyError::positive("v", -3.0).unwrap_err();
        assert_eq!(dom.code(), "buoyancy.out_of_domain");
        assert_eq!(dom.category(), ErrorCategory::Domain);
    }

    #[test]
    fn display_mentions_parameter_and_value() {
        let msg = BuoyancyError::positive("rho_fluid", 0.0)
            .unwrap_err()
            .to_string();
        assert!(msg.contains("rho_fluid"), "got: {msg}");
        assert!(msg.contains('0'), "got: {msg}");

        let msg = BuoyancyError::finite("g", f64::INFINITY)
            .unwrap_err()
            .to_string();
        assert!(msg.contains('g'), "got: {msg}");
        assert!(msg.contains("finite"), "got: {msg}");
    }

    #[test]
    fn error_is_std_error_trait_object() {
        let err: Box<dyn std::error::Error> =
            Box::new(BuoyancyError::positive("x", 0.0).unwrap_err());
        assert!(err.to_string().contains('x'));
    }
}
