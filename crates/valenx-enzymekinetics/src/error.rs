//! Error taxonomy for `valenx-enzymekinetics`.
//!
//! Every fallible constructor and rate function in this crate returns
//! [`Result<_, KineticsError>`](crate::Result). The variants separate the
//! two failure modes an enzyme-kinetics caller cares about:
//!
//! 1. A parameter is outside its physical domain — a negative `Vmax`, a
//!    non-positive `Km`, a Hill coefficient `<= 0`, a negative substrate
//!    or inhibitor concentration ([`KineticsError::NonFinite`] guards the
//!    NaN / infinity case; [`KineticsError::OutOfDomain`] guards the
//!    finite-but-illegal case).
//! 2. The value supplied was not even a finite number
//!    ([`KineticsError::NonFinite`]).
//!
//! The constructors on the parameter structs ([`crate::MichaelisMenten`],
//! [`crate::Hill`], the inhibition types) validate eagerly, so once a
//! parameter object exists every rate evaluation on it is infallible.

use thiserror::Error;

/// Errors produced by `valenx-enzymekinetics`.
///
/// Derives [`thiserror::Error`]; each variant carries a human-readable
/// `Display` message via its `#[error(...)]` attribute.
///
/// Marked `#[non_exhaustive]`: more validation cases may be added as the
/// model surface grows, so downstream matches must include a wildcard
/// arm.
#[derive(Debug, Clone, PartialEq, Error)]
#[non_exhaustive]
pub enum KineticsError {
    /// A parameter is a finite number but lies outside the domain the
    /// model requires (e.g. a non-positive Michaelis constant, a negative
    /// maximal velocity, a Hill coefficient that is zero or negative).
    #[error("parameter `{what}` is out of domain: {value} ({reason})")]
    OutOfDomain {
        /// Logical parameter name (e.g. `"km"`, `"vmax"`, `"n"`).
        what: &'static str,
        /// The offending value, surfaced verbatim for diagnosis.
        value: f64,
        /// Human-readable explanation of the constraint that was broken.
        reason: &'static str,
    },

    /// A parameter was `NaN` or `±∞`. Kept separate from
    /// [`KineticsError::OutOfDomain`] because a non-finite value usually
    /// signals an upstream computation error rather than a merely
    /// out-of-range user input.
    #[error("parameter `{what}` is not finite (was {value})")]
    NonFinite {
        /// Logical parameter name.
        what: &'static str,
        /// The non-finite value (`NaN`, `inf`, or `-inf`).
        value: f64,
    },
}

impl KineticsError {
    /// Convenience constructor for [`KineticsError::OutOfDomain`].
    pub fn out_of_domain(what: &'static str, value: f64, reason: &'static str) -> Self {
        KineticsError::OutOfDomain {
            what,
            value,
            reason,
        }
    }

    /// Convenience constructor for [`KineticsError::NonFinite`].
    pub fn non_finite(what: &'static str, value: f64) -> Self {
        KineticsError::NonFinite { what, value }
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, KineticsError>;

// --- Shared validation helpers ----------------------------------------
//
// These keep the parameter constructors terse and guarantee the same
// rules everywhere (finiteness first, then domain). They are crate-
// internal; callers see only the resulting [`KineticsError`].

/// Require a finite, strictly positive value (`x > 0`).
///
/// Used for the Michaelis constant `Km`, the Hill half-saturation
/// constant `K`, the inhibition constant `Ki`, and the Hill coefficient
/// `n` — none of which has a meaningful zero or negative value.
pub(crate) fn require_positive(what: &'static str, x: f64) -> Result<f64> {
    if !x.is_finite() {
        return Err(KineticsError::non_finite(what, x));
    }
    if x <= 0.0 {
        return Err(KineticsError::out_of_domain(
            what,
            x,
            "must be strictly positive",
        ));
    }
    Ok(x)
}

/// Require a finite, non-negative value (`x >= 0`).
///
/// Used for `Vmax`, substrate concentration `S`, and inhibitor
/// concentration `I`, all of which may legitimately be zero.
pub(crate) fn require_non_negative(what: &'static str, x: f64) -> Result<f64> {
    if !x.is_finite() {
        return Err(KineticsError::non_finite(what, x));
    }
    if x < 0.0 {
        return Err(KineticsError::out_of_domain(
            what,
            x,
            "must be non-negative",
        ));
    }
    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn out_of_domain_display_names_param_and_reason() {
        let e = KineticsError::out_of_domain("km", -1.0, "must be strictly positive");
        let msg = e.to_string();
        assert!(msg.contains("km"), "got: {msg}");
        assert!(msg.contains("strictly positive"), "got: {msg}");
        assert!(msg.contains("-1"), "got: {msg}");
    }

    #[test]
    fn non_finite_display_names_param() {
        let e = KineticsError::non_finite("vmax", f64::NAN);
        let msg = e.to_string();
        assert!(msg.contains("vmax"), "got: {msg}");
        assert!(msg.contains("not finite"), "got: {msg}");
    }

    #[test]
    fn require_positive_accepts_positive() {
        let v = require_positive("km", 2.5).expect("2.5 is positive");
        assert!((v - 2.5).abs() < 1e-12, "got: {v}");
    }

    #[test]
    fn require_positive_rejects_zero_and_negative() {
        assert!(matches!(
            require_positive("km", 0.0),
            Err(KineticsError::OutOfDomain { .. })
        ));
        assert!(matches!(
            require_positive("km", -3.0),
            Err(KineticsError::OutOfDomain { .. })
        ));
    }

    #[test]
    fn require_positive_rejects_non_finite() {
        assert!(matches!(
            require_positive("km", f64::INFINITY),
            Err(KineticsError::NonFinite { .. })
        ));
        assert!(matches!(
            require_positive("km", f64::NAN),
            Err(KineticsError::NonFinite { .. })
        ));
    }

    #[test]
    fn require_non_negative_accepts_zero() {
        let v = require_non_negative("s", 0.0).expect("0 is non-negative");
        assert!(v.abs() < 1e-12, "got: {v}");
    }

    #[test]
    fn require_non_negative_rejects_negative_and_non_finite() {
        assert!(matches!(
            require_non_negative("s", -0.001),
            Err(KineticsError::OutOfDomain { .. })
        ));
        assert!(matches!(
            require_non_negative("s", f64::NEG_INFINITY),
            Err(KineticsError::NonFinite { .. })
        ));
    }

    #[test]
    fn error_is_a_std_error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(KineticsError::non_finite("n", f64::NAN));
        assert!(err.to_string().contains('n'));
    }
}
