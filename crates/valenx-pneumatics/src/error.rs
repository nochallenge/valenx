//! Error taxonomy for `valenx-pneumatics`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, PneumaticsError>`]. The variants are intentionally coarse;
//! a sizing caller usually only needs to distinguish three things:
//!
//! 1. Was a quantity that must be strictly positive instead zero,
//!    negative or non-finite — a bore, a stroke, a pressure, a count
//!    ([`PneumaticsError::NonPositive`])?
//! 2. Was a quantity that must be non-negative instead negative or
//!    non-finite — a gauge pressure that may legitimately be zero
//!    ([`PneumaticsError::Negative`])?
//! 3. Did two related dimensions form an impossible geometry — a piston
//!    rod at least as wide as the bore it runs inside
//!    ([`PneumaticsError::Geometry`])?
//!
//! The validated constructors ([`PneumaticsError::positive`],
//! [`PneumaticsError::non_negative`]) centralise the finite-and-in-range
//! checks so every public model reuses identical validation, and
//! [`PneumaticsError::code`] gives a stable snake-cased tag for logs and
//! telemetry. The shape mirrors the rest of the workspace's `*Error`
//! enums (e.g. `valenx-astro`'s `AstroError`).

use thiserror::Error;

/// Shorthand for `Result<T, PneumaticsError>`.
pub type Result<T> = core::result::Result<T, PneumaticsError>;

/// Anything that can go wrong validating a pneumatic sizing input.
///
/// This enum is `#[non_exhaustive]`: new variants may be added in future
/// releases without it being a breaking change, so downstream `match`
/// arms must include a wildcard.
#[derive(Debug, Error, Clone, PartialEq)]
#[non_exhaustive]
pub enum PneumaticsError {
    /// A quantity that must be strictly positive (a bore diameter, a
    /// stroke length, an absolute pressure, a cycle count, the
    /// ratio-of-specific-heats `k`) was zero, negative, or non-finite.
    /// Carries the parameter name and the offending value.
    #[error("`{what}` must be a finite value greater than zero (got {value})")]
    NonPositive {
        /// Logical parameter name (e.g. `"bore"`, `"stroke"`,
        /// `"supply_pressure"`).
        what: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A quantity that must be non-negative (a *gauge* pressure, which may
    /// legitimately be exactly zero at atmospheric) was negative or
    /// non-finite. Carries the parameter name and the offending value.
    #[error("`{what}` must be a finite value of zero or greater (got {value})")]
    Negative {
        /// Logical parameter name (e.g. `"gauge_pressure"`).
        what: &'static str,
        /// The offending value.
        value: f64,
    },

    /// Two related dimensions describe an impossible geometry — most
    /// commonly a piston rod whose diameter equals or exceeds the cylinder
    /// bore it runs inside, which would leave a zero or negative annular
    /// (retract) area.
    #[error("invalid geometry: {0}")]
    Geometry(&'static str),
}

impl PneumaticsError {
    /// Stable snake-cased error code suitable for log / telemetry tagging.
    ///
    /// Format: `"pneumatics.<sub_id>"`. Codes never change across minor
    /// versions, so they are safe to match in dashboards and alerts.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_pneumatics::error::PneumaticsError;
    /// let err = PneumaticsError::positive("bore", 0.0).unwrap_err();
    /// assert_eq!(err.code(), "pneumatics.non_positive");
    /// ```
    pub fn code(&self) -> &'static str {
        match self {
            PneumaticsError::NonPositive { .. } => "pneumatics.non_positive",
            PneumaticsError::Negative { .. } => "pneumatics.negative",
            PneumaticsError::Geometry(_) => "pneumatics.geometry",
        }
    }

    /// Validate that `value` is finite and strictly greater than zero,
    /// returning it on success or a [`PneumaticsError::NonPositive`]
    /// tagged with `what` on failure.
    ///
    /// This is the single gate every "must be positive" parameter in the
    /// crate passes through, so the finite / sign rules stay identical
    /// across modules.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_pneumatics::error::PneumaticsError;
    /// assert_eq!(PneumaticsError::positive("bore", 0.032).unwrap(), 0.032);
    /// assert!(PneumaticsError::positive("bore", 0.0).is_err());
    /// assert!(PneumaticsError::positive("bore", -1.0).is_err());
    /// assert!(PneumaticsError::positive("bore", f64::NAN).is_err());
    /// ```
    pub fn positive(what: &'static str, value: f64) -> Result<f64> {
        if value.is_finite() && value > 0.0 {
            Ok(value)
        } else {
            Err(PneumaticsError::NonPositive { what, value })
        }
    }

    /// Validate that `value` is finite and greater than or equal to zero,
    /// returning it on success or a [`PneumaticsError::Negative`] tagged
    /// with `what` on failure.
    ///
    /// Used for *gauge* pressures, which equal zero at atmospheric and so
    /// must be allowed through where a strictly-positive check would
    /// wrongly reject them.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_pneumatics::error::PneumaticsError;
    /// assert_eq!(PneumaticsError::non_negative("p", 0.0).unwrap(), 0.0);
    /// assert!(PneumaticsError::non_negative("p", -0.1).is_err());
    /// assert!(PneumaticsError::non_negative("p", f64::INFINITY).is_err());
    /// ```
    pub fn non_negative(what: &'static str, value: f64) -> Result<f64> {
        if value.is_finite() && value >= 0.0 {
            Ok(value)
        } else {
            Err(PneumaticsError::Negative { what, value })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn positive_accepts_positive_finite() {
        assert_eq!(PneumaticsError::positive("x", 1.5).unwrap(), 1.5);
    }

    #[test]
    fn positive_rejects_zero_negative_and_nonfinite() {
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let err = PneumaticsError::positive("x", bad).unwrap_err();
            assert_eq!(err.code(), "pneumatics.non_positive");
            match err {
                PneumaticsError::NonPositive { what, value } => {
                    assert_eq!(what, "x");
                    // NaN never equals itself, so compare bit-patterns for
                    // the NaN case and an ordinary equality otherwise.
                    assert!(value.to_bits() == bad.to_bits() || value == bad);
                }
                other => panic!("wrong variant: {other:?}"),
            }
        }
    }

    #[test]
    fn non_negative_accepts_zero() {
        assert_eq!(PneumaticsError::non_negative("p", 0.0).unwrap(), 0.0);
    }

    #[test]
    fn non_negative_rejects_negative_and_nonfinite() {
        for bad in [-0.001, f64::NAN, f64::INFINITY] {
            let err = PneumaticsError::non_negative("p", bad).unwrap_err();
            assert_eq!(err.code(), "pneumatics.negative");
        }
    }

    #[test]
    fn display_mentions_parameter_and_value() {
        let msg = PneumaticsError::positive("bore", -2.0)
            .unwrap_err()
            .to_string();
        assert!(msg.contains("bore"), "got: {msg}");

        let msg = PneumaticsError::Geometry("rod >= bore").to_string();
        assert!(msg.contains("rod >= bore"), "got: {msg}");
    }

    #[test]
    fn error_is_a_std_error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(PneumaticsError::Geometry("bad"));
        assert!(err.to_string().contains("bad"));
    }

    #[test]
    fn codes_are_distinct_per_variant() {
        assert_eq!(
            PneumaticsError::NonPositive {
                what: "a",
                value: 0.0
            }
            .code(),
            "pneumatics.non_positive"
        );
        assert_eq!(
            PneumaticsError::Negative {
                what: "a",
                value: -1.0
            }
            .code(),
            "pneumatics.negative"
        );
        assert_eq!(PneumaticsError::Geometry("g").code(), "pneumatics.geometry");
    }
}
