//! Error taxonomy for the strain-gauge / Wheatstone-bridge calculator.
//!
//! Every fallible entry point validates its scalar inputs and returns a
//! [`StrainGaugeError`] rather than feeding a non-finite or non-physical
//! value into a division or a `√`/`ln`. The validated constructors
//! ([`StrainGaugeError::positive`], [`StrainGaugeError::finite`]) keep the
//! check sites terse while producing a consistent, machine-stable error.

use thiserror::Error;

/// Shorthand for `Result<T, StrainGaugeError>`.
pub type Result<T> = core::result::Result<T, StrainGaugeError>;

/// Anything that can go wrong validating a gauge, bridge or material
/// input.
///
/// This enum is `#[non_exhaustive]`: new variants may be added in a
/// future release without it being a breaking change, so downstream
/// `match` arms must include a wildcard.
#[derive(Debug, Error, Clone, PartialEq)]
#[non_exhaustive]
pub enum StrainGaugeError {
    /// A parameter that must be strictly positive (gauge factor,
    /// Young's modulus, excitation voltage, nominal resistance) was
    /// zero or negative.
    #[error("parameter `{name}` must be > 0, got {value}")]
    NonPositive {
        /// Which parameter was bad (e.g. `"gauge_factor"`).
        name: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A parameter held a non-finite value (`NaN` or `±∞`), which would
    /// otherwise silently propagate through every downstream formula.
    #[error("parameter `{name}` must be finite, got {value}")]
    NonFinite {
        /// Which parameter was bad.
        name: &'static str,
        /// The offending value.
        value: f64,
    },
}

/// Coarse error category, for callers that route by class rather than by
/// exact variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The input was outside its physical domain (sign / range).
    Domain,
    /// The input was not a finite number.
    NotFinite,
}

impl StrainGaugeError {
    /// Validate that `value` is finite and strictly positive, returning
    /// it on success.
    ///
    /// `NaN`/`±∞` map to [`StrainGaugeError::NonFinite`]; a finite
    /// non-positive value maps to [`StrainGaugeError::NonPositive`].
    ///
    /// ```
    /// use valenx_straingauge::StrainGaugeError;
    ///
    /// assert_eq!(StrainGaugeError::positive("gf", 2.0).unwrap(), 2.0);
    /// assert!(StrainGaugeError::positive("gf", 0.0).is_err());
    /// assert!(StrainGaugeError::positive("gf", f64::NAN).is_err());
    /// ```
    pub fn positive(name: &'static str, value: f64) -> Result<f64> {
        if !value.is_finite() {
            return Err(StrainGaugeError::NonFinite { name, value });
        }
        if value <= 0.0 {
            return Err(StrainGaugeError::NonPositive { name, value });
        }
        Ok(value)
    }

    /// Validate that `value` is finite (any sign, including zero),
    /// returning it on success.
    ///
    /// Used for quantities such as mechanical strain, which is signed
    /// (tension positive, compression negative) and legitimately zero at
    /// the balance point.
    ///
    /// ```
    /// use valenx_straingauge::StrainGaugeError;
    ///
    /// assert_eq!(StrainGaugeError::finite("strain", -1e-3).unwrap(), -1e-3);
    /// assert!(StrainGaugeError::finite("strain", f64::INFINITY).is_err());
    /// ```
    pub fn finite(name: &'static str, value: f64) -> Result<f64> {
        if !value.is_finite() {
            return Err(StrainGaugeError::NonFinite { name, value });
        }
        Ok(value)
    }

    /// Stable kebab-cased identifier, suitable for logs / telemetry.
    pub fn code(&self) -> &'static str {
        match self {
            StrainGaugeError::NonPositive { .. } => "straingauge.non-positive",
            StrainGaugeError::NonFinite { .. } => "straingauge.non-finite",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            StrainGaugeError::NonPositive { .. } => ErrorCategory::Domain,
            StrainGaugeError::NonFinite { .. } => ErrorCategory::NotFinite,
        }
    }
}
