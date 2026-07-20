//! Error taxonomy for `valenx-doseresponse`.
//!
//! Every fallible public function in this crate returns
//! [`Result<T, DoseResponseError>`](crate::Result). The variants are
//! deliberately small: a pharmacology caller almost always wants to know
//! one of three things.
//!
//! 1. Was a parameter non-finite (`NaN` / infinite)?
//! 2. Was a parameter outside the domain physics demands — a
//!    non-positive `EC50` (a concentration), a non-positive Hill slope,
//!    a negative concentration, a non-positive maximal effect?
//! 3. Was a requested response level outside the achievable
//!    `[0, Emax)` window, so the inverse Hill equation has no solution?
//!
//! [`DoseResponseError::code`] gives a stable kebab-cased string for log
//! / telemetry tagging; the constructors ([`non_finite`],
//! [`non_positive`], [`negative`], [`out_of_range`]) do the validation so
//! the model functions never see bad input. The style mirrors
//! `valenx-springs`' `SpringsError`.
//!
//! [`non_finite`]: DoseResponseError::non_finite
//! [`non_positive`]: DoseResponseError::non_positive
//! [`negative`]: DoseResponseError::negative
//! [`out_of_range`]: DoseResponseError::out_of_range

use thiserror::Error;

/// Errors raised when building or evaluating a dose-response model.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum DoseResponseError {
    /// A parameter was `NaN` or infinite. `what` names the parameter.
    #[error("parameter `{what}` must be finite, got {value}")]
    NonFinite {
        /// Logical parameter name (e.g. `"ec50"`, `"hill_slope"`).
        what: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A parameter the physics requires to be strictly positive was
    /// zero or negative — an `EC50` (a concentration), a Hill slope, or
    /// a maximal effect `Emax`.
    #[error("parameter `{what}` must be positive, got {value}")]
    NonPositive {
        /// Logical parameter name.
        what: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A quantity that may be zero but never negative (a drug
    /// concentration) was supplied negative.
    #[error("parameter `{what}` must be non-negative, got {value}")]
    Negative {
        /// Logical parameter name.
        what: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A requested response / fraction was outside the half-open window
    /// the model can reach, so an inverse has no finite solution. For
    /// the Hill equation the achievable absolute response is
    /// `[0, Emax)` and the achievable fraction is `[0, 1)`; the upper
    /// bound is open because `Emax` (fraction `1`) is only approached
    /// as the concentration tends to infinity.
    #[error("requested `{what}` = {value} is outside the achievable range [{lo}, {hi})")]
    OutOfRange {
        /// Logical quantity name (e.g. `"response"`, `"fraction"`).
        what: &'static str,
        /// The offending value.
        value: f64,
        /// Inclusive lower bound.
        lo: f64,
        /// Exclusive upper bound.
        hi: f64,
    },
}

/// Coarse error category for routing / display.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A value was not a finite number.
    NotFinite,
    /// A value was outside its physical domain.
    Domain,
}

impl DoseResponseError {
    /// Stable kebab-cased identifier suitable for logs / telemetry.
    /// Codes never change across minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            DoseResponseError::NonFinite { .. } => "doseresponse.non-finite",
            DoseResponseError::NonPositive { .. } => "doseresponse.non-positive",
            DoseResponseError::Negative { .. } => "doseresponse.negative",
            DoseResponseError::OutOfRange { .. } => "doseresponse.out-of-range",
        }
    }

    /// Coarse category — see [`ErrorCategory`].
    pub fn category(&self) -> ErrorCategory {
        match self {
            DoseResponseError::NonFinite { .. } => ErrorCategory::NotFinite,
            DoseResponseError::NonPositive { .. }
            | DoseResponseError::Negative { .. }
            | DoseResponseError::OutOfRange { .. } => ErrorCategory::Domain,
        }
    }

    /// Validate that `value` is finite, returning it unchanged or a
    /// [`NonFinite`](DoseResponseError::NonFinite) error naming `what`.
    pub fn finite(what: &'static str, value: f64) -> Result<f64, DoseResponseError> {
        if value.is_finite() {
            Ok(value)
        } else {
            Err(DoseResponseError::non_finite(what, value))
        }
    }

    /// Validate that `value` is finite and strictly positive — the rule
    /// for concentrations such as `EC50`, the Hill slope, and `Emax`.
    pub fn positive(what: &'static str, value: f64) -> Result<f64, DoseResponseError> {
        let value = DoseResponseError::finite(what, value)?;
        if value > 0.0 {
            Ok(value)
        } else {
            Err(DoseResponseError::non_positive(what, value))
        }
    }

    /// Validate that `value` is finite and non-negative — the rule for a
    /// drug concentration, which may be exactly zero.
    pub fn non_negative(what: &'static str, value: f64) -> Result<f64, DoseResponseError> {
        let value = DoseResponseError::finite(what, value)?;
        if value >= 0.0 {
            Ok(value)
        } else {
            Err(DoseResponseError::negative(what, value))
        }
    }

    /// Construct a [`NonFinite`](DoseResponseError::NonFinite).
    pub fn non_finite(what: &'static str, value: f64) -> Self {
        DoseResponseError::NonFinite { what, value }
    }

    /// Construct a [`NonPositive`](DoseResponseError::NonPositive).
    pub fn non_positive(what: &'static str, value: f64) -> Self {
        DoseResponseError::NonPositive { what, value }
    }

    /// Construct a [`Negative`](DoseResponseError::Negative).
    pub fn negative(what: &'static str, value: f64) -> Self {
        DoseResponseError::Negative { what, value }
    }

    /// Construct an [`OutOfRange`](DoseResponseError::OutOfRange).
    pub fn out_of_range(what: &'static str, value: f64, lo: f64, hi: f64) -> Self {
        DoseResponseError::OutOfRange {
            what,
            value,
            lo,
            hi,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let e = DoseResponseError::non_finite("ec50", f64::NAN);
        assert_eq!(e.code(), "doseresponse.non-finite");
        assert_eq!(e.category(), ErrorCategory::NotFinite);

        let e = DoseResponseError::non_positive("ec50", 0.0);
        assert_eq!(e.code(), "doseresponse.non-positive");
        assert_eq!(e.category(), ErrorCategory::Domain);

        let e = DoseResponseError::negative("concentration", -1.0);
        assert_eq!(e.code(), "doseresponse.negative");
        assert_eq!(e.category(), ErrorCategory::Domain);

        let e = DoseResponseError::out_of_range("fraction", 1.5, 0.0, 1.0);
        assert_eq!(e.code(), "doseresponse.out-of-range");
        assert_eq!(e.category(), ErrorCategory::Domain);
    }

    #[test]
    fn finite_accepts_and_rejects() {
        assert_eq!(DoseResponseError::finite("x", 2.5).unwrap(), 2.5);
        assert!(DoseResponseError::finite("x", f64::INFINITY).is_err());
        assert!(DoseResponseError::finite("x", f64::NAN).is_err());
    }

    #[test]
    fn positive_rejects_zero_and_negative() {
        assert_eq!(DoseResponseError::positive("x", 1e-9).unwrap(), 1e-9);
        assert!(DoseResponseError::positive("x", 0.0).is_err());
        assert!(DoseResponseError::positive("x", -3.0).is_err());
        // Non-finite is caught first and reported as NonFinite.
        assert_eq!(
            DoseResponseError::positive("x", f64::NAN)
                .unwrap_err()
                .code(),
            "doseresponse.non-finite"
        );
    }

    #[test]
    fn non_negative_allows_zero_rejects_negative() {
        assert_eq!(DoseResponseError::non_negative("c", 0.0).unwrap(), 0.0);
        assert_eq!(DoseResponseError::non_negative("c", 5.0).unwrap(), 5.0);
        assert!(DoseResponseError::non_negative("c", -0.001).is_err());
    }

    #[test]
    fn display_is_informative() {
        let msg = DoseResponseError::non_positive("ec50", -2.0).to_string();
        assert!(msg.contains("ec50"), "got: {msg}");
        assert!(msg.contains("positive"), "got: {msg}");

        let msg = DoseResponseError::out_of_range("response", 200.0, 0.0, 100.0).to_string();
        assert!(msg.contains("200"), "got: {msg}");
        assert!(msg.contains("100"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let e: Box<dyn std::error::Error> = Box::new(DoseResponseError::negative("c", -1.0));
        assert!(e.to_string().contains('c'));
    }
}
