//! Error taxonomy for `valenx-samplingtheory`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, SamplingError>`]. The variants are intentionally coarse:
//! a DSP sampling caller usually only cares about two things.
//!
//! 1. Did the caller pass an argument outside its physical domain — a
//!    non-positive sampling rate, a negative frequency, a bit depth of
//!    zero, a measurement range that is not positive
//!    ([`SamplingError::Invalid`])?
//! 2. Do two arguments together violate a modelling assumption — a
//!    measured value that falls outside the stated full-scale range
//!    ([`SamplingError::OutOfRange`])?
//!
//! Use [`SamplingError::code`] for stable log / telemetry tagging and
//! [`SamplingError::category`] / [`SamplingError::category_enum`] to
//! bucket failures without matching every variant. The shape mirrors
//! the error taxonomies of the other pure Valenx crates (for example
//! `valenx-springs`'s `SpringsError`).

use thiserror::Error;

/// Errors produced by `valenx-samplingtheory`.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum SamplingError {
    /// Caller passed an argument outside its physical domain: a
    /// non-positive sampling rate, a negative frequency, a bit depth of
    /// zero, a non-positive full-scale range, etc. A property of a
    /// single argument's value.
    #[error("invalid `{what}`: {reason}")]
    Invalid {
        /// Logical parameter name (e.g. `"sample_rate_hz"`, `"bits"`).
        what: &'static str,
        /// Human-readable reason the value is rejected.
        reason: String,
    },

    /// A value falls outside an interval implied by another argument:
    /// most commonly a sample whose amplitude exceeds the converter's
    /// full-scale range. A property of two arguments considered
    /// together.
    #[error("{what} value {value} is outside the range [{low}, {high}]")]
    OutOfRange {
        /// Logical name of the quantity that is out of range.
        what: &'static str,
        /// The offending value.
        value: f64,
        /// Inclusive lower bound of the permitted interval.
        low: f64,
        /// Inclusive upper bound of the permitted interval.
        high: f64,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions: switch a single `match` on this rather
/// than on every error variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A single argument is outside its physical domain.
    Input,
    /// Two arguments are mutually inconsistent (a value outside a
    /// range implied by another argument).
    Domain,
}

impl SamplingError {
    /// Construct a [`SamplingError::Invalid`].
    ///
    /// Convenience wrapper so call sites read
    /// `SamplingError::invalid("bits", "must be >= 1")` rather than
    /// spelling out the struct literal.
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        SamplingError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Construct a [`SamplingError::OutOfRange`].
    pub fn out_of_range(what: &'static str, value: f64, low: f64, high: f64) -> Self {
        SamplingError::OutOfRange {
            what,
            value,
            low,
            high,
        }
    }

    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"samplingtheory.<sub_id>"`. Codes never change
    /// across minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            SamplingError::Invalid { .. } => "samplingtheory.invalid",
            SamplingError::OutOfRange { .. } => "samplingtheory.out_of_range",
        }
    }

    /// Coarse category string. See [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            SamplingError::Invalid { .. } => "input",
            SamplingError::OutOfRange { .. } => "domain",
        }
    }

    /// Typed category enum, for callers that prefer to `match` on a
    /// stable enum rather than compare the [`category`](Self::category)
    /// string.
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            SamplingError::Invalid { .. } => ErrorCategory::Input,
            SamplingError::OutOfRange { .. } => ErrorCategory::Domain,
        }
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, SamplingError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = SamplingError::invalid("bits", "must be >= 1");
        assert_eq!(err.code(), "samplingtheory.invalid");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = SamplingError::out_of_range("sample", 2.0, -1.0, 1.0);
        assert_eq!(err.code(), "samplingtheory.out_of_range");
        assert_eq!(err.category(), "domain");
        assert_eq!(err.category_enum(), ErrorCategory::Domain);
    }

    #[test]
    fn display_is_informative() {
        let msg = SamplingError::invalid("sample_rate_hz", "must be positive").to_string();
        assert!(msg.contains("sample_rate_hz"), "got: {msg}");
        assert!(msg.contains("must be positive"), "got: {msg}");

        let msg = SamplingError::out_of_range("sample", 5.0, -1.0, 1.0).to_string();
        assert!(msg.contains('5'), "got: {msg}");
        assert!(msg.contains("sample"), "got: {msg}");
    }

    #[test]
    fn error_is_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(SamplingError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }

    #[test]
    fn error_is_clone_and_eq() {
        let a = SamplingError::invalid("bits", "zero");
        let b = a.clone();
        assert_eq!(a, b);
    }
}
