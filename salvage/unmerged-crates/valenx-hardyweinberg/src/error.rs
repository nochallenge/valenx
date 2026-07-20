//! Error taxonomy for `valenx-hardyweinberg`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, HweError>`](crate::error::Result). The variants are
//! deliberately small — a Hardy-Weinberg caller only ever trips over
//! three things:
//!
//! 1. A numeric argument is not finite ([`HweError::NotFinite`]) — a
//!    `NaN` or `inf` slipped in.
//! 2. A numeric argument is outside its physical domain
//!    ([`HweError::OutOfDomain`]) — a frequency outside `[0, 1]`, a
//!    negative count, a pair of allele frequencies that do not sum to
//!    one, or a sample whose total is zero.
//! 3. The chi-square test has no degrees of freedom left
//!    ([`HweError::NoDegreesOfFreedom`]) — the contingency reduces to a
//!    single fixed cell so the goodness-of-fit statistic is undefined.
//!
//! Use [`HweError::code`] for stable log / telemetry tagging; the codes
//! never change across minor versions. The pattern mirrors
//! `valenx-springs`'s `SpringsError` and `valenx-popgen`'s
//! `PopgenError`.

use thiserror::Error;

/// Errors produced by `valenx-hardyweinberg`.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum HweError {
    /// A numeric argument was `NaN` or infinite. `what` is the logical
    /// parameter name (e.g. `"p"`, `"allele_a_count"`).
    #[error("argument `{what}` is not finite (NaN or infinite)")]
    NotFinite {
        /// Logical parameter name.
        what: &'static str,
    },

    /// A numeric argument fell outside its admissible domain — a
    /// frequency outside `[0, 1]`, a negative genotype count, allele
    /// frequencies that fail to sum to one, or an all-zero sample.
    #[error("argument `{what}` is out of domain: {reason}")]
    OutOfDomain {
        /// Logical parameter name.
        what: &'static str,
        /// Human-readable reason the value is rejected.
        reason: String,
    },

    /// The chi-square goodness-of-fit test has zero (or fewer) degrees
    /// of freedom, so the statistic is undefined. For a single
    /// di-allelic locus this happens when one allele is fixed.
    #[error("chi-square test has no positive degrees of freedom: {reason}")]
    NoDegreesOfFreedom {
        /// Human-readable reason.
        reason: String,
    },
}

/// Coarse category for routing / display. Stable across crate versions:
/// switch a single `match` on this rather than on each error variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A numeric input was malformed (non-finite or out of domain).
    Input,
    /// The requested statistic is undefined for the given data.
    Undefined,
}

impl HweError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"hwe.<sub_id>"`. Codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            HweError::NotFinite { .. } => "hwe.not_finite",
            HweError::OutOfDomain { .. } => "hwe.out_of_domain",
            HweError::NoDegreesOfFreedom { .. } => "hwe.no_degrees_of_freedom",
        }
    }

    /// Coarse [`ErrorCategory`] for the variant.
    pub fn category(&self) -> ErrorCategory {
        match self {
            HweError::NotFinite { .. } | HweError::OutOfDomain { .. } => ErrorCategory::Input,
            HweError::NoDegreesOfFreedom { .. } => ErrorCategory::Undefined,
        }
    }

    /// Convenience constructor for [`HweError::NotFinite`].
    pub fn not_finite(what: &'static str) -> Self {
        HweError::NotFinite { what }
    }

    /// Convenience constructor for [`HweError::OutOfDomain`].
    pub fn out_of_domain(what: &'static str, reason: impl Into<String>) -> Self {
        HweError::OutOfDomain {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`HweError::NoDegreesOfFreedom`].
    pub fn no_degrees_of_freedom(reason: impl Into<String>) -> Self {
        HweError::NoDegreesOfFreedom {
            reason: reason.into(),
        }
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, HweError>;

/// Reject a value that is `NaN` or infinite.
///
/// Internal helper used by the public model functions to validate every
/// floating-point argument before it reaches a formula.
pub(crate) fn finite(what: &'static str, value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(HweError::not_finite(what))
    }
}

/// Reject a frequency that is non-finite or outside the closed unit
/// interval `[0, 1]`.
pub(crate) fn unit_interval(what: &'static str, value: f64) -> Result<f64> {
    let v = finite(what, value)?;
    if (0.0..=1.0).contains(&v) {
        Ok(v)
    } else {
        Err(HweError::out_of_domain(
            what,
            format!("frequency {v} is outside [0, 1]"),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = HweError::not_finite("p");
        assert_eq!(err.code(), "hwe.not_finite");
        assert_eq!(err.category(), ErrorCategory::Input);

        let err = HweError::out_of_domain("q", "must be in [0, 1]");
        assert_eq!(err.code(), "hwe.out_of_domain");
        assert_eq!(err.category(), ErrorCategory::Input);

        let err = HweError::no_degrees_of_freedom("locus is fixed");
        assert_eq!(err.code(), "hwe.no_degrees_of_freedom");
        assert_eq!(err.category(), ErrorCategory::Undefined);
    }

    #[test]
    fn display_is_informative() {
        let msg = HweError::out_of_domain("p", "frequency 1.5 is outside [0, 1]").to_string();
        assert!(msg.contains('p'), "got: {msg}");
        assert!(msg.contains("1.5"), "got: {msg}");
    }

    #[test]
    fn finite_helper_rejects_nan_and_inf() {
        assert_eq!(finite("x", 0.5).unwrap(), 0.5);
        assert_eq!(finite("x", f64::NAN), Err(HweError::not_finite("x")));
        assert_eq!(finite("x", f64::INFINITY), Err(HweError::not_finite("x")));
    }

    #[test]
    fn unit_interval_helper_bounds() {
        assert_eq!(unit_interval("p", 0.0).unwrap(), 0.0);
        assert_eq!(unit_interval("p", 1.0).unwrap(), 1.0);
        assert!(unit_interval("p", -0.001).is_err());
        assert!(unit_interval("p", 1.001).is_err());
        assert!(unit_interval("p", f64::NAN).is_err());
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(HweError::not_finite("p"));
        assert!(err.to_string().contains('p'));
    }
}
