//! Error taxonomy for `valenx-columnsteel`.
//!
//! Every fallible public function returns
//! [`Result<_, ColumnError>`]. The variants are deliberately coarse — a
//! column-design caller usually only cares about two things:
//!
//! 1. Did the caller pass a non-physical argument — a non-positive
//!    modulus, a zero radius of gyration, a negative effective length
//!    ([`ColumnError::Invalid`])?
//! 2. Is the requested slenderness outside the range the chosen method
//!    can represent ([`ColumnError::OutOfRange`])?
//!
//! Use [`ColumnError::code`] for stable log / telemetry tagging and
//! [`ColumnError::category`] to bucket failures without matching every
//! variant. The pattern mirrors the other Valenx engineering crates'
//! error types (`SpringsError`, `GearsError`).

use thiserror::Error;

/// Errors raised by steel-column slenderness and buckling calculations.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum ColumnError {
    /// Caller passed an argument the model cannot accept: a non-positive
    /// Young's modulus or yield stress, a zero / negative radius of
    /// gyration, a negative unbraced length, a non-positive effective-
    /// length factor, or a non-positive factor of safety. A property of
    /// the *call*, not of any parsed data.
    #[error("invalid `{what}`: {reason}")]
    Invalid {
        /// Logical parameter name (e.g. `"youngs_modulus"`, `"radius_of_gyration"`).
        what: &'static str,
        /// Human-readable reason the value was rejected.
        reason: String,
    },

    /// A computed or supplied slenderness ratio fell outside the domain
    /// the requested formula represents — for example asking for the
    /// Johnson (short-column) stress at a slenderness beyond the
    /// transition `Cc`, or the Euler (long-column) stress below it.
    #[error("slenderness {slenderness} out of range for {regime}: {reason}")]
    OutOfRange {
        /// Slenderness ratio `KL/r` that was rejected.
        slenderness: f64,
        /// Regime whose domain was violated (`"euler"` or `"johnson"`).
        regime: &'static str,
        /// Human-readable reason.
        reason: String,
    },
}

/// Coarse error category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this rather
/// than on the individual error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A call argument was non-physical.
    Input,
    /// A slenderness ratio was outside the requested formula's domain.
    Domain,
}

impl ColumnError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"columnsteel.<sub_id>"`. Codes never change
    /// across minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            ColumnError::Invalid { .. } => "columnsteel.invalid",
            ColumnError::OutOfRange { .. } => "columnsteel.out_of_range",
        }
    }

    /// Coarse category — see [`ErrorCategory`].
    pub fn category(&self) -> ErrorCategory {
        match self {
            ColumnError::Invalid { .. } => ErrorCategory::Input,
            ColumnError::OutOfRange { .. } => ErrorCategory::Domain,
        }
    }

    /// Construct a [`ColumnError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        ColumnError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Construct a [`ColumnError::OutOfRange`].
    pub fn out_of_range(slenderness: f64, regime: &'static str, reason: impl Into<String>) -> Self {
        ColumnError::OutOfRange {
            slenderness,
            regime,
            reason: reason.into(),
        }
    }

    /// Validate that a quantity is strictly positive and finite,
    /// returning it on success or a [`ColumnError::Invalid`] naming the
    /// parameter on failure. The shared guard behind every public
    /// constructor in [`crate::column`].
    pub fn require_positive(value: f64, what: &'static str) -> Result<f64> {
        if !value.is_finite() {
            return Err(ColumnError::invalid(what, "must be a finite number"));
        }
        if value <= 0.0 {
            return Err(ColumnError::invalid(what, "must be strictly positive"));
        }
        Ok(value)
    }

    /// Validate that a quantity is finite and non-negative (zero is
    /// allowed — e.g. an unbraced length of zero is a degenerate but
    /// physical input that yields zero slenderness).
    pub fn require_non_negative(value: f64, what: &'static str) -> Result<f64> {
        if !value.is_finite() {
            return Err(ColumnError::invalid(what, "must be a finite number"));
        }
        if value < 0.0 {
            return Err(ColumnError::invalid(what, "must be non-negative"));
        }
        Ok(value)
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, ColumnError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = ColumnError::invalid("youngs_modulus", "must be strictly positive");
        assert_eq!(err.code(), "columnsteel.invalid");
        assert_eq!(err.category(), ErrorCategory::Input);

        let err = ColumnError::out_of_range(150.0, "johnson", "beyond Cc");
        assert_eq!(err.code(), "columnsteel.out_of_range");
        assert_eq!(err.category(), ErrorCategory::Domain);
    }

    #[test]
    fn display_is_informative() {
        let msg =
            ColumnError::invalid("radius_of_gyration", "must be strictly positive").to_string();
        assert!(msg.contains("radius_of_gyration"), "got: {msg}");
        assert!(msg.contains("strictly positive"), "got: {msg}");

        let msg = ColumnError::out_of_range(123.0, "euler", "below Cc").to_string();
        assert!(msg.contains("123"), "got: {msg}");
        assert!(msg.contains("euler"), "got: {msg}");
    }

    #[test]
    fn require_positive_accepts_and_rejects() {
        assert_eq!(ColumnError::require_positive(2.0, "x").unwrap(), 2.0);
        assert!(ColumnError::require_positive(0.0, "x").is_err());
        assert!(ColumnError::require_positive(-1.0, "x").is_err());
        assert!(ColumnError::require_positive(f64::NAN, "x").is_err());
        assert!(ColumnError::require_positive(f64::INFINITY, "x").is_err());
    }

    #[test]
    fn require_non_negative_allows_zero_rejects_negative() {
        assert_eq!(ColumnError::require_non_negative(0.0, "x").unwrap(), 0.0);
        assert_eq!(ColumnError::require_non_negative(5.0, "x").unwrap(), 5.0);
        assert!(ColumnError::require_non_negative(-0.1, "x").is_err());
        assert!(ColumnError::require_non_negative(f64::NAN, "x").is_err());
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(ColumnError::invalid("k", "bad"));
        assert!(err.to_string().contains('k'));
    }
}
