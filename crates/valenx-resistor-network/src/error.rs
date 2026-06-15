//! Error taxonomy for the resistor-network helpers.
//!
//! Every fallible entry point in this crate returns
//! [`ResistorError`]. Construct values through the validated
//! constructors here ([`ResistorError::non_positive`],
//! [`ResistorError::non_finite`], [`ResistorError::empty_network`])
//! rather than building the variants by hand, so the wording and
//! the stable [`ResistorError::code`] string stay consistent.

use thiserror::Error;

/// Errors raised while evaluating a resistor network.
///
/// The models in this crate are purely real-valued (DC, ideal,
/// lumped). Anything that would make a closed-form expression
/// ill-defined — a non-positive resistance, a non-finite input,
/// or an empty combination — surfaces as one of these variants.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum ResistorError {
    /// A resistance (or other strictly-positive quantity) was zero
    /// or negative.
    ///
    /// Series and parallel reductions assume every element is a
    /// real, strictly-positive ohmic resistance; a parallel branch
    /// of `0` would imply an ideal short with infinite conductance,
    /// which the closed-form `1/R = sum(1/Ri)` model does not cover.
    #[error("non-positive value for `{name}`: {value} (must be > 0)")]
    NonPositive {
        /// Name of the offending quantity, e.g. `"resistance"`.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// An input was `NaN` or infinite.
    #[error("non-finite value for `{name}`: {value}")]
    NonFinite {
        /// Name of the offending quantity.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A series or parallel combination was requested with no
    /// resistors.
    ///
    /// The empty series sum (`0`) and the empty parallel sum
    /// (an open circuit, `infinity`) are both degenerate, so the
    /// reductions reject an empty slice rather than return a value
    /// that callers are likely to misuse.
    #[error("empty resistor network: at least one resistor is required")]
    EmptyNetwork,
}

/// Coarse category for a [`ResistorError`], handy for callers that
/// want to branch on the kind of failure without matching every
/// variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The caller supplied an invalid value (out of domain).
    Input,
    /// The requested operation is degenerate for the given network
    /// topology (e.g. empty).
    Topology,
}

impl ResistorError {
    /// Build a [`ResistorError::NonPositive`] for `name`/`value`.
    ///
    /// ```
    /// use valenx_resistor_network::error::ResistorError;
    /// let e = ResistorError::non_positive("resistance", -1.0);
    /// assert_eq!(e.code(), "resistor.non_positive");
    /// ```
    pub fn non_positive(name: &'static str, value: f64) -> Self {
        ResistorError::NonPositive { name, value }
    }

    /// Build a [`ResistorError::NonFinite`] for `name`/`value`.
    pub fn non_finite(name: &'static str, value: f64) -> Self {
        ResistorError::NonFinite { name, value }
    }

    /// Build a [`ResistorError::EmptyNetwork`].
    pub fn empty_network() -> Self {
        ResistorError::EmptyNetwork
    }

    /// Stable, kebab-cased identifier for this error.
    ///
    /// The string is part of the crate's public contract and is
    /// safe to match on or log; the human-readable
    /// [`Display`](std::fmt::Display) text is not.
    pub fn code(&self) -> &'static str {
        match self {
            ResistorError::NonPositive { .. } => "resistor.non_positive",
            ResistorError::NonFinite { .. } => "resistor.non_finite",
            ResistorError::EmptyNetwork => "resistor.empty_network",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            ResistorError::NonPositive { .. } | ResistorError::NonFinite { .. } => {
                ErrorCategory::Input
            }
            ResistorError::EmptyNetwork => ErrorCategory::Topology,
        }
    }
}

/// Validate that `value` is finite and strictly positive, returning
/// it on success.
///
/// This is the shared gate used by every constructor in the crate:
/// a `NaN`/infinite input yields [`ResistorError::NonFinite`] and a
/// zero-or-negative input yields [`ResistorError::NonPositive`].
pub(crate) fn check_positive(name: &'static str, value: f64) -> Result<f64, ResistorError> {
    if !value.is_finite() {
        return Err(ResistorError::non_finite(name, value));
    }
    if value <= 0.0 {
        return Err(ResistorError::non_positive(name, value));
    }
    Ok(value)
}

/// Validate that `value` is finite (it may be zero or negative),
/// returning it on success.
///
/// Used for quantities such as a source voltage, which is allowed
/// to be zero or negative but must still be a real number.
pub(crate) fn check_finite(name: &'static str, value: f64) -> Result<f64, ResistorError> {
    if !value.is_finite() {
        return Err(ResistorError::non_finite(name, value));
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_are_stable() {
        assert_eq!(
            ResistorError::non_positive("resistance", 0.0).code(),
            "resistor.non_positive"
        );
        assert_eq!(
            ResistorError::non_finite("resistance", f64::NAN).code(),
            "resistor.non_finite"
        );
        assert_eq!(
            ResistorError::empty_network().code(),
            "resistor.empty_network"
        );
    }

    #[test]
    fn categories_partition_variants() {
        assert_eq!(
            ResistorError::non_positive("r", -2.0).category(),
            ErrorCategory::Input
        );
        assert_eq!(
            ResistorError::non_finite("r", f64::INFINITY).category(),
            ErrorCategory::Input
        );
        assert_eq!(
            ResistorError::empty_network().category(),
            ErrorCategory::Topology
        );
    }

    #[test]
    fn check_positive_accepts_positive_finite() {
        let v = check_positive("resistance", 47.0).expect("47 ohm is valid");
        assert!((v - 47.0).abs() < 1e-12);
    }

    #[test]
    fn check_positive_rejects_zero_and_negative() {
        assert_eq!(
            check_positive("resistance", 0.0),
            Err(ResistorError::non_positive("resistance", 0.0))
        );
        assert_eq!(
            check_positive("resistance", -5.0),
            Err(ResistorError::non_positive("resistance", -5.0))
        );
    }

    #[test]
    fn check_positive_rejects_non_finite() {
        match check_positive("resistance", f64::NAN) {
            Err(ResistorError::NonFinite { name, .. }) => assert_eq!(name, "resistance"),
            other => panic!("expected NonFinite, got {other:?}"),
        }
        match check_positive("resistance", f64::INFINITY) {
            Err(ResistorError::NonFinite { name, .. }) => assert_eq!(name, "resistance"),
            other => panic!("expected NonFinite, got {other:?}"),
        }
    }

    #[test]
    fn check_finite_allows_zero_and_negative() {
        assert!((check_finite("voltage", 0.0).expect("0 V valid") - 0.0).abs() < 1e-12);
        assert!((check_finite("voltage", -12.0).expect("-12 V valid") - (-12.0)).abs() < 1e-12);
    }

    #[test]
    fn check_finite_rejects_nan() {
        match check_finite("voltage", f64::NAN) {
            Err(ResistorError::NonFinite { name, .. }) => assert_eq!(name, "voltage"),
            other => panic!("expected NonFinite, got {other:?}"),
        }
    }
}
