//! Error taxonomy for the cable-sag models.
//!
//! Every public constructor in this crate validates its inputs through
//! the helpers here, so a [`CableError`] is the single failure channel.
//! Inputs must be finite and physically admissible (positive where the
//! mechanics demands positivity, in-domain where a closed form has a
//! restricted argument range).

use thiserror::Error;

/// Failure modes for cable sag and tension calculations.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum CableError {
    /// A numeric input was `NaN` or infinite.
    #[error("parameter `{name}` must be finite, got {value}")]
    NotFinite {
        /// Offending parameter name.
        name: &'static str,
        /// The non-finite value supplied.
        value: f64,
    },

    /// A quantity that physics requires to be strictly positive was not.
    ///
    /// Examples: span `L`, horizontal tension `H`, unit load `w`.
    #[error("parameter `{name}` must be > 0, got {value}")]
    NonPositive {
        /// Offending parameter name.
        name: &'static str,
        /// The non-positive value supplied.
        value: f64,
    },

    /// An input fell outside the admissible domain for a closed form.
    ///
    /// Carries the violated bound so callers can report it without
    /// re-deriving the constraint.
    #[error("parameter `{name}` out of domain ({value}): {reason}")]
    OutOfDomain {
        /// Offending parameter name.
        name: &'static str,
        /// The value supplied.
        value: f64,
        /// Human-readable statement of the violated constraint.
        reason: &'static str,
    },
}

impl CableError {
    /// Stable, kebab-cased identifier suitable for logs and tests.
    ///
    /// The string is part of the crate's contract and will not change
    /// for a given variant.
    pub fn code(&self) -> &'static str {
        match self {
            CableError::NotFinite { .. } => "cablesag.not-finite",
            CableError::NonPositive { .. } => "cablesag.non-positive",
            CableError::OutOfDomain { .. } => "cablesag.out-of-domain",
        }
    }
}

/// Reject `NaN` and infinities, returning the value unchanged when finite.
pub(crate) fn finite(name: &'static str, value: f64) -> Result<f64, CableError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(CableError::NotFinite { name, value })
    }
}

/// Require a finite, strictly positive value (e.g. a length or a load).
pub(crate) fn positive(name: &'static str, value: f64) -> Result<f64, CableError> {
    let v = finite(name, value)?;
    if v > 0.0 {
        Ok(v)
    } else {
        Err(CableError::NonPositive { name, value: v })
    }
}

/// Require a finite value within the closed interval `[lo, hi]`.
///
/// Used for arguments such as a sample station `x` constrained to lie on
/// the span. `reason` documents the bound for the error message.
pub(crate) fn in_range(
    name: &'static str,
    value: f64,
    lo: f64,
    hi: f64,
    reason: &'static str,
) -> Result<f64, CableError> {
    let v = finite(name, value)?;
    if v >= lo && v <= hi {
        Ok(v)
    } else {
        Err(CableError::OutOfDomain {
            name,
            value: v,
            reason,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finite_accepts_finite_rejects_nan_inf() {
        assert_eq!(finite("x", 3.5).unwrap(), 3.5);
        assert!(matches!(
            finite("x", f64::NAN),
            Err(CableError::NotFinite { name: "x", .. })
        ));
        assert!(matches!(
            finite("x", f64::INFINITY),
            Err(CableError::NotFinite { .. })
        ));
    }

    #[test]
    fn positive_rejects_zero_and_negative() {
        assert_eq!(positive("L", 10.0).unwrap(), 10.0);
        assert!(matches!(
            positive("L", 0.0),
            Err(CableError::NonPositive { name: "L", .. })
        ));
        assert!(matches!(
            positive("L", -2.0),
            Err(CableError::NonPositive { .. })
        ));
    }

    #[test]
    fn in_range_enforces_bounds() {
        assert_eq!(in_range("x", 5.0, 0.0, 10.0, "0..L").unwrap(), 5.0);
        assert!(matches!(
            in_range("x", 11.0, 0.0, 10.0, "0..L"),
            Err(CableError::OutOfDomain { name: "x", .. })
        ));
        assert!(matches!(
            in_range("x", -1.0, 0.0, 10.0, "0..L"),
            Err(CableError::OutOfDomain { .. })
        ));
    }

    #[test]
    fn codes_are_stable() {
        assert_eq!(
            CableError::NotFinite {
                name: "x",
                value: f64::NAN
            }
            .code(),
            "cablesag.not-finite"
        );
        assert_eq!(
            CableError::NonPositive {
                name: "L",
                value: 0.0
            }
            .code(),
            "cablesag.non-positive"
        );
        assert_eq!(
            CableError::OutOfDomain {
                name: "x",
                value: 9.0,
                reason: "r"
            }
            .code(),
            "cablesag.out-of-domain"
        );
    }
}
