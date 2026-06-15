//! Error taxonomy for open-channel hydraulics.
//!
//! Every fallible constructor in this crate validates its inputs and
//! returns [`OpenChannelError`] on a domain violation (a non-positive
//! geometry dimension, a negative roughness or slope, a non-finite
//! value, or a numerical method that failed to bracket / converge).

use thiserror::Error;

/// Errors raised while building channel geometry or solving flow.
///
/// Construct these through the validated helpers
/// ([`OpenChannelError::non_positive`], [`OpenChannelError::negative`],
/// [`OpenChannelError::not_finite`]) rather than the variants directly,
/// so the message wording stays consistent across the crate.
#[derive(Debug, Error, Clone, PartialEq)]
#[non_exhaustive]
pub enum OpenChannelError {
    /// A quantity that must be strictly positive was zero or negative.
    #[error("parameter `{name}` must be > 0, got {value}")]
    NonPositive {
        /// Offending parameter name (static, stable identifier).
        name: &'static str,
        /// The value that violated the constraint.
        value: f64,
    },

    /// A quantity that must be non-negative was negative.
    #[error("parameter `{name}` must be >= 0, got {value}")]
    Negative {
        /// Offending parameter name (static, stable identifier).
        name: &'static str,
        /// The value that violated the constraint.
        value: f64,
    },

    /// A value was `NaN` or infinite where a finite number is required.
    #[error("parameter `{name}` must be finite, got {value}")]
    NotFinite {
        /// Offending parameter name (static, stable identifier).
        name: &'static str,
        /// The non-finite value (`NaN` or `±inf`).
        value: f64,
    },

    /// An iterative solver failed to bracket a root or to converge
    /// within its iteration budget.
    #[error("numerical solver did not converge: {0}")]
    Convergence(String),
}

impl OpenChannelError {
    /// Validate that `value` is finite and strictly positive.
    ///
    /// Returns `Ok(value)` when the constraint holds. Yields
    /// [`OpenChannelError::NotFinite`] for a `NaN` / infinite input and
    /// [`OpenChannelError::NonPositive`] for a zero / negative input.
    pub fn non_positive(name: &'static str, value: f64) -> Result<f64, Self> {
        if !value.is_finite() {
            return Err(Self::NotFinite { name, value });
        }
        if value <= 0.0 {
            return Err(Self::NonPositive { name, value });
        }
        Ok(value)
    }

    /// Validate that `value` is finite and non-negative.
    ///
    /// Returns `Ok(value)` when the constraint holds. Yields
    /// [`OpenChannelError::NotFinite`] for a `NaN` / infinite input and
    /// [`OpenChannelError::Negative`] for a negative input. Zero is
    /// accepted.
    pub fn negative(name: &'static str, value: f64) -> Result<f64, Self> {
        if !value.is_finite() {
            return Err(Self::NotFinite { name, value });
        }
        if value < 0.0 {
            return Err(Self::Negative { name, value });
        }
        Ok(value)
    }

    /// Validate that `value` is finite (rejects `NaN` and `±inf`).
    ///
    /// Returns `Ok(value)` when finite, otherwise
    /// [`OpenChannelError::NotFinite`].
    pub fn not_finite(name: &'static str, value: f64) -> Result<f64, Self> {
        if value.is_finite() {
            Ok(value)
        } else {
            Err(Self::NotFinite { name, value })
        }
    }

    /// Stable kebab-cased identifier for the variant, suitable for logs
    /// and machine matching.
    pub fn code(&self) -> &'static str {
        match self {
            Self::NonPositive { .. } => "openchannel.non-positive",
            Self::Negative { .. } => "openchannel.negative",
            Self::NotFinite { .. } => "openchannel.not-finite",
            Self::Convergence(_) => "openchannel.convergence",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn non_positive_accepts_positive() {
        assert_eq!(OpenChannelError::non_positive("x", 2.5).unwrap(), 2.5);
    }

    #[test]
    fn non_positive_rejects_zero_and_negative() {
        assert!(matches!(
            OpenChannelError::non_positive("x", 0.0),
            Err(OpenChannelError::NonPositive { .. })
        ));
        assert!(matches!(
            OpenChannelError::non_positive("x", -1.0),
            Err(OpenChannelError::NonPositive { .. })
        ));
    }

    #[test]
    fn negative_accepts_zero_rejects_negative() {
        assert_eq!(OpenChannelError::negative("s", 0.0).unwrap(), 0.0);
        assert!(matches!(
            OpenChannelError::negative("s", -0.5),
            Err(OpenChannelError::Negative { .. })
        ));
    }

    #[test]
    fn non_finite_inputs_are_rejected() {
        assert!(matches!(
            OpenChannelError::non_positive("x", f64::NAN),
            Err(OpenChannelError::NotFinite { .. })
        ));
        assert!(matches!(
            OpenChannelError::negative("x", f64::INFINITY),
            Err(OpenChannelError::NotFinite { .. })
        ));
        assert!(matches!(
            OpenChannelError::not_finite("x", f64::NEG_INFINITY),
            Err(OpenChannelError::NotFinite { .. })
        ));
    }

    #[test]
    fn codes_are_stable() {
        assert_eq!(
            OpenChannelError::NonPositive {
                name: "x",
                value: -1.0
            }
            .code(),
            "openchannel.non-positive"
        );
        assert_eq!(
            OpenChannelError::Convergence("boom".into()).code(),
            "openchannel.convergence"
        );
    }
}
