//! Error taxonomy for the AC power-triangle calculator.
//!
//! All public constructors funnel their validation failures through
//! [`PowerError`]. Every variant carries enough context (the offending
//! parameter name plus the value or reason) to be actionable without a
//! debugger.

use thiserror::Error;

/// Errors raised when validating AC-power inputs.
///
/// These are returned by the validated constructors in
/// [`crate::triangle`] and [`crate::correction`] rather than panicking,
/// so a caller feeding user-supplied numbers can recover gracefully.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum PowerError {
    /// A quantity that must be strictly positive was zero or negative.
    ///
    /// Voltage and current magnitudes, and apparent power, all fall in
    /// this bucket: a meaningful power triangle needs a non-degenerate
    /// source.
    #[error("`{name}` must be strictly positive, got {value}")]
    NonPositive {
        /// Name of the offending parameter (e.g. `"voltage_v"`).
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A quantity that must be non-negative was negative.
    ///
    /// Used for magnitudes that may legitimately be zero (e.g. a real
    /// power of `0 W`) but never negative.
    #[error("`{name}` must be non-negative, got {value}")]
    Negative {
        /// Name of the offending parameter.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A power-factor value fell outside the closed unit interval
    /// `[0, 1]`.
    ///
    /// Power factor is defined as `cos(phi)` for a phase angle `phi` in
    /// `[0, pi/2]`, so by construction it is bounded by `[0, 1]`. A
    /// value outside that range indicates an inconsistent input.
    #[error("power factor must lie in [0, 1], got {value}")]
    PowerFactorOutOfRange {
        /// The rejected value.
        value: f64,
    },

    /// A supplied value was not a finite number (`NaN` or infinity).
    ///
    /// Non-finite inputs propagate silently through floating-point
    /// arithmetic, so they are rejected up front.
    #[error("`{name}` must be finite, got {value}")]
    NotFinite {
        /// Name of the offending parameter.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A power-factor correction targeted a power factor no better than
    /// the starting one.
    ///
    /// Shunt-capacitor correction only ever *raises* the power factor
    /// (it cancels lagging reactive power). Requesting a target that is
    /// the same as, or worse than, the present power factor is rejected
    /// because the resulting capacitor rating would be zero or negative.
    #[error(
        "target power factor {target} must exceed the present power factor {present} for correction"
    )]
    NoCorrectionNeeded {
        /// The present (pre-correction) power factor.
        present: f64,
        /// The requested target power factor.
        target: f64,
    },
}

impl PowerError {
    /// Reject a value that is not a finite number.
    ///
    /// Returns `Ok(value)` when `value` is finite, otherwise a
    /// [`PowerError::NotFinite`]. Used as the first gate in every
    /// validated constructor.
    pub(crate) fn finite(name: &'static str, value: f64) -> Result<f64, PowerError> {
        if value.is_finite() {
            Ok(value)
        } else {
            Err(PowerError::NotFinite { name, value })
        }
    }

    /// Reject a value that is not strictly positive.
    ///
    /// Returns `Ok(value)` when `value` is finite and `> 0`, otherwise
    /// [`PowerError::NotFinite`] or [`PowerError::NonPositive`].
    pub(crate) fn positive(name: &'static str, value: f64) -> Result<f64, PowerError> {
        let value = Self::finite(name, value)?;
        if value > 0.0 {
            Ok(value)
        } else {
            Err(PowerError::NonPositive { name, value })
        }
    }

    /// Reject a value that is negative (zero is allowed).
    ///
    /// Returns `Ok(value)` when `value` is finite and `>= 0`, otherwise
    /// [`PowerError::NotFinite`] or [`PowerError::Negative`].
    pub(crate) fn non_negative(name: &'static str, value: f64) -> Result<f64, PowerError> {
        let value = Self::finite(name, value)?;
        if value >= 0.0 {
            Ok(value)
        } else {
            Err(PowerError::Negative { name, value })
        }
    }

    /// Reject a power factor outside the closed unit interval `[0, 1]`.
    ///
    /// Returns `Ok(value)` when `value` is finite and within `[0, 1]`,
    /// otherwise [`PowerError::NotFinite`] or
    /// [`PowerError::PowerFactorOutOfRange`].
    pub(crate) fn power_factor(name: &'static str, value: f64) -> Result<f64, PowerError> {
        let value = Self::finite(name, value)?;
        if (0.0..=1.0).contains(&value) {
            Ok(value)
        } else {
            Err(PowerError::PowerFactorOutOfRange { value })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn positive_accepts_and_rejects() {
        assert_eq!(PowerError::positive("x", 5.0), Ok(5.0));
        assert_eq!(
            PowerError::positive("x", 0.0),
            Err(PowerError::NonPositive {
                name: "x",
                value: 0.0
            })
        );
        assert_eq!(
            PowerError::positive("x", -1.0),
            Err(PowerError::NonPositive {
                name: "x",
                value: -1.0
            })
        );
    }

    #[test]
    fn non_negative_allows_zero() {
        assert_eq!(PowerError::non_negative("x", 0.0), Ok(0.0));
        assert!(matches!(
            PowerError::non_negative("x", -0.5),
            Err(PowerError::Negative { .. })
        ));
    }

    #[test]
    fn non_finite_is_rejected_first() {
        assert!(matches!(
            PowerError::positive("x", f64::NAN),
            Err(PowerError::NotFinite { .. })
        ));
        assert!(matches!(
            PowerError::power_factor("pf", f64::INFINITY),
            Err(PowerError::NotFinite { .. })
        ));
    }

    #[test]
    fn power_factor_bounds() {
        assert_eq!(PowerError::power_factor("pf", 0.0), Ok(0.0));
        assert_eq!(PowerError::power_factor("pf", 1.0), Ok(1.0));
        assert!(matches!(
            PowerError::power_factor("pf", 1.2),
            Err(PowerError::PowerFactorOutOfRange { .. })
        ));
        assert!(matches!(
            PowerError::power_factor("pf", -0.1),
            Err(PowerError::PowerFactorOutOfRange { .. })
        ));
    }

    #[test]
    fn display_is_informative() {
        let e = PowerError::NonPositive {
            name: "voltage_v",
            value: -3.0,
        };
        let msg = e.to_string();
        assert!(msg.contains("voltage_v"));
        assert!(msg.contains("-3"));
    }
}
