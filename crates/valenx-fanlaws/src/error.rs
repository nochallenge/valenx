//! Error taxonomy for the fan-laws crate.
//!
//! Every fallible constructor in this crate validates its inputs and
//! returns a [`FanLawError`] on bad data rather than panicking, so the
//! affinity / power routines downstream can assume finite, physically
//! admissible operating points.

use thiserror::Error;

/// Errors raised when constructing or scaling fan operating points.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum FanLawError {
    /// A quantity that must be strictly positive was zero or negative
    /// (e.g. a rotational speed, density, or impeller diameter).
    #[error("`{name}` must be strictly positive, got {value}")]
    NonPositive {
        /// Name of the offending quantity.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A quantity that must be non-negative was negative (e.g. a
    /// volumetric flow rate or a pressure rise, either of which may
    /// legitimately be zero at shut-off / no-flow).
    #[error("`{name}` must be non-negative, got {value}")]
    Negative {
        /// Name of the offending quantity.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// An efficiency outside the physically admissible half-open
    /// interval `(0, 1]`. A fan can never have zero efficiency (the
    /// power identity would divide by zero) nor exceed unity (that
    /// would extract more flow work than the shaft supplies).
    #[error("efficiency must lie in (0, 1], got {value}")]
    EfficiencyOutOfRange {
        /// The rejected efficiency.
        value: f64,
    },

    /// A supplied value was not finite (NaN or infinity), which can
    /// never describe a real operating point.
    #[error("`{name}` must be finite, got {value}")]
    NotFinite {
        /// Name of the offending quantity.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },
}

/// Coarse category for a [`FanLawError`], mirroring the convention used
/// by the sibling parametric crates so a host application can route all
/// crate errors through one banner.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The error originates in user-supplied input.
    Input,
    /// The error originates in an out-of-range tunable knob.
    Config,
}

impl FanLawError {
    /// Stable kebab-cased identifier for this error, suitable for logs
    /// and test assertions that should not break when the human-facing
    /// message is reworded.
    pub fn code(&self) -> &'static str {
        match self {
            FanLawError::NonPositive { .. } => "fanlaws.non_positive",
            FanLawError::Negative { .. } => "fanlaws.negative",
            FanLawError::EfficiencyOutOfRange { .. } => "fanlaws.efficiency_out_of_range",
            FanLawError::NotFinite { .. } => "fanlaws.not_finite",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            FanLawError::EfficiencyOutOfRange { .. } => ErrorCategory::Config,
            FanLawError::NonPositive { .. }
            | FanLawError::Negative { .. }
            | FanLawError::NotFinite { .. } => ErrorCategory::Input,
        }
    }
}

/// Validate that `value` is finite, returning [`FanLawError::NotFinite`]
/// otherwise. Shared by the typed constructors below.
pub(crate) fn require_finite(name: &'static str, value: f64) -> Result<f64, FanLawError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(FanLawError::NotFinite { name, value })
    }
}

/// Validate that `value` is finite and strictly positive.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64, FanLawError> {
    let value = require_finite(name, value)?;
    if value > 0.0 {
        Ok(value)
    } else {
        Err(FanLawError::NonPositive { name, value })
    }
}

/// Validate that `value` is finite and non-negative.
pub(crate) fn require_non_negative(name: &'static str, value: f64) -> Result<f64, FanLawError> {
    let value = require_finite(name, value)?;
    if value >= 0.0 {
        Ok(value)
    } else {
        Err(FanLawError::Negative { name, value })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn require_positive_rejects_zero_and_negative() {
        assert!(matches!(
            require_positive("speed", 0.0),
            Err(FanLawError::NonPositive { name: "speed", .. })
        ));
        assert!(matches!(
            require_positive("speed", -3.0),
            Err(FanLawError::NonPositive { .. })
        ));
        assert_eq!(require_positive("speed", 1500.0).unwrap(), 1500.0);
    }

    #[test]
    fn require_non_negative_allows_zero_but_not_negative() {
        assert_eq!(require_non_negative("flow", 0.0).unwrap(), 0.0);
        assert!(matches!(
            require_non_negative("flow", -1.0),
            Err(FanLawError::Negative { name: "flow", .. })
        ));
    }

    #[test]
    fn non_finite_inputs_are_rejected_first() {
        assert!(matches!(
            require_positive("density", f64::NAN),
            Err(FanLawError::NotFinite {
                name: "density",
                ..
            })
        ));
        assert!(matches!(
            require_non_negative("pressure", f64::INFINITY),
            Err(FanLawError::NotFinite { .. })
        ));
    }

    #[test]
    fn codes_and_categories_are_stable() {
        assert_eq!(
            FanLawError::NonPositive {
                name: "n",
                value: 0.0
            }
            .code(),
            "fanlaws.non_positive"
        );
        assert_eq!(
            FanLawError::EfficiencyOutOfRange { value: 1.5 }.category(),
            ErrorCategory::Config
        );
        assert_eq!(
            FanLawError::Negative {
                name: "q",
                value: -1.0
            }
            .category(),
            ErrorCategory::Input
        );
    }
}
