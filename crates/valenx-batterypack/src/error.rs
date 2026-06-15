//! Battery-pack configuration error taxonomy.
//!
//! Every fallible constructor in this crate funnels through
//! [`BatteryPackError`]. The variants carry the offending parameter
//! name and a human-readable reason, plus stable
//! [`code`](BatteryPackError::code) /
//! [`category`](BatteryPackError::category) accessors suitable for
//! telemetry or UI grouping.

use thiserror::Error;

/// Errors raised while building or sizing a battery pack.
#[derive(Debug, Error)]
pub enum BatteryPackError {
    /// A continuous physical quantity was non-finite (`NaN` / `inf`) or
    /// outside its required range (e.g. a non-positive voltage or
    /// capacity).
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Name of the offending parameter.
        name: &'static str,
        /// Why the value was rejected.
        reason: String,
    },

    /// An integer count (series or parallel multiplicity) was zero. A
    /// pack must contain at least one cell in each dimension.
    #[error("bad count `{name}`: {reason}")]
    BadCount {
        /// Name of the offending count.
        name: &'static str,
        /// Why the count was rejected.
        reason: String,
    },
}

/// Coarse error category, useful for grouping in a UI or in telemetry.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The caller supplied an invalid value.
    Input,
}

impl BatteryPackError {
    /// Construct a [`BadParameter`](BatteryPackError::BadParameter)
    /// error with the given parameter name and reason.
    pub fn bad_parameter(name: &'static str, reason: impl Into<String>) -> Self {
        BatteryPackError::BadParameter {
            name,
            reason: reason.into(),
        }
    }

    /// Construct a [`BadCount`](BatteryPackError::BadCount) error with
    /// the given count name and reason.
    pub fn bad_count(name: &'static str, reason: impl Into<String>) -> Self {
        BatteryPackError::BadCount {
            name,
            reason: reason.into(),
        }
    }

    /// Stable kebab-cased identifier for this error, suitable for logs
    /// and metrics. Stable across releases.
    pub fn code(&self) -> &'static str {
        match self {
            BatteryPackError::BadParameter { .. } => "batterypack.bad_parameter",
            BatteryPackError::BadCount { .. } => "batterypack.bad_count",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            BatteryPackError::BadParameter { .. } | BatteryPackError::BadCount { .. } => {
                ErrorCategory::Input
            }
        }
    }
}

/// Validate that `value` is finite and strictly greater than zero,
/// returning a [`BadParameter`](BatteryPackError::BadParameter) error
/// naming `name` otherwise. Shared by the cell and rate constructors.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64, BatteryPackError> {
    if !value.is_finite() {
        return Err(BatteryPackError::bad_parameter(
            name,
            format!("must be finite, got {value}"),
        ));
    }
    if value <= 0.0 {
        return Err(BatteryPackError::bad_parameter(
            name,
            format!("must be > 0, got {value}"),
        ));
    }
    Ok(value)
}

/// Validate that `value` is finite and not negative (zero is allowed),
/// returning a [`BadParameter`](BatteryPackError::BadParameter) error
/// naming `name` otherwise.
pub(crate) fn require_non_negative(
    name: &'static str,
    value: f64,
) -> Result<f64, BatteryPackError> {
    if !value.is_finite() {
        return Err(BatteryPackError::bad_parameter(
            name,
            format!("must be finite, got {value}"),
        ));
    }
    if value < 0.0 {
        return Err(BatteryPackError::bad_parameter(
            name,
            format!("must be >= 0, got {value}"),
        ));
    }
    Ok(value)
}

/// Validate that an integer count is at least one, returning a
/// [`BadCount`](BatteryPackError::BadCount) error naming `name`
/// otherwise.
pub(crate) fn require_at_least_one(
    name: &'static str,
    value: u32,
) -> Result<u32, BatteryPackError> {
    if value == 0 {
        return Err(BatteryPackError::bad_count(
            name,
            "must be >= 1, got 0".to_string(),
        ));
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn require_positive_accepts_positive() {
        assert!((require_positive("x", 3.7).unwrap() - 3.7).abs() < 1e-12);
    }

    #[test]
    fn require_positive_rejects_zero_and_negative() {
        assert!(matches!(
            require_positive("x", 0.0),
            Err(BatteryPackError::BadParameter { .. })
        ));
        assert!(matches!(
            require_positive("x", -1.0),
            Err(BatteryPackError::BadParameter { .. })
        ));
    }

    #[test]
    fn require_positive_rejects_non_finite() {
        assert!(matches!(
            require_positive("x", f64::NAN),
            Err(BatteryPackError::BadParameter { .. })
        ));
        assert!(matches!(
            require_positive("x", f64::INFINITY),
            Err(BatteryPackError::BadParameter { .. })
        ));
    }

    #[test]
    fn require_non_negative_accepts_zero() {
        assert!((require_non_negative("x", 0.0).unwrap()).abs() < 1e-12);
    }

    #[test]
    fn require_at_least_one_rejects_zero() {
        assert!(matches!(
            require_at_least_one("series", 0),
            Err(BatteryPackError::BadCount { .. })
        ));
        assert_eq!(require_at_least_one("series", 1).unwrap(), 1);
    }

    #[test]
    fn code_and_category_are_stable() {
        let p = BatteryPackError::bad_parameter("v", "bad");
        assert_eq!(p.code(), "batterypack.bad_parameter");
        assert_eq!(p.category(), ErrorCategory::Input);

        let c = BatteryPackError::bad_count("n", "bad");
        assert_eq!(c.code(), "batterypack.bad_count");
        assert_eq!(c.category(), ErrorCategory::Input);
    }
}
