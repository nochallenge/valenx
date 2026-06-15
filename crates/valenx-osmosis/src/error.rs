//! Error taxonomy for the osmosis / fluid-balance models.
//!
//! Every fallible constructor and computation in this crate returns
//! [`Result<_, OsmosisError>`](OsmosisError). The error carries stable
//! [`code`](OsmosisError::code) and [`category`](OsmosisError::category)
//! accessors so callers (telemetry, a GUI surface) can branch on the
//! failure class without string-matching the human-readable message.

use thiserror::Error;

/// Errors raised when building or evaluating an osmosis model.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum OsmosisError {
    /// A scalar input fell outside its physically meaningful domain.
    ///
    /// Carries the offending parameter `name`, the `value` that was
    /// supplied, and a short `reason` describing the constraint.
    #[error("invalid parameter `{name}` = {value}: {reason}")]
    InvalidParameter {
        /// Parameter identifier (e.g. `"concentration"`, `"temperature_k"`).
        name: &'static str,
        /// The numeric value that violated the constraint.
        value: f64,
        /// Why the value is rejected (e.g. `"must be > 0"`).
        reason: &'static str,
    },
}

impl OsmosisError {
    /// Construct an [`OsmosisError::InvalidParameter`].
    ///
    /// Small helper so the validated constructors elsewhere in the crate
    /// stay terse.
    pub fn invalid(name: &'static str, value: f64, reason: &'static str) -> Self {
        OsmosisError::InvalidParameter {
            name,
            value,
            reason,
        }
    }

    /// Stable, kebab/dot-cased identifier for this error.
    ///
    /// Intended for logs and machine consumers; never localized.
    pub fn code(&self) -> &'static str {
        match self {
            OsmosisError::InvalidParameter { .. } => "osmosis.invalid_parameter",
        }
    }

    /// Coarse category for routing / metrics.
    pub fn category(&self) -> ErrorCategory {
        match self {
            OsmosisError::InvalidParameter { .. } => ErrorCategory::Input,
        }
    }
}

/// Coarse classification of an [`OsmosisError`].
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// A user-supplied value was out of range or malformed.
    Input,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn invalid_helper_populates_fields() {
        let e = OsmosisError::invalid("concentration", -1.0, "must be > 0");
        match e {
            OsmosisError::InvalidParameter {
                name,
                value,
                reason,
            } => {
                assert_eq!(name, "concentration");
                assert!((value - (-1.0)).abs() < 1e-12);
                assert_eq!(reason, "must be > 0");
            }
        }
    }

    #[test]
    fn code_and_category_are_stable() {
        let e = OsmosisError::invalid("temperature_k", 0.0, "must be > 0");
        assert_eq!(e.code(), "osmosis.invalid_parameter");
        assert_eq!(e.category(), ErrorCategory::Input);
    }

    #[test]
    fn display_includes_name_value_and_reason() {
        let e = OsmosisError::invalid("sigma", 2.0, "must be in [0, 1]");
        let msg = format!("{e}");
        assert!(msg.contains("sigma"));
        assert!(msg.contains("must be in [0, 1]"));
    }
}
