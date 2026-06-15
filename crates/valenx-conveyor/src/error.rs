//! Belt-conveyor error taxonomy.
//!
//! A single [`ConveyorError`] enum covers every failure raised by the
//! crate. Constructors that perform input validation (for example
//! [`crate::belt::Belt::new`]) return `Result<_, ConveyorError>` so an
//! invalid physical parameter is rejected eagerly rather than producing
//! a silently wrong number downstream.

use thiserror::Error;

/// Errors raised while building belt-conveyor inputs or evaluating the
/// closed-form models.
///
/// Every variant carries enough context to report which quantity was
/// rejected and why. Use [`ConveyorError::code`] for a stable,
/// machine-readable identifier and [`ConveyorError::category`] for a
/// coarse grouping suitable for UI surfaces.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum ConveyorError {
    /// A parameter that must be strictly positive was zero or negative.
    ///
    /// Examples: belt width, belt speed, bulk density, load area.
    #[error("parameter `{name}` must be strictly positive, got {value}")]
    NonPositive {
        /// Offending parameter name (stable, `snake_case`).
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A parameter that must be finite (not NaN / not infinite) was not.
    #[error("parameter `{name}` must be finite, got {value}")]
    NotFinite {
        /// Offending parameter name (stable, `snake_case`).
        name: &'static str,
        /// The rejected value.
        value: f64,
    },

    /// A parameter fell outside an inclusive `[min, max]` range.
    #[error("parameter `{name}` must be in [{min}, {max}], got {value}")]
    OutOfRange {
        /// Offending parameter name (stable, `snake_case`).
        name: &'static str,
        /// The rejected value.
        value: f64,
        /// Inclusive lower bound.
        min: f64,
        /// Inclusive upper bound.
        max: f64,
    },
}

impl ConveyorError {
    /// Validate that `value` is finite and strictly positive.
    ///
    /// Returns `value` unchanged on success. Used by the validated
    /// constructors throughout the crate so the positivity rule lives in
    /// exactly one place.
    ///
    /// # Errors
    ///
    /// Returns [`ConveyorError::NotFinite`] if `value` is NaN or
    /// infinite, or [`ConveyorError::NonPositive`] if `value <= 0`.
    pub fn require_positive(name: &'static str, value: f64) -> Result<f64, ConveyorError> {
        if !value.is_finite() {
            return Err(ConveyorError::NotFinite { name, value });
        }
        if value <= 0.0 {
            return Err(ConveyorError::NonPositive { name, value });
        }
        Ok(value)
    }

    /// Validate that `value` is finite and within the inclusive range
    /// `[min, max]`.
    ///
    /// Returns `value` unchanged on success.
    ///
    /// # Errors
    ///
    /// Returns [`ConveyorError::NotFinite`] if `value` is NaN or
    /// infinite, or [`ConveyorError::OutOfRange`] if it falls outside
    /// `[min, max]`.
    pub fn require_range(
        name: &'static str,
        value: f64,
        min: f64,
        max: f64,
    ) -> Result<f64, ConveyorError> {
        if !value.is_finite() {
            return Err(ConveyorError::NotFinite { name, value });
        }
        if value < min || value > max {
            return Err(ConveyorError::OutOfRange {
                name,
                value,
                min,
                max,
            });
        }
        Ok(value)
    }

    /// Stable kebab-cased identifier for the variant.
    ///
    /// Suitable for logs, telemetry keys, or localization lookups; the
    /// string is part of the crate's public contract and will not change
    /// for an existing variant.
    pub fn code(&self) -> &'static str {
        match self {
            ConveyorError::NonPositive { .. } => "conveyor.non-positive",
            ConveyorError::NotFinite { .. } => "conveyor.not-finite",
            ConveyorError::OutOfRange { .. } => "conveyor.out-of-range",
        }
    }

    /// Coarse category for grouping in a UI.
    pub fn category(&self) -> ErrorCategory {
        match self {
            ConveyorError::NonPositive { .. }
            | ConveyorError::NotFinite { .. }
            | ConveyorError::OutOfRange { .. } => ErrorCategory::Input,
        }
    }
}

/// Coarse grouping of [`ConveyorError`] variants for presentation.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
#[non_exhaustive]
pub enum ErrorCategory {
    /// The caller supplied an invalid physical input.
    Input,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn require_positive_accepts_positive() {
        let v = ConveyorError::require_positive("x", 2.5).expect("2.5 is positive");
        assert!((v - 2.5).abs() < 1e-12);
    }

    #[test]
    fn require_positive_rejects_zero_and_negative() {
        let zero = ConveyorError::require_positive("x", 0.0);
        assert!(matches!(
            zero,
            Err(ConveyorError::NonPositive { name: "x", .. })
        ));
        let neg = ConveyorError::require_positive("x", -1.0);
        assert!(matches!(neg, Err(ConveyorError::NonPositive { .. })));
    }

    #[test]
    fn require_positive_rejects_nan_and_inf() {
        assert!(matches!(
            ConveyorError::require_positive("x", f64::NAN),
            Err(ConveyorError::NotFinite { .. })
        ));
        assert!(matches!(
            ConveyorError::require_positive("x", f64::INFINITY),
            Err(ConveyorError::NotFinite { .. })
        ));
    }

    #[test]
    fn require_range_enforces_bounds() {
        assert!(ConveyorError::require_range("a", 0.5, 0.0, 1.0).is_ok());
        assert!(matches!(
            ConveyorError::require_range("a", 1.5, 0.0, 1.0),
            Err(ConveyorError::OutOfRange { .. })
        ));
        assert!(matches!(
            ConveyorError::require_range("a", -0.1, 0.0, 1.0),
            Err(ConveyorError::OutOfRange { .. })
        ));
    }

    #[test]
    fn codes_are_stable_and_distinct() {
        let a = ConveyorError::NonPositive {
            name: "x",
            value: 0.0,
        };
        let b = ConveyorError::NotFinite {
            name: "x",
            value: f64::NAN,
        };
        let c = ConveyorError::OutOfRange {
            name: "x",
            value: 2.0,
            min: 0.0,
            max: 1.0,
        };
        assert_eq!(a.code(), "conveyor.non-positive");
        assert_eq!(b.code(), "conveyor.not-finite");
        assert_eq!(c.code(), "conveyor.out-of-range");
        assert_eq!(a.category(), ErrorCategory::Input);
        assert_eq!(c.category(), ErrorCategory::Input);
    }

    #[test]
    fn display_includes_parameter_name() {
        let e = ConveyorError::NonPositive {
            name: "belt_speed",
            value: -3.0,
        };
        let msg = format!("{e}");
        assert!(msg.contains("belt_speed"));
    }
}
