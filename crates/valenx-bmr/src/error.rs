//! Error taxonomy for `valenx-bmr`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, BmrError>`]. The variants are intentionally small — a
//! BMR / energy-balance caller only ever passes anthropometric numbers
//! and a couple of multipliers, so the failure modes are:
//!
//! 1. A value is outside its physical domain — a non-positive or
//!    NaN mass, height or age, an activity factor below `1.0`
//!    ([`BmrError::OutOfRange`]).
//! 2. A value is non-finite where a finite number is required
//!    ([`BmrError::NotFinite`]).
//!
//! Use [`BmrError::code`] for stable log / telemetry tagging. The
//! pattern mirrors `valenx-springs`'s `SpringsError`.

use thiserror::Error;

/// Errors raised by the BMR / energy-balance calculators.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum BmrError {
    /// A parameter fell outside its accepted range. `name` is the
    /// logical parameter (`"mass_kg"`, `"activity_factor"`, …),
    /// `value` is what was supplied, and `reason` explains the bound.
    #[error("parameter `{name}` = {value} is out of range: {reason}")]
    OutOfRange {
        /// Logical parameter name.
        name: &'static str,
        /// The offending value.
        value: f64,
        /// Human-readable description of the violated bound.
        reason: &'static str,
    },

    /// A parameter was `NaN` or `±∞` where a finite number is required.
    #[error("parameter `{name}` must be a finite number, got {value}")]
    NotFinite {
        /// Logical parameter name.
        name: &'static str,
        /// The offending value (`NaN` or infinite).
        value: f64,
    },
}

/// Coarse error category for routing / display.
///
/// Stable across crate versions — `match` on this rather than on the
/// individual variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// User-supplied input is invalid (bad value or non-finite).
    Input,
}

impl BmrError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"bmr.<sub_id>"`. Codes never change across
    /// minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            BmrError::OutOfRange { .. } => "bmr.out_of_range",
            BmrError::NotFinite { .. } => "bmr.not_finite",
        }
    }

    /// Coarse category — every variant in this crate is caller input.
    pub fn category(&self) -> ErrorCategory {
        match self {
            BmrError::OutOfRange { .. } | BmrError::NotFinite { .. } => ErrorCategory::Input,
        }
    }

    /// Validate that `value` is finite and strictly greater than zero.
    ///
    /// Used for masses, heights and ages, which are all strictly
    /// positive quantities. Returns the value unchanged on success so
    /// it can be used inline.
    ///
    /// # Errors
    ///
    /// [`BmrError::NotFinite`] if `value` is `NaN` or infinite;
    /// [`BmrError::OutOfRange`] if `value <= 0`.
    pub fn require_positive(name: &'static str, value: f64) -> Result<f64> {
        if !value.is_finite() {
            return Err(BmrError::NotFinite { name, value });
        }
        if value <= 0.0 {
            return Err(BmrError::OutOfRange {
                name,
                value,
                reason: "must be strictly positive",
            });
        }
        Ok(value)
    }

    /// Validate a physical-activity multiplier: finite and `>= 1.0`.
    ///
    /// A factor of `1.0` is pure basal expenditure (complete bed rest);
    /// activity can only add to it, so values below `1.0` are rejected.
    ///
    /// # Errors
    ///
    /// [`BmrError::NotFinite`] if `value` is `NaN` or infinite;
    /// [`BmrError::OutOfRange`] if `value < 1.0`.
    pub fn require_activity_factor(value: f64) -> Result<f64> {
        let name = "activity_factor";
        if !value.is_finite() {
            return Err(BmrError::NotFinite { name, value });
        }
        if value < 1.0 {
            return Err(BmrError::OutOfRange {
                name,
                value,
                reason: "must be at least 1.0 (basal expenditure)",
            });
        }
        Ok(value)
    }

    /// Validate any value that merely has to be finite (it may be
    /// negative, e.g. an energy balance that represents a deficit).
    ///
    /// # Errors
    ///
    /// [`BmrError::NotFinite`] if `value` is `NaN` or infinite.
    pub fn require_finite(name: &'static str, value: f64) -> Result<f64> {
        if !value.is_finite() {
            return Err(BmrError::NotFinite { name, value });
        }
        Ok(value)
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, BmrError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = BmrError::OutOfRange {
            name: "mass_kg",
            value: -1.0,
            reason: "must be strictly positive",
        };
        assert_eq!(err.code(), "bmr.out_of_range");
        assert_eq!(err.category(), ErrorCategory::Input);

        let err = BmrError::NotFinite {
            name: "height_cm",
            value: f64::NAN,
        };
        assert_eq!(err.code(), "bmr.not_finite");
        assert_eq!(err.category(), ErrorCategory::Input);
    }

    #[test]
    fn require_positive_accepts_and_rejects() {
        assert_eq!(BmrError::require_positive("mass_kg", 80.0), Ok(80.0));

        let too_small = BmrError::require_positive("mass_kg", 0.0);
        assert!(matches!(
            too_small,
            Err(BmrError::OutOfRange {
                name: "mass_kg",
                ..
            })
        ));

        let nan = BmrError::require_positive("mass_kg", f64::NAN);
        assert!(matches!(nan, Err(BmrError::NotFinite { .. })));

        let inf = BmrError::require_positive("mass_kg", f64::INFINITY);
        assert!(matches!(inf, Err(BmrError::NotFinite { .. })));
    }

    #[test]
    fn require_activity_factor_floor_is_one() {
        assert_eq!(BmrError::require_activity_factor(1.0), Ok(1.0));
        assert_eq!(BmrError::require_activity_factor(1.55), Ok(1.55));

        let below = BmrError::require_activity_factor(0.9);
        assert!(matches!(below, Err(BmrError::OutOfRange { .. })));
    }

    #[test]
    fn require_finite_allows_negative() {
        // A negative energy balance (deficit) is a valid finite input.
        assert_eq!(
            BmrError::require_finite("energy_balance_kcal", -500.0),
            Ok(-500.0)
        );
        let nan = BmrError::require_finite("energy_balance_kcal", f64::NAN);
        assert!(matches!(nan, Err(BmrError::NotFinite { .. })));
    }

    #[test]
    fn display_is_informative() {
        let msg = BmrError::OutOfRange {
            name: "activity_factor",
            value: 0.5,
            reason: "must be at least 1.0 (basal expenditure)",
        }
        .to_string();
        assert!(msg.contains("activity_factor"), "got: {msg}");
        assert!(msg.contains("0.5"), "got: {msg}");

        let msg = BmrError::NotFinite {
            name: "mass_kg",
            value: f64::INFINITY,
        }
        .to_string();
        assert!(msg.contains("mass_kg"), "got: {msg}");
        assert!(msg.contains("finite"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(BmrError::NotFinite {
            name: "age_years",
            value: f64::NAN,
        });
        assert!(err.to_string().contains("age_years"));
    }
}
