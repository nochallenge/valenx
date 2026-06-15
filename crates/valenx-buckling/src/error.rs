//! Error taxonomy for the column-buckling calculator.
//!
//! Every fallible constructor in this crate funnels its rejection
//! through [`BucklingError`]. The variants distinguish *which* physical
//! quantity was invalid so a caller (CLI, GUI field validator, test)
//! can map the failure back to the offending input.

use thiserror::Error;

/// Errors raised while validating column-buckling inputs.
///
/// All physical inputs to the Euler model must be strictly positive and
/// finite: a non-positive Young's modulus, second moment of area,
/// length, cross-sectional area, or effective-length factor has no
/// physical meaning and would otherwise produce a `NaN`, an infinity,
/// or a sign-flipped "critical load".
#[derive(Debug, Error, Clone, PartialEq)]
pub enum BucklingError {
    /// A quantity that must be strictly positive and finite was not.
    ///
    /// Covers Young's modulus `E`, second moment of area `I`,
    /// unsupported length `L`, cross-sectional area `A`, and the
    /// effective-length factor `K`.
    #[error("`{name}` must be a positive, finite number (got {value})")]
    NonPositive {
        /// Name of the offending quantity (e.g. `"E"`, `"I"`, `"L"`).
        name: &'static str,
        /// The rejected value, echoed back for diagnostics.
        value: f64,
    },
}

impl BucklingError {
    /// Validate that `value` is strictly positive and finite.
    ///
    /// Returns `Ok(value)` on success so it can be used inline inside a
    /// constructor; otherwise yields [`BucklingError::NonPositive`]
    /// tagged with `name`.
    ///
    /// ```
    /// use valenx_buckling::error::BucklingError;
    ///
    /// assert!(BucklingError::require_positive("E", 200.0e9).is_ok());
    /// assert!(BucklingError::require_positive("L", 0.0).is_err());
    /// assert!(BucklingError::require_positive("L", -1.0).is_err());
    /// assert!(BucklingError::require_positive("I", f64::NAN).is_err());
    /// assert!(BucklingError::require_positive("I", f64::INFINITY).is_err());
    /// ```
    pub fn require_positive(name: &'static str, value: f64) -> Result<f64, Self> {
        if value.is_finite() && value > 0.0 {
            Ok(value)
        } else {
            Err(BucklingError::NonPositive { name, value })
        }
    }

    /// Stable kebab-cased identifier for this error, for logs / tests.
    pub fn code(&self) -> &'static str {
        match self {
            BucklingError::NonPositive { .. } => "buckling.non-positive",
        }
    }
}
