//! Error taxonomy for the linear-elastic fracture-mechanics calculators.

use thiserror::Error;

/// Shorthand for `Result<T, FractureError>`.
pub type Result<T> = core::result::Result<T, FractureError>;

/// Anything that can go wrong constructing a fracture input or evaluating
/// a closed-form LEFM expression.
///
/// Every variant corresponds to a non-physical input that would otherwise
/// feed a silent `NaN` / `Inf` into a `√`, a division, or a power — for
/// example a negative crack length under the square root, a zero applied
/// stress in the critical-crack denominator, or a non-finite material
/// property. The constructors in [`crate::material`] and the free
/// functions in [`crate::crack`] validate up front and return one of these
/// rather than propagating garbage.
///
/// This enum is `#[non_exhaustive]`: new variants may be added in a future
/// release without it being a breaking change, so downstream `match` arms
/// must include a wildcard.
#[derive(Debug, Error, Clone, PartialEq)]
#[non_exhaustive]
pub enum FractureError {
    /// A length quantity (crack size, plate dimension) was negative,
    /// or non-finite. Carries the offending field name and value.
    ///
    /// Crack length `a` is permitted to be exactly zero (an uncracked
    /// body has `K = 0`); it is only rejected when it appears in a
    /// denominator. Other length-like inputs are rejected at zero too.
    #[error("invalid length `{name}` = {value} (must be finite and non-negative)")]
    InvalidLength {
        /// Which length field was bad (e.g. `"crack_length"`).
        name: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A quantity that must be strictly positive (fracture toughness,
    /// yield strength, geometry factor, or an applied stress that sits in
    /// a denominator) was zero, negative, or non-finite. Carries the
    /// offending field name and value.
    #[error("invalid positive quantity `{name}` = {value} (must be finite and > 0)")]
    NonPositive {
        /// Which field was bad (e.g. `"fracture_toughness"`).
        name: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A stress value was supplied that is finite but non-physical for the
    /// requested operation — currently only a negative applied stress fed
    /// to a stress-intensity evaluation, where the Mode-I formula assumes a
    /// non-negative (opening) far-field tension. Carries field and value.
    #[error("invalid stress `{name}` = {value} (must be finite and non-negative)")]
    InvalidStress {
        /// Which stress field was bad (e.g. `"applied_stress"`).
        name: &'static str,
        /// The offending value.
        value: f64,
    },
}

impl FractureError {
    /// Stable kebab-cased identifier for the variant, suitable for logging
    /// or mapping to a UI message catalogue.
    ///
    /// ```
    /// use valenx_fracture::FractureError;
    /// let e = FractureError::NonPositive { name: "geometry_factor", value: 0.0 };
    /// assert_eq!(e.code(), "fracture.non_positive");
    /// ```
    pub fn code(&self) -> &'static str {
        match self {
            FractureError::InvalidLength { .. } => "fracture.invalid_length",
            FractureError::NonPositive { .. } => "fracture.non_positive",
            FractureError::InvalidStress { .. } => "fracture.invalid_stress",
        }
    }
}
