//! Error taxonomy for the Rankine retaining-wall model.
//!
//! Every fallible public entry point validates its inputs through one of
//! the constructors here, so a [`RetainingWallError`] is always the
//! result of an out-of-domain argument rather than an internal panic.

use thiserror::Error;

/// Errors raised when constructing or evaluating a Rankine earth-pressure
/// model.
///
/// The Rankine closed form is only physically meaningful for a soil
/// friction angle strictly inside `0 <= phi < 90 degrees`, a strictly
/// positive unit weight, and a non-negative depth/height. Each variant
/// names the offending quantity and the value that violated the domain.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum RetainingWallError {
    /// Soil friction angle `phi` was outside the half-open interval
    /// `[0, 90)` degrees.
    ///
    /// At `phi = 90 deg` the active coefficient collapses to zero and the
    /// passive coefficient diverges, so the upper bound is exclusive.
    #[error("friction angle phi = {phi_deg} deg is out of range; expected 0 <= phi < 90")]
    FrictionAngleOutOfRange {
        /// The offending friction angle, in degrees.
        phi_deg: f64,
    },

    /// Soil unit weight `gamma` was not strictly positive.
    #[error("unit weight gamma = {gamma} must be strictly positive")]
    NonPositiveUnitWeight {
        /// The offending unit weight (caller's units, e.g. kN/m^3).
        gamma: f64,
    },

    /// A wall height or depth was negative.
    #[error("depth/height value = {value} must be non-negative")]
    NegativeDepth {
        /// The offending depth or height (caller's length units, e.g. m).
        value: f64,
    },

    /// A target lateral thrust was negative.
    ///
    /// Thrust magnitudes are non-negative — a zero thrust corresponds to a
    /// zero-height wall — so a negative target has no real wall height.
    #[error("thrust value = {value} must be non-negative")]
    NegativeThrust {
        /// The offending thrust (caller's force-per-length units, e.g.
        /// kN/m).
        value: f64,
    },

    /// A supplied quantity was not a finite number (NaN or infinite).
    #[error("value `{name}` must be finite, got {value}")]
    NonFinite {
        /// Name of the offending parameter.
        name: &'static str,
        /// The non-finite value that was supplied.
        value: f64,
    },
}

/// Coarse classification of a [`RetainingWallError`], useful for callers
/// that want to react to a class of problem without matching every
/// variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// The caller supplied an out-of-domain or non-finite input.
    Input,
}

impl RetainingWallError {
    /// Stable, kebab/dot-cased identifier for this error, suitable for
    /// logging or matching in tests without depending on the human
    /// message text.
    pub fn code(&self) -> &'static str {
        match self {
            RetainingWallError::FrictionAngleOutOfRange { .. } => {
                "retainingwall.friction_angle_out_of_range"
            }
            RetainingWallError::NonPositiveUnitWeight { .. } => {
                "retainingwall.non_positive_unit_weight"
            }
            RetainingWallError::NegativeDepth { .. } => "retainingwall.negative_depth",
            RetainingWallError::NegativeThrust { .. } => "retainingwall.negative_thrust",
            RetainingWallError::NonFinite { .. } => "retainingwall.non_finite",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    ///
    /// Every current variant is an input-domain failure, but the method
    /// exists so callers can switch on category rather than the exact
    /// variant and stay forward-compatible if new categories are added.
    pub fn category(&self) -> ErrorCategory {
        match self {
            RetainingWallError::FrictionAngleOutOfRange { .. }
            | RetainingWallError::NonPositiveUnitWeight { .. }
            | RetainingWallError::NegativeDepth { .. }
            | RetainingWallError::NegativeThrust { .. }
            | RetainingWallError::NonFinite { .. } => ErrorCategory::Input,
        }
    }
}
