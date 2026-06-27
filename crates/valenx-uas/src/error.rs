//! Fail-loud error taxonomy for `valenx-uas`.
//!
//! Every constructor and analysis returns [`Result`] with one of these
//! variants rather than producing a silent `NaN`, panic, or a physically
//! impossible number. Errors from the composed crates ([`valenx_drone`],
//! [`valenx_rotor`], [`valenx_fixedwing`]) are wrapped so a caller sees a
//! single error type.

use thiserror::Error;

/// An out-of-domain or otherwise rejected `valenx-uas` input.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum UasError {
    /// A quantity that must be finite and strictly positive was not.
    #[error("{quantity} must be finite and positive, got {value}")]
    NonPositive {
        /// The offending quantity's name.
        quantity: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A quantity that must be finite (any sign) was not.
    #[error("{quantity} must be finite, got {value}")]
    NotFinite {
        /// The offending quantity's name.
        quantity: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A fraction (e.g. usable-battery fraction, efficiency) fell outside the
    /// physical half-open range `(0, 1]`.
    #[error("{quantity} must be in (0, 1], got {value}")]
    OutOfUnitRange {
        /// The offending quantity's name.
        quantity: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A trade study was set up with no design points or no objectives.
    #[error("trade study is empty: {what}")]
    EmptyStudy {
        /// What was empty ("design points" / "objectives").
        what: &'static str,
    },

    /// An error bubbled up from the composed multirotor (momentum-theory)
    /// crate.
    #[error("multirotor model: {0}")]
    Drone(#[from] valenx_drone::DroneError),

    /// An error bubbled up from the composed BEMT rotor crate.
    #[error("rotor (BEMT) model: {0}")]
    Rotor(#[from] valenx_rotor::RotorError),

    /// An error bubbled up from the composed fixed-wing point-performance
    /// crate.
    #[error("fixed-wing model: {0}")]
    FixedWing(#[from] valenx_fixedwing::FixedWingError),
}

/// Return `value` when finite and strictly positive, else [`UasError::NonPositive`].
pub(crate) fn require_positive(quantity: &'static str, value: f64) -> Result<f64, UasError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(UasError::NonPositive { quantity, value })
    }
}

/// Return `value` when finite, else [`UasError::NotFinite`].
pub(crate) fn require_finite(quantity: &'static str, value: f64) -> Result<f64, UasError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(UasError::NotFinite { quantity, value })
    }
}

/// Return `value` when finite and in `(0, 1]`, else [`UasError::OutOfUnitRange`].
pub(crate) fn require_unit_fraction(quantity: &'static str, value: f64) -> Result<f64, UasError> {
    if value.is_finite() && value > 0.0 && value <= 1.0 {
        Ok(value)
    } else {
        Err(UasError::OutOfUnitRange { quantity, value })
    }
}
