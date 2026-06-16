//! Error taxonomy for the induction-motor kinematics crate.

use thiserror::Error;

/// Errors raised when constructing or evaluating an induction-motor
/// model with out-of-domain parameters.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum InductionMotorError {
    /// The supplied pole count was zero, odd, or otherwise invalid.
    ///
    /// A physical poly-phase winding always produces an **even**,
    /// strictly positive number of magnetic poles (2, 4, 6, 8, ...),
    /// so an odd or zero value is rejected.
    #[error("invalid pole count {poles}: must be a positive even integer (2, 4, 6, ...)")]
    InvalidPoles {
        /// The rejected pole count.
        poles: u32,
    },

    /// The supplied supply frequency was not strictly positive and finite.
    #[error("invalid supply frequency {hz} Hz: must be finite and > 0")]
    InvalidFrequency {
        /// The rejected line frequency in hertz.
        hz: f64,
    },

    /// The supplied rotor mechanical speed was not finite.
    #[error("invalid rotor speed {rpm} rev/min: must be finite")]
    InvalidRotorSpeed {
        /// The rejected rotor speed in revolutions per minute.
        rpm: f64,
    },

    /// The implied slip fell outside the physical motoring range `[0, 1]`.
    ///
    /// Slip is `0` exactly at synchronous speed and `1` at standstill;
    /// a value `< 0` means the rotor is running **faster** than the
    /// rotating field (generating / hyper-synchronous), and a value
    /// `> 1` means the rotor is being driven backwards against the
    /// field (plugging / counter-rotation). Both lie outside ordinary
    /// motoring, so the validated constructor rejects them.
    #[error("slip {slip} is outside the motoring range [0, 1] for Ns = {sync_rpm} rev/min, Nr = {rotor_rpm} rev/min")]
    SlipOutOfRange {
        /// The computed fractional slip.
        slip: f64,
        /// The synchronous speed in revolutions per minute.
        sync_rpm: f64,
        /// The rotor mechanical speed in revolutions per minute.
        rotor_rpm: f64,
    },

    /// The supplied slip value was not finite or not within `[0, 1]`.
    #[error("invalid slip {slip}: must be finite and within [0, 1]")]
    InvalidSlip {
        /// The rejected fractional slip.
        slip: f64,
    },

    /// The supplied air-gap power was not finite or was negative.
    ///
    /// Air-gap power (the power crossing into the rotor) is non-negative
    /// for ordinary motoring; `NaN`, infinity, and negative values are
    /// rejected by the power-split helpers.
    #[error("invalid air-gap power {watts} W: must be finite and >= 0")]
    InvalidPower {
        /// The rejected air-gap power in watts.
        watts: f64,
    },
}

/// Coarse classification of an [`InductionMotorError`], useful for
/// callers that want to branch on the kind of failure without
/// matching every variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// A user-supplied input value was structurally invalid
    /// (non-finite, non-positive, or an invalid pole count).
    Input,
    /// The inputs were individually valid but together implied a
    /// physically out-of-domain operating point (slip outside `[0, 1]`).
    Domain,
}

impl InductionMotorError {
    /// Stable kebab-cased identifier for this error, suitable for
    /// logging or matching across versions.
    pub fn code(&self) -> &'static str {
        match self {
            InductionMotorError::InvalidPoles { .. } => "inductionmotor.invalid_poles",
            InductionMotorError::InvalidFrequency { .. } => "inductionmotor.invalid_frequency",
            InductionMotorError::InvalidRotorSpeed { .. } => "inductionmotor.invalid_rotor_speed",
            InductionMotorError::SlipOutOfRange { .. } => "inductionmotor.slip_out_of_range",
            InductionMotorError::InvalidSlip { .. } => "inductionmotor.invalid_slip",
            InductionMotorError::InvalidPower { .. } => "inductionmotor.invalid_power",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            InductionMotorError::InvalidPoles { .. }
            | InductionMotorError::InvalidFrequency { .. }
            | InductionMotorError::InvalidRotorSpeed { .. }
            | InductionMotorError::InvalidSlip { .. }
            | InductionMotorError::InvalidPower { .. } => ErrorCategory::Input,
            InductionMotorError::SlipOutOfRange { .. } => ErrorCategory::Domain,
        }
    }
}
