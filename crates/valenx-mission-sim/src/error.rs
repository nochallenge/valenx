//! Fail-loud error taxonomy for `valenx-mission-sim`.
//!
//! Every constructor and analysis returns [`Result`] with one of these variants
//! rather than producing a silent `NaN`, a panic, or a physically impossible
//! number. Errors from the reused [`valenx_uas`] detection geometry are wrapped
//! so a caller sees a single error type.

use thiserror::Error;

/// An out-of-domain or otherwise rejected `valenx-mission-sim` input.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum MissionError {
    /// A quantity that must be finite and strictly positive was not.
    #[error("{quantity} must be finite and positive, got {value}")]
    NonPositive {
        /// The offending quantity's name.
        quantity: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A quantity that must be finite and non-negative (`>= 0`) was not.
    #[error("{quantity} must be finite and non-negative, got {value}")]
    Negative {
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

    /// A probability (e.g. a probability-of-kill `Pk`) fell outside `[0, 1]`.
    #[error("{quantity} must be a probability in [0, 1], got {value}")]
    OutOfProbabilityRange {
        /// The offending quantity's name.
        quantity: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A scheduled event time lay at or before the current simulation clock,
    /// which would break the monotonic-time invariant of the scheduler.
    #[error("event time {event_time} is not in the future of clock {now}")]
    NonMonotonicEvent {
        /// The clock value at scheduling time (s).
        now: f64,
        /// The rejected event time (s).
        event_time: f64,
    },

    /// An entity id referenced a slot that does not exist in the scenario.
    #[error("unknown entity id {0}")]
    UnknownEntity(usize),

    /// A waypoint route was supplied with no waypoints.
    #[error("waypoint mover has an empty route")]
    EmptyRoute,

    /// An error bubbled up from the reused counter-UAS detection geometry.
    #[error("detection geometry: {0}")]
    Uas(#[from] valenx_uas::UasError),
}

/// Return `value` when finite and strictly positive, else
/// [`MissionError::NonPositive`].
pub(crate) fn require_positive(quantity: &'static str, value: f64) -> Result<f64, MissionError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(MissionError::NonPositive { quantity, value })
    }
}

/// Return `value` when finite and non-negative, else [`MissionError::Negative`].
pub(crate) fn require_non_negative(
    quantity: &'static str,
    value: f64,
) -> Result<f64, MissionError> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(MissionError::Negative { quantity, value })
    }
}

/// Return `value` when finite, else [`MissionError::NotFinite`].
pub(crate) fn require_finite(quantity: &'static str, value: f64) -> Result<f64, MissionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(MissionError::NotFinite { quantity, value })
    }
}

/// Return `value` when finite and in `[0, 1]`, else
/// [`MissionError::OutOfProbabilityRange`].
pub(crate) fn require_probability(quantity: &'static str, value: f64) -> Result<f64, MissionError> {
    if value.is_finite() && (0.0..=1.0).contains(&value) {
        Ok(value)
    } else {
        Err(MissionError::OutOfProbabilityRange { quantity, value })
    }
}
