//! Frequency-response data point shared by every filter model.
//!
//! A [`Response`] bundles the linear magnitude (gain) and the phase
//! shift of a filter's transfer function evaluated at one frequency,
//! plus the convenience conversion to decibels.

use serde::{Deserialize, Serialize};

/// The complex frequency-response of a filter at a single frequency,
/// stored as a (magnitude, phase) pair.
///
/// `magnitude` is the **linear** voltage gain `|H(f)|` (dimensionless,
/// `>= 0`); `phase_rad` is the argument `arg H(f)` in radians. For the
/// passive first-order sections modelled here the magnitude never
/// exceeds 1 (they only attenuate) and the phase lies in
/// `[-pi/2, +pi/2]`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Response {
    /// Linear magnitude `|H(f)|` (voltage gain, dimensionless).
    pub magnitude: f64,
    /// Phase shift `arg H(f)` in radians.
    pub phase_rad: f64,
}

impl Response {
    /// Build a response from a linear magnitude and a phase in radians.
    #[must_use]
    pub fn new(magnitude: f64, phase_rad: f64) -> Self {
        Self {
            magnitude,
            phase_rad,
        }
    }

    /// The magnitude expressed in decibels, `20 * log10(|H|)`.
    ///
    /// Returns [`f64::NEG_INFINITY`] at a magnitude of exactly zero (a
    /// perfect null), matching the mathematical limit of the dB scale.
    #[must_use]
    pub fn magnitude_db(&self) -> f64 {
        20.0 * self.magnitude.log10()
    }

    /// The phase shift converted from radians to degrees.
    #[must_use]
    pub fn phase_deg(&self) -> f64 {
        self.phase_rad.to_degrees()
    }
}
