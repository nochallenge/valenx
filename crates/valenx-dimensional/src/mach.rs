//! Mach number — flow speed relative to the local speed of sound.
//!
//! The Mach number
//!
//! ```text
//! Ma = v / c
//! ```
//!
//! is the ratio of a flow (or body) speed `v` to the local speed of
//! sound `c` in the medium. It governs the importance of compressibility
//! effects. It is dimensionless in any consistent unit system.
//!
//! The conventional speed regimes keyed on the Mach number are: subsonic
//! below 1, sonic exactly at 1, and supersonic above 1. (The narrow band
//! around `Ma = 1`, roughly 0.8 to 1.2, is often called transonic; this
//! crate uses the simple three-way split — see [`SpeedRegime`].)

use crate::error::{require_non_negative, require_positive, DimensionlessError};
use serde::{Deserialize, Serialize};

/// Mach number `Ma = v / c`, dimensionless.
///
/// Construct with [`Mach::new`]. The inner value is always finite and
/// non-negative.
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct Mach(f64);

impl Mach {
    /// Build a Mach number from `Ma = v / c`.
    ///
    /// - `speed` (`v`) must be finite and non-negative (it may be zero).
    /// - `speed_of_sound` (`c`) must be strictly positive (it is the
    ///   denominator).
    ///
    /// # Errors
    ///
    /// Returns [`DimensionlessError`] if either input is non-finite or
    /// violates the domain above.
    pub fn new(speed: f64, speed_of_sound: f64) -> Result<Self, DimensionlessError> {
        let v = require_non_negative("speed", speed)?;
        let c = require_positive("speed_of_sound", speed_of_sound)?;
        Ok(Mach(v / c))
    }

    /// The raw dimensionless value.
    pub fn value(&self) -> f64 {
        self.0
    }

    /// Classify the flow speed into a [`SpeedRegime`] using the sonic
    /// boundary `Ma = 1`. The equality case (`Ma == 1`) is reported as
    /// [`SpeedRegime::Sonic`].
    pub fn speed_regime(&self) -> SpeedRegime {
        if self.0 < 1.0 {
            SpeedRegime::Subsonic
        } else if self.0 > 1.0 {
            SpeedRegime::Supersonic
        } else {
            SpeedRegime::Sonic
        }
    }
}

/// Coarse flow-speed regime as a function of the Mach number.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum SpeedRegime {
    /// `Ma < 1`: slower than the local speed of sound.
    Subsonic,
    /// `Ma == 1`: exactly at the local speed of sound.
    Sonic,
    /// `Ma > 1`: faster than the local speed of sound.
    Supersonic,
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn mach_matches_definition() {
        // v=340, c=340  ->  Ma = 1.
        let m = Mach::new(340.0, 340.0).unwrap();
        assert!((m.value() - 1.0).abs() < EPS);
    }

    #[test]
    fn airliner_cruise_is_subsonic() {
        // v=250 m/s at c=295 m/s (cruise altitude) -> Ma ~ 0.847.
        let m = Mach::new(250.0, 295.0).unwrap();
        assert!((m.value() - (250.0 / 295.0)).abs() < 1e-9);
        assert_eq!(m.speed_regime(), SpeedRegime::Subsonic);
    }

    #[test]
    fn fast_jet_is_supersonic() {
        // v=680, c=340 -> Ma = 2.
        let m = Mach::new(680.0, 340.0).unwrap();
        assert!((m.value() - 2.0).abs() < EPS);
        assert_eq!(m.speed_regime(), SpeedRegime::Supersonic);
    }

    #[test]
    fn zero_speed_is_subsonic() {
        let m = Mach::new(0.0, 340.0).unwrap();
        assert!(m.value().abs() < EPS);
        assert_eq!(m.speed_regime(), SpeedRegime::Subsonic);
    }

    #[test]
    fn exactly_sonic() {
        let m = Mach::new(300.0, 300.0).unwrap();
        assert_eq!(m.speed_regime(), SpeedRegime::Sonic);
    }

    #[test]
    fn rejects_non_positive_speed_of_sound() {
        let err = Mach::new(100.0, 0.0).unwrap_err();
        assert_eq!(err.parameter(), "speed_of_sound");
        assert_eq!(err.code(), "dimensionless.out-of-domain");
    }

    #[test]
    fn rejects_negative_speed() {
        let err = Mach::new(-1.0, 340.0).unwrap_err();
        assert_eq!(err.parameter(), "speed");
    }

    #[test]
    fn rejects_non_finite_speed() {
        let err = Mach::new(f64::INFINITY, 340.0).unwrap_err();
        assert_eq!(err.code(), "dimensionless.not-finite");
        assert_eq!(err.parameter(), "speed");
    }
}
