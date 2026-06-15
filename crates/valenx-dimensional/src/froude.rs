//! Froude number — ratio of inertial to gravitational forces.
//!
//! The Froude number
//!
//! ```text
//! Fr = v / sqrt(g L)
//! ```
//!
//! compares a flow's inertia to gravity, where `v` is a characteristic
//! velocity, `g` the gravitational acceleration, and `L` a
//! characteristic length (for open-channel flow, the hydraulic depth).
//! The denominator `sqrt(g L)` is the speed of a shallow-water gravity
//! wave. It is dimensionless in any consistent unit system.
//!
//! In open-channel hydraulics the Froude number sets the flow regime:
//! subcritical (tranquil) below 1, critical at 1, and supercritical
//! (rapid / shooting) above 1. See [`ChannelRegime`].

use crate::error::{require_non_negative, require_positive, DimensionlessError};
use serde::{Deserialize, Serialize};

/// Froude number `Fr = v / sqrt(g L)`, dimensionless.
///
/// Construct with [`Froude::new`]. The inner value is always finite and
/// non-negative.
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct Froude(f64);

impl Froude {
    /// Build a Froude number from `Fr = v / sqrt(g L)`.
    ///
    /// - `velocity` (`v`) must be finite and non-negative (it may be
    ///   zero).
    /// - `gravity` (`g`) and `length` (`L`) must be strictly positive;
    ///   their product appears under the square root in the denominator.
    ///
    /// # Errors
    ///
    /// Returns [`DimensionlessError`] if any input is non-finite or
    /// violates the domain above.
    pub fn new(velocity: f64, gravity: f64, length: f64) -> Result<Self, DimensionlessError> {
        let v = require_non_negative("velocity", velocity)?;
        let g = require_positive("gravity", gravity)?;
        let length = require_positive("length", length)?;
        Ok(Froude(v / (g * length).sqrt()))
    }

    /// The raw dimensionless value.
    pub fn value(&self) -> f64 {
        self.0
    }

    /// Classify open-channel flow into a [`ChannelRegime`] using the
    /// critical boundary `Fr = 1`. The equality case (`Fr == 1`) is
    /// reported as [`ChannelRegime::Critical`].
    pub fn channel_regime(&self) -> ChannelRegime {
        if self.0 < 1.0 {
            ChannelRegime::Subcritical
        } else if self.0 > 1.0 {
            ChannelRegime::Supercritical
        } else {
            ChannelRegime::Critical
        }
    }
}

/// Coarse open-channel flow regime as a function of the Froude number.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum ChannelRegime {
    /// `Fr < 1`: tranquil / subcritical flow; surface waves can travel
    /// upstream.
    Subcritical,
    /// `Fr == 1`: critical flow; the flow speed equals the gravity-wave
    /// speed.
    Critical,
    /// `Fr > 1`: rapid / shooting / supercritical flow; disturbances
    /// cannot propagate upstream.
    Supercritical,
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn froude_matches_definition() {
        // v=2, g=9.81, L=1  ->  Fr = 2/sqrt(9.81) = 0.6386
        let fr = Froude::new(2.0, 9.81, 1.0).unwrap();
        assert!((fr.value() - (2.0 / 9.81_f64.sqrt())).abs() < 1e-9);
    }

    #[test]
    fn subcritical_below_one() {
        // The definition test above: 0.6386 < 1.
        let fr = Froude::new(2.0, 9.81, 1.0).unwrap();
        assert!(fr.value() < 1.0);
        assert_eq!(fr.channel_regime(), ChannelRegime::Subcritical);
    }

    #[test]
    fn supercritical_above_one() {
        // v=5, g=9.81, L=0.5  ->  Fr = 5/sqrt(4.905) = 2.258.
        let fr = Froude::new(5.0, 9.81, 0.5).unwrap();
        assert!(fr.value() > 1.0);
        assert_eq!(fr.channel_regime(), ChannelRegime::Supercritical);
    }

    #[test]
    fn critical_at_one() {
        // Choose v = sqrt(g L) exactly so Fr = 1.
        let g = 9.81_f64;
        let l = 2.0_f64;
        let v = (g * l).sqrt();
        let fr = Froude::new(v, g, l).unwrap();
        assert!((fr.value() - 1.0).abs() < EPS);
        assert_eq!(fr.channel_regime(), ChannelRegime::Critical);
    }

    #[test]
    fn unit_inputs_give_unit_value() {
        let fr = Froude::new(1.0, 1.0, 1.0).unwrap();
        assert!((fr.value() - 1.0).abs() < EPS);
    }

    #[test]
    fn zero_velocity_is_subcritical() {
        let fr = Froude::new(0.0, 9.81, 1.0).unwrap();
        assert!(fr.value().abs() < EPS);
        assert_eq!(fr.channel_regime(), ChannelRegime::Subcritical);
    }

    #[test]
    fn rejects_non_positive_gravity() {
        let err = Froude::new(1.0, 0.0, 1.0).unwrap_err();
        assert_eq!(err.parameter(), "gravity");
        assert_eq!(err.code(), "dimensionless.out-of-domain");
    }

    #[test]
    fn rejects_non_positive_length() {
        let err = Froude::new(1.0, 9.81, -1.0).unwrap_err();
        assert_eq!(err.parameter(), "length");
    }

    #[test]
    fn rejects_non_finite_velocity() {
        let err = Froude::new(f64::NAN, 9.81, 1.0).unwrap_err();
        assert_eq!(err.code(), "dimensionless.not-finite");
        assert_eq!(err.parameter(), "velocity");
    }
}
