//! Horizontal wind models for the ascent drag calculation.
//!
//! A [`WindModel`] returns the eastward (downrange-positive) wind speed
//! at a given altitude. The ascent sim subtracts it — along with the
//! co-rotating atmosphere velocity — from the vehicle velocity before
//! computing drag, so a head- or tail-wind changes the air-relative
//! speed, the dynamic pressure, and the trajectory.
//!
//! The model is a `Copy` enum (so [`crate::config::AscentConfig`] stays
//! `Copy`): no wind, a uniform wind, or a Gaussian "jet" profile peaking
//! at a chosen altitude — a reasonable stand-in for a tropospheric jet
//! stream. Vertical winds and full wind fields are out of scope here.

use serde::{Deserialize, Serialize};

/// Eastward (downrange-positive) wind as a function of altitude.
#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
pub enum WindModel {
    /// Still air.
    #[default]
    None,
    /// Uniform eastward wind (m/s); negative is a westward head-wind.
    Constant(f64),
    /// Gaussian jet: `peak_speed · exp(−((h − peak_altitude)/width)²)`.
    Jet {
        /// Peak eastward speed (m/s).
        peak_speed: f64,
        /// Altitude of the peak (m).
        peak_altitude: f64,
        /// Gaussian half-width (m).
        width: f64,
    },
}

impl WindModel {
    /// Eastward wind speed (m/s) at a geometric `altitude` (m).
    pub fn speed_at(&self, altitude: f64) -> f64 {
        match *self {
            WindModel::None => 0.0,
            WindModel::Constant(s) => s,
            WindModel::Jet {
                peak_speed,
                peak_altitude,
                width,
            } => {
                if width <= 0.0 {
                    return 0.0;
                }
                let z = (altitude - peak_altitude) / width;
                peak_speed * (-z * z).exp()
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn none_is_always_zero() {
        assert_eq!(WindModel::None.speed_at(0.0), 0.0);
        assert_eq!(WindModel::None.speed_at(12_000.0), 0.0);
    }

    #[test]
    fn constant_is_uniform() {
        let w = WindModel::Constant(-40.0);
        assert_eq!(w.speed_at(0.0), -40.0);
        assert_eq!(w.speed_at(50_000.0), -40.0);
    }

    #[test]
    fn jet_peaks_at_its_altitude() {
        let w = WindModel::Jet {
            peak_speed: 60.0,
            peak_altitude: 12_000.0,
            width: 4_000.0,
        };
        // Maximum at the peak altitude.
        assert!((w.speed_at(12_000.0) - 60.0).abs() < 1e-9);
        // Falls off away from the peak.
        assert!(w.speed_at(12_000.0) > w.speed_at(0.0));
        assert!(w.speed_at(12_000.0) > w.speed_at(30_000.0));
        // Effectively zero far away.
        assert!(w.speed_at(40_000.0).abs() < 1.0);
    }
}
