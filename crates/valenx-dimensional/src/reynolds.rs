//! Reynolds number — ratio of inertial to viscous forces.
//!
//! The Reynolds number
//!
//! ```text
//! Re = rho v L / mu = v L / nu
//! ```
//!
//! compares inertial forces (`rho v^2`) to viscous shear (`mu v / L`).
//! Here `rho` is density, `v` a characteristic velocity, `L` a
//! characteristic length, `mu` the dynamic viscosity, and
//! `nu = mu / rho` the kinematic viscosity. It is dimensionless in any
//! consistent unit system.
//!
//! For flow in a circular pipe the textbook regime boundaries on the
//! Reynolds number (based on the pipe diameter) are approximately:
//! laminar below about 2300, transitional from about 2300 to 4000, and
//! turbulent above about 4000. See [`PipeRegime`].

use crate::error::{require_non_negative, require_positive, DimensionlessError};
use serde::{Deserialize, Serialize};

/// Lower transition threshold for fully developed pipe flow: below this
/// Reynolds number the flow is treated as laminar. The classical value
/// is approximately 2300.
pub const PIPE_LAMINAR_LIMIT: f64 = 2300.0;

/// Upper transition threshold for fully developed pipe flow: above this
/// Reynolds number the flow is treated as fully turbulent. The classical
/// value is approximately 4000.
pub const PIPE_TURBULENT_LIMIT: f64 = 4000.0;

/// Reynolds number `Re = rho v L / mu`, dimensionless.
///
/// Construct with [`Reynolds::new`] (density / viscosity form) or
/// [`Reynolds::from_kinematic`] (kinematic-viscosity form). The inner
/// value is always finite and non-negative.
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct Reynolds(f64);

impl Reynolds {
    /// Build a Reynolds number from `Re = rho v L / mu`.
    ///
    /// - `density` (`rho`) and `length` (`L`) must be strictly positive.
    /// - `dynamic_viscosity` (`mu`) must be strictly positive (it is the
    ///   denominator).
    /// - `velocity` (`v`) must be finite and non-negative.
    ///
    /// # Errors
    ///
    /// Returns [`DimensionlessError`] if any input is non-finite or
    /// violates the domain above.
    pub fn new(
        density: f64,
        velocity: f64,
        length: f64,
        dynamic_viscosity: f64,
    ) -> Result<Self, DimensionlessError> {
        let density = require_positive("density", density)?;
        let velocity = require_non_negative("velocity", velocity)?;
        let length = require_positive("length", length)?;
        let mu = require_positive("dynamic_viscosity", dynamic_viscosity)?;
        Ok(Reynolds(density * velocity * length / mu))
    }

    /// Build a Reynolds number from the kinematic-viscosity form
    /// `Re = v L / nu`, where `nu = mu / rho`.
    ///
    /// - `length` (`L`) and `kinematic_viscosity` (`nu`) must be strictly
    ///   positive.
    /// - `velocity` (`v`) must be finite and non-negative.
    ///
    /// # Errors
    ///
    /// Returns [`DimensionlessError`] if any input is non-finite or
    /// violates the domain above.
    pub fn from_kinematic(
        velocity: f64,
        length: f64,
        kinematic_viscosity: f64,
    ) -> Result<Self, DimensionlessError> {
        let velocity = require_non_negative("velocity", velocity)?;
        let length = require_positive("length", length)?;
        let nu = require_positive("kinematic_viscosity", kinematic_viscosity)?;
        Ok(Reynolds(velocity * length / nu))
    }

    /// The raw dimensionless value.
    pub fn value(&self) -> f64 {
        self.0
    }

    /// Classify fully developed circular-pipe flow into a [`PipeRegime`]
    /// using the classical thresholds at [`PIPE_LAMINAR_LIMIT`] and
    /// [`PIPE_TURBULENT_LIMIT`]. These boundaries are approximate rules
    /// of thumb, not sharp physical transitions.
    pub fn pipe_regime(&self) -> PipeRegime {
        if self.0 < PIPE_LAMINAR_LIMIT {
            PipeRegime::Laminar
        } else if self.0 < PIPE_TURBULENT_LIMIT {
            PipeRegime::Transitional
        } else {
            PipeRegime::Turbulent
        }
    }
}

/// Coarse flow regime for fully developed flow in a circular pipe, as a
/// function of the diameter-based Reynolds number.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum PipeRegime {
    /// Smooth, layered flow (`Re < ~2300`).
    Laminar,
    /// Intermittent / transitional flow (`~2300 <= Re < ~4000`).
    Transitional,
    /// Chaotic, eddying flow (`Re >= ~4000`).
    Turbulent,
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn reynolds_matches_definition() {
        // rho=1000, v=2, L=0.05, mu=1e-3  ->  Re = 1000*2*0.05/1e-3 = 1e5
        let re = Reynolds::new(1000.0, 2.0, 0.05, 1.0e-3).unwrap();
        assert!((re.value() - 1.0e5).abs() < 1e-3);
    }

    #[test]
    fn reynolds_unit_inputs_give_unit_value() {
        let re = Reynolds::new(1.0, 1.0, 1.0, 1.0).unwrap();
        assert!((re.value() - 1.0).abs() < EPS);
    }

    #[test]
    fn kinematic_form_agrees_with_density_form() {
        // nu = mu / rho, so both constructors must agree.
        let rho = 998.0;
        let mu = 1.002e-3;
        let nu = mu / rho;
        let v = 1.5;
        let l = 0.1;
        let a = Reynolds::new(rho, v, l, mu).unwrap();
        let b = Reynolds::from_kinematic(v, l, nu).unwrap();
        assert!((a.value() - b.value()).abs() < 1e-6 * a.value());
    }

    #[test]
    fn zero_velocity_is_zero_reynolds() {
        let re = Reynolds::new(1.0, 0.0, 1.0, 1.0).unwrap();
        assert!(re.value().abs() < EPS);
        assert_eq!(re.pipe_regime(), PipeRegime::Laminar);
    }

    #[test]
    fn pipe_regime_thresholds() {
        // Just below the laminar limit -> laminar.
        let lam = Reynolds::new(1.0, 2299.0, 1.0, 1.0).unwrap();
        assert_eq!(lam.pipe_regime(), PipeRegime::Laminar);

        // Exactly at the laminar limit -> transitional (>= boundary).
        let at = Reynolds::new(1.0, PIPE_LAMINAR_LIMIT, 1.0, 1.0).unwrap();
        assert_eq!(at.pipe_regime(), PipeRegime::Transitional);

        // Between the two limits -> transitional.
        let tr = Reynolds::new(1.0, 3000.0, 1.0, 1.0).unwrap();
        assert_eq!(tr.pipe_regime(), PipeRegime::Transitional);

        // Above the turbulent limit -> turbulent.
        let turb = Reynolds::new(1.0, 5000.0, 1.0, 1.0).unwrap();
        assert_eq!(turb.pipe_regime(), PipeRegime::Turbulent);
    }

    #[test]
    fn rejects_non_positive_viscosity() {
        let err = Reynolds::new(1.0, 1.0, 1.0, 0.0).unwrap_err();
        assert_eq!(err.code(), "dimensionless.out-of-domain");
        assert_eq!(err.parameter(), "dynamic_viscosity");
    }

    #[test]
    fn rejects_negative_velocity() {
        let err = Reynolds::new(1.0, -1.0, 1.0, 1.0).unwrap_err();
        assert_eq!(err.parameter(), "velocity");
    }

    #[test]
    fn rejects_non_finite_density() {
        let err = Reynolds::new(f64::NAN, 1.0, 1.0, 1.0).unwrap_err();
        assert_eq!(err.code(), "dimensionless.not-finite");
        assert_eq!(err.parameter(), "density");
    }

    #[test]
    fn serde_round_trip() {
        let re = Reynolds::new(1000.0, 2.0, 0.05, 1.0e-3).unwrap();
        let json = serde_json::to_string(&re).unwrap();
        let back: Reynolds = serde_json::from_str(&json).unwrap();
        assert!((re.value() - back.value()).abs() < EPS);
    }
}
