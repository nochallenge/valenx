//! Fluids, gravity, and shared physical constants.
//!
//! A [`Fluid`] is just a validated mass density in kilograms per cubic
//! metre. Standard-gravity and a handful of reference densities live
//! here as named constants so every example and test uses the same
//! textbook numbers.

use crate::error::{require_positive, Result};
use serde::{Deserialize, Serialize};

/// Standard gravitational acceleration at the Earth's surface, in
/// metres per second squared (the CGPM-defined value `g_0`).
pub const STANDARD_GRAVITY: f64 = 9.806_65;

/// Density of fresh water at 4 °C, in kilograms per cubic metre — the
/// density maximum, and the historical definition of the kilogram.
pub const DENSITY_WATER_4C: f64 = 1000.0;

/// Representative density of sea water, in kilograms per cubic metre
/// (about 2.5 % denser than fresh water).
pub const DENSITY_SEAWATER: f64 = 1025.0;

/// Density of mercury at 20 °C, in kilograms per cubic metre — the
/// classic manometer fluid.
pub const DENSITY_MERCURY: f64 = 13_534.0;

/// Density of dry air at 15 °C and one atmosphere, in kilograms per
/// cubic metre.
pub const DENSITY_AIR_15C: f64 = 1.225;

/// One standard atmosphere of absolute pressure, in pascals.
pub const STANDARD_ATMOSPHERE_PA: f64 = 101_325.0;

/// A homogeneous fluid characterised by its (constant) mass density.
///
/// Fluid statics treats the fluid as incompressible and of uniform
/// density, so a single positive number fully describes it for the
/// purposes of this crate.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Fluid {
    /// Mass density, in kilograms per cubic metre. Always strictly
    /// positive (enforced by [`Fluid::new`]).
    density: f64,
}

impl Fluid {
    /// Construct a fluid from a strictly-positive mass density in
    /// kilograms per cubic metre.
    ///
    /// # Errors
    ///
    /// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
    /// if `density` is not finite or not strictly positive.
    pub fn new(density: f64) -> Result<Self> {
        let density = require_positive("density", density)?;
        Ok(Fluid { density })
    }

    /// Fresh water at 4 °C (1000 kg/m³).
    pub fn water() -> Self {
        Fluid {
            density: DENSITY_WATER_4C,
        }
    }

    /// Sea water (1025 kg/m³).
    pub fn seawater() -> Self {
        Fluid {
            density: DENSITY_SEAWATER,
        }
    }

    /// Mercury at 20 °C (13 534 kg/m³).
    pub fn mercury() -> Self {
        Fluid {
            density: DENSITY_MERCURY,
        }
    }

    /// Dry air at 15 °C (1.225 kg/m³).
    pub fn air() -> Self {
        Fluid {
            density: DENSITY_AIR_15C,
        }
    }

    /// The fluid's mass density, in kilograms per cubic metre.
    pub fn density(&self) -> f64 {
        self.density
    }

    /// The fluid's specific weight `gamma = rho * g`, in newtons per
    /// cubic metre, under the supplied gravitational acceleration.
    ///
    /// # Errors
    ///
    /// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
    /// if `gravity` is not finite or not strictly positive.
    pub fn specific_weight(&self, gravity: f64) -> Result<f64> {
        let gravity = require_positive("gravity", gravity)?;
        Ok(self.density * gravity)
    }

    /// Relative density (specific gravity) of this fluid with respect to
    /// a `reference` fluid — the dimensionless ratio of mass densities.
    pub fn relative_density(&self, reference: &Fluid) -> f64 {
        self.density / reference.density
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn presets_have_expected_densities() {
        assert!((Fluid::water().density() - 1000.0).abs() < EPS);
        assert!((Fluid::seawater().density() - 1025.0).abs() < EPS);
        assert!((Fluid::mercury().density() - 13_534.0).abs() < EPS);
        assert!((Fluid::air().density() - 1.225).abs() < EPS);
    }

    #[test]
    fn new_rejects_non_positive_density() {
        assert!(Fluid::new(0.0).is_err());
        assert!(Fluid::new(-5.0).is_err());
        assert!(Fluid::new(f64::NAN).is_err());
        assert!(Fluid::new(998.0).is_ok());
    }

    #[test]
    fn specific_weight_is_rho_g() {
        // gamma for water at standard gravity = 1000 * 9.80665 = 9806.65 N/m^3.
        let gamma = Fluid::water().specific_weight(STANDARD_GRAVITY).unwrap();
        assert!((gamma - 9806.65).abs() < 1e-6, "got {gamma}");
    }

    #[test]
    fn specific_weight_rejects_bad_gravity() {
        assert!(Fluid::water().specific_weight(0.0).is_err());
        assert!(Fluid::water().specific_weight(-9.81).is_err());
    }

    #[test]
    fn relative_density_of_mercury_to_water() {
        // Specific gravity of mercury ~ 13.534.
        let sg = Fluid::mercury().relative_density(&Fluid::water());
        assert!((sg - 13.534).abs() < EPS, "got {sg}");
        // A fluid relative to itself is exactly 1.
        let one = Fluid::water().relative_density(&Fluid::water());
        assert!((one - 1.0).abs() < EPS, "got {one}");
    }
}
