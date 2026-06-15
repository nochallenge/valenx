//! Prandtl number — ratio of momentum to thermal diffusivity.
//!
//! The Prandtl number
//!
//! ```text
//! Pr = cp mu / k = nu / alpha
//! ```
//!
//! relates the rate of viscous (momentum) diffusion to the rate of
//! thermal diffusion. Here `cp` is the specific heat at constant
//! pressure, `mu` the dynamic viscosity, `k` the thermal conductivity,
//! `nu = mu / rho` the kinematic viscosity, and
//! `alpha = k / (rho cp)` the thermal diffusivity. It is a property of
//! the fluid alone (it contains no length or velocity scale) and is
//! dimensionless in any consistent unit system.
//!
//! Order-of-magnitude values: gases sit near `Pr ~ 0.7`, water near
//! `Pr ~ 7`, and viscous oils reach the hundreds or thousands.

use crate::error::{require_positive, DimensionlessError};
use serde::{Deserialize, Serialize};

/// Prandtl number `Pr = cp mu / k`, dimensionless.
///
/// Construct with [`Prandtl::new`] (property form) or
/// [`Prandtl::from_diffusivities`] (`Pr = nu / alpha`). The inner value
/// is always finite and strictly positive.
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct Prandtl(f64);

impl Prandtl {
    /// Build a Prandtl number from `Pr = cp mu / k`.
    ///
    /// All three inputs are physical transport / thermodynamic
    /// properties and must be strictly positive.
    ///
    /// - `specific_heat` (`cp`)
    /// - `dynamic_viscosity` (`mu`)
    /// - `thermal_conductivity` (`k`) — the denominator
    ///
    /// # Errors
    ///
    /// Returns [`DimensionlessError`] if any input is non-finite or not
    /// strictly positive.
    pub fn new(
        specific_heat: f64,
        dynamic_viscosity: f64,
        thermal_conductivity: f64,
    ) -> Result<Self, DimensionlessError> {
        let cp = require_positive("specific_heat", specific_heat)?;
        let mu = require_positive("dynamic_viscosity", dynamic_viscosity)?;
        let k = require_positive("thermal_conductivity", thermal_conductivity)?;
        Ok(Prandtl(cp * mu / k))
    }

    /// Build a Prandtl number from the diffusivity form
    /// `Pr = nu / alpha`, where `nu` is the kinematic viscosity and
    /// `alpha` the thermal diffusivity. Both must be strictly positive.
    ///
    /// # Errors
    ///
    /// Returns [`DimensionlessError`] if either input is non-finite or
    /// not strictly positive.
    pub fn from_diffusivities(
        kinematic_viscosity: f64,
        thermal_diffusivity: f64,
    ) -> Result<Self, DimensionlessError> {
        let nu = require_positive("kinematic_viscosity", kinematic_viscosity)?;
        let alpha = require_positive("thermal_diffusivity", thermal_diffusivity)?;
        Ok(Prandtl(nu / alpha))
    }

    /// The raw dimensionless value.
    pub fn value(&self) -> f64 {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn prandtl_matches_definition() {
        // cp=4180, mu=1e-3, k=0.6  ->  Pr = 4180*1e-3/0.6 = 6.9667
        let pr = Prandtl::new(4180.0, 1.0e-3, 0.6).unwrap();
        assert!((pr.value() - (4180.0 * 1.0e-3 / 0.6)).abs() < 1e-9);
    }

    #[test]
    fn water_is_about_seven() {
        let pr = Prandtl::new(4180.0, 1.0e-3, 0.6).unwrap();
        assert!((pr.value() - 6.9667).abs() < 1e-3);
    }

    #[test]
    fn air_is_below_one() {
        // cp=1005, mu=1.81e-5, k=0.0257  ->  Pr ~ 0.708
        let pr = Prandtl::new(1005.0, 1.81e-5, 0.0257).unwrap();
        assert!((pr.value() - 0.7077).abs() < 1e-2);
        assert!(pr.value() < 1.0);
    }

    #[test]
    fn diffusivity_form_agrees_with_property_form() {
        // alpha = k / (rho cp); nu = mu / rho; Pr = nu/alpha = cp mu / k.
        let rho = 1000.0;
        let cp = 4180.0;
        let mu = 1.0e-3;
        let k = 0.6;
        let nu = mu / rho;
        let alpha = k / (rho * cp);
        let a = Prandtl::new(cp, mu, k).unwrap();
        let b = Prandtl::from_diffusivities(nu, alpha).unwrap();
        assert!((a.value() - b.value()).abs() < 1e-9 * a.value());
    }

    #[test]
    fn unit_inputs_give_unit_value() {
        let pr = Prandtl::new(1.0, 1.0, 1.0).unwrap();
        assert!((pr.value() - 1.0).abs() < EPS);
    }

    #[test]
    fn rejects_non_positive_conductivity() {
        let err = Prandtl::new(1.0, 1.0, 0.0).unwrap_err();
        assert_eq!(err.parameter(), "thermal_conductivity");
        assert_eq!(err.code(), "dimensionless.out-of-domain");
    }

    #[test]
    fn rejects_non_finite_input() {
        let err = Prandtl::new(f64::INFINITY, 1.0, 1.0).unwrap_err();
        assert_eq!(err.code(), "dimensionless.not-finite");
        assert_eq!(err.parameter(), "specific_heat");
    }
}
