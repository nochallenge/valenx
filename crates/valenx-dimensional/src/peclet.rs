//! Peclet number — ratio of advective to diffusive transport.
//!
//! The thermal Peclet number measures how strongly bulk fluid motion
//! (advection) carries heat compared with molecular diffusion. It has
//! two equivalent definitions:
//!
//! ```text
//! Pe = rho cp v L / k        (direct, from properties)
//! Pe = Re * Pr               (product of Reynolds and Prandtl)
//! ```
//!
//! The two agree because
//! `Re * Pr = (rho v L / mu) * (cp mu / k) = rho cp v L / k`; the
//! viscosity `mu` cancels. A large Peclet number means advection
//! dominates; a small one means diffusion dominates. It is dimensionless
//! in any consistent unit system.
//!
//! Equivalently `Pe = v L / alpha` with thermal diffusivity
//! `alpha = k / (rho cp)`.

use crate::error::{require_non_negative, require_positive, DimensionlessError};
use crate::prandtl::Prandtl;
use crate::reynolds::Reynolds;
use serde::{Deserialize, Serialize};

/// Peclet number, dimensionless.
///
/// Construct directly with [`Peclet::new`], from a thermal diffusivity
/// with [`Peclet::from_diffusivity`], or as the product of a Reynolds
/// and Prandtl number with [`Peclet::from_reynolds_prandtl`]. The inner
/// value is always finite and non-negative.
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct Peclet(f64);

impl Peclet {
    /// Build a Peclet number directly from
    /// `Pe = rho cp v L / k`.
    ///
    /// - `density` (`rho`), `specific_heat` (`cp`), `length` (`L`), and
    ///   `thermal_conductivity` (`k`) must be strictly positive.
    /// - `velocity` (`v`) must be finite and non-negative.
    ///
    /// # Errors
    ///
    /// Returns [`DimensionlessError`] if any input is non-finite or
    /// violates the domain above.
    pub fn new(
        density: f64,
        specific_heat: f64,
        velocity: f64,
        length: f64,
        thermal_conductivity: f64,
    ) -> Result<Self, DimensionlessError> {
        let rho = require_positive("density", density)?;
        let cp = require_positive("specific_heat", specific_heat)?;
        let v = require_non_negative("velocity", velocity)?;
        let length = require_positive("length", length)?;
        let k = require_positive("thermal_conductivity", thermal_conductivity)?;
        Ok(Peclet(rho * cp * v * length / k))
    }

    /// Build a Peclet number from the diffusivity form
    /// `Pe = v L / alpha`, where `alpha` is the thermal diffusivity.
    ///
    /// - `length` (`L`) and `thermal_diffusivity` (`alpha`) must be
    ///   strictly positive.
    /// - `velocity` (`v`) must be finite and non-negative.
    ///
    /// # Errors
    ///
    /// Returns [`DimensionlessError`] if any input is non-finite or
    /// violates the domain above.
    pub fn from_diffusivity(
        velocity: f64,
        length: f64,
        thermal_diffusivity: f64,
    ) -> Result<Self, DimensionlessError> {
        let v = require_non_negative("velocity", velocity)?;
        let length = require_positive("length", length)?;
        let alpha = require_positive("thermal_diffusivity", thermal_diffusivity)?;
        Ok(Peclet(v * length / alpha))
    }

    /// Build a Peclet number as the product `Pe = Re * Pr` of an
    /// already-constructed [`Reynolds`] and [`Prandtl`] number.
    ///
    /// Because both inputs were validated at their own construction, this
    /// is infallible — the product of two finite, non-negative values is
    /// finite and non-negative.
    pub fn from_reynolds_prandtl(reynolds: Reynolds, prandtl: Prandtl) -> Self {
        Peclet(reynolds.value() * prandtl.value())
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
    fn peclet_matches_definition() {
        // rho=1000, cp=4180, v=1, L=0.1, k=0.6
        // Pe = 1000*4180*1*0.1/0.6 = 696666.67
        let pe = Peclet::new(1000.0, 4180.0, 1.0, 0.1, 0.6).unwrap();
        let expected = 1000.0 * 4180.0 * 1.0 * 0.1 / 0.6;
        assert!((pe.value() - expected).abs() < 1e-3);
    }

    #[test]
    fn product_form_equals_direct_form() {
        // The central validation: Pe = Re * Pr must equal rho cp v L / k.
        let rho = 998.0;
        let cp = 4182.0;
        let mu = 1.0e-3;
        let k = 0.598;
        let v = 1.2;
        let l = 0.05;

        let direct = Peclet::new(rho, cp, v, l, k).unwrap();

        let re = Reynolds::new(rho, v, l, mu).unwrap();
        let pr = Prandtl::new(cp, mu, k).unwrap();
        let product = Peclet::from_reynolds_prandtl(re, pr);

        // mu cancels analytically, so the two must agree to round-off.
        assert!((direct.value() - product.value()).abs() < 1e-6 * direct.value());
    }

    #[test]
    fn product_form_re_times_pr() {
        // Direct numeric check of the product wrapper.
        let re = Reynolds::new(1.0, 100.0, 1.0, 1.0).unwrap(); // Re = 100
        let pr = Prandtl::new(7.0, 1.0, 1.0).unwrap(); // Pr = 7
        let pe = Peclet::from_reynolds_prandtl(re, pr);
        assert!((pe.value() - 700.0).abs() < 1e-6);
    }

    #[test]
    fn diffusivity_form_agrees_with_direct() {
        // alpha = k / (rho cp); Pe = v L / alpha = rho cp v L / k.
        let rho = 1000.0;
        let cp = 4180.0;
        let k = 0.6;
        let v = 1.0;
        let l = 0.1;
        let alpha = k / (rho * cp);
        let a = Peclet::new(rho, cp, v, l, k).unwrap();
        let b = Peclet::from_diffusivity(v, l, alpha).unwrap();
        assert!((a.value() - b.value()).abs() < 1e-6 * a.value());
    }

    #[test]
    fn unit_inputs_give_unit_value() {
        let pe = Peclet::new(1.0, 1.0, 1.0, 1.0, 1.0).unwrap();
        assert!((pe.value() - 1.0).abs() < EPS);
    }

    #[test]
    fn zero_velocity_gives_zero() {
        let pe = Peclet::new(1000.0, 4180.0, 0.0, 0.1, 0.6).unwrap();
        assert!(pe.value().abs() < EPS);
    }

    #[test]
    fn rejects_non_positive_conductivity() {
        let err = Peclet::new(1.0, 1.0, 1.0, 1.0, 0.0).unwrap_err();
        assert_eq!(err.parameter(), "thermal_conductivity");
        assert_eq!(err.code(), "dimensionless.out-of-domain");
    }

    #[test]
    fn rejects_negative_velocity() {
        let err = Peclet::new(1.0, 1.0, -1.0, 1.0, 1.0).unwrap_err();
        assert_eq!(err.parameter(), "velocity");
    }

    #[test]
    fn rejects_non_finite_density() {
        let err = Peclet::new(f64::NAN, 1.0, 1.0, 1.0, 1.0).unwrap_err();
        assert_eq!(err.code(), "dimensionless.not-finite");
        assert_eq!(err.parameter(), "density");
    }
}
