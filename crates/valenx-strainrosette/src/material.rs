//! Linear-elastic material and the 2D (plane-stress) Hooke's law that
//! maps a Cartesian strain state to stress.
//!
//! # Constitutive model
//!
//! For an isotropic linear-elastic material under *plane stress* the
//! in-plane stresses follow from the strains by
//!
//! ```text
//! sigma_x  = E/(1 - nu^2) (eps_x + nu eps_y)
//! sigma_y  = E/(1 - nu^2) (eps_y + nu eps_x)
//! tau_xy   = G gamma_xy,   G = E / (2 (1 + nu))
//! ```
//!
//! where `E` is Young's modulus and `nu` is Poisson's ratio. In matrix
//! form the constitutive (stiffness) matrix is
//!
//! ```text
//!            E      [ 1   nu  0          ]
//! C = ------------- [ nu  1   0          ]
//!     (1 - nu^2)    [ 0   0   (1 - nu)/2 ]
//! ```
//!
//! acting on the engineering-strain vector
//! `[eps_x, eps_y, gamma_xy]^T`. [`ElasticMaterial::plane_stress`]
//! builds exactly this `C` (as an `nalgebra` matrix) and multiplies.

use nalgebra::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};

use crate::error::{ensure_finite, RosetteError};
use crate::rosette::CartesianStrain;

/// An isotropic linear-elastic material defined by Young's modulus and
/// Poisson's ratio.
///
/// Units for `youngs_modulus` are the caller's choice (for example MPa
/// or psi); the resulting stresses come out in the same unit, since the
/// strains are dimensionless. `poisson_ratio` is dimensionless.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ElasticMaterial {
    youngs_modulus: f64,
    poisson_ratio: f64,
}

/// In-plane stress state produced by the constitutive law.
///
/// Same unit as the material's Young's modulus.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PlaneStress {
    /// Direct stress along x.
    pub sigma_x: f64,
    /// Direct stress along y.
    pub sigma_y: f64,
    /// Shear stress in the xy-plane.
    pub tau_xy: f64,
}

impl ElasticMaterial {
    /// Validated constructor.
    ///
    /// # Errors
    ///
    /// Returns [`RosetteError::NonFinite`] if either argument is `NaN`
    /// or infinite, and [`RosetteError::InvalidMaterial`] if
    /// `youngs_modulus` is not strictly positive or `poisson_ratio` is
    /// outside the thermodynamically admissible isotropic range
    /// `(-1, 0.5)`.
    pub fn new(youngs_modulus: f64, poisson_ratio: f64) -> Result<Self, RosetteError> {
        ensure_finite("youngs_modulus", youngs_modulus)?;
        ensure_finite("poisson_ratio", poisson_ratio)?;

        if youngs_modulus <= 0.0 {
            return Err(RosetteError::InvalidMaterial {
                name: "youngs_modulus",
                reason: "must be strictly positive",
                value: youngs_modulus,
            });
        }
        // For isotropic elasticity the bulk and shear moduli are positive
        // only when -1 < nu < 0.5.
        if poisson_ratio <= -1.0 || poisson_ratio >= 0.5 {
            return Err(RosetteError::InvalidMaterial {
                name: "poisson_ratio",
                reason: "must lie in the open interval (-1, 0.5)",
                value: poisson_ratio,
            });
        }

        Ok(Self {
            youngs_modulus,
            poisson_ratio,
        })
    }

    /// Young's modulus `E`.
    pub fn youngs_modulus(&self) -> f64 {
        self.youngs_modulus
    }

    /// Poisson's ratio `nu`.
    pub fn poisson_ratio(&self) -> f64 {
        self.poisson_ratio
    }

    /// Shear modulus `G = E / (2 (1 + nu))`.
    pub fn shear_modulus(&self) -> f64 {
        self.youngs_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }

    /// The 3x3 plane-stress constitutive (stiffness) matrix `C` such
    /// that `[sigma_x, sigma_y, tau_xy]^T = C [eps_x, eps_y, gamma_xy]^T`.
    pub fn plane_stress_matrix(&self) -> Matrix3<f64> {
        let e = self.youngs_modulus;
        let nu = self.poisson_ratio;
        let factor = e / (1.0 - nu * nu);
        Matrix3::new(
            factor,
            factor * nu,
            0.0,
            factor * nu,
            factor,
            0.0,
            0.0,
            0.0,
            factor * (1.0 - nu) / 2.0,
        )
    }

    /// Map a Cartesian strain state to plane stress via 2D Hooke's law.
    ///
    /// Builds the constitutive matrix [`plane_stress_matrix`] and
    /// multiplies the engineering-strain vector
    /// `[eps_x, eps_y, gamma_xy]^T`. Note `C[2][2] = G` exactly, so the
    /// shear row reduces to `tau_xy = G gamma_xy`.
    ///
    /// [`plane_stress_matrix`]: ElasticMaterial::plane_stress_matrix
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_strainrosette::{CartesianStrain, ElasticMaterial};
    ///
    /// // Steel-like: E = 200 GPa (in MPa), nu = 0.3.
    /// let mat = ElasticMaterial::new(200_000.0, 0.3).unwrap();
    /// // Uniaxial stress: apply eps_x and the matching lateral
    /// // contraction eps_y = -nu * eps_x, expect sigma_y = 0.
    /// let eps_x = 0.001;
    /// let s = CartesianStrain::new(eps_x, -0.3 * eps_x, 0.0);
    /// let stress = mat.plane_stress(s);
    /// assert!((stress.sigma_x - 200_000.0 * eps_x).abs() < 1e-6);
    /// assert!(stress.sigma_y.abs() < 1e-9);
    /// ```
    pub fn plane_stress(&self, strain: CartesianStrain) -> PlaneStress {
        let c = self.plane_stress_matrix();
        let eps = Vector3::new(strain.eps_x, strain.eps_y, strain.gamma_xy);
        let sigma = c * eps;
        PlaneStress {
            sigma_x: sigma[0],
            sigma_y: sigma[1],
            tau_xy: sigma[2],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn rejects_non_positive_modulus() {
        let err = ElasticMaterial::new(0.0, 0.3).unwrap_err();
        assert_eq!(err.code(), "rosette.invalid-material");
        assert!(ElasticMaterial::new(-1.0, 0.3).is_err());
    }

    #[test]
    fn rejects_out_of_range_poisson() {
        assert!(ElasticMaterial::new(1.0, 0.5).is_err());
        assert!(ElasticMaterial::new(1.0, -1.0).is_err());
        assert!(ElasticMaterial::new(1.0, 0.6).is_err());
        // Interior values are fine.
        assert!(ElasticMaterial::new(1.0, 0.0).is_ok());
        assert!(ElasticMaterial::new(1.0, 0.49).is_ok());
    }

    #[test]
    fn rejects_non_finite_inputs() {
        assert_eq!(
            ElasticMaterial::new(f64::NAN, 0.3).unwrap_err().code(),
            "rosette.non-finite"
        );
        assert_eq!(
            ElasticMaterial::new(1.0, f64::INFINITY).unwrap_err().code(),
            "rosette.non-finite"
        );
    }

    #[test]
    fn shear_modulus_matches_closed_form() {
        let mat = ElasticMaterial::new(200_000.0, 0.3).unwrap();
        let expected = 200_000.0 / (2.0 * 1.3);
        assert!(
            (mat.shear_modulus() - expected).abs() < 1e-6,
            "G = {g}",
            g = mat.shear_modulus()
        );
    }

    #[test]
    fn uniaxial_stress_recovers_hookes_law_sigma_equals_e_eps() {
        // A bar in uniaxial stress along x has lateral strain
        // eps_y = -nu * eps_x and eps_z out of plane; the plane-stress
        // law must then give sigma_x = E * eps_x exactly and sigma_y = 0.
        let e = 70_000.0; // aluminium-like, MPa
        let nu = 0.33;
        let mat = ElasticMaterial::new(e, nu).unwrap();
        let eps_x = 1.5e-3;
        let strain = CartesianStrain::new(eps_x, -nu * eps_x, 0.0);
        let s = mat.plane_stress(strain);
        assert!(
            (s.sigma_x - e * eps_x).abs() < 1e-6,
            "sigma_x = {sx}, expected {ex}",
            sx = s.sigma_x,
            ex = e * eps_x
        );
        assert!(s.sigma_y.abs() < EPS, "sigma_y = {sy}", sy = s.sigma_y);
        assert!(s.tau_xy.abs() < EPS, "tau_xy = {t}", t = s.tau_xy);
    }

    #[test]
    fn equibiaxial_strain_matches_closed_form() {
        // eps_x = eps_y = eps, gamma = 0 => sigma_x = sigma_y =
        // E*eps/(1-nu) by the plane-stress law; tau_xy = 0.
        let e = 100_000.0;
        let nu = 0.25;
        let mat = ElasticMaterial::new(e, nu).unwrap();
        let eps = 8.0e-4;
        let s = mat.plane_stress(CartesianStrain::new(eps, eps, 0.0));
        let expected = e * eps / (1.0 - nu);
        assert!(
            (s.sigma_x - expected).abs() < 1e-6,
            "sigma_x = {sx}",
            sx = s.sigma_x
        );
        assert!(
            (s.sigma_y - expected).abs() < 1e-6,
            "sigma_y = {sy}",
            sy = s.sigma_y
        );
        assert!(s.tau_xy.abs() < EPS);
    }

    #[test]
    fn shear_row_reduces_to_g_times_gamma() {
        // Pure shear strain: tau_xy must equal G * gamma_xy and the
        // direct stresses vanish.
        let mat = ElasticMaterial::new(210_000.0, 0.3).unwrap();
        let gamma = 4.0e-4;
        let s = mat.plane_stress(CartesianStrain::new(0.0, 0.0, gamma));
        assert!(s.sigma_x.abs() < EPS);
        assert!(s.sigma_y.abs() < EPS);
        assert!(
            (s.tau_xy - mat.shear_modulus() * gamma).abs() < 1e-6,
            "tau_xy = {t}",
            t = s.tau_xy
        );
    }

    #[test]
    fn constitutive_matrix_is_symmetric_with_g_in_shear_slot() {
        let mat = ElasticMaterial::new(123_456.0, 0.31).unwrap();
        let c = mat.plane_stress_matrix();
        // Symmetry of the upper-left 2x2 coupling block.
        assert!((c[(0, 1)] - c[(1, 0)]).abs() < 1e-9);
        // Off-axis coupling to shear is zero (isotropy).
        assert!(c[(0, 2)].abs() < EPS && c[(1, 2)].abs() < EPS);
        assert!(c[(2, 0)].abs() < EPS && c[(2, 1)].abs() < EPS);
        // The (2,2) entry equals the shear modulus.
        assert!(
            (c[(2, 2)] - mat.shear_modulus()).abs() < 1e-6,
            "C22 = {c22}, G = {g}",
            c22 = c[(2, 2)],
            g = mat.shear_modulus()
        );
    }
}
