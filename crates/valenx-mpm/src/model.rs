//! Elastic constitutive models.
//!
//! Both supplied models map a deformation gradient `F` to the first
//! Piola–Kirchhoff (PK1) stress `P`, which the MLS-MPM step converts into the
//! grid force contribution. Parameters are the Lamé constants `(mu, lambda)`,
//! derivable from an engineering Young's modulus `E` and Poisson's ratio `nu`.

use crate::math::Mat2;

/// Lamé parameters for an isotropic linear-elastic-at-small-strain material.
///
/// Frame: solid mechanics of engineering materials (metals, polymers, soils),
/// large-deformation regimes where mesh-based FEM elements would invert.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElasticParams {
    /// Shear modulus `μ` (first Lamé parameter), in pascals.
    pub mu: f64,
    /// Lamé's `λ` (second Lamé parameter), in pascals.
    pub lambda: f64,
}

impl ElasticParams {
    /// Builds Lamé parameters from Young's modulus `e` (Pa) and Poisson's
    /// ratio `nu` (dimensionless, `-1 < nu < 0.5`).
    ///
    /// `μ = E / (2(1+ν))`, `λ = Eν / ((1+ν)(1-2ν))`.
    ///
    /// # Panics
    /// Panics if `e <= 0`, or `nu` is outside `(-1, 0.5)`, or non-finite —
    /// these produce a non-physical (non-positive-definite) stiffness.
    #[must_use]
    pub fn from_youngs(e: f64, nu: f64) -> Self {
        assert!(
            e.is_finite() && e > 0.0,
            "Young's modulus must be finite and positive"
        );
        assert!(
            nu.is_finite() && nu > -1.0 && nu < 0.5,
            "Poisson's ratio must lie in (-1, 0.5)"
        );
        let mu = e / (2.0 * (1.0 + nu));
        let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        Self { mu, lambda }
    }

    /// Returns `true` iff both parameters are finite.
    #[must_use]
    pub fn is_finite(self) -> bool {
        self.mu.is_finite() && self.lambda.is_finite()
    }
}

/// The elastic energy density / stress law applied to each particle.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConstitutiveModel {
    /// Neo-Hookean solid (compressible). PK1:
    /// `P = μ (F − F⁻ᵀ) + λ ln(J) F⁻ᵀ`, with `J = det F`.
    NeoHookean,
    /// Fixed-corotated solid (Stomakhin et al. 2012). PK1:
    /// `P = 2μ (F − R) + λ (J − 1) J F⁻ᵀ`, with `F = R S` the polar
    /// decomposition and `J = det F`. Robust under large rotations.
    FixedCorotated,
}

impl ConstitutiveModel {
    /// Computes the first Piola–Kirchhoff stress `P(F)` for this model.
    ///
    /// Returns `None` when `F` is (numerically) singular or non-invertible, so
    /// callers fail loud rather than propagating `NaN`/`inf` through the grid.
    #[must_use]
    pub fn pk1_stress(self, f: Mat2, p: ElasticParams) -> Option<Mat2> {
        if !f.is_finite() {
            return None;
        }
        let j = f.determinant();
        if !j.is_finite() || j.abs() <= 1e-12 {
            return None;
        }
        // F⁻ᵀ = (1/J) · [[m11, -m10], [-m01, m00]]  (inverse-transpose, 2x2).
        let inv = 1.0 / j;
        let f_inv_t = Mat2::new(f.m11 * inv, -f.m10 * inv, -f.m01 * inv, f.m00 * inv);

        let stress = match self {
            ConstitutiveModel::NeoHookean => {
                let ln_j = j.ln();
                if !ln_j.is_finite() {
                    return None;
                }
                // μ(F − F⁻ᵀ) + λ ln(J) F⁻ᵀ
                f.minus(f_inv_t)
                    .scale(p.mu)
                    .plus(f_inv_t.scale(p.lambda * ln_j))
            }
            ConstitutiveModel::FixedCorotated => {
                let r = f.polar_rotation();
                // 2μ(F − R) + λ(J − 1)J F⁻ᵀ
                f.minus(r)
                    .scale(2.0 * p.mu)
                    .plus(f_inv_t.scale(p.lambda * (j - 1.0) * j))
            }
        };
        if stress.is_finite() {
            Some(stress)
        } else {
            None
        }
    }
}
