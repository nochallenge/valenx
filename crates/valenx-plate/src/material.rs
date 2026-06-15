//! Linear-elastic isotropic plate material and its flexural rigidity.
//!
//! ## Model
//!
//! A thin flat plate of uniform thickness `t` made of a linear-elastic,
//! isotropic material with Young's modulus `E` and Poisson's ratio `nu`
//! resists transverse bending with a *flexural rigidity*
//!
//! ```text
//! D = E t^3 / (12 (1 - nu^2))
//! ```
//!
//! `D` plays the role for a plate that the bending stiffness `E I` plays
//! for a beam: it is the constant of proportionality between curvature and
//! bending moment per unit width. It carries SI units of N·m (force times
//! length) when `E` is in pascals and `t` in metres.
//!
//! ## Honest scope
//!
//! Isotropic linear elasticity only — no orthotropy, plasticity,
//! temperature dependence, or large-deflection membrane stiffening. The
//! `1 / (1 - nu^2)` factor is exactly the plane-stress-to-bending
//! correction of Kirchhoff-Love thin-plate theory and is undefined as
//! `nu -> 1`, which is why [`PlateMaterial::new`] rejects `nu` outside the
//! open interval `(-1, 0.5)`.

use serde::{Deserialize, Serialize};

use crate::error::{require_positive, PlateError};

/// A linear-elastic, isotropic plate material of uniform thickness.
///
/// Construct with [`PlateMaterial::new`], which validates every field, then
/// read the derived [`flexural_rigidity`](PlateMaterial::flexural_rigidity).
///
/// Units are not fixed by the type, but must be *consistent*: with `E` in
/// pascals (N/m^2) and `thickness` in metres, the flexural rigidity comes
/// out in newton-metres (N·m).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PlateMaterial {
    /// Young's modulus `E` (must be finite and strictly positive).
    pub youngs_modulus: f64,
    /// Poisson's ratio `nu` (must lie in the open interval `(-1, 0.5)`).
    pub poisson_ratio: f64,
    /// Plate thickness `t` (must be finite and strictly positive).
    pub thickness: f64,
}

impl PlateMaterial {
    /// The exclusive lower bound on Poisson's ratio (thermodynamic limit).
    pub const POISSON_MIN: f64 = -1.0;
    /// The exclusive upper bound on Poisson's ratio for an isotropic solid.
    pub const POISSON_MAX: f64 = 0.5;

    /// Build a validated [`PlateMaterial`].
    ///
    /// # Errors
    ///
    /// Returns [`PlateError::InvalidParameter`] if `youngs_modulus` or
    /// `thickness` is not finite and strictly positive, or if
    /// `poisson_ratio` lies outside the open interval
    /// (`POISSON_MIN`, `POISSON_MAX`). The open interval is required because
    /// the `1 - nu^2` denominator of the flexural rigidity vanishes at the
    /// boundaries and changes sign beyond them.
    pub fn new(
        youngs_modulus: f64,
        poisson_ratio: f64,
        thickness: f64,
    ) -> Result<Self, PlateError> {
        let youngs_modulus = require_positive("youngs_modulus", youngs_modulus)?;
        let thickness = require_positive("thickness", thickness)?;

        if !poisson_ratio.is_finite() {
            return Err(PlateError::invalid(
                "poisson_ratio",
                "must be a finite number",
                poisson_ratio,
            ));
        }
        if poisson_ratio <= Self::POISSON_MIN || poisson_ratio >= Self::POISSON_MAX {
            return Err(PlateError::invalid(
                "poisson_ratio",
                format!(
                    "must lie in the open interval ({lo}, {hi})",
                    lo = Self::POISSON_MIN,
                    hi = Self::POISSON_MAX
                ),
                poisson_ratio,
            ));
        }

        Ok(Self {
            youngs_modulus,
            poisson_ratio,
            thickness,
        })
    }

    /// Flexural rigidity `D = E t^3 / (12 (1 - nu^2))`.
    ///
    /// Always finite and strictly positive for a validly-constructed
    /// material, since `E, t > 0` and `0 < 1 - nu^2 <= 1` on the admissible
    /// `nu` interval.
    pub fn flexural_rigidity(&self) -> f64 {
        let e = self.youngs_modulus;
        let t = self.thickness;
        let nu = self.poisson_ratio;
        e * t * t * t / (12.0 * (1.0 - nu * nu))
    }
}
