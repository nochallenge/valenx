//! Uniformly-loaded circular plate: centre deflection and maximum stress.
//!
//! ## Model
//!
//! A flat circular plate of radius `a` and uniform thickness `t`, carrying
//! a uniform transverse pressure `p` over its whole face, bends according
//! to small-deflection Kirchhoff-Love thin-plate theory. The maximum
//! transverse deflection occurs at the centre and is written
//!
//! ```text
//! w = k * p * a^4 / D
//! ```
//!
//! where `D` is the [flexural rigidity](crate::PlateMaterial::flexural_rigidity)
//! and `k` is a dimensionless coefficient fixed by the edge support. The
//! coefficients below are the classical results (Timoshenko & Woinowsky-
//! Krieger, *Theory of Plates and Shells*, 2nd ed., Art. 16-17):
//!
//! ```text
//! clamped edge:            k = 1 / 64                      ~= 0.015625
//! simply-supported edge:   k = (5 + nu) / (64 (1 + nu))
//! ```
//!
//! The extreme fibre bending stress is `sigma = 6 M / t^2`, evaluated where
//! the bending moment per unit width `M` is largest:
//!
//! ```text
//! clamped:           M_max = p a^2 / 8          (radial moment, at the edge)
//!                    sigma_max = 3 p a^2 / (4 t^2)         = 0.75  * p (a/t)^2
//! simply-supported:  M_max = (3 + nu) p a^2 / 16   (at the centre)
//!                    sigma_max = 3 (3 + nu) p a^2 / (8 t^2)
//! ```
//!
//! Both stresses scale linearly with `p` and quadratically with the
//! radius-to-thickness ratio `a / t`, and the clamped plate is the stiffer
//! of the two (smaller centre deflection, since for any admissible `nu` the
//! simply-supported coefficient exceeds `1/64`).
//!
//! ## Honest scope
//!
//! Small-deflection linear theory only: deflections are assumed small
//! relative to the thickness, so there is no membrane (stretching)
//! stiffening, no contact, no plasticity, and no transverse-shear
//! correction. Results are a textbook first estimate, not a substitute for
//! a validated finite-element analysis.

use serde::{Deserialize, Serialize};

use crate::error::{require_positive, PlateError};
use crate::material::PlateMaterial;

/// How the rim of a circular plate is supported.
///
/// The choice fixes the dimensionless deflection coefficient and the
/// location / magnitude of the peak bending moment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum EdgeSupport {
    /// Built-in / encastre edge: zero deflection *and* zero slope at the rim.
    Clamped,
    /// Simply-supported edge: zero deflection but free rotation at the rim.
    SimplySupported,
}

impl EdgeSupport {
    /// Dimensionless centre-deflection coefficient `k` for `w = k p a^4 / D`.
    ///
    /// For a [`Clamped`](EdgeSupport::Clamped) edge this is the constant
    /// `1/64`; for a [`SimplySupported`](EdgeSupport::SimplySupported) edge
    /// it depends on Poisson's ratio, `(5 + nu) / (64 (1 + nu))`.
    pub fn deflection_coefficient(self, poisson_ratio: f64) -> f64 {
        match self {
            EdgeSupport::Clamped => 1.0 / 64.0,
            EdgeSupport::SimplySupported => (5.0 + poisson_ratio) / (64.0 * (1.0 + poisson_ratio)),
        }
    }
}

/// A uniformly-loaded, simply-shaped circular plate problem.
///
/// Bundles the [`PlateMaterial`] (which carries `E`, `nu` and `t`) with the
/// plate radius, the applied uniform pressure, and the edge support, and
/// exposes the closed-form centre deflection and maximum bending stress.
///
/// Construct via [`CircularPlate::new`]; all inputs are validated.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CircularPlate {
    /// Plate material and thickness.
    pub material: PlateMaterial,
    /// Plate radius `a` (must be finite and strictly positive).
    pub radius: f64,
    /// Uniform transverse pressure `p` over the face (finite, positive).
    pub pressure: f64,
    /// Edge-support condition.
    pub support: EdgeSupport,
}

impl CircularPlate {
    /// The minimum radius-to-thickness ratio `a / t` this crate accepts.
    ///
    /// Below roughly this value the plate is no longer "thin": transverse
    /// shear and thickness effects that Kirchhoff-Love theory ignores stop
    /// being negligible, so the closed-form results would mislead. The
    /// threshold of `5` is deliberately permissive — engineering practice
    /// usually takes `a / t >= 10..20` as "thin" — and is enforced so the
    /// API never silently returns a result for a manifestly thick plate.
    pub const MIN_RADIUS_THICKNESS_RATIO: f64 = 5.0;

    /// Build a validated circular-plate problem.
    ///
    /// # Errors
    ///
    /// Returns [`PlateError::InvalidParameter`] if `radius` or `pressure` is
    /// not finite and strictly positive, and [`PlateError::NotThin`] if the
    /// ratio `radius / material.thickness` is below
    /// [`MIN_RADIUS_THICKNESS_RATIO`](Self::MIN_RADIUS_THICKNESS_RATIO).
    pub fn new(
        material: PlateMaterial,
        radius: f64,
        pressure: f64,
        support: EdgeSupport,
    ) -> Result<Self, PlateError> {
        let radius = require_positive("radius", radius)?;
        let pressure = require_positive("pressure", pressure)?;

        let ratio = radius / material.thickness;
        if ratio < Self::MIN_RADIUS_THICKNESS_RATIO {
            return Err(PlateError::NotThin {
                ratio,
                min: Self::MIN_RADIUS_THICKNESS_RATIO,
            });
        }

        Ok(Self {
            material,
            radius,
            pressure,
            support,
        })
    }

    /// Flexural rigidity `D` of the plate (delegates to the material).
    pub fn flexural_rigidity(&self) -> f64 {
        self.material.flexural_rigidity()
    }

    /// Dimensionless centre-deflection coefficient `k` for this problem.
    pub fn deflection_coefficient(&self) -> f64 {
        self.support
            .deflection_coefficient(self.material.poisson_ratio)
    }

    /// Maximum (centre) transverse deflection `w = k p a^4 / D`.
    ///
    /// Returned in the same length unit as the inputs (metres if `E` is in
    /// pascals, `p` in pascals, and lengths in metres). Always positive for
    /// a validly-constructed plate.
    pub fn center_deflection(&self) -> f64 {
        let k = self.deflection_coefficient();
        let a = self.radius;
        let a4 = a * a * a * a;
        k * self.pressure * a4 / self.flexural_rigidity()
    }

    /// Largest extreme-fibre bending stress `sigma = 6 M_max / t^2`.
    ///
    /// For a [`Clamped`](EdgeSupport::Clamped) plate the peak is the radial
    /// bending stress at the rim, `3 p a^2 / (4 t^2)`. For a
    /// [`SimplySupported`](EdgeSupport::SimplySupported) plate the peak is at
    /// the centre, `3 (3 + nu) p a^2 / (8 t^2)`. Both scale with `p` and
    /// with `(a / t)^2`. Returned as a positive magnitude in the pressure's
    /// stress unit.
    pub fn max_bending_stress(&self) -> f64 {
        let p = self.pressure;
        let a = self.radius;
        let t = self.material.thickness;
        let nu = self.material.poisson_ratio;
        let a_over_t = a / t;
        let scale = p * a_over_t * a_over_t; // p (a/t)^2

        match self.support {
            EdgeSupport::Clamped => 0.75 * scale,
            EdgeSupport::SimplySupported => 3.0 * (3.0 + nu) / 8.0 * scale,
        }
    }
}
