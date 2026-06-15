//! Circular cross-sections and their polar second moment of area.
//!
//! ## Model
//!
//! For a **solid** round bar of diameter `d` the polar second moment of
//! area about the centroidal axis is
//!
//! ```text
//! J = pi * d^4 / 32
//! ```
//!
//! For a **hollow** round bar (outer diameter `D`, bore `d`) the annulus
//! is the difference of two solid sections:
//!
//! ```text
//! J = pi * (D^4 - d^4) / 32
//! ```
//!
//! `J` carries units of length to the fourth power; whatever length unit
//! the diameter is given in, `J` is in that unit to the fourth power.
//!
//! These are the exact elastic-theory results for a prismatic bar; for a
//! circular section the polar second moment of area equals the torsion
//! constant, so the same `J` drives both the stress and twist formulas in
//! [`crate::response`].

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

use crate::error::{require_positive, TorsionError};

/// A circular shaft cross-section: either a solid bar or a hollow tube.
///
/// Construct one with [`Shaft::solid`] or [`Shaft::hollow`]; both validate
/// their inputs and return a [`TorsionError`] on a non-positive diameter
/// or an inverted annulus. Once built, a `Shaft` is guaranteed to have a
/// strictly positive [`Shaft::polar_moment`].
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Shaft {
    /// A solid round bar.
    Solid {
        /// Bar diameter (`> 0`), in the caller's length unit.
        diameter: f64,
    },
    /// A hollow round tube (a circular annulus).
    Hollow {
        /// Outer diameter (`> 0`).
        outer_diameter: f64,
        /// Inner / bore diameter (`> 0` and strictly less than the outer).
        inner_diameter: f64,
    },
}

impl Shaft {
    /// Construct a solid round shaft of the given `diameter`.
    ///
    /// # Errors
    ///
    /// Returns [`TorsionError::NonPositive`] if `diameter` is not finite
    /// and strictly positive.
    pub fn solid(diameter: f64) -> Result<Self, TorsionError> {
        let diameter = require_positive("diameter", diameter)?;
        Ok(Shaft::Solid { diameter })
    }

    /// Construct a hollow round shaft from its outer and inner diameters.
    ///
    /// # Errors
    ///
    /// Returns [`TorsionError::NonPositive`] if either diameter is not
    /// finite and strictly positive, or [`TorsionError::InvertedAnnulus`]
    /// if `inner_diameter` is not strictly less than `outer_diameter`.
    pub fn hollow(outer_diameter: f64, inner_diameter: f64) -> Result<Self, TorsionError> {
        let outer_diameter = require_positive("outer_diameter", outer_diameter)?;
        let inner_diameter = require_positive("inner_diameter", inner_diameter)?;
        if inner_diameter >= outer_diameter {
            return Err(TorsionError::InvertedAnnulus {
                inner: inner_diameter,
                outer: outer_diameter,
            });
        }
        Ok(Shaft::Hollow {
            outer_diameter,
            inner_diameter,
        })
    }

    /// Outer diameter of the section.
    ///
    /// For a solid bar this is simply its diameter; for a tube it is the
    /// outer diameter.
    pub fn outer_diameter(&self) -> f64 {
        match *self {
            Shaft::Solid { diameter } => diameter,
            Shaft::Hollow { outer_diameter, .. } => outer_diameter,
        }
    }

    /// Inner (bore) diameter of the section: zero for a solid bar.
    pub fn inner_diameter(&self) -> f64 {
        match *self {
            Shaft::Solid { .. } => 0.0,
            Shaft::Hollow { inner_diameter, .. } => inner_diameter,
        }
    }

    /// Outer radius of the section (half the outer diameter).
    ///
    /// This is the radius at which the shear stress is largest, and the
    /// value used by [`crate::response::max_shear_stress`].
    pub fn outer_radius(&self) -> f64 {
        self.outer_diameter() / 2.0
    }

    /// Inner (bore) radius of the section: zero for a solid bar.
    pub fn inner_radius(&self) -> f64 {
        self.inner_diameter() / 2.0
    }

    /// Polar second moment of area `J` of the section.
    ///
    /// `J = pi d^4 / 32` for a solid bar and
    /// `J = pi (D^4 - d^4) / 32` for a tube. The returned value is always
    /// strictly positive because the constructors reject degenerate
    /// sections.
    pub fn polar_moment(&self) -> f64 {
        match *self {
            Shaft::Solid { diameter } => PI * diameter.powi(4) / 32.0,
            Shaft::Hollow {
                outer_diameter,
                inner_diameter,
            } => PI * (outer_diameter.powi(4) - inner_diameter.powi(4)) / 32.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tight tolerance for comparing two exact closed-form expressions.
    const EPS: f64 = 1e-9;

    #[test]
    fn solid_polar_moment_matches_pi_d4_over_32() {
        let d = 50.0_f64;
        let shaft = Shaft::solid(d).unwrap();
        let expected = PI * d.powi(4) / 32.0;
        assert!((shaft.polar_moment() - expected).abs() < EPS * expected);
    }

    #[test]
    fn hollow_polar_moment_matches_pi_big4_minus_small4_over_32() {
        let big = 60.0_f64;
        let small = 40.0_f64;
        let shaft = Shaft::hollow(big, small).unwrap();
        let expected = PI * (big.powi(4) - small.powi(4)) / 32.0;
        assert!((shaft.polar_moment() - expected).abs() < EPS * expected);
    }

    #[test]
    fn hollow_with_zero_bore_limits_to_solid() {
        // A vanishingly small bore should approach the solid result.
        let d = 30.0_f64;
        let solid = Shaft::solid(d).unwrap();
        let nearly_solid = Shaft::hollow(d, 1e-6).unwrap();
        assert!((solid.polar_moment() - nearly_solid.polar_moment()).abs() < 1e-3);
    }

    #[test]
    fn doubling_solid_diameter_multiplies_polar_moment_by_sixteen() {
        // J scales as d^4, so 2x diameter gives 2^4 = 16x J.
        let base = Shaft::solid(10.0).unwrap();
        let doubled = Shaft::solid(20.0).unwrap();
        let ratio = doubled.polar_moment() / base.polar_moment();
        assert!((ratio - 16.0).abs() < EPS);
    }

    #[test]
    fn radii_are_half_the_diameters() {
        let solid = Shaft::solid(8.0).unwrap();
        assert!((solid.outer_radius() - 4.0).abs() < EPS);
        assert!((solid.inner_radius() - 0.0).abs() < EPS);

        let hollow = Shaft::hollow(8.0, 6.0).unwrap();
        assert!((hollow.outer_radius() - 4.0).abs() < EPS);
        assert!((hollow.inner_radius() - 3.0).abs() < EPS);
    }

    #[test]
    fn solid_rejects_non_positive_diameter() {
        assert!(matches!(
            Shaft::solid(0.0),
            Err(TorsionError::NonPositive {
                name: "diameter",
                ..
            })
        ));
        assert!(matches!(
            Shaft::solid(-3.0),
            Err(TorsionError::NonPositive { .. })
        ));
    }

    #[test]
    fn hollow_rejects_inverted_or_equal_annulus() {
        assert!(matches!(
            Shaft::hollow(10.0, 10.0),
            Err(TorsionError::InvertedAnnulus { .. })
        ));
        assert!(matches!(
            Shaft::hollow(10.0, 12.0),
            Err(TorsionError::InvertedAnnulus { .. })
        ));
    }

    #[test]
    fn hollow_rejects_non_positive_diameters() {
        assert!(matches!(
            Shaft::hollow(10.0, 0.0),
            Err(TorsionError::NonPositive {
                name: "inner_diameter",
                ..
            })
        ));
        assert!(matches!(
            Shaft::hollow(-1.0, 0.5),
            Err(TorsionError::NonPositive {
                name: "outer_diameter",
                ..
            })
        ));
    }
}
