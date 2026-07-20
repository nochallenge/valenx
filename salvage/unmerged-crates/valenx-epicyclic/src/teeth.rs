//! Tooth counts and the standard meshing constraint for a simple
//! single-stage sun / planet / ring (annulus) epicyclic set.
//!
//! ## Geometry
//!
//! All three central wheels are coaxial and share one diametral pitch,
//! so a tooth count is proportional to a pitch radius. For the sun to
//! mesh with a planet and the planet to mesh with the internal ring on
//! the same centre line, the pitch radii must satisfy
//! `r_ring = r_sun + 2 r_planet`, i.e. in teeth:
//!
//! ```text
//! R = S + 2 P
//! ```
//!
//! where `S`, `P`, `R` are the sun, planet and ring tooth counts.
//!
//! ## Basic train ratio
//!
//! Holding the carrier stationary turns the set into an ordinary
//! (reverted) gear train. Sun to ring then runs through the planet as
//! an idler, giving the signed basic ratio
//!
//! ```text
//! e = omega_ring / omega_sun  |_(carrier fixed)  =  -S / R
//! ```
//!
//! The minus sign records that, with the carrier held, the ring turns
//! opposite to the sun (one external mesh sun-planet, one internal mesh
//! planet-ring). This `e` is the quantity the Willis equation uses.

use crate::error::{require_positive_teeth, EpicyclicError};

/// Validated tooth counts for a simple planetary set.
///
/// Construct with [`TeethSet::new`], which enforces positivity and the
/// `R = S + 2 P` meshing constraint, or with [`TeethSet::from_sun_planet`]
/// which derives the ring from the constraint.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TeethSet {
    sun: u32,
    planet: u32,
    ring: u32,
}

impl TeethSet {
    /// Build a tooth set from explicit sun, planet and ring counts.
    ///
    /// All three are validated as strictly positive, then checked
    /// against the meshing constraint `ring == sun + 2*planet`.
    ///
    /// # Errors
    ///
    /// Returns [`EpicyclicError::NonPositiveTeeth`] if any count is below
    /// one, or [`EpicyclicError::MeshingConstraint`] if the ring count is
    /// inconsistent with the sun and planet counts.
    ///
    /// # Example
    ///
    /// ```
    /// use valenx_epicyclic::TeethSet;
    /// let t = TeethSet::new(24, 24, 72).unwrap();
    /// assert_eq!(t.ring(), 72);
    /// ```
    pub fn new(sun: i64, planet: i64, ring: i64) -> Result<Self, EpicyclicError> {
        let sun = require_positive_teeth("sun", sun)?;
        let planet = require_positive_teeth("planet", planet)?;
        let ring = require_positive_teeth("ring", ring)?;
        let expected = sun + 2 * planet;
        if ring != expected {
            return Err(EpicyclicError::MeshingConstraint {
                sun,
                planet,
                ring,
                expected,
            });
        }
        Ok(Self { sun, planet, ring })
    }

    /// Build a tooth set from the sun and planet counts, deriving the
    /// ring count from the meshing constraint `ring = sun + 2*planet`.
    ///
    /// # Errors
    ///
    /// Returns [`EpicyclicError::NonPositiveTeeth`] if either supplied
    /// count is below one. The derived ring count is always consistent
    /// by construction.
    ///
    /// # Example
    ///
    /// ```
    /// use valenx_epicyclic::TeethSet;
    /// let t = TeethSet::from_sun_planet(24, 24).unwrap();
    /// assert_eq!(t.ring(), 72);
    /// ```
    pub fn from_sun_planet(sun: i64, planet: i64) -> Result<Self, EpicyclicError> {
        let sun = require_positive_teeth("sun", sun)?;
        let planet = require_positive_teeth("planet", planet)?;
        let ring = sun + 2 * planet;
        Ok(Self { sun, planet, ring })
    }

    /// Sun (central) tooth count `S`.
    #[inline]
    pub fn sun(&self) -> u32 {
        self.sun
    }

    /// Planet (orbiting pinion) tooth count `P`.
    #[inline]
    pub fn planet(&self) -> u32 {
        self.planet
    }

    /// Ring / annulus (internal) tooth count `R`.
    #[inline]
    pub fn ring(&self) -> u32 {
        self.ring
    }

    /// Signed basic train ratio `e = omega_ring / omega_sun` with the
    /// carrier held fixed, equal to `-S / R`.
    ///
    /// This is the value fed to the Willis equation. It is always
    /// negative for a sun-planet-ring set because the single planet
    /// reverses the sense between the external sun mesh and the internal
    /// ring mesh.
    ///
    /// # Example
    ///
    /// ```
    /// use valenx_epicyclic::TeethSet;
    /// let t = TeethSet::new(24, 24, 72).unwrap();
    /// assert!((t.basic_train_ratio() - (-1.0 / 3.0)).abs() < 1e-12);
    /// ```
    #[inline]
    pub fn basic_train_ratio(&self) -> f64 {
        -(self.sun as f64) / (self.ring as f64)
    }

    /// Signed ratio of planet spin to sun spin with the carrier held
    /// fixed, `omega_planet / omega_sun = -S / P`.
    ///
    /// The sun drives the planet through one external mesh, so the
    /// planet turns opposite to the sun (hence the minus sign) at a
    /// magnitude set by the inverse tooth ratio.
    ///
    /// # Example
    ///
    /// ```
    /// use valenx_epicyclic::TeethSet;
    /// let t = TeethSet::new(24, 24, 72).unwrap();
    /// assert!((t.planet_sun_ratio_carrier_fixed() - (-1.0)).abs() < 1e-12);
    /// ```
    #[inline]
    pub fn planet_sun_ratio_carrier_fixed(&self) -> f64 {
        -(self.sun as f64) / (self.planet as f64)
    }
}
