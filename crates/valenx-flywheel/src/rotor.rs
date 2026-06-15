//! Rotor geometries and their mass moments of inertia.
//!
//! The [`Rotor`] enum captures the canonical axisymmetric rotor shapes
//! used in flywheel sizing, each spinning about its geometric axis of
//! symmetry. Every variant is built through a validated constructor so a
//! constructed [`Rotor`] is always physically meaningful (positive mass,
//! positive radii, an annulus with `r_in < r_out`).
//!
//! All inertias are the textbook rigid-body results for a body of
//! uniform density. See the crate-level docs for the closed forms and
//! their idealising assumptions.

use serde::{Deserialize, Serialize};

use crate::error::FlywheelError;

/// A rigid rotor of uniform density spinning about its symmetry axis.
///
/// Construct via the validated builders ([`Rotor::solid_disk`],
/// [`Rotor::thin_ring`], [`Rotor::annular_disk`], [`Rotor::point_mass`],
/// [`Rotor::thin_rod`]) rather than the variants directly, so invariants
/// (positive mass / radii, `r_in < r_out`) are guaranteed.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Rotor {
    /// A solid uniform disk (or cylinder) of mass `m` and radius `r`
    /// about its central axis. Inertia `I = 1/2 m r^2`.
    SolidDisk {
        /// Rotor mass (kg).
        mass: f64,
        /// Outer radius (m).
        radius: f64,
    },

    /// A thin ring / rim with effectively all of its mass concentrated
    /// at radius `r`. Inertia `I = m r^2`.
    ThinRing {
        /// Rotor mass (kg).
        mass: f64,
        /// Ring radius (m).
        radius: f64,
    },

    /// A uniform annular (hollow) disk of mass `m` between an inner
    /// radius `r_in` and an outer radius `r_out`. Inertia
    /// `I = 1/2 m (r_in^2 + r_out^2)`.
    AnnularDisk {
        /// Rotor mass (kg).
        mass: f64,
        /// Inner radius (m), with `r_in < r_out`.
        r_in: f64,
        /// Outer radius (m).
        r_out: f64,
    },

    /// A point mass `m` orbiting at radius `r`. Inertia `I = m r^2`
    /// (identical form to a thin ring; kept distinct for modelling
    /// clarity).
    PointMass {
        /// Point mass (kg).
        mass: f64,
        /// Orbit radius (m).
        radius: f64,
    },

    /// A thin uniform rod of length `L` rotating about a transverse axis
    /// through its centre. Inertia `I = 1/12 m L^2`.
    ThinRod {
        /// Rod mass (kg).
        mass: f64,
        /// Rod length (m).
        length: f64,
    },
}

impl Rotor {
    /// Build a [`Rotor::SolidDisk`] after checking `mass > 0` and
    /// `radius > 0`.
    ///
    /// # Errors
    ///
    /// Returns [`FlywheelError::InvalidParameter`] if either argument is
    /// non-finite or non-positive.
    pub fn solid_disk(mass: f64, radius: f64) -> Result<Self, FlywheelError> {
        let mass = FlywheelError::require_positive("mass", mass)?;
        let radius = FlywheelError::require_positive("radius", radius)?;
        Ok(Self::SolidDisk { mass, radius })
    }

    /// Build a [`Rotor::ThinRing`] after checking `mass > 0` and
    /// `radius > 0`.
    ///
    /// # Errors
    ///
    /// Returns [`FlywheelError::InvalidParameter`] if either argument is
    /// non-finite or non-positive.
    pub fn thin_ring(mass: f64, radius: f64) -> Result<Self, FlywheelError> {
        let mass = FlywheelError::require_positive("mass", mass)?;
        let radius = FlywheelError::require_positive("radius", radius)?;
        Ok(Self::ThinRing { mass, radius })
    }

    /// Build a [`Rotor::AnnularDisk`] after checking `mass > 0`,
    /// `r_in > 0`, `r_out > 0`, and `r_in < r_out`.
    ///
    /// # Errors
    ///
    /// Returns [`FlywheelError::InvalidParameter`] for a non-finite or
    /// non-positive scalar, or [`FlywheelError::Inconsistent`] if the
    /// inner radius is not strictly less than the outer radius.
    pub fn annular_disk(mass: f64, r_in: f64, r_out: f64) -> Result<Self, FlywheelError> {
        let mass = FlywheelError::require_positive("mass", mass)?;
        let r_in = FlywheelError::require_positive("r_in", r_in)?;
        let r_out = FlywheelError::require_positive("r_out", r_out)?;
        if r_in >= r_out {
            return Err(FlywheelError::Inconsistent("annulus requires r_in < r_out"));
        }
        Ok(Self::AnnularDisk { mass, r_in, r_out })
    }

    /// Build a [`Rotor::PointMass`] after checking `mass > 0` and
    /// `radius > 0`.
    ///
    /// # Errors
    ///
    /// Returns [`FlywheelError::InvalidParameter`] if either argument is
    /// non-finite or non-positive.
    pub fn point_mass(mass: f64, radius: f64) -> Result<Self, FlywheelError> {
        let mass = FlywheelError::require_positive("mass", mass)?;
        let radius = FlywheelError::require_positive("radius", radius)?;
        Ok(Self::PointMass { mass, radius })
    }

    /// Build a [`Rotor::ThinRod`] after checking `mass > 0` and
    /// `length > 0`.
    ///
    /// # Errors
    ///
    /// Returns [`FlywheelError::InvalidParameter`] if either argument is
    /// non-finite or non-positive.
    pub fn thin_rod(mass: f64, length: f64) -> Result<Self, FlywheelError> {
        let mass = FlywheelError::require_positive("mass", mass)?;
        let length = FlywheelError::require_positive("length", length)?;
        Ok(Self::ThinRod { mass, length })
    }

    /// The rotor mass (kg).
    #[must_use]
    pub fn mass(&self) -> f64 {
        match *self {
            Self::SolidDisk { mass, .. }
            | Self::ThinRing { mass, .. }
            | Self::AnnularDisk { mass, .. }
            | Self::PointMass { mass, .. }
            | Self::ThinRod { mass, .. } => mass,
        }
    }

    /// The outermost radius of the rotor (m), used by the rim-stress
    /// model as the radius of the fastest-moving material.
    ///
    /// For a [`Rotor::ThinRod`] this is half the rod length (the tip
    /// radius about the centre).
    #[must_use]
    pub fn outer_radius(&self) -> f64 {
        match *self {
            Self::SolidDisk { radius, .. }
            | Self::ThinRing { radius, .. }
            | Self::PointMass { radius, .. } => radius,
            Self::AnnularDisk { r_out, .. } => r_out,
            Self::ThinRod { length, .. } => 0.5 * length,
        }
    }

    /// The mass moment of inertia about the spin axis (kg.m^2).
    ///
    /// Uses the closed-form rigid-body result for each shape: a solid
    /// disk gives `1/2 m r^2`; a thin ring or point mass gives `m r^2`;
    /// an annular disk gives `1/2 m (r_in^2 + r_out^2)`; and a thin rod
    /// about its centre gives `1/12 m L^2`.
    #[must_use]
    pub fn moment_of_inertia(&self) -> f64 {
        match *self {
            Self::SolidDisk { mass, radius } => 0.5 * mass * radius * radius,
            Self::ThinRing { mass, radius } | Self::PointMass { mass, radius } => {
                mass * radius * radius
            }
            Self::AnnularDisk { mass, r_in, r_out } => 0.5 * mass * (r_in * r_in + r_out * r_out),
            Self::ThinRod { mass, length } => mass * length * length / 12.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    #[test]
    fn solid_disk_inertia_is_half_m_r_squared() {
        // m = 4 kg, r = 0.5 m -> I = 0.5 * 4 * 0.25 = 0.5
        let rotor = Rotor::solid_disk(4.0, 0.5).unwrap();
        assert!((rotor.moment_of_inertia() - 0.5).abs() < EPS);
    }

    #[test]
    fn thin_ring_inertia_is_m_r_squared() {
        // m = 3 kg, r = 2 m -> I = 3 * 4 = 12
        let rotor = Rotor::thin_ring(3.0, 2.0).unwrap();
        assert!((rotor.moment_of_inertia() - 12.0).abs() < EPS);
    }

    #[test]
    fn ring_inertia_is_double_disk_for_same_mass_and_radius() {
        // Classic result: a ring stores twice the inertia of a solid
        // disk of equal mass and radius.
        let m = 5.0;
        let r = 1.3;
        let disk = Rotor::solid_disk(m, r).unwrap();
        let ring = Rotor::thin_ring(m, r).unwrap();
        assert!((ring.moment_of_inertia() - 2.0 * disk.moment_of_inertia()).abs() < EPS);
    }

    #[test]
    fn point_mass_matches_thin_ring() {
        let pm = Rotor::point_mass(2.5, 0.8).unwrap();
        let ring = Rotor::thin_ring(2.5, 0.8).unwrap();
        assert!((pm.moment_of_inertia() - ring.moment_of_inertia()).abs() < EPS);
    }

    #[test]
    fn annular_disk_inertia_matches_closed_form() {
        // m = 2, r_in = 1, r_out = 2 -> I = 0.5 * 2 * (1 + 4) = 5
        let rotor = Rotor::annular_disk(2.0, 1.0, 2.0).unwrap();
        assert!((rotor.moment_of_inertia() - 5.0).abs() < EPS);
    }

    #[test]
    fn annular_disk_degenerates_to_thin_ring_as_inner_approaches_outer() {
        // With r_in -> r_out = r, I -> 0.5 m (r^2 + r^2) = m r^2 (ring).
        let m = 7.0;
        let r = 1.5;
        let nearly = Rotor::annular_disk(m, r - 1e-7, r).unwrap();
        let ring = Rotor::thin_ring(m, r).unwrap();
        assert!((nearly.moment_of_inertia() - ring.moment_of_inertia()).abs() < 1e-4);
    }

    #[test]
    fn thin_rod_inertia_is_m_l_squared_over_twelve() {
        // m = 12 kg, L = 1 m -> I = 12 * 1 / 12 = 1
        let rotor = Rotor::thin_rod(12.0, 1.0).unwrap();
        assert!((rotor.moment_of_inertia() - 1.0).abs() < EPS);
    }

    #[test]
    fn outer_radius_of_rod_is_half_length() {
        let rotor = Rotor::thin_rod(2.0, 3.0).unwrap();
        assert!((rotor.outer_radius() - 1.5).abs() < EPS);
    }

    #[test]
    fn constructors_reject_non_positive() {
        assert!(Rotor::solid_disk(0.0, 1.0).is_err());
        assert!(Rotor::solid_disk(1.0, -1.0).is_err());
        assert!(Rotor::thin_ring(-2.0, 1.0).is_err());
        assert!(Rotor::thin_rod(1.0, 0.0).is_err());
    }

    #[test]
    fn constructors_reject_non_finite() {
        assert!(Rotor::solid_disk(f64::NAN, 1.0).is_err());
        assert!(Rotor::solid_disk(1.0, f64::INFINITY).is_err());
    }

    #[test]
    fn annular_rejects_inner_ge_outer() {
        let err = Rotor::annular_disk(1.0, 2.0, 1.0).unwrap_err();
        assert_eq!(
            err,
            FlywheelError::Inconsistent("annulus requires r_in < r_out")
        );
        assert!(Rotor::annular_disk(1.0, 2.0, 2.0).is_err());
    }

    #[test]
    fn mass_accessor_round_trips() {
        assert!((Rotor::solid_disk(4.0, 0.5).unwrap().mass() - 4.0).abs() < EPS);
        assert!((Rotor::annular_disk(2.0, 1.0, 2.0).unwrap().mass() - 2.0).abs() < EPS);
    }
}
