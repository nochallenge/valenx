//! Thin-wall (membrane) sphere under internal pressure.
//!
//! ## Model
//!
//! A thin spherical shell of internal radius `r` and wall thickness `t`
//! under gauge internal pressure `p` carries the pressure as a single
//! membrane stress that is, by symmetry, **identical in every in-plane
//! direction**. Balancing the bursting thrust `p * pi * r^2` across any
//! great-circle plane against the ring of wall material `2 * pi * r * t`
//! gives
//!
//! ```text
//! sigma = p * r / (2 t)
//! ```
//!
//! This is exactly the *longitudinal* stress of a cylinder of the same
//! `r`, `t`, `p` — and therefore exactly **half** that cylinder's hoop
//! stress. For a fixed volume the sphere is the most weight-efficient
//! pressure vessel: its single membrane stress is half the governing hoop
//! stress of the equivalent cylinder, which is why high-pressure gas is
//! stored in spheres (see Roark, ch. 13; Boresi & Schmidt, *Advanced
//! Mechanics of Materials*).
//!
//! ## Honest scope
//!
//! Membrane theory only — valid roughly for `r / t > 10`. The radial
//! stress (`-p` at the bore, `0` outside) is neglected against `sigma`,
//! and the model assigns the third principal stress as zero. Bending,
//! nozzle / penetration discontinuities, thermal and residual stress, and
//! all code allowables are out of scope. For thick walls use the exact
//! [`crate::thick`] Lame sphere counterpart of this membrane limit.

use crate::error::{Result, VesselError};
use serde::{Deserialize, Serialize};

/// A thin-wall spherical pressure vessel and its membrane stress.
///
/// Construct with [`ThinSphere::new`]. Stresses are returned in the same
/// unit as the supplied `pressure`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ThinSphere {
    pressure: f64,
    radius: f64,
    thickness: f64,
}

impl ThinSphere {
    /// Build a thin-wall sphere from internal `pressure`, internal
    /// `radius` and wall `thickness` (consistent units).
    ///
    /// # Errors
    ///
    /// Returns [`VesselError::BadParameter`] if any argument is
    /// non-finite, zero or negative.
    pub fn new(pressure: f64, radius: f64, thickness: f64) -> Result<Self> {
        let pressure = VesselError::require_positive("pressure", pressure)?;
        let radius = VesselError::require_positive("radius", radius)?;
        let thickness = VesselError::require_positive("thickness", thickness)?;
        Ok(Self {
            pressure,
            radius,
            thickness,
        })
    }

    /// Internal gauge pressure `p`.
    pub fn pressure(&self) -> f64 {
        self.pressure
    }

    /// Internal radius `r`.
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// Wall thickness `t`.
    pub fn thickness(&self) -> f64 {
        self.thickness
    }

    /// Radius-to-thickness ratio `r / t`.
    pub fn radius_to_thickness(&self) -> f64 {
        self.radius / self.thickness
    }

    /// Whether the thin-wall assumption is reasonable (`r / t >= 10`).
    /// Advisory only.
    pub fn is_thin_wall(&self) -> bool {
        self.radius_to_thickness() >= 10.0
    }

    /// Membrane stress `sigma = p r / (2 t)`, equal in all in-plane
    /// directions. Positive (tensile) for internal pressure.
    pub fn membrane_stress(&self) -> f64 {
        self.pressure * self.radius / (2.0 * self.thickness)
    }

    /// Maximum shear stress. The two in-plane principal stresses are
    /// equal (`sigma`), so the governing pair involves the (near-zero)
    /// radial stress: `(sigma - 0) / 2 = sigma / 2 = p r / (4 t)`.
    pub fn max_shear_stress(&self) -> f64 {
        self.membrane_stress() / 2.0
    }

    /// Von Mises equivalent stress. With equal biaxial principal
    /// stresses `sigma_1 = sigma_2 = sigma` and `sigma_3 = 0`, the von
    /// Mises stress reduces to exactly `sigma` — an equibiaxial state is
    /// as severe (by this criterion) as a uniaxial one of the same
    /// magnitude.
    pub fn von_mises_stress(&self) -> f64 {
        self.membrane_stress()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cylinder::ThinCylinder;

    /// Loose absolute tolerance for f64 stress comparisons.
    const EPS: f64 = 1e-9;

    #[test]
    fn rejects_non_positive_inputs() {
        assert!(ThinSphere::new(0.0, 1.0, 0.1).is_err());
        assert!(ThinSphere::new(1.0, 0.0, 0.1).is_err());
        assert!(ThinSphere::new(1.0, 1.0, -0.1).is_err());
        assert!(ThinSphere::new(1.0, f64::INFINITY, 0.1).is_err());
    }

    #[test]
    fn membrane_matches_closed_form() {
        // p = 4, r = 0.5, t = 0.01 -> sigma = 4*0.5/(2*0.01) = 100.
        let s = ThinSphere::new(4.0, 0.5, 0.01).unwrap();
        assert!((s.membrane_stress() - 100.0).abs() < EPS);
    }

    #[test]
    fn sphere_is_half_cylinder_hoop() {
        // The defining cross-check: for the SAME p, r, t the sphere
        // membrane stress is half the cylinder hoop stress.
        for &(p, r, t) in &[(2.0, 0.5, 0.01), (9.0, 1.3, 0.04), (0.5, 7.0, 0.2)] {
            let sphere = ThinSphere::new(p, r, t).unwrap();
            let cyl = ThinCylinder::new(p, r, t).unwrap();
            assert!(
                (sphere.membrane_stress() - 0.5 * cyl.hoop_stress()).abs() < EPS,
                "sphere != half cylinder hoop for p={p}, r={r}, t={t}"
            );
            // ...and equal to that cylinder's longitudinal stress.
            assert!((sphere.membrane_stress() - cyl.longitudinal_stress()).abs() < EPS);
        }
    }

    #[test]
    fn stress_scales_with_pressure_and_radius_over_thickness() {
        let a = ThinSphere::new(3.0, 1.0, 0.05).unwrap(); // r/t = 20
        let b = ThinSphere::new(6.0, 1.0, 0.05).unwrap(); // 2x pressure
        assert!((b.membrane_stress() - 2.0 * a.membrane_stress()).abs() < EPS);

        // Same r/t, same p -> same stress (scale-invariant in r/t).
        let c = ThinSphere::new(3.0, 2.0, 0.10).unwrap(); // r/t = 20
        assert!((c.membrane_stress() - a.membrane_stress()).abs() < EPS);
    }

    #[test]
    fn shear_is_quarter_pr_over_t() {
        let s = ThinSphere::new(4.0, 0.5, 0.01).unwrap();
        let expected = s.pressure() * s.radius() / (4.0 * s.thickness());
        assert!((s.max_shear_stress() - expected).abs() < EPS);
    }

    #[test]
    fn von_mises_equals_membrane_stress() {
        let s = ThinSphere::new(4.0, 0.5, 0.01).unwrap();
        assert!((s.von_mises_stress() - s.membrane_stress()).abs() < EPS);
    }
}
