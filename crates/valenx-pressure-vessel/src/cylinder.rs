//! Thin-wall (membrane) cylinder under internal pressure.
//!
//! ## Model
//!
//! A closed cylindrical vessel of internal radius `r` and wall thickness
//! `t`, carrying gauge internal pressure `p`, is treated as a *membrane*:
//! when `t ≪ r` the through-thickness stress gradient is negligible and
//! the wall carries pressure as two near-uniform direct stresses (see
//! Timoshenko, *Theory of Plates and Shells*; Roark, ch. 13).
//!
//! Force balance on a longitudinal cut (the wall resisting the bursting
//! force across a diametral plane) gives the **hoop** (circumferential)
//! stress
//!
//! ```text
//! sigma_h = p * r / t
//! ```
//!
//! Force balance on a transverse cut (the wall resisting the end-cap
//! thrust `p * pi * r^2` over the ring area `2 * pi * r * t`) gives the
//! **longitudinal** (axial) stress
//!
//! ```text
//! sigma_l = p * r / (2 t)
//! ```
//!
//! so the hoop stress is exactly twice the longitudinal — the reason
//! pressurised pipes split along their length, not across. The radial
//! stress varies from `-p` at the bore to `0` at the free surface; both
//! are negligible beside `sigma_h` for a thin wall and the membrane model
//! takes the third principal stress as zero.
//!
//! ## Honest scope
//!
//! Membrane theory only — valid roughly for `r / t > 10` (equivalently
//! `t / r < 0.1`). It ignores the through-wall radial stress, bending and
//! discontinuity stresses at end caps / nozzles / supports, residual and
//! thermal stress, and every code allowable. For thick walls use the
//! exact [`crate::thick`] Lame solution, which this model is the
//! thin-wall limit of.

use crate::error::{Result, VesselError};
use serde::{Deserialize, Serialize};

/// A thin-wall cylindrical pressure vessel and its membrane stresses.
///
/// Construct with [`ThinCylinder::new`], which validates the inputs and
/// pre-computes the stress field. All stress accessors return values in
/// the same pressure unit as the supplied `pressure` (the formulae are
/// dimensionless in `r / t`).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ThinCylinder {
    pressure: f64,
    radius: f64,
    thickness: f64,
}

impl ThinCylinder {
    /// Build a thin-wall cylinder from internal `pressure`, internal
    /// `radius` and wall `thickness` (consistent units; pressure may be
    /// any positive unit and the stresses come out in that same unit).
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

    /// Radius-to-thickness ratio `r / t`. Membrane theory is accurate for
    /// roughly `r / t > 10`; see [`Self::is_thin_wall`].
    pub fn radius_to_thickness(&self) -> f64 {
        self.radius / self.thickness
    }

    /// Whether the thin-wall assumption is reasonable, i.e.
    /// `r / t >= 10`. Purely advisory — the formulae still evaluate
    /// below this, with growing error versus the exact [`crate::thick`]
    /// solution.
    pub fn is_thin_wall(&self) -> bool {
        self.radius_to_thickness() >= 10.0
    }

    /// Hoop (circumferential) stress `sigma_h = p r / t`. The largest of
    /// the membrane stresses; positive (tensile) for internal pressure.
    pub fn hoop_stress(&self) -> f64 {
        self.pressure * self.radius / self.thickness
    }

    /// Longitudinal (axial) stress `sigma_l = p r / (2 t)` for a closed
    /// cylinder with end caps. Exactly half the hoop stress.
    pub fn longitudinal_stress(&self) -> f64 {
        self.pressure * self.radius / (2.0 * self.thickness)
    }

    /// Maximum **in-plane** shear stress, from the hoop and longitudinal
    /// principal stresses only:
    /// `(sigma_h - sigma_l) / 2 = p r / (4 t)`.
    ///
    /// This is the value a 2-D (plane-stress) Mohr's-circle reading of
    /// the surface gives. The true 3-D maximum shear additionally
    /// involves the (near-zero) radial stress — use
    /// [`Self::max_shear_stress`] for that.
    pub fn max_in_plane_shear(&self) -> f64 {
        (self.hoop_stress() - self.longitudinal_stress()) / 2.0
    }

    /// Maximum (3-D, absolute) shear stress. With principal stresses
    /// `sigma_h > sigma_l > sigma_r ~= 0`, the governing pair is
    /// `(sigma_h - sigma_r) / 2 ~= sigma_h / 2 = p r / (2 t)`.
    ///
    /// Taking the outer-surface radial stress as zero (the membrane
    /// idealisation), this equals half the hoop stress and exceeds the
    /// in-plane value of [`Self::max_in_plane_shear`].
    pub fn max_shear_stress(&self) -> f64 {
        self.hoop_stress() / 2.0
    }

    /// Von Mises equivalent stress for the biaxial membrane state
    /// (`sigma_r = 0`):
    /// `sqrt(sigma_h^2 - sigma_h sigma_l + sigma_l^2)`.
    ///
    /// With `sigma_l = sigma_h / 2` this collapses to
    /// `sigma_h * sqrt(3) / 2`, the standard equivalent stress used to
    /// compare against a uniaxial yield strength.
    pub fn von_mises_stress(&self) -> f64 {
        let sh = self.hoop_stress();
        let sl = self.longitudinal_stress();
        (sh * sh - sh * sl + sl * sl).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Loose absolute tolerance for f64 stress comparisons.
    const EPS: f64 = 1e-9;

    #[test]
    fn rejects_non_positive_inputs() {
        assert!(ThinCylinder::new(0.0, 1.0, 0.1).is_err());
        assert!(ThinCylinder::new(1.0, -1.0, 0.1).is_err());
        assert!(ThinCylinder::new(1.0, 1.0, 0.0).is_err());
        assert!(ThinCylinder::new(f64::NAN, 1.0, 0.1).is_err());
    }

    #[test]
    fn hoop_matches_closed_form() {
        // p = 2, r = 0.5, t = 0.01 -> sigma_h = 2*0.5/0.01 = 100.
        let c = ThinCylinder::new(2.0, 0.5, 0.01).unwrap();
        assert!((c.hoop_stress() - 100.0).abs() < EPS);
        assert!((c.longitudinal_stress() - 50.0).abs() < EPS);
    }

    #[test]
    fn hoop_is_exactly_twice_longitudinal() {
        // The defining validation: sigma_h = 2 * sigma_l for ANY inputs.
        for &(p, r, t) in &[(1.0, 1.0, 0.05), (7.3, 0.123, 0.004), (50.0, 9.0, 0.5)] {
            let c = ThinCylinder::new(p, r, t).unwrap();
            assert!(
                (c.hoop_stress() - 2.0 * c.longitudinal_stress()).abs() < EPS,
                "hoop != 2*long for p={p}, r={r}, t={t}"
            );
        }
    }

    #[test]
    fn stress_scales_linearly_with_pressure() {
        // Doubling p doubles every stress.
        let a = ThinCylinder::new(3.0, 0.4, 0.02).unwrap();
        let b = ThinCylinder::new(6.0, 0.4, 0.02).unwrap();
        assert!((b.hoop_stress() - 2.0 * a.hoop_stress()).abs() < EPS);
        assert!((b.longitudinal_stress() - 2.0 * a.longitudinal_stress()).abs() < EPS);
    }

    #[test]
    fn stress_scales_with_radius_over_thickness() {
        // Hoop depends only on p and r/t: same r/t -> same hoop stress.
        let a = ThinCylinder::new(5.0, 1.0, 0.05).unwrap(); // r/t = 20
        let b = ThinCylinder::new(5.0, 2.0, 0.10).unwrap(); // r/t = 20
        assert!((a.radius_to_thickness() - b.radius_to_thickness()).abs() < EPS);
        assert!((a.hoop_stress() - b.hoop_stress()).abs() < EPS);

        // Halving the thickness (doubling r/t) doubles the hoop stress.
        let c = ThinCylinder::new(5.0, 1.0, 0.025).unwrap(); // r/t = 40
        assert!((c.hoop_stress() - 2.0 * a.hoop_stress()).abs() < EPS);
    }

    #[test]
    fn shear_relationships_hold() {
        let c = ThinCylinder::new(2.0, 0.5, 0.01).unwrap();
        // In-plane max shear = (sigma_h - sigma_l)/2 = p r /(4 t).
        let expected_in_plane = c.pressure() * c.radius() / (4.0 * c.thickness());
        assert!((c.max_in_plane_shear() - expected_in_plane).abs() < EPS);
        // 3-D max shear = sigma_h / 2, and exceeds the in-plane value.
        assert!((c.max_shear_stress() - c.hoop_stress() / 2.0).abs() < EPS);
        assert!(c.max_shear_stress() > c.max_in_plane_shear());
    }

    #[test]
    fn von_mises_collapses_to_root3_over_2_times_hoop() {
        // With sigma_l = sigma_h/2, von Mises = sigma_h * sqrt(3)/2.
        let c = ThinCylinder::new(4.0, 0.3, 0.006).unwrap();
        let expected = c.hoop_stress() * (3.0f64).sqrt() / 2.0;
        assert!((c.von_mises_stress() - expected).abs() < 1e-6);
        // Equivalent stress sits between the longitudinal and hoop values.
        assert!(c.von_mises_stress() > c.longitudinal_stress());
        assert!(c.von_mises_stress() < c.hoop_stress());
    }

    #[test]
    fn thin_wall_flag() {
        assert!(ThinCylinder::new(1.0, 1.0, 0.05).unwrap().is_thin_wall()); // r/t=20
        assert!(!ThinCylinder::new(1.0, 1.0, 0.2).unwrap().is_thin_wall()); // r/t=5
    }
}
