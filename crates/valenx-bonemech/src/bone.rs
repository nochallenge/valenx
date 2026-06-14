//! The [`Bone`] material/section model and its axial-loading mechanics.
//!
//! A [`Bone`] bundles one cross-section's geometry and one bone-tissue
//! material into a small validated struct, then exposes the textbook
//! relations for **uniaxial** (pure tension/compression) loading:
//!
//! - axial stress from an applied force,
//! - linear-elastic (Hooke's-law) strain from a stress,
//! - the ultimate axial force the section carries before fracture, and
//! - a utilisation / safety-factor read-out against the ultimate force.
//!
//! ## Units
//!
//! | Quantity                | Symbol      | Unit  |
//! |-------------------------|-------------|-------|
//! | Cross-sectional area    | `A`         | mm^2  |
//! | Axial force             | `F`         | N     |
//! | Axial / ultimate stress | `sigma`     | MPa   |
//! | Elastic (Young's) mod.  | `E`         | GPa   |
//! | Strain                  | `epsilon`   | —     |
//!
//! Because `1 MPa = 1 N/mm^2`, a force in newtons divided by an area in
//! mm^2 lands directly in MPa with no conversion factor. Hooke's law
//! divides a stress (MPa) by a modulus (GPa); since `1 GPa = 1000 MPa`
//! the modulus is multiplied by `1000.0` to put both on the MPa scale,
//! leaving strain dimensionless.
//!
//! ## Model
//!
//! Bone tissue is treated as a **homogeneous, isotropic, linear-elastic**
//! solid up to a single ultimate stress — a deliberate simplification.
//! Real cortical and trabecular bone are anisotropic, viscoelastic,
//! tension/compression-asymmetric, and fail by progressive
//! micro-damage; none of that is modelled here. The defaults
//! ([`Bone::cortical`]) use representative textbook cortical-bone values:
//! a Young's modulus of `17.0 GPa` and an ultimate stress of
//! `150.0 MPa`, both round figures within the commonly cited ranges
//! (cortical `E ~ 15–20 GPa`; ultimate strength order `100–200 MPa`,
//! load-mode dependent).

use crate::error::{BoneError, Result};
use serde::{Deserialize, Serialize};

/// Representative cortical-bone Young's modulus, in GPa.
///
/// A round figure inside the commonly cited `15–20 GPa` range for human
/// cortical bone loaded along the osteon (long) axis.
pub const CORTICAL_MODULUS_GPA: f64 = 17.0;

/// Representative cortical-bone ultimate stress, in MPa.
///
/// A round figure on the order of cortical bone's ultimate strength
/// (`~100–200 MPa`, depending on load mode and orientation).
pub const CORTICAL_ULTIMATE_STRESS_MPA: f64 = 150.0;

/// A bone cross-section paired with its tissue material.
///
/// Construct with [`Bone::new`] (explicit material) or [`Bone::cortical`]
/// (cortical-bone defaults), both of which validate their inputs. The
/// fields are public for read access and `serde`-(de)serialisable, but
/// prefer the constructors so the invariants (positive area, positive
/// finite modulus and ultimate stress) are guaranteed to hold.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Bone {
    /// Cross-sectional area carrying the axial load, in mm^2. Strictly
    /// positive.
    pub cross_sectional_area_mm2: f64,
    /// Young's (elastic) modulus of the tissue, in GPa. Strictly
    /// positive and finite.
    pub elastic_modulus_gpa: f64,
    /// Ultimate (failure) stress of the tissue, in MPa. Strictly
    /// positive and finite.
    pub ultimate_stress_mpa: f64,
}

impl Bone {
    /// Build a [`Bone`] from an explicit cross-sectional area, elastic
    /// modulus and ultimate stress.
    ///
    /// # Errors
    ///
    /// Returns [`BoneError::Invalid`] if `area_mm2 <= 0`, if
    /// `elastic_modulus_gpa <= 0`, if `ultimate_stress_mpa <= 0`, or if
    /// any argument is non-finite.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_bonemech::Bone;
    /// let b = Bone::new(300.0, 17.0, 150.0).unwrap();
    /// assert!((b.cross_sectional_area_mm2 - 300.0).abs() < 1e-12);
    /// ```
    pub fn new(area_mm2: f64, elastic_modulus_gpa: f64, ultimate_stress_mpa: f64) -> Result<Self> {
        if !area_mm2.is_finite() || area_mm2 <= 0.0 {
            return Err(BoneError::invalid(
                "area_mm2",
                "cross-sectional area must be positive and finite",
            ));
        }
        if !elastic_modulus_gpa.is_finite() || elastic_modulus_gpa <= 0.0 {
            return Err(BoneError::invalid(
                "elastic_modulus_gpa",
                "elastic modulus must be positive and finite",
            ));
        }
        if !ultimate_stress_mpa.is_finite() || ultimate_stress_mpa <= 0.0 {
            return Err(BoneError::invalid(
                "ultimate_stress_mpa",
                "ultimate stress must be positive and finite",
            ));
        }
        Ok(Bone {
            cross_sectional_area_mm2: area_mm2,
            elastic_modulus_gpa,
            ultimate_stress_mpa,
        })
    }

    /// Build a [`Bone`] of a given cross-sectional area using the
    /// representative **cortical-bone** defaults
    /// ([`CORTICAL_MODULUS_GPA`] = 17.0 GPa,
    /// [`CORTICAL_ULTIMATE_STRESS_MPA`] = 150.0 MPa).
    ///
    /// # Errors
    ///
    /// Returns [`BoneError::Invalid`] if `area_mm2 <= 0` or is
    /// non-finite.
    pub fn cortical(area_mm2: f64) -> Result<Self> {
        Bone::new(area_mm2, CORTICAL_MODULUS_GPA, CORTICAL_ULTIMATE_STRESS_MPA)
    }

    /// Axial (normal) stress produced by an applied axial force, in MPa.
    ///
    /// Computes `sigma = force_n / area_mm2`. Because `1 MPa = 1 N/mm^2`
    /// the newtons-over-mm^2 quotient is already in MPa. The force may be
    /// signed: a positive force is tensile, a negative force compressive,
    /// and the sign carries through to the stress.
    ///
    /// # Errors
    ///
    /// Returns [`BoneError::Invalid`] if `force_n` is non-finite.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_bonemech::Bone;
    /// // 1500 N over 300 mm^2 -> 5 MPa.
    /// let b = Bone::cortical(300.0).unwrap();
    /// let s = b.axial_stress_mpa(1500.0).unwrap();
    /// assert!((s - 5.0).abs() < 1e-12);
    /// ```
    pub fn axial_stress_mpa(&self, force_n: f64) -> Result<f64> {
        if !force_n.is_finite() {
            return Err(BoneError::invalid("force_n", "axial force must be finite"));
        }
        Ok(force_n / self.cross_sectional_area_mm2)
    }

    /// Linear-elastic (Hooke's-law) strain produced by a given axial
    /// stress, dimensionless.
    ///
    /// Computes `epsilon = stress_mpa / (E_gpa * 1000.0)`. The modulus is
    /// multiplied by `1000.0` because `1 GPa = 1000 MPa`, putting stress
    /// and modulus on the same MPa scale so the quotient is a pure
    /// (unitless) strain. Valid only in the elastic regime — this is
    /// `sigma = E * epsilon` rearranged, with no yield or plasticity.
    ///
    /// # Errors
    ///
    /// Returns [`BoneError::Invalid`] if `stress_mpa` is non-finite.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_bonemech::Bone;
    /// // 17 MPa at E = 17 GPa = 17000 MPa -> strain 0.001.
    /// let b = Bone::cortical(300.0).unwrap();
    /// let eps = b.strain(17.0).unwrap();
    /// assert!((eps - 0.001).abs() < 1e-12);
    /// ```
    pub fn strain(&self, stress_mpa: f64) -> Result<f64> {
        if !stress_mpa.is_finite() {
            return Err(BoneError::invalid("stress_mpa", "stress must be finite"));
        }
        Ok(stress_mpa / (self.elastic_modulus_gpa * 1000.0))
    }

    /// The ultimate axial force the section can carry before fracture, in
    /// newtons.
    ///
    /// Computes `F_ult = ultimate_stress_mpa * area_mm2`. With stress in
    /// MPa (= N/mm^2) and area in mm^2 the product is in newtons. This is
    /// simply [`axial_stress_mpa`](Self::axial_stress_mpa) inverted at the
    /// ultimate stress, and is always strictly positive (both factors are
    /// positive by construction).
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_bonemech::Bone;
    /// // 150 MPa * 300 mm^2 = 45 000 N.
    /// let b = Bone::cortical(300.0).unwrap();
    /// assert!((b.fracture_load_n() - 45_000.0).abs() < 1e-9);
    /// ```
    pub fn fracture_load_n(&self) -> f64 {
        self.ultimate_stress_mpa * self.cross_sectional_area_mm2
    }

    /// Utilisation under an axial force: `|F| / F_ult`, dimensionless.
    ///
    /// A value below `1.0` means the section is below its fracture load
    /// (tension or compression, hence the magnitude); `1.0` is exactly at
    /// the ultimate force; above `1.0` the linear model predicts
    /// fracture. The reciprocal is the safety factor — see
    /// [`safety_factor`](Self::safety_factor).
    ///
    /// # Errors
    ///
    /// Returns [`BoneError::Invalid`] if `force_n` is non-finite.
    pub fn utilisation(&self, force_n: f64) -> Result<f64> {
        if !force_n.is_finite() {
            return Err(BoneError::invalid("force_n", "axial force must be finite"));
        }
        Ok(force_n.abs() / self.fracture_load_n())
    }

    /// Safety factor against fracture under an axial force:
    /// `F_ult / |F|`, dimensionless.
    ///
    /// A value above `1.0` is safe in this linear model; `1.0` is exactly
    /// at the fracture load. For a zero force the safety factor is
    /// infinite ([`f64::INFINITY`]).
    ///
    /// # Errors
    ///
    /// Returns [`BoneError::Invalid`] if `force_n` is non-finite.
    pub fn safety_factor(&self, force_n: f64) -> Result<f64> {
        let u = self.utilisation(force_n)?;
        Ok(if u == 0.0 { f64::INFINITY } else { 1.0 / u })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cortical_defaults_are_the_documented_constants() {
        let b = Bone::cortical(300.0).unwrap();
        assert!((b.elastic_modulus_gpa - 17.0).abs() < 1e-12);
        assert!((b.ultimate_stress_mpa - 150.0).abs() < 1e-12);
        assert!((b.cross_sectional_area_mm2 - 300.0).abs() < 1e-12);
    }

    /// sigma = F / A on a worked example: 1500 N over 300 mm^2 = 5 MPa.
    #[test]
    fn axial_stress_is_force_over_area() {
        let b = Bone::cortical(300.0).unwrap();
        let s = b.axial_stress_mpa(1500.0).unwrap();
        assert!((s - 5.0).abs() < 1e-12, "got {s}");
    }

    /// Compression carries the sign through: a negative force gives a
    /// negative (compressive) stress equal and opposite to the tensile
    /// case.
    #[test]
    fn axial_stress_sign_follows_force() {
        let b = Bone::cortical(300.0).unwrap();
        let tens = b.axial_stress_mpa(1500.0).unwrap();
        let comp = b.axial_stress_mpa(-1500.0).unwrap();
        assert!((tens + comp).abs() < 1e-12, "expected opposite signs");
        assert!(comp < 0.0);
    }

    /// Hooke's law: epsilon = sigma / (E * 1000). At sigma = 17 MPa and
    /// E = 17 GPa = 17000 MPa the strain is exactly 0.001.
    #[test]
    fn strain_is_hookes_law() {
        let b = Bone::cortical(300.0).unwrap();
        let eps = b.strain(17.0).unwrap();
        assert!((eps - 0.001).abs() < 1e-12, "got {eps}");
    }

    /// Round-trip: stress -> strain -> back to stress via sigma = E*eps
    /// recovers the original stress. This is the defining identity of the
    /// linear model.
    #[test]
    fn stress_strain_roundtrip() {
        let b = Bone::new(250.0, 18.5, 140.0).unwrap();
        let sigma = 42.0;
        let eps = b.strain(sigma).unwrap();
        let recovered = eps * b.elastic_modulus_gpa * 1000.0;
        assert!((recovered - sigma).abs() < 1e-9, "got {recovered}");
    }

    /// Fracture load = sigma_ult * A: 150 MPa * 300 mm^2 = 45 000 N.
    #[test]
    fn fracture_load_is_ultimate_stress_times_area() {
        let b = Bone::cortical(300.0).unwrap();
        let f = b.fracture_load_n();
        assert!((f - 45_000.0).abs() < 1e-9, "got {f}");
    }

    /// A force exactly at the fracture load gives utilisation 1.0 and
    /// safety factor 1.0; half that gives 0.5 and 2.0 respectively.
    #[test]
    fn utilisation_and_safety_factor_are_reciprocal() {
        let b = Bone::cortical(300.0).unwrap();
        let f_ult = b.fracture_load_n();

        let u_full = b.utilisation(f_ult).unwrap();
        assert!((u_full - 1.0).abs() < 1e-12, "got {u_full}");
        assert!((b.safety_factor(f_ult).unwrap() - 1.0).abs() < 1e-12);

        let u_half = b.utilisation(f_ult / 2.0).unwrap();
        assert!((u_half - 0.5).abs() < 1e-12, "got {u_half}");
        assert!((b.safety_factor(f_ult / 2.0).unwrap() - 2.0).abs() < 1e-12);
    }

    /// Compression utilisation uses the force magnitude, so a compressive
    /// force matches the equal-magnitude tensile force.
    #[test]
    fn utilisation_uses_force_magnitude() {
        let b = Bone::cortical(300.0).unwrap();
        let f_ult = b.fracture_load_n();
        let u = b.utilisation(-f_ult).unwrap();
        assert!((u - 1.0).abs() < 1e-12, "got {u}");
    }

    #[test]
    fn zero_force_has_infinite_safety_factor() {
        let b = Bone::cortical(300.0).unwrap();
        assert!(b.safety_factor(0.0).unwrap().is_infinite());
        assert!((b.utilisation(0.0).unwrap()).abs() < 1e-12);
    }

    #[test]
    fn constructors_reject_bad_inputs() {
        assert!(Bone::new(0.0, 17.0, 150.0).is_err());
        assert!(Bone::new(-1.0, 17.0, 150.0).is_err());
        assert!(Bone::new(300.0, 0.0, 150.0).is_err());
        assert!(Bone::new(300.0, 17.0, 0.0).is_err());
        assert!(Bone::new(f64::NAN, 17.0, 150.0).is_err());
        assert!(Bone::new(300.0, f64::INFINITY, 150.0).is_err());
        assert!(Bone::cortical(-5.0).is_err());
    }

    #[test]
    fn non_finite_loads_rejected() {
        let b = Bone::cortical(300.0).unwrap();
        assert!(b.axial_stress_mpa(f64::NAN).is_err());
        assert!(b.strain(f64::INFINITY).is_err());
        assert!(b.utilisation(f64::NAN).is_err());
    }
}
