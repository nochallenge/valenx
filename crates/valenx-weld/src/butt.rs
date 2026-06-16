//! Complete-joint-penetration (full-penetration) butt-weld strength.
//!
//! A complete-joint-penetration (CJP) groove — or "butt" — weld fuses
//! the full thickness of the joined plates, so the weld is treated as
//! continuous with the base metal: the load-carrying section is the
//! plate thickness times the weld length, and a tensile or compressive
//! force normal to the weld produces a uniform *normal* (direct) stress
//!
//! ```text
//! sigma = F / (thickness * length).
//! ```
//!
//! The effective area uses the thickness of the *thinner* connected part
//! (the thinner plate governs the throat of a CJP weld); pass that
//! thickness as `thickness`.
//!
//! ## Honest scope
//!
//! The same first-pass simplifications apply as for fillets: full
//! penetration (effective throat equals the plate thickness), uniform
//! stress over the section, force normal to the weld, the base metal not
//! governing, and no fatigue, residual stress or eccentricity. Partial
//! penetration grooves (which have a reduced effective throat) and shear
//! on the groove are out of scope for this module. Use a governing code
//! and a qualified engineer for any real joint.

use crate::error::{require_positive, Result};
use serde::{Deserialize, Serialize};

/// A complete-joint-penetration butt weld of a given (thinner-part)
/// thickness and length.
///
/// Construct with [`ButtWeld::new`], which validates both inputs. The
/// methods return the load-carrying [`area`](ButtWeld::area), the direct
/// [`normal_stress`](ButtWeld::normal_stress) under a force, and the
/// allowable [`capacity`](ButtWeld::capacity) for an allowable normal
/// stress.
///
/// Units are caller-defined but must be consistent (e.g. thickness and
/// length in millimetres, force in newtons, stress in `N/mm^2 = MPa`).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ButtWeld {
    /// Effective thickness — the thickness of the thinner connected part
    /// (length units).
    thickness: f64,
    /// Effective length of the weld (length units).
    length: f64,
}

impl ButtWeld {
    /// Build a CJP butt weld of effective `thickness` and `length`.
    ///
    /// # Errors
    ///
    /// Returns [`WeldError::BadParameter`](crate::error::WeldError::BadParameter) when either `thickness` or
    /// `length` is not a strictly positive finite number.
    pub fn new(thickness: f64, length: f64) -> Result<Self> {
        let thickness = require_positive("thickness", thickness)?;
        let length = require_positive("length", length)?;
        Ok(Self { thickness, length })
    }

    /// Effective thickness of the weld (length units).
    pub fn thickness(&self) -> f64 {
        self.thickness
    }

    /// Effective length of the weld (length units).
    pub fn length(&self) -> f64 {
        self.length
    }

    /// Load-carrying area `thickness * length` (area units).
    pub fn area(&self) -> f64 {
        self.thickness * self.length
    }

    /// Direct normal stress `sigma = F / (thickness * length)` under an
    /// applied normal force `force`.
    ///
    /// Units follow the inputs: a force in newtons over an area in
    /// `mm^2` yields a stress in `N/mm^2 = MPa`.
    ///
    /// # Errors
    ///
    /// Returns [`WeldError::BadParameter`](crate::error::WeldError::BadParameter) when `force` is not a strictly
    /// positive finite number.
    pub fn normal_stress(&self, force: f64) -> Result<f64> {
        let force = require_positive("force", force)?;
        Ok(force / self.area())
    }

    /// Allowable load `F_allow = sigma_allow * thickness * length` for an
    /// allowable normal stress `allowable_normal_stress`.
    ///
    /// The inverse of [`normal_stress`](Self::normal_stress).
    ///
    /// # Errors
    ///
    /// Returns [`WeldError::BadParameter`](crate::error::WeldError::BadParameter) when
    /// `allowable_normal_stress` is not a strictly positive finite
    /// number.
    pub fn capacity(&self, allowable_normal_stress: f64) -> Result<f64> {
        let sigma = require_positive("allowable_normal_stress", allowable_normal_stress)?;
        Ok(sigma * self.area())
    }
}

/// Free-function form of the direct butt-weld normal stress
/// `sigma = F / (thickness * length)`.
///
/// A thin convenience over [`ButtWeld::normal_stress`].
///
/// # Errors
///
/// Returns [`WeldError::BadParameter`](crate::error::WeldError::BadParameter) when any of `force`, `thickness`
/// or `length` is not a strictly positive finite number.
pub fn normal_stress(force: f64, thickness: f64, length: f64) -> Result<f64> {
    ButtWeld::new(thickness, length)?.normal_stress(force)
}

/// Free-function form of the allowable butt-weld load
/// `F_allow = sigma_allow * thickness * length`.
///
/// A thin convenience over [`ButtWeld::capacity`].
///
/// # Errors
///
/// Returns [`WeldError::BadParameter`](crate::error::WeldError::BadParameter) when any of
/// `allowable_normal_stress`, `thickness` or `length` is not a strictly
/// positive finite number.
pub fn capacity(allowable_normal_stress: f64, thickness: f64, length: f64) -> Result<f64> {
    ButtWeld::new(thickness, length)?.capacity(allowable_normal_stress)
}

/// Weld **length required** to carry a normal `force` at an allowable
/// normal stress, for a CJP butt weld of effective `thickness`:
///
/// ```text
/// length = force / (sigma_allow * thickness)
/// ```
///
/// This is [`capacity`] solved for the length — the weld-sizing inverse.
/// A [`ButtWeld`] built with this length carries exactly `force` at the
/// allowable stress (its [`normal_stress`](ButtWeld::normal_stress) then
/// equals `sigma_allow`).
///
/// # Errors
///
/// Returns [`WeldError::BadParameter`](crate::error::WeldError::BadParameter)
/// when any of `force`, `thickness` or `allowable_normal_stress` is not a
/// strictly positive finite number.
pub fn required_length(force: f64, thickness: f64, allowable_normal_stress: f64) -> Result<f64> {
    let force = require_positive("force", force)?;
    let thickness = require_positive("thickness", thickness)?;
    let sigma = require_positive("allowable_normal_stress", allowable_normal_stress)?;
    Ok(force / (sigma * thickness))
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // area = thickness * length.
    #[test]
    fn area_is_thickness_times_length() {
        let w = ButtWeld::new(12.0, 200.0).expect("valid");
        assert!((w.area() - 2400.0).abs() < EPS, "got: {}", w.area());
    }

    // butt stress = F / (t * L). Ground-truth hand calculation:
    // t 12, L 200 -> area 2400 mm^2; F = 360_000 N -> sigma = 150 MPa.
    #[test]
    fn normal_stress_equals_force_over_area() {
        let w = ButtWeld::new(12.0, 200.0).expect("valid");
        let sigma = w.normal_stress(360_000.0).expect("positive force");
        assert!((sigma - 150.0).abs() < 1e-6, "got: {sigma}");

        let expected = 360_000.0 / (12.0 * 200.0);
        assert!((sigma - expected).abs() < EPS, "got: {sigma}");
    }

    // Free function agrees with the method form.
    #[test]
    fn free_normal_stress_matches_method() {
        let a = normal_stress(120_000.0, 10.0, 150.0).expect("valid");
        let b = ButtWeld::new(10.0, 150.0)
            .unwrap()
            .normal_stress(120_000.0)
            .unwrap();
        assert!((a - b).abs() < EPS, "free {a} vs method {b}");
    }

    // capacity = allowable * t * L, and it is the exact inverse of
    // normal_stress.
    #[test]
    fn capacity_is_inverse_of_normal_stress() {
        let w = ButtWeld::new(16.0, 250.0).expect("valid");
        let allow = 160.0; // MPa
        let cap = w.capacity(allow).expect("positive allowable");
        let expected = allow * 16.0 * 250.0;
        assert!((cap - expected).abs() < 1e-6, "got: {cap}");

        let sigma_at_cap = w.normal_stress(cap).expect("positive");
        assert!((sigma_at_cap - allow).abs() < 1e-9, "got: {sigma_at_cap}");
    }

    // Capacity scales linearly with both length and thickness.
    #[test]
    fn capacity_scales_with_length_and_thickness() {
        let base = ButtWeld::new(10.0, 100.0).unwrap().capacity(150.0).unwrap();

        let long = ButtWeld::new(10.0, 200.0).unwrap().capacity(150.0).unwrap();
        assert!((long - 2.0 * base).abs() < 1e-6, "base {base} long {long}");

        let thick = ButtWeld::new(20.0, 100.0).unwrap().capacity(150.0).unwrap();
        assert!(
            (thick - 2.0 * base).abs() < 1e-6,
            "base {base} thick {thick}"
        );
    }

    #[test]
    fn rejects_non_positive_inputs() {
        assert!(ButtWeld::new(0.0, 10.0).is_err());
        assert!(ButtWeld::new(10.0, 0.0).is_err());
        assert!(ButtWeld::new(-2.0, 10.0).is_err());

        let w = ButtWeld::new(10.0, 100.0).unwrap();
        assert!(w.normal_stress(-1.0).is_err());
        assert!(w.capacity(f64::INFINITY).is_err());
    }

    // required_length = F / (sigma_allow * thickness). Hand calc:
    // F=360000, t=12, sigma=150 -> 360000/(150*12) = 200 mm.
    #[test]
    fn required_length_matches_hand_calc() {
        let l = required_length(360_000.0, 12.0, 150.0).expect("valid");
        assert!((l - 200.0).abs() < 1e-9, "got: {l}");
    }

    // Inverse of capacity / normal_stress: a weld of that length carries
    // exactly the force at the allowable stress.
    #[test]
    fn required_length_round_trips_capacity_and_stress() {
        let (force, t, allow) = (300_000.0, 16.0, 160.0);
        let l = required_length(force, t, allow).unwrap();
        let w = ButtWeld::new(t, l).unwrap();
        assert!((w.normal_stress(force).unwrap() - allow).abs() < 1e-9);
        assert!((w.capacity(allow).unwrap() - force).abs() < 1e-6);
        assert!((capacity(allow, t, l).unwrap() - force).abs() < 1e-6);
    }

    // Sizing scales with force and inversely with allowable stress.
    #[test]
    fn required_length_scaling() {
        let base = required_length(100_000.0, 10.0, 150.0).unwrap();
        let double_f = required_length(200_000.0, 10.0, 150.0).unwrap();
        let half_allow = required_length(100_000.0, 10.0, 75.0).unwrap();
        assert!((double_f - 2.0 * base).abs() < 1e-9);
        assert!((half_allow - 2.0 * base).abs() < 1e-9);
    }

    #[test]
    fn required_length_rejects_bad_inputs() {
        assert!(required_length(0.0, 10.0, 150.0).is_err());
        assert!(required_length(100_000.0, 10.0, -150.0).is_err());
        assert!(required_length(f64::INFINITY, 10.0, 150.0).is_err());
    }
}
