//! Fillet-weld static strength.
//!
//! A fillet weld deposits a roughly triangular bead in the corner
//! between two members. For the common equal-leg, 45-degree fillet the
//! load-carrying section is the *effective throat* — the shortest line
//! from the root to the face of the triangle — and design codes treat
//! every fillet weld (transverse or longitudinal) as failing in shear on
//! that throat.
//!
//! ## Model
//!
//! For an equal-leg fillet of leg size `w` the geometric throat is the
//! altitude of an isosceles right triangle:
//!
//! ```text
//! throat = w * cos(45 deg) = w / sqrt(2) = 0.7071067... * w
//! ```
//!
//! The AWS D1.1 / AISC textbook convention rounds the factor to the
//! tabulated constant [`FILLET_THROAT_FACTOR`] = `0.707`, and this crate
//! uses that exact tabulated value so results match a hand calculation
//! done from a code table. The throat carries the force as a uniform
//! average shear stress over the area `throat * length`:
//!
//! ```text
//! tau = F / (throat * length)
//! ```
//!
//! and, inverting, the allowable load for a given allowable shear stress
//! `tau_allow` is
//!
//! ```text
//! F_allow = tau_allow * throat * length.
//! ```
//!
//! ## Honest scope
//!
//! These are the standard simplifying assumptions of a first-pass weld
//! check: equal-leg 45-degree fillet, uniform stress over the throat,
//! concentric in-plane loading, the base metal not governing, and no
//! fatigue, residual stress, eccentricity or weld-group geometry. The
//! directional-strength increase that AISC allows for transversely
//! loaded fillets is deliberately *not* applied — the result is the
//! plain, conservative `F / (throat * length)`. Use a governing code and
//! a qualified engineer for any real joint.

use crate::error::{require_positive, Result};
use serde::{Deserialize, Serialize};

/// Effective-throat factor for an equal-leg, 45-degree fillet weld:
/// `throat = FILLET_THROAT_FACTOR * leg`.
///
/// This is the AWS D1.1 / AISC tabulated value `0.707`, the two-decimal
/// rounding of `cos(45 deg) = 1 / sqrt(2) = 0.7071067...`. Using the
/// tabulated constant (rather than the full-precision `1/sqrt(2)`) keeps
/// results identical to a by-hand calculation taken from a code table.
pub const FILLET_THROAT_FACTOR: f64 = 0.707;

/// Effective throat of an equal-leg fillet weld of leg size `leg`.
///
/// Computes `FILLET_THROAT_FACTOR * leg`. Units follow the input: a leg
/// in millimetres yields a throat in millimetres.
///
/// # Errors
///
/// Returns [`WeldError::BadParameter`](crate::error::WeldError::BadParameter) when `leg` is not a strictly
/// positive finite number.
///
/// # Examples
///
/// ```
/// use valenx_weld::fillet::throat;
/// let t = throat(10.0).unwrap();
/// assert!((t - 7.07).abs() < 1e-9);
/// ```
pub fn throat(leg: f64) -> Result<f64> {
    let leg = require_positive("leg", leg)?;
    Ok(FILLET_THROAT_FACTOR * leg)
}

/// An equal-leg fillet weld of a given leg size and length.
///
/// Construct with [`FilletWeld::new`], which validates both inputs. The
/// methods then return the [`throat`](FilletWeld::throat),
/// [`throat_area`](FilletWeld::throat_area), an average
/// [`shear_stress`](FilletWeld::shear_stress) under a force, and the
/// allowable [`capacity`](FilletWeld::capacity) for an allowable shear
/// stress.
///
/// Units are caller-defined but must be consistent (e.g. leg and length
/// in millimetres, force in newtons, stress in `N/mm^2 = MPa`).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FilletWeld {
    /// Leg size `w` of the equal-leg fillet (length units).
    leg: f64,
    /// Effective length of the weld run (length units).
    length: f64,
}

impl FilletWeld {
    /// Build a fillet weld of `leg` size and effective `length`.
    ///
    /// # Errors
    ///
    /// Returns [`WeldError::BadParameter`](crate::error::WeldError::BadParameter) when either `leg` or `length`
    /// is not a strictly positive finite number.
    pub fn new(leg: f64, length: f64) -> Result<Self> {
        let leg = require_positive("leg", leg)?;
        let length = require_positive("length", length)?;
        Ok(Self { leg, length })
    }

    /// Leg size `w` of the fillet (length units).
    pub fn leg(&self) -> f64 {
        self.leg
    }

    /// Effective length of the weld run (length units).
    pub fn length(&self) -> f64 {
        self.length
    }

    /// Effective throat `0.707 * leg` (length units).
    pub fn throat(&self) -> f64 {
        FILLET_THROAT_FACTOR * self.leg
    }

    /// Throat (shear) area `throat * length` (area units).
    pub fn throat_area(&self) -> f64 {
        self.throat() * self.length
    }

    /// Average shear stress `tau = F / (throat * length)` carried by the
    /// throat under applied force `force`.
    ///
    /// Units follow the inputs: a force in newtons over an area in
    /// `mm^2` yields a stress in `N/mm^2 = MPa`.
    ///
    /// # Errors
    ///
    /// Returns [`WeldError::BadParameter`](crate::error::WeldError::BadParameter) when `force` is not a strictly
    /// positive finite number.
    pub fn shear_stress(&self, force: f64) -> Result<f64> {
        let force = require_positive("force", force)?;
        Ok(force / self.throat_area())
    }

    /// Allowable load `F_allow = tau_allow * throat * length` for an
    /// allowable shear stress `allowable_shear_stress`.
    ///
    /// This is the inverse of [`shear_stress`](Self::shear_stress): the
    /// largest concentric force the throat can carry before its average
    /// shear stress reaches the allowable.
    ///
    /// # Errors
    ///
    /// Returns [`WeldError::BadParameter`](crate::error::WeldError::BadParameter) when
    /// `allowable_shear_stress` is not a strictly positive finite
    /// number.
    pub fn capacity(&self, allowable_shear_stress: f64) -> Result<f64> {
        let tau = require_positive("allowable_shear_stress", allowable_shear_stress)?;
        Ok(tau * self.throat_area())
    }

    /// Utilisation `tau / tau_allow` under `force` against an allowable
    /// shear stress.
    ///
    /// A value at or below `1.0` means the weld passes this first-pass
    /// check; above `1.0` it is overstressed. Equivalent to
    /// `force / capacity(allowable_shear_stress)`.
    ///
    /// # Errors
    ///
    /// Returns [`WeldError::BadParameter`](crate::error::WeldError::BadParameter) when `force` or
    /// `allowable_shear_stress` is not a strictly positive finite
    /// number.
    pub fn utilisation(&self, force: f64, allowable_shear_stress: f64) -> Result<f64> {
        let stress = self.shear_stress(force)?;
        let allow = require_positive("allowable_shear_stress", allowable_shear_stress)?;
        Ok(stress / allow)
    }
}

/// Free-function form of the average fillet throat shear stress
/// `tau = F / (0.707 * leg * length)`.
///
/// A thin convenience over [`FilletWeld::shear_stress`] for callers that
/// do not want to hold a [`FilletWeld`] value.
///
/// # Errors
///
/// Returns [`WeldError::BadParameter`](crate::error::WeldError::BadParameter) when any of `force`, `leg` or
/// `length` is not a strictly positive finite number.
pub fn shear_stress(force: f64, leg: f64, length: f64) -> Result<f64> {
    FilletWeld::new(leg, length)?.shear_stress(force)
}

/// Free-function form of the allowable fillet load
/// `F_allow = tau_allow * 0.707 * leg * length`.
///
/// A thin convenience over [`FilletWeld::capacity`].
///
/// # Errors
///
/// Returns [`WeldError::BadParameter`](crate::error::WeldError::BadParameter) when any of
/// `allowable_shear_stress`, `leg` or `length` is not a strictly
/// positive finite number.
pub fn capacity(allowable_shear_stress: f64, leg: f64, length: f64) -> Result<f64> {
    FilletWeld::new(leg, length)?.capacity(allowable_shear_stress)
}

/// Weld run **length required** to carry a concentric `force` at an
/// allowable shear stress, for an equal-leg fillet of leg size `leg`:
///
/// ```text
/// length = force / (tau_allow * throat) = force / (tau_allow * 0.707 * leg)
/// ```
///
/// This is [`capacity`] solved for the length — the weld-sizing inverse
/// that closes the design loop. A [`FilletWeld`] built with this length
/// carries exactly `force` at the allowable stress: its
/// [`shear_stress`](FilletWeld::shear_stress) then equals `tau_allow` and
/// its [`utilisation`](FilletWeld::utilisation) is `1.0`.
///
/// # Errors
///
/// Returns [`WeldError::BadParameter`](crate::error::WeldError::BadParameter)
/// when any of `force`, `leg` or `allowable_shear_stress` is not a
/// strictly positive finite number.
pub fn required_length(force: f64, leg: f64, allowable_shear_stress: f64) -> Result<f64> {
    let force = require_positive("force", force)?;
    let leg = require_positive("leg", leg)?;
    let tau = require_positive("allowable_shear_stress", allowable_shear_stress)?;
    Ok(force / (tau * FILLET_THROAT_FACTOR * leg))
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // throat = 0.707 * leg, exactly the tabulated factor.
    #[test]
    fn throat_is_factor_times_leg() {
        let t = throat(10.0).expect("10 mm leg");
        assert!((t - 7.07).abs() < EPS, "got: {t}");

        let t = throat(8.0).expect("8 mm leg");
        assert!((t - 0.707 * 8.0).abs() < EPS, "got: {t}");
    }

    // The tabulated 0.707 is within 1.5e-4 of the exact 1/sqrt(2) throat.
    #[test]
    fn tabulated_factor_tracks_exact_geometry() {
        let leg = 12.0_f64;
        let exact = leg / 2.0_f64.sqrt();
        let tabulated = throat(leg).expect("12 mm leg");
        assert!(
            (tabulated - exact).abs() < 2.0e-3,
            "tabulated {tabulated} vs exact {exact}"
        );
    }

    // throat_area = throat * length = 0.707 * leg * length.
    #[test]
    fn throat_area_matches_hand_calc() {
        let w = FilletWeld::new(6.0, 100.0).expect("valid");
        let expected = 0.707 * 6.0 * 100.0; // 424.2 mm^2
        assert!(
            (w.throat_area() - expected).abs() < EPS,
            "got: {}",
            w.throat_area()
        );
    }

    // shear = F / throat-area. Ground-truth hand calculation:
    // leg 6, length 100 -> throat 4.242, area 424.2; F = 50_000 N
    // -> tau = 50000 / 424.2 = 117.8689297... MPa.
    #[test]
    fn shear_stress_equals_force_over_throat_area() {
        let w = FilletWeld::new(6.0, 100.0).expect("valid");
        let tau = w.shear_stress(50_000.0).expect("positive force");
        let expected = 50_000.0 / (0.707 * 6.0 * 100.0);
        assert!((tau - expected).abs() < EPS, "got: {tau}");
        assert!((tau - 117.868_929_75).abs() < 1e-6, "got: {tau}");
    }

    // Free function agrees with the method form.
    #[test]
    fn free_shear_stress_matches_method() {
        let a = shear_stress(30_000.0, 5.0, 80.0).expect("valid");
        let b = FilletWeld::new(5.0, 80.0)
            .unwrap()
            .shear_stress(30_000.0)
            .unwrap();
        assert!((a - b).abs() < EPS, "free {a} vs method {b}");
    }

    // capacity = allowable * throat * length, and it is the exact
    // inverse of shear_stress: stressing a weld at its capacity gives
    // tau == allowable.
    #[test]
    fn capacity_is_inverse_of_shear_stress() {
        let w = FilletWeld::new(8.0, 150.0).expect("valid");
        let allow = 95.0; // MPa
        let cap = w.capacity(allow).expect("positive allowable");
        let expected = allow * 0.707 * 8.0 * 150.0;
        assert!((cap - expected).abs() < 1e-6, "got: {cap}");

        let tau_at_cap = w.shear_stress(cap).expect("positive");
        assert!((tau_at_cap - allow).abs() < 1e-9, "got: {tau_at_cap}");
    }

    // Capacity scales linearly with length: doubling the length doubles
    // the capacity, everything else fixed.
    #[test]
    fn capacity_scales_with_length() {
        let short = FilletWeld::new(6.0, 100.0)
            .unwrap()
            .capacity(100.0)
            .unwrap();
        let long = FilletWeld::new(6.0, 200.0)
            .unwrap()
            .capacity(100.0)
            .unwrap();
        assert!(
            (long - 2.0 * short).abs() < 1e-6,
            "short {short} long {long}"
        );
    }

    // Capacity scales linearly with leg, and DOUBLING the leg DOUBLES the
    // capacity (because throat = 0.707 * leg is linear in leg).
    #[test]
    fn doubling_leg_doubles_capacity() {
        let base = FilletWeld::new(5.0, 120.0)
            .unwrap()
            .capacity(110.0)
            .unwrap();
        let doubled = FilletWeld::new(10.0, 120.0)
            .unwrap()
            .capacity(110.0)
            .unwrap();
        assert!(
            (doubled - 2.0 * base).abs() < 1e-6,
            "base {base} doubled {doubled}"
        );
    }

    // Utilisation == force / capacity; at capacity it is exactly 1.0.
    #[test]
    fn utilisation_is_force_over_capacity() {
        let w = FilletWeld::new(7.0, 90.0).expect("valid");
        let allow = 100.0;
        let cap = w.capacity(allow).unwrap();
        let u = w.utilisation(cap, allow).unwrap();
        assert!((u - 1.0).abs() < 1e-9, "got: {u}");

        let half = w.utilisation(cap / 2.0, allow).unwrap();
        assert!((half - 0.5).abs() < 1e-9, "got: {half}");
    }

    #[test]
    fn rejects_non_positive_inputs() {
        assert!(throat(0.0).is_err());
        assert!(throat(-1.0).is_err());
        assert!(FilletWeld::new(0.0, 10.0).is_err());
        assert!(FilletWeld::new(10.0, -5.0).is_err());

        let w = FilletWeld::new(6.0, 100.0).unwrap();
        assert!(w.shear_stress(0.0).is_err());
        assert!(w.capacity(f64::NAN).is_err());
    }

    // required_length = F / (tau_allow * 0.707 * leg). Hand calc:
    // F=50000, leg=6, tau=100 -> 50000/(100*0.707*6) = 117.8689... mm.
    #[test]
    fn required_length_matches_hand_calc() {
        let l = required_length(50_000.0, 6.0, 100.0).expect("valid");
        let expected = 50_000.0 / (100.0 * 0.707 * 6.0);
        assert!((l - expected).abs() < 1e-9, "got: {l}");
    }

    // required_length is the inverse of capacity / shear_stress: a weld
    // of that length carries exactly the force at the allowable stress.
    #[test]
    fn required_length_round_trips_capacity_and_stress() {
        let (force, leg, allow) = (42_000.0, 8.0, 95.0);
        let l = required_length(force, leg, allow).unwrap();
        let w = FilletWeld::new(leg, l).unwrap();
        assert!((w.shear_stress(force).unwrap() - allow).abs() < 1e-9);
        assert!((w.capacity(allow).unwrap() - force).abs() < 1e-6);
        assert!((w.utilisation(force, allow).unwrap() - 1.0).abs() < 1e-9);
        // And the free-function capacity round-trips too.
        assert!((capacity(allow, leg, l).unwrap() - force).abs() < 1e-6);
    }

    // Sizing scales with force and inversely with allowable stress and leg.
    #[test]
    fn required_length_scaling() {
        let base = required_length(20_000.0, 6.0, 100.0).unwrap();
        let double_f = required_length(40_000.0, 6.0, 100.0).unwrap();
        let half_allow = required_length(20_000.0, 6.0, 50.0).unwrap();
        let half_leg = required_length(20_000.0, 3.0, 100.0).unwrap();
        assert!((double_f - 2.0 * base).abs() < 1e-9);
        assert!((half_allow - 2.0 * base).abs() < 1e-9);
        assert!((half_leg - 2.0 * base).abs() < 1e-9);
    }

    #[test]
    fn required_length_rejects_bad_inputs() {
        assert!(required_length(0.0, 6.0, 100.0).is_err());
        assert!(required_length(50_000.0, -6.0, 100.0).is_err());
        assert!(required_length(50_000.0, 6.0, f64::NAN).is_err());
    }
}
