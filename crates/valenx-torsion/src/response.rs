//! Elastic torsion response: shear stress, angle of twist, and power.
//!
//! ## Model
//!
//! For a prismatic circular bar carrying a torque `T`, linear-elastic
//! St. Venant torsion gives a shear stress that grows linearly with the
//! radius `r` from the axis:
//!
//! ```text
//! tau(r) = T * r / J
//! ```
//!
//! so the maximum shear stress is at the outer surface, `r = d / 2`:
//!
//! ```text
//! tau_max = T * (d / 2) / J
//! ```
//!
//! The total angle of twist over a length `L` of shaft with shear modulus
//! `G` is
//!
//! ```text
//! theta = T * L / (G * J)
//! ```
//!
//! and the mechanical power transmitted by a shaft rotating at angular
//! speed `omega` (radians per unit time) under torque `T` is
//!
//! ```text
//! P = T * omega
//! ```
//!
//! All formulas are exact for the circular section; keep the units
//! consistent (SI: `T` in N·m, lengths in m, `G` in Pa, `omega` in rad/s
//! gives `tau` in Pa, `theta` in rad, `P` in W).

use crate::error::{require_positive, TorsionError};
use crate::shaft::Shaft;

/// Shear stress `tau = T r / J` at radius `r` from the shaft axis.
///
/// The stress varies linearly with radius, vanishing on the axis (or on
/// the bore wall, for a tube the bore is stress-bearing too — the linear
/// profile simply continues inward through the hollow region, which has
/// no material). The query radius must lie within the material:
/// `inner_radius <= r <= outer_radius`.
///
/// # Errors
///
/// Returns [`TorsionError::NonPositive`] if `torque` is not finite and
/// strictly positive, or [`TorsionError::RadiusOutOfRange`] if `radius`
/// falls outside the cross-section.
pub fn shear_stress_at(shaft: &Shaft, torque: f64, radius: f64) -> Result<f64, TorsionError> {
    let torque = require_positive("torque", torque)?;
    let min_r = shaft.inner_radius();
    let max_r = shaft.outer_radius();
    if !(radius.is_finite() && (min_r..=max_r).contains(&radius)) {
        return Err(TorsionError::RadiusOutOfRange {
            radius,
            min_radius: min_r,
            max_radius: max_r,
        });
    }
    Ok(torque * radius / shaft.polar_moment())
}

/// Maximum shear stress `tau_max = T (d / 2) / J`, at the outer surface.
///
/// This is the design-driving stress for the section. It is exactly
/// [`shear_stress_at`] evaluated at the outer radius.
///
/// # Errors
///
/// Returns [`TorsionError::NonPositive`] if `torque` is not finite and
/// strictly positive.
pub fn max_shear_stress(shaft: &Shaft, torque: f64) -> Result<f64, TorsionError> {
    let torque = require_positive("torque", torque)?;
    Ok(torque * shaft.outer_radius() / shaft.polar_moment())
}

/// Angle of twist `theta = T L / (G J)` over a length `L` of shaft.
///
/// The result is in radians (when inputs are in consistent units). Twist
/// is inversely proportional to the **torsional rigidity** `G J`, so a
/// stiffer material or a fatter shaft twists less for the same torque.
///
/// # Errors
///
/// Returns [`TorsionError::NonPositive`] if `torque`, `length`, or
/// `shear_modulus` is not finite and strictly positive.
pub fn angle_of_twist(
    shaft: &Shaft,
    torque: f64,
    length: f64,
    shear_modulus: f64,
) -> Result<f64, TorsionError> {
    let torque = require_positive("torque", torque)?;
    let length = require_positive("length", length)?;
    let shear_modulus = require_positive("shear_modulus", shear_modulus)?;
    Ok(torque * length / (shear_modulus * shaft.polar_moment()))
}

/// Torsional rigidity `G J` of the shaft (stiffness per unit twist-rate).
///
/// `theta = T L / (G J)`, so `G J` is the proportionality between applied
/// torque and twist-per-unit-length. Provided as a convenience for
/// callers comparing candidate sections.
///
/// # Errors
///
/// Returns [`TorsionError::NonPositive`] if `shear_modulus` is not finite
/// and strictly positive.
pub fn torsional_rigidity(shaft: &Shaft, shear_modulus: f64) -> Result<f64, TorsionError> {
    let shear_modulus = require_positive("shear_modulus", shear_modulus)?;
    Ok(shear_modulus * shaft.polar_moment())
}

/// Mechanical power `P = T omega` transmitted by a rotating shaft.
///
/// With `T` in N·m and `omega` in rad/s the result is in watts. This is
/// a pure kinematic relation and does not depend on the section geometry.
///
/// # Errors
///
/// Returns [`TorsionError::NonPositive`] if `torque` or `angular_speed`
/// is not finite and strictly positive.
pub fn power(torque: f64, angular_speed: f64) -> Result<f64, TorsionError> {
    let torque = require_positive("torque", torque)?;
    let angular_speed = require_positive("angular_speed", angular_speed)?;
    Ok(torque * angular_speed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Relative tolerance for comparing closed-form expressions.
    const EPS: f64 = 1e-9;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < EPS * b.abs().max(1.0)
    }

    #[test]
    fn surface_stress_equals_torque_times_half_d_over_j() {
        let d = 40.0_f64;
        let shaft = Shaft::solid(d).unwrap();
        let t = 1_234.0_f64;
        let j = PI * d.powi(4) / 32.0;
        let expected = t * (d / 2.0) / j;
        let got = max_shear_stress(&shaft, t).unwrap();
        assert!(close(got, expected), "got {got}, expected {expected}");
    }

    #[test]
    fn max_shear_stress_equals_shear_stress_at_outer_radius() {
        let shaft = Shaft::hollow(50.0, 30.0).unwrap();
        let t = 900.0;
        let at_surface = shear_stress_at(&shaft, t, shaft.outer_radius()).unwrap();
        let max = max_shear_stress(&shaft, t).unwrap();
        assert!(close(at_surface, max), "{at_surface} vs {max}");
    }

    #[test]
    fn shear_stress_is_linear_in_radius() {
        // tau(r) = T r / J, so tau(2r) = 2 * tau(r).
        let shaft = Shaft::solid(20.0).unwrap();
        let t = 500.0;
        let inner = shear_stress_at(&shaft, t, 2.5).unwrap();
        let outer = shear_stress_at(&shaft, t, 5.0).unwrap();
        assert!(close(outer, 2.0 * inner), "{outer} vs {}", 2.0 * inner);
    }

    #[test]
    fn shear_stress_vanishes_on_the_axis_of_a_solid_shaft() {
        let shaft = Shaft::solid(12.0).unwrap();
        let on_axis = shear_stress_at(&shaft, 750.0, 0.0).unwrap();
        assert!(on_axis.abs() < EPS, "expected ~0 on axis, got {on_axis}");
    }

    #[test]
    fn angle_of_twist_matches_tl_over_gj() {
        let d = 25.0_f64;
        let shaft = Shaft::solid(d).unwrap();
        let t = 300.0_f64;
        let length = 2_000.0_f64;
        let g = 79_300.0_f64;
        let j = PI * d.powi(4) / 32.0;
        let expected = t * length / (g * j);
        let got = angle_of_twist(&shaft, t, length, g).unwrap();
        assert!(close(got, expected), "got {got}, expected {expected}");
    }

    #[test]
    fn twist_doubles_when_length_doubles() {
        let shaft = Shaft::solid(18.0).unwrap();
        let (t, g) = (220.0, 79_300.0);
        let short = angle_of_twist(&shaft, t, 1_000.0, g).unwrap();
        let long = angle_of_twist(&shaft, t, 2_000.0, g).unwrap();
        assert!(close(long, 2.0 * short), "{long} vs {}", 2.0 * short);
    }

    #[test]
    fn twist_halves_when_shear_modulus_doubles() {
        let shaft = Shaft::solid(18.0).unwrap();
        let (t, length) = (220.0, 1_500.0);
        let soft = angle_of_twist(&shaft, t, length, 40_000.0).unwrap();
        let stiff = angle_of_twist(&shaft, t, length, 80_000.0).unwrap();
        assert!(close(stiff, soft / 2.0), "{stiff} vs {}", soft / 2.0);
    }

    #[test]
    fn rigidity_equals_g_times_j_and_is_twist_denominator() {
        let shaft = Shaft::solid(22.0).unwrap();
        let g = 79_300.0;
        let gj = torsional_rigidity(&shaft, g).unwrap();
        assert!(close(gj, g * shaft.polar_moment()));

        // theta == T * L / (G J)
        let (t, length) = (450.0, 1_200.0);
        let theta = angle_of_twist(&shaft, t, length, g).unwrap();
        assert!(close(theta, t * length / gj));
    }

    #[test]
    fn power_equals_torque_times_omega() {
        // Worked check: 100 N·m at 30 rad/s = 3000 W.
        let p = power(100.0, 30.0).unwrap();
        assert!(close(p, 3_000.0), "got {p}");
    }

    #[test]
    fn power_at_one_rev_per_second_uses_two_pi_omega() {
        // 50 N·m spinning at 1 rev/s (omega = 2 pi rad/s).
        let t = 50.0;
        let omega = 2.0 * PI;
        let p = power(t, omega).unwrap();
        assert!(close(p, t * omega), "got {p}, expected {}", t * omega);
    }

    #[test]
    fn stress_query_rejects_radius_outside_section() {
        let shaft = Shaft::hollow(50.0, 30.0).unwrap();
        // Below the bore radius (15.0) and above the outer radius (25.0).
        assert!(matches!(
            shear_stress_at(&shaft, 100.0, 10.0),
            Err(TorsionError::RadiusOutOfRange { .. })
        ));
        assert!(matches!(
            shear_stress_at(&shaft, 100.0, 26.0),
            Err(TorsionError::RadiusOutOfRange { .. })
        ));
    }

    #[test]
    fn response_functions_reject_non_positive_inputs() {
        let shaft = Shaft::solid(10.0).unwrap();
        assert!(matches!(
            max_shear_stress(&shaft, 0.0),
            Err(TorsionError::NonPositive { name: "torque", .. })
        ));
        assert!(matches!(
            angle_of_twist(&shaft, 1.0, -1.0, 1.0),
            Err(TorsionError::NonPositive { name: "length", .. })
        ));
        assert!(matches!(
            angle_of_twist(&shaft, 1.0, 1.0, 0.0),
            Err(TorsionError::NonPositive {
                name: "shear_modulus",
                ..
            })
        ));
        assert!(matches!(
            power(1.0, f64::NAN),
            Err(TorsionError::NonPositive {
                name: "angular_speed",
                ..
            })
        ));
    }
}
