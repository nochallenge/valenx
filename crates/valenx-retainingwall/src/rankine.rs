//! Rankine lateral earth-pressure theory for a smooth vertical wall.
//!
//! # Model
//!
//! For a frictionless (smooth) vertical wall retaining dry, cohesionless
//! soil with a horizontal ground surface, Rankine's theory gives the
//! lateral earth-pressure coefficients purely from the soil's internal
//! friction angle `phi`:
//!
//! ```text
//! Ka = tan^2(45 deg - phi/2)        (active)
//! Kp = tan^2(45 deg + phi/2) = 1/Ka (passive)
//! ```
//!
//! The lateral effective pressure grows linearly with depth `z` below the
//! ground surface,
//!
//! ```text
//! sigma_a(z) = Ka * gamma * z
//! sigma_p(z) = Kp * gamma * z
//! ```
//!
//! and, integrated over a wall of height `H`, the resultant thrust per
//! unit length of wall is the area of that triangular pressure diagram,
//!
//! ```text
//! Pa = 1/2 * Ka * gamma * H^2
//! ```
//!
//! acting horizontally at the centroid of the triangle, i.e. a height of
//! `H/3` above the base of the wall.
//!
//! # Honest scope
//!
//! This is the introductory textbook closed form only: a smooth vertical
//! wall, dry cohesionless soil (`c = 0`), and a level backfill. It does
//! **not** model cohesion, wall friction (Coulomb), sloping or surcharged
//! backfill, water tables / pore pressure, seismic (Mononobe-Okabe)
//! effects, layered soils, or wall stability checks (sliding, overturning,
//! bearing). Use it to learn and sanity-check, not to design.

use crate::error::RetainingWallError;
use serde::{Deserialize, Serialize};

/// Quarter turn, in degrees: the `45 deg` constant in Rankine's formulae.
const QUARTER_TURN_DEG: f64 = 45.0;

/// Convert an angle from degrees to radians.
#[inline]
fn to_radians(deg: f64) -> f64 {
    deg * (core::f64::consts::PI / 180.0)
}

/// Reject a value that is not a finite number.
#[inline]
fn require_finite(name: &'static str, value: f64) -> Result<f64, RetainingWallError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(RetainingWallError::NonFinite { name, value })
    }
}

/// A validated dry, cohesionless soil for Rankine earth-pressure
/// calculations.
///
/// Construct one with [`SoilProfile::new`], which enforces the physical
/// domain (`0 <= phi < 90 deg`, `gamma > 0`). Once built, the stored
/// `phi_deg` and `gamma` are guaranteed finite and in range, so the
/// coefficient and pressure accessors below are infallible.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SoilProfile {
    /// Internal angle of friction `phi`, in degrees, with
    /// `0 <= phi < 90`.
    phi_deg: f64,
    /// Soil unit weight `gamma` (e.g. kN/m^3), strictly positive.
    gamma: f64,
}

impl SoilProfile {
    /// Build a validated soil profile from a friction angle (degrees) and
    /// a unit weight.
    ///
    /// # Errors
    ///
    /// Returns [`RetainingWallError::NonFinite`] if either argument is NaN
    /// or infinite, [`RetainingWallError::FrictionAngleOutOfRange`] if
    /// `phi_deg` is not in `[0, 90)`, or
    /// [`RetainingWallError::NonPositiveUnitWeight`] if `gamma <= 0`.
    pub fn new(phi_deg: f64, gamma: f64) -> Result<Self, RetainingWallError> {
        let phi_deg = require_finite("phi_deg", phi_deg)?;
        let gamma = require_finite("gamma", gamma)?;

        if !(0.0..90.0).contains(&phi_deg) {
            return Err(RetainingWallError::FrictionAngleOutOfRange { phi_deg });
        }
        if gamma <= 0.0 {
            return Err(RetainingWallError::NonPositiveUnitWeight { gamma });
        }

        Ok(Self { phi_deg, gamma })
    }

    /// The internal friction angle `phi`, in degrees.
    #[inline]
    pub fn phi_deg(&self) -> f64 {
        self.phi_deg
    }

    /// The internal friction angle `phi`, in radians.
    #[inline]
    pub fn phi_rad(&self) -> f64 {
        to_radians(self.phi_deg)
    }

    /// The soil unit weight `gamma` (caller's units, e.g. kN/m^3).
    #[inline]
    pub fn gamma(&self) -> f64 {
        self.gamma
    }

    /// Rankine **active** earth-pressure coefficient
    /// `Ka = tan^2(45 deg - phi/2)`.
    ///
    /// For `phi > 0` this is strictly less than 1; at `phi = 0` it equals
    /// the at-rest hydrostatic-like value of exactly 1.
    #[inline]
    pub fn ka(&self) -> f64 {
        let t = to_radians(QUARTER_TURN_DEG - self.phi_deg / 2.0).tan();
        t * t
    }

    /// Rankine **passive** earth-pressure coefficient
    /// `Kp = tan^2(45 deg + phi/2)`.
    ///
    /// This is the exact reciprocal of [`SoilProfile::ka`]; for `phi > 0`
    /// it is strictly greater than 1.
    #[inline]
    pub fn kp(&self) -> f64 {
        let t = to_radians(QUARTER_TURN_DEG + self.phi_deg / 2.0).tan();
        t * t
    }

    /// Lateral **active** pressure `sigma_a = Ka * gamma * z` at a depth
    /// `z` below the ground surface.
    ///
    /// # Errors
    ///
    /// Returns [`RetainingWallError::NonFinite`] if `depth` is not finite,
    /// or [`RetainingWallError::NegativeDepth`] if `depth < 0`.
    pub fn active_pressure_at(&self, depth: f64) -> Result<f64, RetainingWallError> {
        let z = self.validated_depth(depth)?;
        Ok(self.ka() * self.gamma * z)
    }

    /// Lateral **passive** pressure `sigma_p = Kp * gamma * z` at a depth
    /// `z` below the ground surface.
    ///
    /// # Errors
    ///
    /// Returns [`RetainingWallError::NonFinite`] if `depth` is not finite,
    /// or [`RetainingWallError::NegativeDepth`] if `depth < 0`.
    pub fn passive_pressure_at(&self, depth: f64) -> Result<f64, RetainingWallError> {
        let z = self.validated_depth(depth)?;
        Ok(self.kp() * self.gamma * z)
    }

    /// Validate a depth/height argument: must be finite and non-negative.
    #[inline]
    fn validated_depth(&self, value: f64) -> Result<f64, RetainingWallError> {
        let v = require_finite("depth", value)?;
        if v < 0.0 {
            Err(RetainingWallError::NegativeDepth { value: v })
        } else {
            Ok(v)
        }
    }
}

/// The integrated lateral-thrust result for a wall of a given height,
/// reported as the resultant force per unit length of wall plus the height
/// of its line of action above the base.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Thrust {
    /// Resultant horizontal force per unit length of wall (e.g. kN/m):
    /// `1/2 * K * gamma * H^2`.
    pub resultant: f64,
    /// Height of the line of action of [`Thrust::resultant`] above the
    /// base of the wall (e.g. m): `H/3` for the triangular Rankine
    /// diagram.
    pub line_of_action: f64,
}

impl SoilProfile {
    /// Resultant **active** thrust on a wall of height `height`.
    ///
    /// Integrates the triangular active-pressure diagram to give
    /// `Pa = 1/2 * Ka * gamma * H^2` per unit length of wall, acting at
    /// `H/3` above the base.
    ///
    /// # Errors
    ///
    /// Returns [`RetainingWallError::NonFinite`] if `height` is not
    /// finite, or [`RetainingWallError::NegativeDepth`] if `height < 0`.
    pub fn active_thrust(&self, height: f64) -> Result<Thrust, RetainingWallError> {
        let h = self.validated_depth(height)?;
        Ok(Thrust {
            resultant: 0.5 * self.ka() * self.gamma * h * h,
            line_of_action: h / 3.0,
        })
    }

    /// Resultant **passive** thrust on a wall of height `height`.
    ///
    /// Integrates the triangular passive-pressure diagram to give
    /// `Pp = 1/2 * Kp * gamma * H^2` per unit length of wall, acting at
    /// `H/3` above the base.
    ///
    /// # Errors
    ///
    /// Returns [`RetainingWallError::NonFinite`] if `height` is not
    /// finite, or [`RetainingWallError::NegativeDepth`] if `height < 0`.
    pub fn passive_thrust(&self, height: f64) -> Result<Thrust, RetainingWallError> {
        let h = self.validated_depth(height)?;
        Ok(Thrust {
            resultant: 0.5 * self.kp() * self.gamma * h * h,
            line_of_action: h / 3.0,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::ErrorCategory;

    /// Absolute tolerance for floating-point ground-truth comparisons.
    const EPS: f64 = 1e-9;

    /// A convenience constructor used throughout the tests.
    fn soil(phi_deg: f64, gamma: f64) -> SoilProfile {
        SoilProfile::new(phi_deg, gamma).expect("valid soil profile")
    }

    #[test]
    fn ka_matches_closed_form_for_thirty_degrees() {
        // phi = 30 deg => Ka = tan^2(45 - 15) = tan^2(30 deg) = 1/3.
        let s = soil(30.0, 18.0);
        assert!((s.ka() - 1.0 / 3.0).abs() < EPS, "Ka = {}", s.ka());
    }

    #[test]
    fn kp_matches_closed_form_for_thirty_degrees() {
        // phi = 30 deg => Kp = tan^2(45 + 15) = tan^2(60 deg) = 3.
        let s = soil(30.0, 18.0);
        assert!((s.kp() - 3.0).abs() < EPS, "Kp = {}", s.kp());
    }

    #[test]
    fn kp_is_reciprocal_of_ka() {
        // Kp = 1/Ka exactly, for several friction angles.
        for &phi in &[1.0_f64, 10.0, 20.0, 25.0, 33.0, 40.0, 45.0, 60.0, 80.0] {
            let s = soil(phi, 20.0);
            let product = s.ka() * s.kp();
            assert!(
                (product - 1.0).abs() < 1e-9,
                "Ka*Kp = {product} at phi = {phi}"
            );
            assert!(
                (s.kp() - 1.0 / s.ka()).abs() < 1e-9,
                "Kp != 1/Ka at phi = {phi}"
            );
        }
    }

    #[test]
    fn active_less_than_one_less_than_passive() {
        // For phi > 0: Ka < 1 < Kp.
        for &phi in &[5.0_f64, 15.0, 30.0, 45.0, 60.0] {
            let s = soil(phi, 18.0);
            assert!(s.ka() < 1.0, "Ka = {} not < 1 at phi = {phi}", s.ka());
            assert!(s.kp() > 1.0, "Kp = {} not > 1 at phi = {phi}", s.kp());
            assert!(s.ka() < s.kp(), "Ka not < Kp at phi = {phi}");
        }
    }

    #[test]
    fn at_zero_friction_both_coefficients_are_one() {
        // phi = 0 => Ka = Kp = tan^2(45 deg) = 1 (the at-rest limit).
        let s = soil(0.0, 18.0);
        assert!((s.ka() - 1.0).abs() < EPS, "Ka = {}", s.ka());
        assert!((s.kp() - 1.0).abs() < EPS, "Kp = {}", s.kp());
    }

    #[test]
    fn pressure_is_linear_in_depth() {
        // sigma_a(z) = Ka*gamma*z is linear: doubling z doubles pressure,
        // and the value at z equals z * (constant slope Ka*gamma).
        let s = soil(35.0, 19.0);
        let slope = s.ka() * s.gamma();

        let p1 = s.active_pressure_at(2.0).expect("valid depth");
        let p2 = s.active_pressure_at(4.0).expect("valid depth");

        assert!((p1 - slope * 2.0).abs() < EPS, "p(2) = {p1}");
        assert!((p2 - slope * 4.0).abs() < EPS, "p(4) = {p2}");
        // Linearity: p(4) is exactly twice p(2).
        assert!((p2 - 2.0 * p1).abs() < EPS, "p(4) != 2*p(2)");
    }

    #[test]
    fn pressure_at_surface_is_zero() {
        // At z = 0 the lateral pressure vanishes.
        let s = soil(28.0, 17.5);
        let p = s.active_pressure_at(0.0).expect("valid depth");
        assert!(p.abs() < EPS, "surface pressure = {p}");
    }

    #[test]
    fn active_thrust_matches_half_ka_gamma_h_squared_at_third_height() {
        // Pa = 1/2 * Ka * gamma * H^2, acting at H/3 above the base.
        // phi = 30 deg, gamma = 18, H = 5 => Pa = 0.5*(1/3)*18*25 = 75.
        let s = soil(30.0, 18.0);
        let t = s.active_thrust(5.0).expect("valid height");

        let expected_force = 0.5 * (1.0 / 3.0) * 18.0 * 25.0;
        assert!(
            (t.resultant - expected_force).abs() < EPS,
            "Pa = {} expected {expected_force}",
            t.resultant
        );
        assert!((t.resultant - 75.0).abs() < EPS, "Pa = {}", t.resultant);
        assert!(
            (t.line_of_action - 5.0 / 3.0).abs() < EPS,
            "line of action = {}",
            t.line_of_action
        );
    }

    #[test]
    fn passive_thrust_matches_half_kp_gamma_h_squared() {
        // Pp = 1/2 * Kp * gamma * H^2 at H/3.
        // phi = 30 deg => Kp = 3; gamma = 18, H = 5 => Pp = 0.5*3*18*25 = 675.
        let s = soil(30.0, 18.0);
        let t = s.passive_thrust(5.0).expect("valid height");
        assert!((t.resultant - 675.0).abs() < EPS, "Pp = {}", t.resultant);
        assert!(
            (t.line_of_action - 5.0 / 3.0).abs() < EPS,
            "line of action = {}",
            t.line_of_action
        );
    }

    #[test]
    fn thrust_is_area_of_pressure_triangle() {
        // The resultant must equal 1/2 * base_pressure * height, where the
        // base pressure is sigma_a(H). Cross-check the two formulations.
        let s = soil(33.0, 20.0);
        let h = 6.0;
        let base_pressure = s.active_pressure_at(h).expect("valid depth");
        let triangle_area = 0.5 * base_pressure * h;
        let t = s.active_thrust(h).expect("valid height");
        assert!(
            (t.resultant - triangle_area).abs() < EPS,
            "resultant {} != triangle area {triangle_area}",
            t.resultant
        );
    }

    #[test]
    fn higher_phi_lowers_active_pressure_and_thrust() {
        // Strictly increasing phi => strictly decreasing Ka, active
        // pressure (at fixed depth) and active thrust (at fixed height).
        let gamma = 18.0;
        let depth = 4.0;
        let height = 4.0;

        let mut last_ka = f64::INFINITY;
        let mut last_p = f64::INFINITY;
        let mut last_force = f64::INFINITY;

        for &phi in &[10.0_f64, 20.0, 30.0, 40.0, 50.0] {
            let s = soil(phi, gamma);
            let ka = s.ka();
            let p = s.active_pressure_at(depth).expect("valid depth");
            let force = s.active_thrust(height).expect("valid height").resultant;

            assert!(ka < last_ka, "Ka not decreasing at phi = {phi}");
            assert!(p < last_p, "active pressure not decreasing at phi = {phi}");
            assert!(force < last_force, "thrust not decreasing at phi = {phi}");

            last_ka = ka;
            last_p = p;
            last_force = force;
        }
    }

    #[test]
    fn higher_phi_raises_passive_pressure() {
        // The passive coefficient (and pressure) move the opposite way:
        // increasing phi strictly increases Kp.
        let mut last_kp = 0.0;
        for &phi in &[10.0_f64, 20.0, 30.0, 40.0, 50.0] {
            let s = soil(phi, 18.0);
            assert!(s.kp() > last_kp, "Kp not increasing at phi = {phi}");
            last_kp = s.kp();
        }
    }

    #[test]
    fn rejects_friction_angle_at_or_above_ninety() {
        let err = SoilProfile::new(90.0, 18.0).expect_err("phi = 90 is out of range");
        assert_eq!(err.code(), "retainingwall.friction_angle_out_of_range");
        assert_eq!(err.category(), ErrorCategory::Input);

        assert!(SoilProfile::new(120.0, 18.0).is_err());
    }

    #[test]
    fn rejects_negative_friction_angle() {
        let err = SoilProfile::new(-1.0, 18.0).expect_err("negative phi is invalid");
        assert_eq!(err.code(), "retainingwall.friction_angle_out_of_range");
    }

    #[test]
    fn rejects_non_positive_unit_weight() {
        let err = SoilProfile::new(30.0, 0.0).expect_err("gamma = 0 is invalid");
        assert_eq!(err.code(), "retainingwall.non_positive_unit_weight");

        assert!(SoilProfile::new(30.0, -5.0).is_err());
    }

    #[test]
    fn rejects_non_finite_inputs() {
        assert_eq!(
            SoilProfile::new(f64::NAN, 18.0)
                .expect_err("NaN phi")
                .code(),
            "retainingwall.non_finite"
        );
        assert_eq!(
            SoilProfile::new(30.0, f64::INFINITY)
                .expect_err("inf gamma")
                .code(),
            "retainingwall.non_finite"
        );
    }

    #[test]
    fn rejects_negative_depth_and_height() {
        let s = soil(30.0, 18.0);
        assert_eq!(
            s.active_pressure_at(-0.5)
                .expect_err("negative depth")
                .code(),
            "retainingwall.negative_depth"
        );
        assert_eq!(
            s.active_thrust(-2.0).expect_err("negative height").code(),
            "retainingwall.negative_depth"
        );
        assert_eq!(
            s.active_pressure_at(f64::NAN)
                .expect_err("nan depth")
                .code(),
            "retainingwall.non_finite"
        );
    }

    #[test]
    fn soil_profile_serde_round_trips() {
        let s = soil(32.0, 18.5);
        let json = serde_json::to_string(&s).expect("serialize");
        let back: SoilProfile = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(s, back);
    }

    #[test]
    fn thrust_serde_round_trips() {
        let s = soil(30.0, 18.0);
        let t = s.active_thrust(5.0).expect("valid height");
        let json = serde_json::to_string(&t).expect("serialize");
        let back: Thrust = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(t, back);
    }
}
