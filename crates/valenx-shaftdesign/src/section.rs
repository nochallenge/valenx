//! The [`ShaftSection`] input model and its section properties.
//!
//! A section bundles the solid circular shaft diameter with the steady
//! bending moment and torque applied at that station. Construct it with
//! the validated [`ShaftSection::new`] so the diameter is guaranteed
//! positive and the loads finite before any stress is computed.

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use crate::error::ShaftError;

/// A solid circular shaft cross-section carrying a steady bending
/// moment and torque.
///
/// Use SI units consistently — diameter in metres, moment / torque in
/// newton-metres — and the stresses come out in pascals. Any coherent
/// unit set works (e.g. mm and N·mm giving MPa); the formulae are
/// unit-agnostic.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct ShaftSection {
    /// Outer diameter `d` of the solid shaft (> 0).
    pub diameter: f64,
    /// Applied bending moment `M` at the section. May be either sign;
    /// only its magnitude affects the fibre stresses.
    pub bending_moment: f64,
    /// Applied torque `T` at the section. May be either sign; only its
    /// magnitude affects the fibre stresses.
    pub torque: f64,
}

impl ShaftSection {
    /// Build a validated section.
    ///
    /// # Errors
    ///
    /// [`ShaftError::NonPositive`] / [`ShaftError::NotFinite`] if
    /// `diameter` is not strictly positive and finite, or
    /// [`ShaftError::NotFinite`] if either load is not finite.
    pub fn new(diameter: f64, bending_moment: f64, torque: f64) -> Result<Self, ShaftError> {
        let diameter = ShaftError::require_positive("diameter", diameter)?;
        let bending_moment = ShaftError::require_finite("bending_moment", bending_moment)?;
        let torque = ShaftError::require_finite("torque", torque)?;
        Ok(Self {
            diameter,
            bending_moment,
            torque,
        })
    }

    /// Pure-torsion section: a shaft carrying torque only (`M = 0`).
    ///
    /// # Errors
    ///
    /// As [`ShaftSection::new`].
    pub fn pure_torsion(diameter: f64, torque: f64) -> Result<Self, ShaftError> {
        Self::new(diameter, 0.0, torque)
    }

    /// Pure-bending section: a shaft carrying a bending moment only
    /// (`T = 0`).
    ///
    /// # Errors
    ///
    /// As [`ShaftSection::new`].
    pub fn pure_bending(diameter: f64, bending_moment: f64) -> Result<Self, ShaftError> {
        Self::new(diameter, bending_moment, 0.0)
    }

    /// Polar second moment of area `J = pi d^4 / 32` (about the axis).
    ///
    /// Governs torsional shear: `tau = T (d/2) / J`.
    pub fn polar_second_moment(&self) -> f64 {
        PI * self.diameter.powi(4) / 32.0
    }

    /// Second moment of area `I = pi d^4 / 64` (about a diameter).
    ///
    /// Governs bending stress: `sigma = M (d/2) / I`.
    pub fn second_moment(&self) -> f64 {
        PI * self.diameter.powi(4) / 64.0
    }

    /// Polar section modulus `Z_p = J / (d/2) = pi d^3 / 16`.
    ///
    /// The torsional shear is simply `tau = T / Z_p`.
    pub fn polar_section_modulus(&self) -> f64 {
        PI * self.diameter.powi(3) / 16.0
    }

    /// Section modulus `Z = I / (d/2) = pi d^3 / 32`.
    ///
    /// The bending stress is simply `sigma = M / Z`.
    pub fn section_modulus(&self) -> f64 {
        PI * self.diameter.powi(3) / 32.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    #[test]
    fn polar_second_moment_is_pi_d4_over_32() {
        let s = ShaftSection::new(0.05, 0.0, 0.0).unwrap();
        let want = PI * 0.05_f64.powi(4) / 32.0;
        assert!((s.polar_second_moment() - want).abs() < EPS);
    }

    #[test]
    fn second_moment_is_pi_d4_over_64() {
        let s = ShaftSection::new(0.05, 0.0, 0.0).unwrap();
        let want = PI * 0.05_f64.powi(4) / 64.0;
        assert!((s.second_moment() - want).abs() < EPS);
    }

    #[test]
    fn polar_is_exactly_twice_the_diametral_second_moment() {
        let s = ShaftSection::new(0.037, 0.0, 0.0).unwrap();
        assert!((s.polar_second_moment() - 2.0 * s.second_moment()).abs() < EPS);
    }

    #[test]
    fn polar_section_modulus_is_pi_d3_over_16() {
        let s = ShaftSection::new(0.05, 0.0, 0.0).unwrap();
        let want = PI * 0.05_f64.powi(3) / 16.0;
        assert!((s.polar_section_modulus() - want).abs() < EPS);
    }

    #[test]
    fn section_modulus_is_pi_d3_over_32() {
        let s = ShaftSection::new(0.05, 0.0, 0.0).unwrap();
        let want = PI * 0.05_f64.powi(3) / 32.0;
        assert!((s.section_modulus() - want).abs() < EPS);
    }

    #[test]
    fn modulus_equals_second_moment_over_outer_radius() {
        let s = ShaftSection::new(0.05, 0.0, 0.0).unwrap();
        let radius = s.diameter / 2.0;
        assert!((s.section_modulus() - s.second_moment() / radius).abs() < EPS);
        assert!((s.polar_section_modulus() - s.polar_second_moment() / radius).abs() < EPS);
    }

    #[test]
    fn constructor_stores_fields() {
        let s = ShaftSection::new(0.04, 100.0, 200.0).unwrap();
        assert!((s.diameter - 0.04).abs() < EPS);
        assert!((s.bending_moment - 100.0).abs() < EPS);
        assert!((s.torque - 200.0).abs() < EPS);
    }

    #[test]
    fn pure_torsion_zeroes_moment() {
        let s = ShaftSection::pure_torsion(0.04, 200.0).unwrap();
        assert!(s.bending_moment.abs() < EPS);
        assert!((s.torque - 200.0).abs() < EPS);
    }

    #[test]
    fn pure_bending_zeroes_torque() {
        let s = ShaftSection::pure_bending(0.04, 200.0).unwrap();
        assert!(s.torque.abs() < EPS);
        assert!((s.bending_moment - 200.0).abs() < EPS);
    }

    #[test]
    fn rejects_non_positive_diameter() {
        assert!(ShaftSection::new(0.0, 1.0, 1.0).is_err());
        assert!(ShaftSection::new(-1.0, 1.0, 1.0).is_err());
    }

    #[test]
    fn rejects_non_finite_inputs() {
        assert!(ShaftSection::new(f64::NAN, 1.0, 1.0).is_err());
        assert!(ShaftSection::new(0.05, f64::INFINITY, 1.0).is_err());
        assert!(ShaftSection::new(0.05, 1.0, f64::NAN).is_err());
    }

    #[test]
    fn serde_round_trips() {
        let s = ShaftSection::new(0.05, 600.0, 800.0).unwrap();
        let json = serde_json::to_string(&s).unwrap();
        let back: ShaftSection = serde_json::from_str(&json).unwrap();
        assert_eq!(s, back);
    }
}
