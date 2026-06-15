//! Bending normal stress of a solid circular shaft.
//!
//! The outer-fibre normal stress from Euler-Bernoulli bending is
//! `sigma = M (d/2) / I` with `I = pi d^4 / 64`, which simplifies to
//! the closed form `sigma = 32 M / (pi d^3)`. The sign of `M` is
//! irrelevant to the peak magnitude, so the magnitude `|M|` is used.

use std::f64::consts::PI;

use crate::error::ShaftError;
use crate::section::ShaftSection;

/// Maximum bending normal stress `sigma = 32 |M| / (pi d^3)` at the
/// outer fibre of a solid circular shaft.
///
/// Equivalent to `|M| / Z` with the section modulus
/// `Z = pi d^3 / 32` from [`ShaftSection::section_modulus`].
///
/// # Errors
///
/// Propagates the validation error from [`ShaftSection::new`] (the
/// diameter must be positive and finite, the moment finite).
pub fn bending_stress(diameter: f64, bending_moment: f64) -> Result<f64, ShaftError> {
    let section = ShaftSection::pure_bending(diameter, bending_moment)?;
    Ok(bending_stress_section(&section))
}

/// Maximum bending normal stress for an already-validated section.
///
/// Uses the magnitude of [`ShaftSection::bending_moment`]; the torque
/// is ignored.
pub fn bending_stress_section(section: &ShaftSection) -> f64 {
    32.0 * section.bending_moment.abs() / (PI * section.diameter.powi(3))
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn matches_closed_form_32m_over_pi_d_cubed() {
        let d = 0.05;
        let m = 800.0;
        let got = bending_stress(d, m).unwrap();
        let want = 32.0 * m / (PI * d.powi(3));
        assert!((got - want).abs() < EPS, "got {got}, want {want}");
    }

    #[test]
    fn equals_moment_over_section_modulus() {
        let section = ShaftSection::pure_bending(0.04, 250.0).unwrap();
        let via_formula = bending_stress_section(&section);
        let via_modulus = section.bending_moment.abs() / section.section_modulus();
        assert!((via_formula - via_modulus).abs() < EPS);
    }

    #[test]
    fn hand_checked_numeric_value() {
        // d = 0.1 m -> pi d^3 = pi * 1e-3. sigma = 32 * 500 / (pi * 1e-3)
        //            = 16000 / (pi * 1e-3) = 1.6e7 / pi Pa.
        let got = bending_stress(0.1, 500.0).unwrap();
        let want = 1.6e7 / PI;
        assert!((got - want).abs() < 1e-2, "got {got}, want {want}");
    }

    #[test]
    fn bending_is_exactly_twice_torsion_for_equal_load() {
        // sigma = 32 X / (pi d^3) is exactly 2x tau = 16 X / (pi d^3)
        // when M and T have the same magnitude X.
        let d = 0.03;
        let load = 350.0;
        let sigma = bending_stress(d, load).unwrap();
        let tau = crate::torsion::shear_stress(d, load).unwrap();
        assert!((sigma - 2.0 * tau).abs() < EPS, "sigma {sigma}, tau {tau}");
    }

    #[test]
    fn sign_of_moment_does_not_change_magnitude() {
        let pos = bending_stress(0.03, 420.0).unwrap();
        let neg = bending_stress(0.03, -420.0).unwrap();
        assert!((pos - neg).abs() < EPS);
    }

    #[test]
    fn scales_as_inverse_diameter_cubed() {
        let m = 900.0;
        let small = bending_stress(0.02, m).unwrap();
        let large = bending_stress(0.04, m).unwrap();
        assert!(
            (small / large - 8.0).abs() < 1e-6,
            "ratio {}",
            small / large
        );
    }

    #[test]
    fn larger_diameter_lowers_stress() {
        let m = 600.0;
        let a = bending_stress(0.03, m).unwrap();
        let b = bending_stress(0.06, m).unwrap();
        assert!(b < a);
    }

    #[test]
    fn zero_moment_gives_zero_bending() {
        let got = bending_stress(0.05, 0.0).unwrap();
        assert!(got.abs() < EPS);
    }

    #[test]
    fn rejects_non_positive_diameter() {
        assert!(bending_stress(0.0, 100.0).is_err());
        assert!(bending_stress(-0.01, 100.0).is_err());
    }

    #[test]
    fn rejects_non_finite_inputs() {
        assert!(bending_stress(f64::NAN, 100.0).is_err());
        assert!(bending_stress(0.05, f64::NEG_INFINITY).is_err());
    }
}
