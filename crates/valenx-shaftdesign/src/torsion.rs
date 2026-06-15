//! Torsional shear stress of a solid circular shaft.
//!
//! The outer-fibre shear from Saint-Venant torsion is
//! `tau = T (d/2) / J` with `J = pi d^4 / 32`, which simplifies to the
//! closed form `tau = 16 T / (pi d^3)`. The sign of `T` is irrelevant
//! to the peak magnitude, so the magnitude `|T|` is used.

use std::f64::consts::PI;

use crate::error::ShaftError;
use crate::section::ShaftSection;

/// Maximum torsional shear stress `tau = 16 |T| / (pi d^3)` at the
/// outer fibre of a solid circular shaft.
///
/// Equivalent to `|T| / Z_p` with the polar section modulus
/// `Z_p = pi d^3 / 16` from
/// [`ShaftSection::polar_section_modulus`].
///
/// # Errors
///
/// Propagates the validation error from [`ShaftSection::new`] (the
/// diameter must be positive and finite, the torque finite).
pub fn shear_stress(diameter: f64, torque: f64) -> Result<f64, ShaftError> {
    let section = ShaftSection::pure_torsion(diameter, torque)?;
    Ok(shear_stress_section(&section))
}

/// Maximum torsional shear stress for an already-validated section.
///
/// Uses the magnitude of [`ShaftSection::torque`]; the bending moment
/// is ignored.
pub fn shear_stress_section(section: &ShaftSection) -> f64 {
    16.0 * section.torque.abs() / (PI * section.diameter.powi(3))
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn matches_closed_form_16t_over_pi_d_cubed() {
        // d = 0.05 m, T = 1000 N·m.
        let d = 0.05;
        let t = 1000.0;
        let got = shear_stress(d, t).unwrap();
        let want = 16.0 * t / (PI * d.powi(3));
        assert!((got - want).abs() < EPS, "got {got}, want {want}");
    }

    #[test]
    fn equals_torque_over_polar_section_modulus() {
        let section = ShaftSection::pure_torsion(0.04, 250.0).unwrap();
        let via_formula = shear_stress_section(&section);
        let via_modulus = section.torque.abs() / section.polar_section_modulus();
        assert!((via_formula - via_modulus).abs() < EPS);
    }

    #[test]
    fn hand_checked_numeric_value() {
        // d = 0.1 m -> pi d^3 = pi * 1e-3. tau = 16 * 500 / (pi * 1e-3)
        //            = 8000 / (pi * 1e-3) = 8e6 / pi Pa.
        let got = shear_stress(0.1, 500.0).unwrap();
        let want = 8.0e6 / PI;
        assert!((got - want).abs() < 1e-3, "got {got}, want {want}");
    }

    #[test]
    fn sign_of_torque_does_not_change_magnitude() {
        let pos = shear_stress(0.03, 420.0).unwrap();
        let neg = shear_stress(0.03, -420.0).unwrap();
        assert!((pos - neg).abs() < EPS);
    }

    #[test]
    fn scales_as_inverse_diameter_cubed() {
        // Doubling d at fixed torque must cut tau to one eighth.
        let t = 750.0;
        let small = shear_stress(0.02, t).unwrap();
        let large = shear_stress(0.04, t).unwrap();
        assert!(
            (small / large - 8.0).abs() < 1e-6,
            "ratio {}",
            small / large
        );
    }

    #[test]
    fn larger_diameter_lowers_stress() {
        let t = 600.0;
        let a = shear_stress(0.03, t).unwrap();
        let b = shear_stress(0.06, t).unwrap();
        assert!(b < a);
    }

    #[test]
    fn zero_torque_gives_zero_shear() {
        let got = shear_stress(0.05, 0.0).unwrap();
        assert!(got.abs() < EPS);
    }

    #[test]
    fn rejects_non_positive_diameter() {
        assert!(shear_stress(0.0, 100.0).is_err());
        assert!(shear_stress(-0.01, 100.0).is_err());
    }

    #[test]
    fn rejects_non_finite_inputs() {
        assert!(shear_stress(f64::NAN, 100.0).is_err());
        assert!(shear_stress(0.05, f64::INFINITY).is_err());
    }
}
