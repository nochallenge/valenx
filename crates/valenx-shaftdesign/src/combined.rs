//! Combined bending-and-torsion stress on a solid circular shaft.
//!
//! A point on the surface of a shaft loaded by bending moment `M` and
//! torque `T` sees a uniaxial normal stress `sigma = 32 |M| / (pi d^3)`
//! (from bending) together with a shear stress
//! `tau = 16 |T| / (pi d^3)` (from torsion). Mohr's circle for that
//! plane-stress state has centre `sigma / 2` and radius
//! `R = sqrt((sigma / 2)^2 + tau^2)`, so:
//!
//! - maximum shear stress: `tau_max = R`;
//! - maximum (principal) normal stress: `sigma_1 = sigma / 2 + R`.
//!
//! Substituting the fibre formulae, `tau_max = 16 T_e / (pi d^3)` with
//! the **equivalent torque** `T_e = sqrt(M^2 + T^2)`, and
//! `sigma_1 = 32 M_e / (pi d^3)` with the **equivalent bending moment**
//! `M_e = (|M| + sqrt(M^2 + T^2)) / 2` — the standard shaft-design
//! shortcuts that turn a combined load into a single equivalent load
//! reusing the pure-torsion / pure-bending modulus.

use crate::bending::bending_stress_section;
use crate::error::ShaftError;
use crate::section::ShaftSection;
use crate::torsion::shear_stress_section;

/// Equivalent torque `T_e = sqrt(M^2 + T^2)`.
///
/// The single torque that, acting alone, produces the same maximum
/// shear stress as the combined bending moment `M` and torque `T`.
/// Independent of diameter — it is a pure load combination.
///
/// # Errors
///
/// [`ShaftError::NotFinite`] if either load is not finite.
pub fn equivalent_torque(bending_moment: f64, torque: f64) -> Result<f64, ShaftError> {
    let m = ShaftError::require_finite("bending_moment", bending_moment)?;
    let t = ShaftError::require_finite("torque", torque)?;
    Ok(m.hypot(t))
}

/// Equivalent bending moment `M_e = (|M| + sqrt(M^2 + T^2)) / 2`.
///
/// The single bending moment that, acting alone, produces the same
/// maximum (principal) normal stress as the combined load. Independent
/// of diameter.
///
/// # Errors
///
/// [`ShaftError::NotFinite`] if either load is not finite.
pub fn equivalent_bending_moment(bending_moment: f64, torque: f64) -> Result<f64, ShaftError> {
    let m = ShaftError::require_finite("bending_moment", bending_moment)?;
    let t = ShaftError::require_finite("torque", torque)?;
    Ok(0.5 * (m.abs() + m.hypot(t)))
}

/// Combined maximum shear stress
/// `tau_max = sqrt((sigma / 2)^2 + tau^2)` for a validated section.
///
/// Equal to `16 sqrt(M^2 + T^2) / (pi d^3)`, i.e. the shear an
/// [`equivalent_torque`] would produce alone, and to
/// `equivalent_torque / Z_p` with the polar section modulus
/// [`ShaftSection::polar_section_modulus`]. This is the Tresca
/// maximum-shear-stress design quantity.
pub fn max_shear_stress(section: &ShaftSection) -> f64 {
    let sigma = bending_stress_section(section);
    let tau = shear_stress_section(section);
    (0.5 * sigma).hypot(tau)
}

/// Combined maximum (principal) normal stress
/// `sigma_1 = sigma / 2 + sqrt((sigma / 2)^2 + tau^2)` for a validated
/// section.
///
/// Equal to `32 M_e / (pi d^3)` with the
/// [`equivalent_bending_moment`], i.e. `equivalent_bending_moment / Z`
/// with the section modulus [`ShaftSection::section_modulus`]. This is
/// the maximum-principal-stress (Rankine) design quantity.
pub fn max_normal_stress(section: &ShaftSection) -> f64 {
    let sigma = bending_stress_section(section);
    let tau = shear_stress_section(section);
    0.5 * sigma + (0.5 * sigma).hypot(tau)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const EPS: f64 = 1e-9;

    #[test]
    fn equivalent_torque_345_pythagorean() {
        // M = 3, T = 4 -> T_e = 5 (the 3-4-5 triangle).
        let got = equivalent_torque(3.0, 4.0).unwrap();
        assert!((got - 5.0).abs() < EPS, "got {got}");
    }

    #[test]
    fn equivalent_torque_general_formula() {
        let m = 600.0;
        let t = 800.0;
        let got = equivalent_torque(m, t).unwrap();
        let want = (m * m + t * t).sqrt();
        assert!((got - want).abs() < 1e-6, "got {got}, want {want}");
        // 600-800-1000 triangle.
        assert!((got - 1000.0).abs() < 1e-6);
    }

    #[test]
    fn equivalent_torque_pure_torsion_recovers_torque() {
        let got = equivalent_torque(0.0, 900.0).unwrap();
        assert!((got - 900.0).abs() < EPS);
    }

    #[test]
    fn equivalent_torque_pure_bending_recovers_moment() {
        let got = equivalent_torque(900.0, 0.0).unwrap();
        assert!((got - 900.0).abs() < EPS);
    }

    #[test]
    fn equivalent_torque_sign_independent() {
        let a = equivalent_torque(-3.0, -4.0).unwrap();
        assert!((a - 5.0).abs() < EPS);
    }

    #[test]
    fn equivalent_bending_moment_formula() {
        // M = 3, T = 4 -> M_e = (3 + 5) / 2 = 4.
        let got = equivalent_bending_moment(3.0, 4.0).unwrap();
        assert!((got - 4.0).abs() < EPS, "got {got}");
    }

    #[test]
    fn equivalent_bending_moment_pure_bending_recovers_moment() {
        // T = 0 -> M_e = (|M| + |M|) / 2 = |M|.
        let got = equivalent_bending_moment(700.0, 0.0).unwrap();
        assert!((got - 700.0).abs() < EPS);
    }

    #[test]
    fn equivalent_bending_moment_pure_torsion_halves_torque() {
        // M = 0 -> M_e = (0 + |T|) / 2 = |T| / 2.
        let got = equivalent_bending_moment(0.0, 1000.0).unwrap();
        assert!((got - 500.0).abs() < EPS);
    }

    #[test]
    fn max_shear_equals_equivalent_torque_over_polar_modulus() {
        let section = ShaftSection::new(0.05, 600.0, 800.0).unwrap();
        let tau_max = max_shear_stress(&section);
        let te = equivalent_torque(section.bending_moment, section.torque).unwrap();
        let want = te / section.polar_section_modulus();
        assert!(
            (tau_max - want).abs() < 1e-9,
            "tau_max {tau_max}, want {want}"
        );
        // And the explicit 16 T_e / (pi d^3) closed form.
        let closed = 16.0 * te / (PI * section.diameter.powi(3));
        assert!((tau_max - closed).abs() < 1e-9);
    }

    #[test]
    fn max_normal_equals_equivalent_moment_over_section_modulus() {
        let section = ShaftSection::new(0.05, 600.0, 800.0).unwrap();
        let sigma_1 = max_normal_stress(&section);
        let me = equivalent_bending_moment(section.bending_moment, section.torque).unwrap();
        let want = me / section.section_modulus();
        assert!(
            (sigma_1 - want).abs() < 1e-9,
            "sigma_1 {sigma_1}, want {want}"
        );
        // And the explicit 32 M_e / (pi d^3) closed form.
        let closed = 32.0 * me / (PI * section.diameter.powi(3));
        assert!((sigma_1 - closed).abs() < 1e-9);
    }

    #[test]
    fn max_shear_pure_torsion_matches_torsion_module() {
        // With M = 0 the combined max shear must equal the plain
        // torsional shear 16 T / (pi d^3).
        let section = ShaftSection::pure_torsion(0.04, 500.0).unwrap();
        let got = max_shear_stress(&section);
        let want = crate::torsion::shear_stress(0.04, 500.0).unwrap();
        assert!((got - want).abs() < 1e-9, "got {got}, want {want}");
    }

    #[test]
    fn max_normal_pure_bending_matches_bending_module() {
        // With T = 0 the combined max normal must equal the plain
        // bending stress 32 M / (pi d^3): sigma/2 + sqrt((sigma/2)^2)
        // = sigma/2 + sigma/2 = sigma.
        let section = ShaftSection::pure_bending(0.04, 500.0).unwrap();
        let got = max_normal_stress(&section);
        let want = crate::bending::bending_stress(0.04, 500.0).unwrap();
        assert!((got - want).abs() < 1e-9, "got {got}, want {want}");
    }

    #[test]
    fn mohr_circle_hand_check() {
        // Construct a section whose fibre stresses are exactly
        // sigma = 6, tau = 4 in some coherent unit set, by choosing the
        // loads to hit those values: pick d so pi d^3 = 32, then
        // sigma = M, and tau = 0.5 * T. Easiest: assert the Mohr
        // relation directly on the two fibre stresses.
        let section = ShaftSection::new(0.05, 600.0, 800.0).unwrap();
        let sigma = bending_stress_section(&section);
        let tau = shear_stress_section(&section);
        let r = (0.5 * sigma).hypot(tau);
        assert!((max_shear_stress(&section) - r).abs() < 1e-12);
        assert!((max_normal_stress(&section) - (0.5 * sigma + r)).abs() < 1e-12);
        // tau_max >= each individual stress contribution's half-range.
        assert!(max_shear_stress(&section) >= tau);
    }

    #[test]
    fn larger_diameter_lowers_combined_stresses() {
        let small = ShaftSection::new(0.03, 600.0, 800.0).unwrap();
        let large = ShaftSection::new(0.06, 600.0, 800.0).unwrap();
        assert!(max_shear_stress(&large) < max_shear_stress(&small));
        assert!(max_normal_stress(&large) < max_normal_stress(&small));
        // Inverse-cube scaling: 8x reduction for a 2x diameter.
        let ratio = max_shear_stress(&small) / max_shear_stress(&large);
        assert!((ratio - 8.0).abs() < 1e-6, "ratio {ratio}");
    }

    #[test]
    fn rejects_non_finite_loads() {
        assert!(equivalent_torque(f64::NAN, 1.0).is_err());
        assert!(equivalent_torque(1.0, f64::INFINITY).is_err());
        assert!(equivalent_bending_moment(f64::NAN, 1.0).is_err());
    }
}
