//! Cross-section geometry: second moment of area and Euler-Bernoulli
//! bending stress.
//!
//! A long bone shaft (a diaphysis) is idealised here as a **hollow
//! circular tube** — a dense cortical wall around a hollow medullary
//! canal. Its resistance to bending is governed by the **second moment
//! of area** `I` of that annular cross-section. This module computes `I`
//! and the bending (flexural) stress it produces under a bending moment.
//!
//! ## Units
//!
//! | Quantity                  | Symbol | Unit  |
//! |---------------------------|--------|-------|
//! | Diameter                  | `D`    | mm    |
//! | Second moment of area     | `I`    | mm^4  |
//! | Bending moment            | `M`    | N mm  |
//! | Extreme-fibre distance    | `c`    | mm    |
//! | Bending (flexural) stress | `sigma`| MPa   |
//!
//! With those units the flexure formula `sigma = M c / I` comes out
//! directly in MPa, because `1 MPa = 1 N/mm^2` and
//! `(N mm)(mm) / mm^4 = N/mm^2`.
//!
//! ## Model
//!
//! The standard closed form for the second moment of area of an annulus
//! (a hollow circle) about a centroidal diameter is
//!
//! ```text
//! I = (pi / 64) * (D_o^4 - D_i^4)
//! ```
//!
//! where `D_o` is the outer diameter and `D_i` the inner (bore)
//! diameter. Setting `D_i = 0` recovers the solid-circle value
//! `I = (pi / 64) * D^4`. Euler-Bernoulli beam theory then gives the
//! peak normal stress at a fibre a distance `c` from the neutral axis as
//! `sigma = M c / I`; for a symmetric tube the maximum is at the outer
//! surface, `c = D_o / 2`. The relation inverts to `M = sigma I / c`, so
//! feeding the tissue's ultimate stress gives the bending fracture
//! moment of the section.

use crate::error::{BoneError, Result};
use std::f64::consts::PI;

/// Second moment of area of a **hollow circular** cross-section about a
/// centroidal diameter, in mm^4.
///
/// Computes `I = (pi / 64) * (outer_d^4 - inner_d^4)` where both
/// diameters are in millimetres. Pass `inner_d = 0.0` for a solid circle.
///
/// # Errors
///
/// Returns [`BoneError::Invalid`] if `outer_d <= 0`, if `inner_d < 0`, or
/// if either diameter is non-finite; returns [`BoneError::Geometry`] if
/// `inner_d >= outer_d` (the bore would meet or exceed the wall).
///
/// # Examples
///
/// ```
/// use valenx_bonemech::second_moment_hollow_circle_mm4;
///
/// // A solid circle of diameter 20 mm: I = pi/64 * 20^4.
/// let i = second_moment_hollow_circle_mm4(20.0, 0.0).unwrap();
/// assert!((i - std::f64::consts::PI / 64.0 * 20f64.powi(4)).abs() < 1e-9);
/// ```
pub fn second_moment_hollow_circle_mm4(outer_d: f64, inner_d: f64) -> Result<f64> {
    if !outer_d.is_finite() || outer_d <= 0.0 {
        return Err(BoneError::invalid(
            "outer_d",
            "outer diameter must be a positive, finite length in mm",
        ));
    }
    if !inner_d.is_finite() || inner_d < 0.0 {
        return Err(BoneError::invalid(
            "inner_d",
            "inner diameter must be a non-negative, finite length in mm",
        ));
    }
    if inner_d >= outer_d {
        return Err(BoneError::geometry(format!(
            "inner diameter ({inner_d} mm) must be strictly less than outer diameter ({outer_d} mm)"
        )));
    }
    Ok(PI / 64.0 * (outer_d.powi(4) - inner_d.powi(4)))
}

/// Section modulus `S = I / c` of a cross-section, in mm^3.
///
/// The section modulus packages the geometry so that peak bending stress
/// is simply `sigma = M / S`. Here `i_mm4` is the second moment of area
/// (mm^4) and `c_mm` is the extreme-fibre distance from the neutral axis
/// (mm).
///
/// # Errors
///
/// Returns [`BoneError::Invalid`] if `i_mm4 <= 0` or `c_mm <= 0`, or if
/// either value is non-finite.
pub fn section_modulus_mm3(i_mm4: f64, c_mm: f64) -> Result<f64> {
    if !i_mm4.is_finite() || i_mm4 <= 0.0 {
        return Err(BoneError::invalid(
            "i_mm4",
            "second moment of area must be positive and finite",
        ));
    }
    if !c_mm.is_finite() || c_mm <= 0.0 {
        return Err(BoneError::invalid(
            "c_mm",
            "extreme-fibre distance must be positive and finite",
        ));
    }
    Ok(i_mm4 / c_mm)
}

/// Euler-Bernoulli bending (flexural) stress, in MPa.
///
/// Computes `sigma = moment_nmm * c_mm / i_mm4`, the peak normal stress
/// at a fibre a distance `c_mm` from the neutral axis under a bending
/// moment `moment_nmm`. With moment in N mm, `c` in mm and `I` in mm^4
/// the result is in N/mm^2 = MPa.
///
/// The moment may be signed (sagging vs. hogging); the sign carries
/// through, so a negative moment yields a negative (sign-flipped) stress.
///
/// # Errors
///
/// Returns [`BoneError::Invalid`] if `i_mm4 <= 0`, if `c_mm < 0`, or if
/// any argument is non-finite. (A zero moment is permitted and yields a
/// zero stress.)
///
/// # Examples
///
/// ```
/// use valenx_bonemech::{second_moment_hollow_circle_mm4, bending_stress_mpa};
///
/// let i = second_moment_hollow_circle_mm4(20.0, 12.0).unwrap();
/// // Outer fibre is at c = D_o / 2 = 10 mm under a 50 000 N mm moment.
/// let sigma = bending_stress_mpa(50_000.0, 10.0, i).unwrap();
/// assert!(sigma > 0.0);
/// ```
pub fn bending_stress_mpa(moment_nmm: f64, c_mm: f64, i_mm4: f64) -> Result<f64> {
    if !moment_nmm.is_finite() {
        return Err(BoneError::invalid(
            "moment_nmm",
            "bending moment must be finite",
        ));
    }
    if !c_mm.is_finite() || c_mm < 0.0 {
        return Err(BoneError::invalid(
            "c_mm",
            "extreme-fibre distance must be non-negative and finite",
        ));
    }
    if !i_mm4.is_finite() || i_mm4 <= 0.0 {
        return Err(BoneError::invalid(
            "i_mm4",
            "second moment of area must be positive and finite",
        ));
    }
    Ok(moment_nmm * c_mm / i_mm4)
}

/// The bending moment that produces a given extreme-fibre stress, in
/// N mm — the inverse of [`bending_stress_mpa`]:
///
/// ```text
/// M = sigma * I / c = sigma * S
/// ```
///
/// Feeding the tissue's ultimate stress yields the **bending fracture
/// moment** of the section — the flexural counterpart of an axial
/// fracture load. The stress may be signed; the sign carries through to
/// the moment.
///
/// # Errors
///
/// Returns [`BoneError::Invalid`] if `stress_mpa` is non-finite, if
/// `c_mm <= 0`, or if `i_mm4 <= 0` (or either is non-finite).
///
/// # Examples
///
/// ```
/// use valenx_bonemech::{bending_stress_mpa, bending_moment_for_stress};
/// // Round-trip: moment -> stress -> moment.
/// let (c, i) = (5.0, 1000.0);
/// let sigma = bending_stress_mpa(20_000.0, c, i).unwrap();
/// let m = bending_moment_for_stress(sigma, c, i).unwrap();
/// assert!((m - 20_000.0).abs() < 1e-6);
/// ```
pub fn bending_moment_for_stress(stress_mpa: f64, c_mm: f64, i_mm4: f64) -> Result<f64> {
    if !stress_mpa.is_finite() {
        return Err(BoneError::invalid("stress_mpa", "stress must be finite"));
    }
    if !c_mm.is_finite() || c_mm <= 0.0 {
        return Err(BoneError::invalid(
            "c_mm",
            "extreme-fibre distance must be positive and finite",
        ));
    }
    if !i_mm4.is_finite() || i_mm4 <= 0.0 {
        return Err(BoneError::invalid(
            "i_mm4",
            "second moment of area must be positive and finite",
        ));
    }
    Ok(stress_mpa * i_mm4 / c_mm)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Solid-circle limit: with `inner_d = 0` the annulus formula must
    /// reduce to the textbook `I = pi/64 * D^4`. For `D = 10 mm` that is
    /// `pi/64 * 10000 ~= 490.8739 mm^4`.
    #[test]
    fn solid_circle_matches_analytic() {
        let i = second_moment_hollow_circle_mm4(10.0, 0.0).unwrap();
        let expected = PI / 64.0 * 10f64.powi(4);
        assert!((i - expected).abs() < 1e-9, "got {i}, expected {expected}");
        // Cross-check against an independently computed constant.
        assert!((i - 490.873_852_123_405_2).abs() < 1e-6, "got {i}");
    }

    /// Hollow-circle on a known geometry: D_o = 30 mm, D_i = 20 mm gives
    /// I = pi/64 * (30^4 - 20^4) = pi/64 * (810000 - 160000)
    ///   = pi/64 * 650000 ~= 31906.80 mm^4.
    #[test]
    fn hollow_circle_known_geometry() {
        let i = second_moment_hollow_circle_mm4(30.0, 20.0).unwrap();
        let expected = PI / 64.0 * (30f64.powi(4) - 20f64.powi(4));
        assert!((i - expected).abs() < 1e-9, "got {i}, expected {expected}");
        assert!((i - 31_906.800_388_021_336).abs() < 1e-3, "got {i}");
    }

    /// A hollow tube has strictly less `I` than the solid bar of the same
    /// outer diameter — removing central material can only reduce it.
    #[test]
    fn hollow_is_weaker_than_solid_same_outer() {
        let solid = second_moment_hollow_circle_mm4(25.0, 0.0).unwrap();
        let hollow = second_moment_hollow_circle_mm4(25.0, 15.0).unwrap();
        assert!(hollow < solid, "hollow {hollow} should be < solid {solid}");
    }

    /// Section modulus packages I and c so that `sigma = M / S` equals the
    /// direct flexure formula `sigma = M c / I`.
    #[test]
    fn section_modulus_consistent_with_flexure() {
        let i = second_moment_hollow_circle_mm4(20.0, 12.0).unwrap();
        let c = 10.0;
        let m = 75_000.0;
        let s = section_modulus_mm3(i, c).unwrap();
        let via_modulus = m / s;
        let direct = bending_stress_mpa(m, c, i).unwrap();
        assert!(
            (via_modulus - direct).abs() < 1e-9,
            "modulus path {via_modulus} != direct {direct}"
        );
    }

    /// Bending stress equals `M c / I` exactly for a worked example.
    #[test]
    fn bending_stress_worked_example() {
        // I = 1000 mm^4, c = 5 mm, M = 20000 N mm -> sigma = 100 MPa.
        let sigma = bending_stress_mpa(20_000.0, 5.0, 1000.0).unwrap();
        assert!((sigma - 100.0).abs() < 1e-9, "got {sigma}");
    }

    /// A negative (hogging) moment flips the sign of the stress.
    #[test]
    fn bending_stress_sign_follows_moment() {
        let pos = bending_stress_mpa(20_000.0, 5.0, 1000.0).unwrap();
        let neg = bending_stress_mpa(-20_000.0, 5.0, 1000.0).unwrap();
        assert!(
            (pos + neg).abs() < 1e-9,
            "expected opposite signs: {pos}, {neg}"
        );
    }

    #[test]
    fn invalid_diameters_rejected() {
        assert!(second_moment_hollow_circle_mm4(0.0, 0.0).is_err());
        assert!(second_moment_hollow_circle_mm4(-5.0, 0.0).is_err());
        assert!(second_moment_hollow_circle_mm4(10.0, -1.0).is_err());
        // Inner >= outer is a geometry inconsistency, not a scalar domain error.
        let e = second_moment_hollow_circle_mm4(10.0, 10.0).unwrap_err();
        assert_eq!(e.category(), crate::ErrorCategory::Geometry);
        assert!(second_moment_hollow_circle_mm4(10.0, 12.0).is_err());
    }

    #[test]
    fn invalid_bending_inputs_rejected() {
        assert!(bending_stress_mpa(f64::NAN, 5.0, 1000.0).is_err());
        assert!(bending_stress_mpa(100.0, -1.0, 1000.0).is_err());
        assert!(bending_stress_mpa(100.0, 5.0, 0.0).is_err());
        // Zero moment is allowed and yields zero stress.
        let sigma = bending_stress_mpa(0.0, 5.0, 1000.0).unwrap();
        assert!(sigma.abs() < 1e-12, "got {sigma}");
    }

    /// The moment-for-stress inverse round-trips with the flexure formula.
    #[test]
    fn moment_for_stress_inverts_bending_stress() {
        let i = second_moment_hollow_circle_mm4(20.0, 12.0).unwrap();
        let c = 10.0;
        let m0 = 75_000.0;
        let sigma = bending_stress_mpa(m0, c, i).unwrap();
        let m = bending_moment_for_stress(sigma, c, i).unwrap();
        assert!((m - m0).abs() / m0 < 1e-12, "got {m}");
    }

    /// Closed form: sigma=100 MPa, c=5 mm, I=1000 mm^4 -> M=20000 N mm.
    #[test]
    fn moment_for_stress_closed_form() {
        let m = bending_moment_for_stress(100.0, 5.0, 1000.0).unwrap();
        assert!((m - 20_000.0).abs() < 1e-9, "got {m}");
    }

    /// M = sigma * S, tying the inverse to the section modulus.
    #[test]
    fn moment_for_stress_equals_stress_times_section_modulus() {
        let i = second_moment_hollow_circle_mm4(27.0, 17.0).unwrap();
        let c = 13.5;
        let s = section_modulus_mm3(i, c).unwrap();
        let sigma = 150.0;
        let m = bending_moment_for_stress(sigma, c, i).unwrap();
        assert!(
            (m - sigma * s).abs() < 1e-6,
            "M {m} vs sigma*S {}",
            sigma * s
        );
    }

    /// Feeding the ultimate stress gives a fracture moment that, fed back
    /// through the flexure formula, reproduces the ultimate stress exactly.
    #[test]
    fn fracture_moment_brings_outer_fibre_to_ultimate() {
        let i = second_moment_hollow_circle_mm4(27.0, 17.0).unwrap();
        let c = 13.5;
        let sigma_ult = crate::CORTICAL_ULTIMATE_STRESS_MPA;
        let m_frac = bending_moment_for_stress(sigma_ult, c, i).unwrap();
        let sigma_back = bending_stress_mpa(m_frac, c, i).unwrap();
        assert!((sigma_back - sigma_ult).abs() < 1e-9, "got {sigma_back}");
    }

    /// A negative stress yields a sign-flipped moment.
    #[test]
    fn moment_for_stress_sign_carries() {
        let pos = bending_moment_for_stress(120.0, 5.0, 1000.0).unwrap();
        let neg = bending_moment_for_stress(-120.0, 5.0, 1000.0).unwrap();
        assert!(
            (pos + neg).abs() < 1e-9,
            "expected opposite signs: {pos}, {neg}"
        );
    }

    #[test]
    fn moment_for_stress_rejects_bad_inputs() {
        assert!(bending_moment_for_stress(f64::NAN, 5.0, 1000.0).is_err());
        assert!(bending_moment_for_stress(100.0, 0.0, 1000.0).is_err());
        assert!(bending_moment_for_stress(100.0, -1.0, 1000.0).is_err());
        assert!(bending_moment_for_stress(100.0, 5.0, 0.0).is_err());
    }
}
