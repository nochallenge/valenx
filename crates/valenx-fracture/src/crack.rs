//! Core linear-elastic fracture-mechanics relations.
//!
//! All quantities are Mode-I (opening mode). The four closed forms are:
//!
//! - **Stress-intensity factor** `K = Y·σ·√(πa)`
//!   ([`stress_intensity`]).
//! - **Critical crack length** `a_c = (1/π)·(K_Ic / (Y·σ))²`
//!   ([`critical_crack_length`]).
//! - **Fracture stress** `σ_f = K_Ic / (Y·√(πa))`
//!   ([`fracture_stress`]).
//! - **Irwin plastic-zone size** `r_p = (1/2π)·(K / σ_y)²`
//!   ([`plastic_zone_radius`]).
//!
//! `Y` is the dimensionless geometry / configuration factor (1.0 for a
//! central through-crack in an infinite plate; ≈1.12 for an edge crack).
//! `σ` is the remote applied tensile stress, `a` the crack length (half
//! the total length of a central crack, the full depth of an edge crack),
//! `K_Ic` the material fracture toughness and `σ_y` the yield strength.

use crate::error::Result;
use crate::material::{require_non_negative_length, require_non_negative_stress, require_positive};
use crate::Material;
use core::f64::consts::PI;

/// Mode-I stress-intensity factor `K = Y·σ·√(πa)`.
///
/// `K` measures the strength of the `1/√r` crack-tip stress singularity.
/// It scales **linearly** with the applied stress `σ` and with the
/// **square root** of the crack length `a`; doubling `a` multiplies `K`
/// by `√2`, doubling `σ` doubles `K`.
///
/// # Arguments
///
/// - `geometry_factor` — dimensionless `Y`, must be finite and `> 0`.
/// - `applied_stress` — remote tension `σ`, must be finite and `>= 0`.
/// - `crack_length` — crack size `a`, must be finite and `>= 0`.
///
/// Returns `K` in the consistent stress·√length unit (e.g. `MPa·√m`).
///
/// # Errors
///
/// [`crate::FractureError::NonPositive`] for a bad `Y`,
/// [`crate::FractureError::InvalidStress`] for a negative / non-finite
/// `σ`, or [`crate::FractureError::InvalidLength`] for a negative /
/// non-finite `a`.
///
/// ```
/// use valenx_fracture::stress_intensity;
/// // Central crack, Y = 1, sigma = 100 MPa, a = 1 mm = 0.001 m.
/// let k = stress_intensity(1.0, 100.0, 0.001).unwrap();
/// // K = 100 * sqrt(pi * 0.001) ~= 5.605 MPa*sqrt(m).
/// assert!((k - 5.604991216_f64).abs() < 1e-6);
/// ```
pub fn stress_intensity(
    geometry_factor: f64,
    applied_stress: f64,
    crack_length: f64,
) -> Result<f64> {
    require_positive("geometry_factor", geometry_factor)?;
    require_non_negative_stress("applied_stress", applied_stress)?;
    require_non_negative_length("crack_length", crack_length)?;
    Ok(geometry_factor * applied_stress * (PI * crack_length).sqrt())
}

/// Critical crack length `a_c = (1/π)·(K_Ic / (Y·σ))²`.
///
/// The crack size at which the stress-intensity factor reaches the
/// material's fracture toughness under the given applied stress — i.e. the
/// `a` that solves `Y·σ·√(πa) = K_Ic`. A larger applied stress shrinks the
/// tolerable flaw (`a_c ∝ 1/σ²`).
///
/// # Arguments
///
/// - `material` — supplies `K_Ic` (already validated `> 0`).
/// - `geometry_factor` — `Y`, must be finite and `> 0`.
/// - `applied_stress` — `σ`, must be finite and `> 0` (it sits in a
///   denominator, so zero is rejected here).
///
/// Returns `a_c` in the consistent length unit.
///
/// # Errors
///
/// [`crate::FractureError::NonPositive`] if `Y` or `σ` is not finite and
/// strictly positive.
///
/// ```
/// use valenx_fracture::{critical_crack_length, Material};
/// let m = Material::new(24.0, 470.0).unwrap();      // K_Ic = 24 MPa*sqrt(m)
/// let a_c = critical_crack_length(&m, 1.0, 200.0).unwrap();
/// // a_c = (24 / 200)^2 / pi ~= 0.004584 m.
/// assert!((a_c - 0.004583662_f64).abs() < 1e-6);
/// ```
pub fn critical_crack_length(
    material: &Material,
    geometry_factor: f64,
    applied_stress: f64,
) -> Result<f64> {
    require_positive("geometry_factor", geometry_factor)?;
    require_positive("applied_stress", applied_stress)?;
    let ratio = material.fracture_toughness / (geometry_factor * applied_stress);
    Ok(ratio * ratio / PI)
}

/// Fracture stress `σ_f = K_Ic / (Y·√(πa))`.
///
/// The remote tension at which a crack of length `a` becomes unstable —
/// the applied stress that drives `K` up to `K_Ic`. A **larger crack
/// lowers** the fracture stress (`σ_f ∝ 1/√a`).
///
/// # Arguments
///
/// - `material` — supplies `K_Ic`.
/// - `geometry_factor` — `Y`, must be finite and `> 0`.
/// - `crack_length` — `a`, must be finite and `> 0` (it sits under a
///   square root in a denominator, so zero is rejected — an uncracked body
///   has no finite fracture stress in this model).
///
/// Returns `σ_f` in the consistent stress unit.
///
/// # Errors
///
/// [`crate::FractureError::NonPositive`] if `Y` or `a` is not finite and
/// strictly positive.
///
/// ```
/// use valenx_fracture::{fracture_stress, Material};
/// let m = Material::new(24.0, 470.0).unwrap();
/// let s = fracture_stress(&m, 1.0, 0.004583662_f64).unwrap();
/// // Round-trips the a_c example above: ~200 MPa.
/// assert!((s - 200.0).abs() < 1e-3);
/// ```
pub fn fracture_stress(
    material: &Material,
    geometry_factor: f64,
    crack_length: f64,
) -> Result<f64> {
    require_positive("geometry_factor", geometry_factor)?;
    require_positive("crack_length", crack_length)?;
    Ok(material.fracture_toughness / (geometry_factor * (PI * crack_length).sqrt()))
}

/// Irwin first-order plastic-zone radius `r_p = (1/2π)·(K / σ_y)²`.
///
/// The size of the crack-tip region in which the elastic stresses would
/// exceed the yield strength (plane-stress, first-order Irwin estimate).
/// It **grows with the square of `K`**: doubling the stress-intensity
/// factor quadruples the plastic-zone size.
///
/// # Arguments
///
/// - `stress_intensity` — `K`, must be finite and `>= 0`.
/// - `material` — supplies the yield strength `σ_y` (already `> 0`).
///
/// Returns `r_p` in the consistent length unit.
///
/// # Errors
///
/// [`crate::FractureError::NonPositive`] if `K` is negative or non-finite.
/// (`K` is treated as a magnitude here, so it must be `>= 0`; a negative
/// value is a programming error and is rejected.)
///
/// ```
/// use valenx_fracture::{plastic_zone_radius, Material};
/// let m = Material::new(24.0, 470.0).unwrap();   // sigma_y = 470 MPa
/// let r = plastic_zone_radius(24.0, &m).unwrap(); // K = K_Ic
/// // r_p = (24/470)^2 / (2 pi) ~= 4.150e-4 m.
/// assert!((r - 4.149988556e-4_f64).abs() < 1e-9);
/// ```
pub fn plastic_zone_radius(stress_intensity: f64, material: &Material) -> Result<f64> {
    // K is used as a magnitude; reuse the positive-quantity contract but
    // allow exactly zero (an unloaded crack has no plastic zone).
    if !(stress_intensity.is_finite() && stress_intensity >= 0.0) {
        return Err(crate::FractureError::NonPositive {
            name: "stress_intensity",
            value: stress_intensity,
        });
    }
    let ratio = stress_intensity / material.yield_strength;
    Ok(ratio * ratio / (2.0 * PI))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn alu() -> Material {
        Material::new(24.0, 470.0).expect("valid alloy")
    }

    // --- K = Y * sigma * sqrt(pi * a) -------------------------------------

    #[test]
    fn k_matches_closed_form() {
        let k = stress_intensity(1.0, 100.0, 0.001).unwrap();
        let expected = 100.0 * (PI * 0.001_f64).sqrt();
        assert!((k - expected).abs() < 1e-12, "k = {k}, expected {expected}");
    }

    #[test]
    fn k_scales_linearly_with_stress() {
        let k1 = stress_intensity(1.0, 100.0, 0.002).unwrap();
        let k2 = stress_intensity(1.0, 300.0, 0.002).unwrap();
        // Tripling the stress triples K.
        assert!(
            (k2 - 3.0 * k1).abs() < 1e-9,
            "k2 = {k2}, 3*k1 = {}",
            3.0 * k1
        );
    }

    #[test]
    fn k_scales_with_sqrt_of_crack_length() {
        let a = 0.0015_f64;
        let k1 = stress_intensity(1.0, 120.0, a).unwrap();
        let k2 = stress_intensity(1.0, 120.0, 4.0 * a).unwrap();
        // Quadrupling a doubles K (sqrt scaling).
        assert!(
            (k2 - 2.0 * k1).abs() < 1e-9,
            "k2 = {k2}, 2*k1 = {}",
            2.0 * k1
        );
    }

    #[test]
    fn k_geometry_factor_is_linear() {
        let base = stress_intensity(1.0, 150.0, 0.003).unwrap();
        let edge = stress_intensity(1.12, 150.0, 0.003).unwrap();
        assert!((edge - 1.12 * base).abs() < 1e-9, "edge = {edge}");
    }

    #[test]
    fn k_is_zero_for_zero_crack() {
        let k = stress_intensity(1.0, 500.0, 0.0).unwrap();
        assert!(k.abs() < 1e-12, "k = {k}");
    }

    #[test]
    fn k_rejects_negative_stress_and_length() {
        assert!(stress_intensity(1.0, -1.0, 0.001).is_err());
        assert!(stress_intensity(1.0, 100.0, -0.001).is_err());
        assert!(stress_intensity(0.0, 100.0, 0.001).is_err());
    }

    // --- fracture when K = K_Ic ------------------------------------------

    #[test]
    fn fracture_occurs_when_k_reaches_toughness() {
        let m = alu();
        let y = 1.0;
        // Drive a crack to its critical length, then evaluate K there: it
        // must equal K_Ic (that is the definition of a_c).
        let sigma = 180.0;
        let a_c = critical_crack_length(&m, y, sigma).unwrap();
        let k_at_ac = stress_intensity(y, sigma, a_c).unwrap();
        assert!(
            (k_at_ac - m.fracture_toughness).abs() < 1e-9,
            "K(a_c) = {k_at_ac}, K_Ic = {}",
            m.fracture_toughness
        );
    }

    #[test]
    fn fracture_stress_drives_k_to_toughness() {
        let m = alu();
        let y = 1.0;
        let a = 0.0025_f64;
        let sigma_f = fracture_stress(&m, y, a).unwrap();
        let k = stress_intensity(y, sigma_f, a).unwrap();
        assert!(
            (k - m.fracture_toughness).abs() < 1e-9,
            "K = {k}, K_Ic = {}",
            m.fracture_toughness
        );
    }

    // --- larger crack lowers fracture stress -----------------------------

    #[test]
    fn larger_crack_lowers_fracture_stress() {
        let m = alu();
        let small = fracture_stress(&m, 1.0, 0.001).unwrap();
        let large = fracture_stress(&m, 1.0, 0.010).unwrap();
        assert!(large < small, "large = {large} should be < small = {small}");
        // sqrt(1/a) scaling: 10x crack -> sqrt(10) lower stress.
        assert!((small / large - 10.0_f64.sqrt()).abs() < 1e-9);
    }

    #[test]
    fn fracture_stress_matches_closed_form() {
        let m = alu();
        let a = 0.005_f64;
        let s = fracture_stress(&m, 1.0, a).unwrap();
        let expected = m.fracture_toughness / (PI * a).sqrt();
        assert!((s - expected).abs() < 1e-9, "s = {s}, expected {expected}");
    }

    // --- a_c round-trips with fracture_stress ----------------------------

    #[test]
    fn critical_crack_length_round_trips() {
        let m = alu();
        let y = 1.0;
        let sigma = 210.0;
        // a_c at this stress, then the fracture stress at that a_c must
        // come back to the original sigma.
        let a_c = critical_crack_length(&m, y, sigma).unwrap();
        let sigma_back = fracture_stress(&m, y, a_c).unwrap();
        assert!(
            (sigma_back - sigma).abs() < 1e-6,
            "sigma_back = {sigma_back}, sigma = {sigma}"
        );
    }

    #[test]
    fn critical_crack_length_matches_closed_form() {
        let m = alu();
        let sigma = 200.0;
        let a_c = critical_crack_length(&m, 1.0, sigma).unwrap();
        let ratio = m.fracture_toughness / sigma;
        let expected = ratio * ratio / PI;
        assert!(
            (a_c - expected).abs() < 1e-12,
            "a_c = {a_c}, expected {expected}"
        );
    }

    #[test]
    fn critical_crack_length_falls_with_stress_squared() {
        let m = alu();
        let a1 = critical_crack_length(&m, 1.0, 100.0).unwrap();
        let a2 = critical_crack_length(&m, 1.0, 200.0).unwrap();
        // Doubling stress quarters the critical crack length.
        assert!((a1 / a2 - 4.0).abs() < 1e-9, "ratio = {}", a1 / a2);
    }

    #[test]
    fn critical_crack_length_rejects_zero_stress() {
        let m = alu();
        assert!(critical_crack_length(&m, 1.0, 0.0).is_err());
        assert!(critical_crack_length(&m, 0.0, 100.0).is_err());
    }

    // --- plastic zone grows with K^2 -------------------------------------

    #[test]
    fn plastic_zone_matches_closed_form() {
        let m = alu();
        let k = 20.0;
        let r = plastic_zone_radius(k, &m).unwrap();
        let ratio = k / m.yield_strength;
        let expected = ratio * ratio / (2.0 * PI);
        assert!((r - expected).abs() < 1e-15, "r = {r}, expected {expected}");
    }

    #[test]
    fn plastic_zone_grows_with_k_squared() {
        let m = alu();
        let r1 = plastic_zone_radius(10.0, &m).unwrap();
        let r2 = plastic_zone_radius(20.0, &m).unwrap();
        // Doubling K quadruples the plastic-zone size.
        assert!((r2 / r1 - 4.0).abs() < 1e-9, "ratio = {}", r2 / r1);
    }

    #[test]
    fn plastic_zone_is_zero_at_zero_k() {
        let m = alu();
        let r = plastic_zone_radius(0.0, &m).unwrap();
        assert!(r.abs() < 1e-15, "r = {r}");
    }

    #[test]
    fn plastic_zone_rejects_negative_k() {
        let m = alu();
        assert!(plastic_zone_radius(-1.0, &m).is_err());
        assert!(plastic_zone_radius(f64::NAN, &m).is_err());
    }

    // --- a small end-to-end consistency check ----------------------------

    #[test]
    fn end_to_end_unit_consistency() {
        // Worked example: 7075-T6, central crack Y = 1, applied 150 MPa.
        let m = alu();
        let y = 1.0;
        let sigma = 150.0;
        let a_c = critical_crack_length(&m, y, sigma).unwrap();
        // K at half the critical length must be below toughness.
        let k_half = stress_intensity(y, sigma, 0.5 * a_c).unwrap();
        assert!(k_half < m.fracture_toughness);
        // K at twice the critical length must exceed toughness.
        let k_double = stress_intensity(y, sigma, 2.0 * a_c).unwrap();
        assert!(k_double > m.fracture_toughness);
    }
}
