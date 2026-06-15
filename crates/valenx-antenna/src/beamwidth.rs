//! First-order beamwidth estimates for aperture antennas.
//!
//! ## Model
//!
//! For an aperture of width `D` (metres) at wavelength `lambda`, the
//! half-power (3 dB) beamwidth scales as
//!
//! ```text
//! theta ~ k * lambda / D    (radians)
//! ```
//!
//! The dimensionless constant `k` depends on the aperture illumination
//! taper. Common textbook values are `k ~ 0.886` for a uniformly
//! illuminated line source (first-null-based `lambda/D` gives `k = 1`)
//! and `k ~ 1.22` for the first null of a uniformly illuminated circular
//! aperture. This module exposes both the bare `lambda / D` estimate and
//! a `k`-weighted form, plus a coarse directivity-from-beamwidth
//! estimate.
//!
//! These are deliberately **first-order** estimates; a real pattern
//! depends on the full aperture distribution. The honest scope note in
//! the crate root applies.

use crate::error::{require_positive, AntennaError};

/// Beamwidth taper constant `k ~ 0.886` for the half-power beamwidth of
/// a uniformly illuminated line / rectangular aperture
/// (`theta_3dB ~ 0.886 * lambda / D`).
pub const K_UNIFORM_LINE_HPBW: f64 = 0.886;

/// Beamwidth constant `k ~ 1.22` for the first-null beamwidth of a
/// uniformly illuminated circular aperture
/// (`theta_null ~ 1.22 * lambda / D`).
pub const K_UNIFORM_CIRCULAR_NULL: f64 = 1.22;

/// First-order beamwidth `theta ~ lambda / D` (radians).
///
/// This is the bare `lambda/D` rule of thumb (`k = 1`).
///
/// # Errors
///
/// Returns an error if `wavelength_m` or `aperture_dim_m` is not finite
/// and strictly positive.
///
/// # Examples
///
/// ```
/// use valenx_antenna::beamwidth::beamwidth_rad;
/// // lambda = 0.1 m over a 1 m aperture -> 0.1 rad.
/// let t = beamwidth_rad(0.1, 1.0).unwrap();
/// assert!((t - 0.1).abs() < 1e-12);
/// ```
pub fn beamwidth_rad(wavelength_m: f64, aperture_dim_m: f64) -> Result<f64, AntennaError> {
    let lambda = require_positive("wavelength_m", wavelength_m)?;
    let d = require_positive("aperture_dim_m", aperture_dim_m)?;
    Ok(lambda / d)
}

/// `k`-weighted beamwidth `theta ~ k * lambda / D` (radians).
///
/// Pass [`K_UNIFORM_LINE_HPBW`] for a uniform line-source half-power
/// beamwidth, [`K_UNIFORM_CIRCULAR_NULL`] for a circular-aperture first
/// null, or any taper constant of your choosing.
///
/// # Errors
///
/// Returns an error if `k_factor`, `wavelength_m` or `aperture_dim_m` is
/// not finite and strictly positive.
pub fn beamwidth_k_rad(
    k_factor: f64,
    wavelength_m: f64,
    aperture_dim_m: f64,
) -> Result<f64, AntennaError> {
    let k = require_positive("k_factor", k_factor)?;
    let theta = beamwidth_rad(wavelength_m, aperture_dim_m)?;
    Ok(k * theta)
}

/// Convenience wrapper for the radians-to-degrees conversion of a
/// beamwidth value.
///
/// # Errors
///
/// Returns an error if `radians` is not finite.
pub fn radians_to_degrees(radians: f64) -> Result<f64, AntennaError> {
    if !radians.is_finite() {
        return Err(AntennaError::NonFinite { name: "radians" });
    }
    Ok(radians.to_degrees())
}

/// Coarse directivity estimate from two orthogonal half-power beamwidths
/// (radians), using the standard approximation
///
/// ```text
/// D0 ~ 4 * pi / (theta_E * theta_H)
/// ```
///
/// where the solid beam angle is approximated by the product of the two
/// principal-plane beamwidths. Returns a **linear** directivity ratio.
///
/// # Errors
///
/// Returns an error if either beamwidth is not finite and strictly
/// positive.
pub fn directivity_from_beamwidths(
    theta_e_rad: f64,
    theta_h_rad: f64,
) -> Result<f64, AntennaError> {
    let te = require_positive("theta_e_rad", theta_e_rad)?;
    let th = require_positive("theta_h_rad", theta_h_rad)?;
    Ok(4.0 * core::f64::consts::PI / (te * th))
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::PI;

    const EPS: f64 = 1e-12;

    #[test]
    fn bare_lambda_over_d() {
        let t = beamwidth_rad(0.1, 1.0).unwrap();
        assert!((t - 0.1).abs() < EPS);
    }

    #[test]
    fn larger_aperture_gives_narrower_beam() {
        let small_ap = beamwidth_rad(0.1, 1.0).unwrap();
        let large_ap = beamwidth_rad(0.1, 10.0).unwrap();
        assert!(large_ap < small_ap);
        // 10x aperture -> 1/10 the beamwidth.
        assert!((large_ap - small_ap / 10.0).abs() < EPS);
    }

    #[test]
    fn beamwidth_scales_with_wavelength() {
        // Longer wavelength -> wider beam, linearly.
        let t1 = beamwidth_rad(0.1, 2.0).unwrap();
        let t2 = beamwidth_rad(0.2, 2.0).unwrap();
        assert!((t2 - 2.0 * t1).abs() < EPS);
    }

    #[test]
    fn k_factor_scales_beamwidth() {
        let bare = beamwidth_rad(0.05, 0.5).unwrap();
        let tapered = beamwidth_k_rad(K_UNIFORM_LINE_HPBW, 0.05, 0.5).unwrap();
        assert!((tapered - K_UNIFORM_LINE_HPBW * bare).abs() < EPS);
    }

    #[test]
    fn circular_null_k_value() {
        // First-null beamwidth of a circular aperture: 1.22 lambda/D.
        let theta = beamwidth_k_rad(K_UNIFORM_CIRCULAR_NULL, 0.03, 1.0).unwrap();
        assert!((theta - 1.22 * 0.03).abs() < EPS);
    }

    #[test]
    fn degrees_conversion() {
        // pi radians -> 180 degrees.
        let deg = radians_to_degrees(PI).unwrap();
        assert!((deg - 180.0).abs() < 1e-9);
    }

    #[test]
    fn beamwidth_in_degrees_example() {
        // 0.1 rad ~ 5.7296 degrees.
        let t = beamwidth_rad(0.1, 1.0).unwrap();
        let deg = radians_to_degrees(t).unwrap();
        assert!((deg - 5.729_577_951).abs() < 1e-6);
    }

    #[test]
    fn directivity_unit_case() {
        // theta_E = theta_H = 1 rad -> D0 = 4*pi.
        let d0 = directivity_from_beamwidths(1.0, 1.0).unwrap();
        assert!((d0 - 4.0 * PI).abs() < EPS);
    }

    #[test]
    fn narrower_beams_give_higher_directivity() {
        let wide = directivity_from_beamwidths(0.5, 0.5).unwrap();
        let narrow = directivity_from_beamwidths(0.1, 0.1).unwrap();
        assert!(narrow > wide);
        // Beamwidths /5 in each plane -> directivity x25.
        assert!((narrow - 25.0 * wide).abs() < 1e-9);
    }

    #[test]
    fn rejects_bad_inputs() {
        assert!(beamwidth_rad(0.0, 1.0).is_err());
        assert!(beamwidth_rad(0.1, 0.0).is_err());
        assert!(beamwidth_k_rad(0.0, 0.1, 1.0).is_err());
        assert!(directivity_from_beamwidths(-1.0, 1.0).is_err());
        assert!(radians_to_degrees(f64::NAN).is_err());
    }
}
