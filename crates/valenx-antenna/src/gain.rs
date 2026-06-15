//! Antenna gain: aperture gain, effective aperture and decibel
//! (dBi) conversions.
//!
//! ## Model
//!
//! For an aperture antenna with effective aperture area `Ae` (square
//! metres) operating at wavelength `lambda`, the maximum gain is the
//! standard aperture relation
//!
//! ```text
//! G = 4 * pi * Ae / lambda^2
//! ```
//!
//! `G` here is a **dimensionless linear power ratio** relative to an
//! isotropic radiator. Expressed in decibels relative to isotropic,
//!
//! ```text
//! G_dBi = 10 * log10(G)
//! ```
//!
//! The inverse `Ae = G * lambda^2 / (4 * pi)` recovers the effective
//! aperture that a given gain implies. For a physical aperture of area
//! `A` with aperture efficiency `eta` (0 < eta <= 1), the effective
//! aperture is `Ae = eta * A`.

use crate::error::{require_non_negative_gain, require_positive, AntennaError};

/// Maximum aperture gain (linear, dimensionless power ratio) from an
/// effective aperture area and wavelength:
/// `G = 4*pi*Ae / lambda^2`.
///
/// # Errors
///
/// Returns an error if `eff_aperture_m2` or `wavelength_m` is not finite
/// and strictly positive.
///
/// # Examples
///
/// ```
/// use valenx_antenna::gain::gain_from_aperture;
/// // 1 m^2 effective aperture at lambda = 1 m -> G = 4*pi.
/// let g = gain_from_aperture(1.0, 1.0).unwrap();
/// assert!((g - 4.0 * std::f64::consts::PI).abs() < 1e-9);
/// ```
pub fn gain_from_aperture(eff_aperture_m2: f64, wavelength_m: f64) -> Result<f64, AntennaError> {
    let ae = require_positive("eff_aperture_m2", eff_aperture_m2)?;
    let lambda = require_positive("wavelength_m", wavelength_m)?;
    Ok(4.0 * core::f64::consts::PI * ae / (lambda * lambda))
}

/// Effective aperture area (square metres) implied by a linear gain and
/// wavelength: `Ae = G * lambda^2 / (4*pi)`.
///
/// Exact inverse of [`gain_from_aperture`].
///
/// # Errors
///
/// Returns an error if `gain_linear` is not finite and non-negative, or
/// if `wavelength_m` is not finite and strictly positive.
pub fn aperture_from_gain(gain_linear: f64, wavelength_m: f64) -> Result<f64, AntennaError> {
    let g = require_non_negative_gain("gain_linear", gain_linear)?;
    let lambda = require_positive("wavelength_m", wavelength_m)?;
    Ok(g * lambda * lambda / (4.0 * core::f64::consts::PI))
}

/// Effective aperture `Ae = eta * A` for a physical aperture area `A`
/// and aperture efficiency `eta`.
///
/// # Errors
///
/// Returns an error if `physical_area_m2` is not finite and strictly
/// positive, or if `efficiency` is outside `(0, 1]` (or non-finite).
pub fn effective_aperture(physical_area_m2: f64, efficiency: f64) -> Result<f64, AntennaError> {
    let a = require_positive("physical_area_m2", physical_area_m2)?;
    let eta = require_positive("efficiency", efficiency)?;
    if eta > 1.0 {
        return Err(AntennaError::NonPositive {
            name: "efficiency",
            value: efficiency,
        });
    }
    Ok(eta * a)
}

/// Convert a linear power ratio to decibels relative to isotropic:
/// `G_dBi = 10 * log10(G)`.
///
/// # Errors
///
/// Returns an error if `gain_linear` is not finite and strictly
/// positive (the logarithm of zero / negative is undefined).
///
/// # Examples
///
/// ```
/// use valenx_antenna::gain::to_dbi;
/// // A linear gain of 100 is 20 dBi.
/// assert!((to_dbi(100.0).unwrap() - 20.0).abs() < 1e-9);
/// ```
pub fn to_dbi(gain_linear: f64) -> Result<f64, AntennaError> {
    let g = require_positive("gain_linear", gain_linear)?;
    Ok(10.0 * g.log10())
}

/// Convert decibels relative to isotropic to a linear power ratio:
/// `G = 10^(G_dBi / 10)`.
///
/// Exact inverse of [`to_dbi`]. Any finite input is valid (negative
/// dBi simply maps to a gain below unity).
///
/// # Errors
///
/// Returns an error only if `gain_dbi` is `NaN` or infinite.
pub fn from_dbi(gain_dbi: f64) -> Result<f64, AntennaError> {
    if !gain_dbi.is_finite() {
        return Err(AntennaError::NonFinite { name: "gain_dbi" });
    }
    Ok(10.0_f64.powf(gain_dbi / 10.0))
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::PI;

    const EPS: f64 = 1e-9;

    #[test]
    fn aperture_gain_unit_case() {
        // Ae = 1, lambda = 1 -> G = 4*pi.
        let g = gain_from_aperture(1.0, 1.0).unwrap();
        assert!((g - 4.0 * PI).abs() < EPS);
    }

    #[test]
    fn gain_scales_linearly_with_aperture() {
        let g1 = gain_from_aperture(1.0, 0.5).unwrap();
        let g2 = gain_from_aperture(2.0, 0.5).unwrap();
        // Doubling the aperture doubles the gain.
        assert!((g2 - 2.0 * g1).abs() < EPS);
    }

    #[test]
    fn larger_aperture_gives_higher_gain() {
        let small = gain_from_aperture(0.5, 0.1).unwrap();
        let large = gain_from_aperture(5.0, 0.1).unwrap();
        assert!(large > small);
    }

    #[test]
    fn gain_scales_with_inverse_wavelength_squared() {
        // Halving lambda quadruples the gain (1/lambda^2).
        let g_long = gain_from_aperture(1.0, 0.2).unwrap();
        let g_short = gain_from_aperture(1.0, 0.1).unwrap();
        assert!((g_short - 4.0 * g_long).abs() < 1e-6);
    }

    #[test]
    fn gain_aperture_roundtrip() {
        let lambda = 0.125;
        for &ae in &[0.01, 0.5, 3.0, 25.0] {
            let g = gain_from_aperture(ae, lambda).unwrap();
            let ae_back = aperture_from_gain(g, lambda).unwrap();
            assert!((ae_back - ae).abs() < 1e-9, "roundtrip failed at Ae={ae}");
        }
    }

    #[test]
    fn dish_gain_sanity() {
        // A 1 m diameter dish, ~55% efficiency, at 10 GHz.
        // lambda = c/f ~ 0.02998 m. Physical area = pi*(0.5)^2.
        let lambda = crate::wave::wavelength_from_frequency(10.0e9).unwrap();
        let phys = PI * 0.5 * 0.5;
        let ae = effective_aperture(phys, 0.55).unwrap();
        let g = gain_from_aperture(ae, lambda).unwrap();
        let g_dbi = to_dbi(g).unwrap();
        // Hand-computed ground truth: Ae = 0.55*pi*0.25 = 0.43197 m^2,
        // lambda = 0.0299792 m, G = 4*pi*Ae/lambda^2 = 6039.78,
        // 10*log10(6039.78) = 37.810 dBi.
        assert!(
            (g_dbi - 37.810_210).abs() < 1e-3,
            "expected ~37.810 dBi, got {g_dbi}"
        );
    }

    #[test]
    fn effective_aperture_applies_efficiency() {
        let ae = effective_aperture(4.0, 0.5).unwrap();
        assert!((ae - 2.0).abs() < EPS);
    }

    #[test]
    fn efficiency_above_one_rejected() {
        assert!(effective_aperture(1.0, 1.5).is_err());
    }

    #[test]
    fn efficiency_of_one_allowed() {
        let ae = effective_aperture(3.0, 1.0).unwrap();
        assert!((ae - 3.0).abs() < EPS);
    }

    #[test]
    fn dbi_of_known_ratios() {
        // 1 -> 0 dBi, 2 -> ~3.0103 dBi, 10 -> 10 dBi, 100 -> 20 dBi.
        assert!((to_dbi(1.0).unwrap() - 0.0).abs() < EPS);
        assert!((to_dbi(2.0).unwrap() - 3.010_299_956).abs() < 1e-6);
        assert!((to_dbi(10.0).unwrap() - 10.0).abs() < EPS);
        assert!((to_dbi(100.0).unwrap() - 20.0).abs() < EPS);
    }

    #[test]
    fn dbi_roundtrip() {
        for &g in &[0.25, 1.0, 4.0, 50.0, 1234.0] {
            let db = to_dbi(g).unwrap();
            let g_back = from_dbi(db).unwrap();
            assert!((g_back - g).abs() / g < 1e-12, "roundtrip failed at {g}");
        }
    }

    #[test]
    fn negative_dbi_is_loss() {
        // -3 dBi -> ~0.5 linear.
        let g = from_dbi(-3.010_299_956).unwrap();
        assert!((g - 0.5).abs() < 1e-6);
    }

    #[test]
    fn to_dbi_rejects_non_positive() {
        assert!(to_dbi(0.0).is_err());
        assert!(to_dbi(-1.0).is_err());
    }
}
