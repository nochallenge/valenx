//! Stefan-Boltzmann law: total hemispherical emissive power of a black
//! body and its inverse.
//!
//! The total emissive power radiated per unit area by an ideal black
//! body at absolute temperature `T` (kelvin) is
//!
//! ```text
//! Eb = sigma * T^4          [W / m^2]
//! ```
//!
//! where `sigma` is the Stefan-Boltzmann constant. A real (gray)
//! surface of emissivity `eps` radiates `E = eps * Eb`; that scaling
//! lives in the [`graybody`](crate::graybody) module.
//!
//! Ground-truth checks in the tests use the CODATA value
//! `sigma = 5.670374419e-8 W m^-2 K^-4` and the textbook round number
//! `Eb(1000 K) = 56.7 kW/m^2` (Incropera, *Fundamentals of Heat and
//! Mass Transfer*).

use crate::error::{check_temperature, finite, RadiationError, Result};

/// Stefan-Boltzmann constant, `sigma`, in `W m^-2 K^-4`.
///
/// CODATA 2018 recommended value. It is exact in SI by definition of
/// the kelvin since the 2019 redefinition:
/// `sigma = 2 pi^5 k^4 / (15 h^3 c^2)`.
pub const SIGMA: f64 = 5.670_374_419e-8;

/// Total hemispherical black-body emissive power `Eb = sigma T^4`.
///
/// # Arguments
///
/// `temperature_k` — absolute temperature in kelvin (`>= 0`).
///
/// # Returns
///
/// Emissive power in `W / m^2`.
///
/// # Errors
///
/// Returns [`RadiationError::OutOfDomain`] for a negative temperature
/// and [`RadiationError::NotFinite`] for a non-finite one.
///
/// # Example
///
/// ```
/// use valenx_radiation::stefan_boltzmann::black_body_emissive_power;
/// let eb = black_body_emissive_power(1000.0).unwrap();
/// // Textbook value ~56.7 kW/m^2.
/// assert!((eb - 56_703.74).abs() < 1.0);
/// ```
pub fn black_body_emissive_power(temperature_k: f64) -> Result<f64> {
    let t = check_temperature("temperature_k", temperature_k)?;
    Ok(SIGMA * t.powi(4))
}

/// Gray-surface total emissive power `E = eps * sigma T^4`.
///
/// A convenience over [`black_body_emissive_power`] that applies the
/// total hemispherical emissivity `eps` of a diffuse gray surface.
///
/// # Arguments
///
/// `emissivity` — total hemispherical emissivity in `[0, 1]`;
/// `temperature_k` — absolute temperature in kelvin (`>= 0`).
///
/// # Errors
///
/// Propagates the domain checks on both arguments.
pub fn gray_emissive_power(emissivity: f64, temperature_k: f64) -> Result<f64> {
    let eps = crate::error::check_unit_interval("emissivity", emissivity)?;
    let eb = black_body_emissive_power(temperature_k)?;
    Ok(eps * eb)
}

/// Inverse Stefan-Boltzmann law: the black-body temperature whose
/// emissive power equals `eb`, i.e. `T = (eb / sigma)^(1/4)`.
///
/// This is the radiative (effective) temperature corresponding to a
/// measured hemispherical emissive flux.
///
/// # Arguments
///
/// `eb` — emissive power in `W / m^2` (`>= 0`).
///
/// # Returns
///
/// Absolute temperature in kelvin.
///
/// # Errors
///
/// Returns [`RadiationError::OutOfDomain`] for a negative flux and
/// [`RadiationError::NotFinite`] for a non-finite one.
pub fn temperature_from_emissive_power(eb: f64) -> Result<f64> {
    let e = finite("eb", eb)?;
    if e < 0.0 {
        return Err(RadiationError::OutOfDomain {
            what: "eb",
            value: e,
            reason: "emissive power must be >= 0",
        });
    }
    Ok((e / SIGMA).powf(0.25))
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn zero_kelvin_radiates_nothing() {
        assert!(black_body_emissive_power(0.0).unwrap().abs() < EPS);
    }

    #[test]
    fn matches_textbook_1000k() {
        // Incropera round number: Eb(1000 K) ~= 56.7 kW/m^2.
        // Exact analytic ground truth: sigma * 1000^4 = sigma * 1e12.
        let eb = black_body_emissive_power(1000.0).unwrap();
        let truth = SIGMA * 1.0e12;
        assert!((eb - truth).abs() < 1e-6, "got {eb}");
        // And within 1 W/m^2 of the published rounded figure.
        assert!((eb - 56_703.74).abs() < 1.0, "got {eb}");
    }

    #[test]
    fn quartic_scaling_doubling_temperature() {
        // Eb(2T) / Eb(T) == 16 exactly for the T^4 law.
        let lo = black_body_emissive_power(300.0).unwrap();
        let hi = black_body_emissive_power(600.0).unwrap();
        assert!((hi / lo - 16.0).abs() < 1e-9, "ratio {}", hi / lo);
    }

    #[test]
    fn hand_worked_300k() {
        // 300^4 = 8.1e9; Eb = 5.670374419e-8 * 8.1e9 = 459.3003... W/m^2.
        let eb = black_body_emissive_power(300.0).unwrap();
        let truth = 5.670_374_419e-8 * 8.1e9;
        assert!((eb - truth).abs() < 1e-9, "got {eb}");
        assert!((eb - 459.300_3).abs() < 1e-3, "got {eb}");
    }

    #[test]
    fn inverse_round_trips() {
        for &t in &[1.0, 50.0, 300.0, 1234.5, 5778.0] {
            let eb = black_body_emissive_power(t).unwrap();
            let back = temperature_from_emissive_power(eb).unwrap();
            assert!((back - t).abs() < 1e-6, "t={t} back={back}");
        }
    }

    #[test]
    fn solar_effective_temperature() {
        // The Sun's surface emits ~6.32e7 W/m^2; effective T ~5778 K.
        let t = temperature_from_emissive_power(6.32e7).unwrap();
        assert!((t - 5778.0).abs() < 5.0, "got {t}");
    }

    #[test]
    fn gray_surface_scales_by_emissivity() {
        let eb = black_body_emissive_power(500.0).unwrap();
        let e = gray_emissive_power(0.4, 500.0).unwrap();
        assert!((e - 0.4 * eb).abs() < 1e-9, "got {e}");
        // Black body (eps = 1) equals Eb.
        let bb = gray_emissive_power(1.0, 500.0).unwrap();
        assert!((bb - eb).abs() < EPS, "got {bb}");
    }

    #[test]
    fn rejects_bad_inputs() {
        assert!(black_body_emissive_power(-1.0).is_err());
        assert!(black_body_emissive_power(f64::NAN).is_err());
        assert!(temperature_from_emissive_power(-1.0).is_err());
        assert!(gray_emissive_power(1.5, 300.0).is_err());
        assert!(gray_emissive_power(-0.1, 300.0).is_err());
    }
}
