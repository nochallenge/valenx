//! Wave relations: frequency, wavelength and the free-space wavenumber.
//!
//! These are the elementary relations that tie an antenna's operating
//! frequency to the wavelength used throughout the link-budget formulas.
//!
//! ## Model
//!
//! For a plane wave in free space (or any non-dispersive medium with
//! phase velocity `v`), frequency `f` and wavelength `lambda` satisfy
//!
//! ```text
//! v = f * lambda
//! ```
//!
//! so `lambda = v / f` and `f = v / lambda`. In free space `v = c`, the
//! speed of light. The angular wavenumber is `k = 2*pi / lambda`.

use crate::error::{require_positive, AntennaError};

/// Speed of light in vacuum, in metres per second.
///
/// Exact by the 2019 SI definition of the metre / second:
/// `c = 299_792_458 m/s`.
pub const SPEED_OF_LIGHT_M_S: f64 = 299_792_458.0;

/// Free-space wavelength `lambda = c / f` for a frequency `freq_hz`
/// (hertz), returned in metres.
///
/// # Errors
///
/// Returns an error if `freq_hz` is not finite and strictly positive.
///
/// # Examples
///
/// ```
/// use valenx_antenna::wave::wavelength_from_frequency;
/// // 300 MHz <-> 1 metre (a standard rule-of-thumb).
/// let lambda = wavelength_from_frequency(299_792_458.0).unwrap();
/// assert!((lambda - 1.0).abs() < 1e-9);
/// ```
pub fn wavelength_from_frequency(freq_hz: f64) -> Result<f64, AntennaError> {
    let f = require_positive("freq_hz", freq_hz)?;
    Ok(SPEED_OF_LIGHT_M_S / f)
}

/// Wavelength `lambda = v / f` for an arbitrary (non-dispersive) phase
/// velocity `velocity_m_s`, in metres.
///
/// Use this when the antenna is embedded in a medium whose phase
/// velocity differs from `c` (a dielectric, a waveguide, ...).
///
/// # Errors
///
/// Returns an error if either argument is not finite and strictly
/// positive.
pub fn wavelength_in_medium(velocity_m_s: f64, freq_hz: f64) -> Result<f64, AntennaError> {
    let v = require_positive("velocity_m_s", velocity_m_s)?;
    let f = require_positive("freq_hz", freq_hz)?;
    Ok(v / f)
}

/// Free-space frequency `f = c / lambda` for a wavelength
/// `wavelength_m` (metres), returned in hertz.
///
/// This is the exact inverse of [`wavelength_from_frequency`].
///
/// # Errors
///
/// Returns an error if `wavelength_m` is not finite and strictly
/// positive.
pub fn frequency_from_wavelength(wavelength_m: f64) -> Result<f64, AntennaError> {
    let lambda = require_positive("wavelength_m", wavelength_m)?;
    Ok(SPEED_OF_LIGHT_M_S / lambda)
}

/// Angular wavenumber `k = 2*pi / lambda` (radians per metre) for a
/// wavelength `wavelength_m`.
///
/// # Errors
///
/// Returns an error if `wavelength_m` is not finite and strictly
/// positive.
pub fn wavenumber(wavelength_m: f64) -> Result<f64, AntennaError> {
    let lambda = require_positive("wavelength_m", wavelength_m)?;
    Ok(2.0 * core::f64::consts::PI / lambda)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn three_hundred_mhz_is_one_metre() {
        // 300 MHz -> ~1 m. Using exact c, f = c gives exactly 1 m only
        // at f = c Hz; check the familiar engineering value at 3e8 Hz.
        let lambda = wavelength_from_frequency(3.0e8).unwrap();
        // c / 3e8 = 0.99930819... m
        assert!((lambda - SPEED_OF_LIGHT_M_S / 3.0e8).abs() < EPS);
        assert!((lambda - 0.999_308_193).abs() < 1e-6);
    }

    #[test]
    fn one_ghz_wavelength() {
        // 1 GHz -> ~0.2998 m.
        let lambda = wavelength_from_frequency(1.0e9).unwrap();
        assert!((lambda - 0.299_792_458).abs() < EPS);
    }

    #[test]
    fn two_point_four_ghz_wifi() {
        // 2.4 GHz Wi-Fi -> ~12.49 cm.
        let lambda = wavelength_from_frequency(2.4e9).unwrap();
        assert!((lambda - 0.124_913_524).abs() < 1e-6);
    }

    #[test]
    fn frequency_wavelength_roundtrip() {
        for &f in &[1.0e6, 100.0e6, 2.4e9, 28.0e9] {
            let lambda = wavelength_from_frequency(f).unwrap();
            let f_back = frequency_from_wavelength(lambda).unwrap();
            assert!((f_back - f).abs() / f < 1e-12, "roundtrip failed at {f}");
        }
    }

    #[test]
    fn medium_velocity_matches_free_space_at_c() {
        let a = wavelength_in_medium(SPEED_OF_LIGHT_M_S, 1.0e9).unwrap();
        let b = wavelength_from_frequency(1.0e9).unwrap();
        assert!((a - b).abs() < EPS);
    }

    #[test]
    fn slower_medium_gives_shorter_wavelength() {
        // Phase velocity halved (e.g. relative permittivity 4) ->
        // wavelength halved at the same frequency.
        let free = wavelength_from_frequency(1.0e9).unwrap();
        let slow = wavelength_in_medium(SPEED_OF_LIGHT_M_S / 2.0, 1.0e9).unwrap();
        assert!((slow - free / 2.0).abs() < EPS);
    }

    #[test]
    fn wavenumber_of_one_metre() {
        let k = wavenumber(1.0).unwrap();
        assert!((k - 2.0 * core::f64::consts::PI).abs() < EPS);
    }

    #[test]
    fn wavenumber_scales_inversely_with_wavelength() {
        let k1 = wavenumber(2.0).unwrap();
        let k2 = wavenumber(1.0).unwrap();
        // Halving the wavelength doubles k.
        assert!((k2 - 2.0 * k1).abs() < EPS);
    }

    #[test]
    fn rejects_non_positive_frequency() {
        assert!(wavelength_from_frequency(0.0).is_err());
        assert!(wavelength_from_frequency(-1.0).is_err());
        assert!(wavelength_from_frequency(f64::NAN).is_err());
    }
}
