//! Physical constants and the thermal-voltage helper.
//!
//! The single-diode equation is parameterised by the *thermal voltage*
//! `Vt = k*T/q`, where `k` is the Boltzmann constant, `T` the absolute
//! cell temperature in kelvin, and `q` the elementary charge. The
//! constants below are the 2019-redefinition SI exact values, so `Vt` is
//! reproducible to the last bit across platforms.

use crate::error::{Result, SolarPvError};

/// Elementary charge `q`, in coulombs (exact, SI 2019 redefinition).
pub const ELEMENTARY_CHARGE_C: f64 = 1.602_176_634e-19;

/// Boltzmann constant `k`, in joules per kelvin (exact, SI 2019
/// redefinition).
pub const BOLTZMANN_J_PER_K: f64 = 1.380_649e-23;

/// Standard Test Conditions (STC) cell temperature, 25 degrees Celsius,
/// expressed in kelvin. The de-facto datasheet reference point.
pub const STC_TEMPERATURE_K: f64 = 298.15;

/// Standard Test Conditions (STC) plane-of-array irradiance, in watts per
/// square metre. The "one sun" reference used on PV datasheets.
pub const STC_IRRADIANCE_W_PER_M2: f64 = 1000.0;

/// Zero degrees Celsius in kelvin. Convenience for converting datasheet
/// temperatures.
pub const CELSIUS_ZERO_K: f64 = 273.15;

/// Thermal voltage `Vt = k*T/q` in volts for absolute temperature
/// `temperature_k` (kelvin).
///
/// At the STC temperature of 298.15 K this is approximately
/// `0.025_693` V (about 25.7 mV), the textbook room-temperature value.
///
/// # Errors
///
/// Returns [`SolarPvError::Invalid`] if `temperature_k` is not strictly
/// positive (zero or negative absolute temperature is non-physical).
pub fn thermal_voltage(temperature_k: f64) -> Result<f64> {
    if !(temperature_k.is_finite()) || temperature_k <= 0.0 {
        return Err(SolarPvError::invalid(
            "temperature_k",
            format!("absolute temperature must be finite and > 0, got {temperature_k}"),
        ));
    }
    Ok(BOLTZMANN_J_PER_K * temperature_k / ELEMENTARY_CHARGE_C)
}

/// Convert a temperature in degrees Celsius to kelvin.
///
/// This is a pure unit conversion and never fails; the result may be
/// non-physical (<= 0 K) if the input is below absolute zero, which the
/// solver-facing entry points reject downstream.
pub fn celsius_to_kelvin(celsius: f64) -> f64 {
    celsius + CELSIUS_ZERO_K
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Ground truth: at STC the thermal voltage is ~25.7 mV. Computed
    /// directly from the SI exact constants:
    /// 1.380649e-23 * 298.15 / 1.602176634e-19 = 0.0256926... V.
    #[test]
    fn thermal_voltage_at_stc_is_textbook() {
        let vt = thermal_voltage(STC_TEMPERATURE_K).unwrap();
        let expected = BOLTZMANN_J_PER_K * STC_TEMPERATURE_K / ELEMENTARY_CHARGE_C;
        assert!((vt - expected).abs() < 1e-18, "vt = {vt}");
        // And it lands in the well-known ~25.7 mV band.
        assert!((vt - 0.025_692_6).abs() < 1e-6, "vt = {vt}");
    }

    /// Vt is linear in T: doubling the absolute temperature doubles Vt.
    #[test]
    fn thermal_voltage_is_linear_in_temperature() {
        let v1 = thermal_voltage(300.0).unwrap();
        let v2 = thermal_voltage(600.0).unwrap();
        assert!((v2 - 2.0 * v1).abs() < 1e-15, "v1 = {v1}, v2 = {v2}");
    }

    #[test]
    fn thermal_voltage_rejects_non_positive_temperature() {
        assert!(thermal_voltage(0.0).is_err());
        assert!(thermal_voltage(-10.0).is_err());
        assert!(thermal_voltage(f64::NAN).is_err());
        assert!(thermal_voltage(f64::INFINITY).is_err());
    }

    #[test]
    fn celsius_round_trips() {
        assert!((celsius_to_kelvin(25.0) - STC_TEMPERATURE_K).abs() < 1e-12);
        assert!((celsius_to_kelvin(0.0) - CELSIUS_ZERO_K).abs() < 1e-12);
        assert!((celsius_to_kelvin(-273.15)).abs() < 1e-9);
    }
}
