//! Physical constants used throughout the crate.
//!
//! Values are the 2018 CODATA recommended constants. The molar gas constant
//! `R` and the Faraday constant `F` are now exact (they derive from the
//! exact `2019` SI definitions of the Boltzmann constant, the Avogadro
//! constant, and the elementary charge), so no measurement uncertainty
//! enters the textbook formulas that consume them.

/// Molar gas constant `R`, in joules per mole-kelvin (`J / (mol K)`).
///
/// Exact under the 2019 SI: `R = N_A k_B`.
pub const GAS_CONSTANT_J_PER_MOL_K: f64 = 8.314_462_618_153_24;

/// Faraday constant `F`, in coulombs per mole (`C / mol`).
///
/// The charge of one mole of elementary charges, `F = N_A e`. Exact under
/// the 2019 SI.
pub const FARADAY_C_PER_MOL: f64 = 96_485.332_123_310_02;

/// Standard reference temperature, in kelvin (25 degrees Celsius).
///
/// The conventional temperature at which standard electrode potentials are
/// tabulated and at which the "0.0592 V per decade" Nernst slope is quoted.
pub const STANDARD_TEMPERATURE_K: f64 = 298.15;

/// Offset between the Celsius and Kelvin scales, in kelvin.
///
/// `T_kelvin = T_celsius + ZERO_CELSIUS_IN_KELVIN`.
pub const ZERO_CELSIUS_IN_KELVIN: f64 = 273.15;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_temperature_is_25_celsius() {
        // 298.15 K must be exactly 25 C above the Celsius zero.
        assert!(
            (STANDARD_TEMPERATURE_K - (25.0 + ZERO_CELSIUS_IN_KELVIN)).abs() < 1e-12,
            "standard T should be 25 C, got {STANDARD_TEMPERATURE_K} K"
        );
    }

    #[test]
    fn rt_over_f_recovers_textbook_thermal_voltage() {
        // R T / F at the standard temperature is the ~0.0257 V scale.
        let vt = GAS_CONSTANT_J_PER_MOL_K * STANDARD_TEMPERATURE_K / FARADAY_C_PER_MOL;
        assert!(
            (vt - 0.025_693).abs() < 1e-5,
            "RT/F should be ~0.0257 V, got {vt}"
        );
    }

    #[test]
    fn constants_are_in_expected_ballpark() {
        // Guard against accidental edits / typos in the literals.
        assert!((GAS_CONSTANT_J_PER_MOL_K - 8.314_462_618).abs() < 1e-6);
        assert!((FARADAY_C_PER_MOL - 96_485.332_12).abs() < 1e-2);
    }
}
