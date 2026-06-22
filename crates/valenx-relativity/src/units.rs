//! Conversions between geometrized (`G = c = 1`) units and SI, with named
//! physical constants, so the dimensionless geometric results can be expressed
//! for real black holes (kilometres, kelvin, joules-per-kelvin, years).
//!
//! Constants are CODATA/IAU values.

use std::f64::consts::PI;

/// Newton's gravitational constant `G` (m³·kg⁻¹·s⁻²).
pub const G: f64 = 6.674_30e-11;
/// Speed of light `c` (m·s⁻¹).
pub const C: f64 = 2.997_924_58e8;
/// Reduced Planck constant `ħ` (J·s).
pub const HBAR: f64 = 1.054_571_817e-34;
/// Boltzmann constant `k_B` (J·K⁻¹).
pub const K_B: f64 = 1.380_649e-23;
/// Solar mass `M_☉` (kg).
pub const SOLAR_MASS: f64 = 1.988_92e30;
/// Seconds in a Julian year.
pub const YEAR_SECONDS: f64 = 3.155_76e7;

/// Gravitational length `GM/c²` (metres) of a mass given in kilograms.
pub fn mass_to_length_m(mass_kg: f64) -> f64 {
    G * mass_kg / (C * C)
}

/// Schwarzschild radius `2GM/c²` in kilometres, for a mass in solar masses
/// (≈ 2.95 km for the Sun).
pub fn schwarzschild_radius_km(mass_solar: f64) -> f64 {
    2.0 * mass_to_length_m(mass_solar * SOLAR_MASS) / 1000.0
}

/// Hawking temperature of a Schwarzschild hole in kelvin, mass in solar masses:
/// `T = ħc³ / (8π G M k_B)` (≈ 6.17 × 10⁻⁸ K for one solar mass).
pub fn hawking_temperature_kelvin(mass_solar: f64) -> f64 {
    HBAR * C.powi(3) / (8.0 * PI * G * mass_solar * SOLAR_MASS * K_B)
}

/// Bekenstein–Hawking entropy of a Schwarzschild hole in J/K, mass in solar
/// masses: `S = 4π G M² k_B / (ħ c)`.
pub fn entropy_si(mass_solar: f64) -> f64 {
    let m = mass_solar * SOLAR_MASS;
    4.0 * PI * G * m * m * K_B / (HBAR * C)
}

/// Evaporation time of a Schwarzschild hole in years, mass in solar masses:
/// the standard Page estimate `t = 5120 π G² M³ / (ħ c⁴)` (≈ 2.1 × 10⁶⁷ years
/// for one solar mass). Photon-dominated; greybody detail folded into the
/// prefactor.
pub fn evaporation_time_years(mass_solar: f64) -> f64 {
    let m = mass_solar * SOLAR_MASS;
    let t_sec = 5120.0 * PI * G * G * m * m * m / (HBAR * C.powi(4));
    t_sec / YEAR_SECONDS
}
