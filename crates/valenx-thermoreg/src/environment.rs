//! The thermal environment surrounding the body.
//!
//! An [`Environment`] bundles the air temperature, the mean radiant
//! temperature of the surroundings, and the convective heat-transfer
//! coefficient. These feed the sensible-heat (convection + radiation)
//! terms of the balance in [`crate::balance`].

use crate::error::{non_negative, positive, ThermoregError};
use serde::{Deserialize, Serialize};

/// Stefan-Boltzmann constant, `5.670374419e-8 W m^-2 K^-4`.
pub const STEFAN_BOLTZMANN: f64 = 5.670_374_419e-8;

/// Conversion offset between Celsius and Kelvin, `273.15`.
pub const KELVIN_OFFSET: f64 = 273.15;

/// Typical long-wave emissivity of human skin / clothing, `0.95`
/// (dimensionless, 0–1). Skin is very nearly a black body in the
/// infrared.
pub const SKIN_EMISSIVITY: f64 = 0.95;

/// The surrounding thermal environment.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Environment {
    /// Ambient (dry-bulb) air temperature (°C).
    pub air_temp_c: f64,
    /// Mean radiant temperature of the surrounding surfaces (°C) — the
    /// area-weighted temperature the body "sees" radiatively. Often
    /// equal to `air_temp_c` indoors, but a cold window or a hot stove
    /// pulls it away.
    pub mean_radiant_temp_c: f64,
    /// Convective heat-transfer coefficient `h` (W m^-2 K^-1). Must be
    /// ≥ 0. Rises with air speed; ~3–4 for still air, tens in wind.
    pub convective_coeff: f64,
    /// Long-wave emissivity of the body surface (0–1).
    pub emissivity: f64,
}

impl Environment {
    /// Construct a validated [`Environment`].
    ///
    /// # Errors
    ///
    /// [`ThermoregError::NotFinite`] for non-finite inputs;
    /// [`ThermoregError::OutOfRange`] if `convective_coeff` is negative
    /// or `emissivity` is outside `(0, 1]`.
    pub fn new(
        air_temp_c: f64,
        mean_radiant_temp_c: f64,
        convective_coeff: f64,
        emissivity: f64,
    ) -> Result<Self, ThermoregError> {
        let emissivity = positive("emissivity", emissivity)?;
        if emissivity > 1.0 {
            return Err(ThermoregError::OutOfRange {
                name: "emissivity",
                value: emissivity,
                reason: "must be within (0, 1]",
            });
        }
        Ok(Self {
            air_temp_c: crate::error::finite("air_temp_c", air_temp_c)?,
            mean_radiant_temp_c: crate::error::finite("mean_radiant_temp_c", mean_radiant_temp_c)?,
            convective_coeff: non_negative("convective_coeff", convective_coeff)?,
            emissivity,
        })
    }

    /// A "still indoor air" environment: air and radiant temperatures
    /// equal, `h = 3.5` (free convection over a clothed person), skin
    /// emissivity [`SKIN_EMISSIVITY`].
    ///
    /// # Errors
    ///
    /// Propagates validation from [`Environment::new`].
    pub fn still_indoor(air_temp_c: f64) -> Result<Self, ThermoregError> {
        Self::new(air_temp_c, air_temp_c, 3.5, SKIN_EMISSIVITY)
    }

    /// Convective heat loss per unit area (W m^-2) from a surface at
    /// `skin_temp_c`, by Newton's law of cooling
    /// `q'' = h * (T_skin - T_air)`.
    ///
    /// Positive when the skin is warmer than the air.
    pub fn convective_flux(&self, skin_temp_c: f64) -> f64 {
        self.convective_coeff * (skin_temp_c - self.air_temp_c)
    }

    /// Total convective heat loss (W) from a surface of area
    /// `area_m2`, `Q = h * A * (T_skin - T_air)`.
    pub fn convective_power(&self, skin_temp_c: f64, area_m2: f64) -> f64 {
        self.convective_flux(skin_temp_c) * area_m2
    }

    /// Net long-wave radiative heat loss per unit area (W m^-2),
    /// `q'' = ε * σ * (T_skin^4 - T_radiant^4)` with both temperatures
    /// in kelvin.
    ///
    /// Positive when the skin is radiatively warmer than its
    /// surroundings.
    pub fn radiative_flux(&self, skin_temp_c: f64) -> f64 {
        let ts = skin_temp_c + KELVIN_OFFSET;
        let tr = self.mean_radiant_temp_c + KELVIN_OFFSET;
        self.emissivity * STEFAN_BOLTZMANN * (ts.powi(4) - tr.powi(4))
    }

    /// Total net long-wave radiative heat loss (W) from a surface of
    /// area `area_m2`.
    pub fn radiative_power(&self, skin_temp_c: f64, area_m2: f64) -> f64 {
        self.radiative_flux(skin_temp_c) * area_m2
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn env_rejects_bad_emissivity_and_negative_h() {
        assert!(Environment::new(20.0, 20.0, -1.0, 0.95).is_err());
        assert!(Environment::new(20.0, 20.0, 3.5, 0.0).is_err());
        assert!(Environment::new(20.0, 20.0, 3.5, 1.5).is_err());
        assert!(Environment::new(20.0, 20.0, 3.5, 1.0).is_ok());
    }

    #[test]
    fn convective_flux_is_h_times_delta_t() {
        let env = Environment::new(20.0, 20.0, 4.0, 0.95).unwrap();
        // 4 * (33 - 20) = 52 W/m^2.
        assert!(approx(env.convective_flux(33.0), 52.0, 1e-9));
    }

    #[test]
    fn convective_power_scales_with_area() {
        let env = Environment::new(20.0, 20.0, 4.0, 0.95).unwrap();
        // 52 W/m^2 * 1.8 m^2 = 93.6 W.
        assert!(approx(env.convective_power(33.0, 1.8), 93.6, 1e-9));
    }

    #[test]
    fn convection_reverses_sign_when_air_hotter_than_skin() {
        let env = Environment::new(40.0, 40.0, 4.0, 0.95).unwrap();
        assert!(env.convective_flux(33.0) < 0.0);
    }

    #[test]
    fn no_convective_loss_at_thermal_equilibrium() {
        let env = Environment::new(33.0, 33.0, 4.0, 0.95).unwrap();
        assert!(approx(env.convective_flux(33.0), 0.0, 1e-12));
    }

    #[test]
    fn radiative_flux_matches_stefan_boltzmann_hand_calc() {
        // Skin 33 °C = 306.15 K, surroundings 20 °C = 293.15 K, ε = 1.
        // q'' = σ (306.15^4 - 293.15^4).
        let env = Environment::new(20.0, 20.0, 3.5, 1.0).unwrap();
        let ts = 306.15_f64;
        let tr = 293.15_f64;
        let expected = STEFAN_BOLTZMANN * (ts.powi(4) - tr.powi(4));
        assert!(approx(env.radiative_flux(33.0), expected, 1e-9));
        // Sanity: a ~13 K gradient near room temperature radiates on
        // the order of 70–80 W/m^2.
        assert!(env.radiative_flux(33.0) > 60.0 && env.radiative_flux(33.0) < 90.0);
    }

    #[test]
    fn radiative_flux_zero_when_skin_equals_radiant() {
        let env = Environment::new(20.0, 30.0, 3.5, 0.95).unwrap();
        assert!(approx(env.radiative_flux(30.0), 0.0, 1e-12));
    }

    #[test]
    fn emissivity_scales_radiation_linearly() {
        let full = Environment::new(20.0, 20.0, 3.5, 1.0).unwrap();
        let half = Environment::new(20.0, 20.0, 3.5, 0.5).unwrap();
        assert!(approx(
            half.radiative_flux(33.0),
            0.5 * full.radiative_flux(33.0),
            1e-9
        ));
    }

    #[test]
    fn still_indoor_sets_radiant_equal_to_air() {
        let env = Environment::still_indoor(22.0).unwrap();
        assert!(approx(env.air_temp_c, env.mean_radiant_temp_c, 1e-12));
        assert!(approx(env.emissivity, SKIN_EMISSIVITY, 1e-12));
    }
}
