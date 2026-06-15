//! The thermal description of a human body.
//!
//! A [`Body`] is the lumped-capacitance object that stores heat: a
//! single mass at a single (core) temperature, with a specific heat
//! and a skin-surface area through which it exchanges heat with the
//! environment. This is the classic "one-node" thermoregulation model
//! — the body is treated as one well-mixed thermal mass, exactly as in
//! an introductory heat-transfer or physiology text. It is *not* a
//! multi-segment (e.g. Stolwijk / Fiala) model.

use crate::error::{non_negative, positive, ThermoregError};
use serde::{Deserialize, Serialize};

/// Approximate specific heat capacity of the human body,
/// `3492 J kg^-1 K^-1` (~0.83 kcal kg^-1 °C^-1).
///
/// The whole-body average is dominated by the high water content of
/// tissue; this is the standard textbook value.
pub const BODY_SPECIFIC_HEAT: f64 = 3492.0;

/// The lumped thermal description of a body.
///
/// All fields are SI. The struct is validated on construction via
/// [`Body::new`]; the public fields are read-only by convention (use a
/// constructor to make a new, re-validated `Body`).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Body {
    /// Body mass (kg). Must be strictly positive.
    pub mass_kg: f64,
    /// Specific heat capacity (J kg^-1 K^-1). Must be strictly
    /// positive; see [`BODY_SPECIFIC_HEAT`].
    pub specific_heat: f64,
    /// Skin / DuBois surface area (m^2). Must be strictly positive.
    pub surface_area_m2: f64,
    /// Core (deep-body) temperature (°C).
    pub core_temp_c: f64,
    /// Mean skin temperature (°C).
    pub skin_temp_c: f64,
}

impl Body {
    /// Construct a validated [`Body`].
    ///
    /// # Errors
    ///
    /// Returns [`ThermoregError::NotFinite`] if any argument is `NaN`
    /// or infinite, and [`ThermoregError::OutOfRange`] if `mass_kg`,
    /// `specific_heat`, or `surface_area_m2` is not strictly positive.
    /// Temperatures may be any finite value.
    pub fn new(
        mass_kg: f64,
        specific_heat: f64,
        surface_area_m2: f64,
        core_temp_c: f64,
        skin_temp_c: f64,
    ) -> Result<Self, ThermoregError> {
        Ok(Self {
            mass_kg: positive("mass_kg", mass_kg)?,
            specific_heat: positive("specific_heat", specific_heat)?,
            surface_area_m2: positive("surface_area_m2", surface_area_m2)?,
            core_temp_c: crate::error::finite("core_temp_c", core_temp_c)?,
            skin_temp_c: crate::error::finite("skin_temp_c", skin_temp_c)?,
        })
    }

    /// The DuBois body-surface-area estimate (m^2),
    /// `A = 0.007184 * mass^0.425 * height^0.725`, with `mass` in kg
    /// and `height` in cm.
    ///
    /// This is the 1916 DuBois & DuBois formula, still the most-cited
    /// BSA estimator. For a 70 kg, 170 cm adult it gives ≈ 1.81 m^2.
    ///
    /// # Errors
    ///
    /// Returns an error if `mass_kg` or `height_cm` is not strictly
    /// positive (or is non-finite).
    pub fn dubois_surface_area(mass_kg: f64, height_cm: f64) -> Result<f64, ThermoregError> {
        let m = positive("mass_kg", mass_kg)?;
        let h = positive("height_cm", height_cm)?;
        Ok(0.007184 * m.powf(0.425) * h.powf(0.725))
    }

    /// Build a body for a "standard" adult from mass and height,
    /// using the [`Body::dubois_surface_area`] BSA and the textbook
    /// [`BODY_SPECIFIC_HEAT`].
    ///
    /// A convenience over [`Body::new`]; core and skin temperatures
    /// are supplied by the caller.
    ///
    /// # Errors
    ///
    /// Propagates any validation error from
    /// [`Body::dubois_surface_area`] / [`Body::new`].
    pub fn standard_adult(
        mass_kg: f64,
        height_cm: f64,
        core_temp_c: f64,
        skin_temp_c: f64,
    ) -> Result<Self, ThermoregError> {
        let area = Self::dubois_surface_area(mass_kg, height_cm)?;
        Self::new(mass_kg, BODY_SPECIFIC_HEAT, area, core_temp_c, skin_temp_c)
    }

    /// Whole-body heat capacity (J K^-1), `C = m * c`.
    ///
    /// The amount of heat that raises the body temperature by 1 K.
    pub fn heat_capacity(&self) -> f64 {
        self.mass_kg * self.specific_heat
    }

    /// Return a copy of this body with its core temperature shifted by
    /// `delta_c` kelvin (equivalently °C).
    ///
    /// Used to advance the lumped-capacitance state by an integration
    /// step; see [`crate::balance`].
    pub fn with_core_shifted(&self, delta_c: f64) -> Self {
        Self {
            core_temp_c: self.core_temp_c + delta_c,
            ..*self
        }
    }

    /// The temperature gradient driving sensible heat loss,
    /// `T_skin - T_air`, given an air temperature (°C).
    ///
    /// Positive when the skin is warmer than the air (the usual case,
    /// heat flows *out* of the body).
    pub fn skin_minus_air(&self, air_temp_c: f64) -> f64 {
        self.skin_temp_c - air_temp_c
    }
}

/// A sweat (evaporative) state: how fast sweat evaporates and how much
/// energy that carries away per unit mass.
///
/// Kept separate from [`Body`] because the sweat rate is a *control
/// action*, not a fixed body property — the thermoregulatory system
/// modulates it.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Sweat {
    /// Mass rate of sweat *evaporated* (kg s^-1). Must be ≥ 0.
    ///
    /// Note: this is evaporated sweat, not secreted sweat — sweat that
    /// drips off the body without evaporating removes no heat.
    pub evap_rate_kg_s: f64,
    /// Latent heat of vaporisation of sweat (J kg^-1). Must be > 0.
    ///
    /// See [`LATENT_HEAT_SWEAT`] for the standard skin-temperature
    /// value.
    pub latent_heat_j_kg: f64,
}

/// Latent heat of vaporisation of water/sweat at skin temperature
/// (~30 °C), `2_426_000 J kg^-1` (≈ 2426 kJ kg^-1 ≈ 0.58 kcal g^-1).
///
/// Slightly higher than the 100 °C value (2257 kJ kg^-1) because
/// evaporation happens at skin temperature.
pub const LATENT_HEAT_SWEAT: f64 = 2_426_000.0;

impl Sweat {
    /// Construct a validated [`Sweat`].
    ///
    /// # Errors
    ///
    /// [`ThermoregError::NotFinite`] for non-finite inputs;
    /// [`ThermoregError::OutOfRange`] if `evap_rate_kg_s` is negative
    /// or `latent_heat_j_kg` is not strictly positive.
    pub fn new(evap_rate_kg_s: f64, latent_heat_j_kg: f64) -> Result<Self, ThermoregError> {
        Ok(Self {
            evap_rate_kg_s: non_negative("evap_rate_kg_s", evap_rate_kg_s)?,
            latent_heat_j_kg: positive("latent_heat_j_kg", latent_heat_j_kg)?,
        })
    }

    /// A sweat state from an evaporated rate alone, using the textbook
    /// [`LATENT_HEAT_SWEAT`].
    ///
    /// # Errors
    ///
    /// Propagates validation from [`Sweat::new`].
    pub fn from_rate(evap_rate_kg_s: f64) -> Result<Self, ThermoregError> {
        Self::new(evap_rate_kg_s, LATENT_HEAT_SWEAT)
    }

    /// No sweating — a zero evaporative rate.
    pub fn none() -> Self {
        Self {
            evap_rate_kg_s: 0.0,
            latent_heat_j_kg: LATENT_HEAT_SWEAT,
        }
    }

    /// Evaporative heat-loss power (W), `E = m_dot * L`.
    ///
    /// The energy carried away per second by evaporating sweat. Larger
    /// sweat rates remove proportionally more heat.
    pub fn evaporative_power(&self) -> f64 {
        self.evap_rate_kg_s * self.latent_heat_j_kg
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    fn approx(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn body_rejects_non_positive_mass_and_area() {
        assert!(Body::new(0.0, BODY_SPECIFIC_HEAT, 1.8, 37.0, 33.0).is_err());
        assert!(Body::new(-1.0, BODY_SPECIFIC_HEAT, 1.8, 37.0, 33.0).is_err());
        assert!(Body::new(70.0, BODY_SPECIFIC_HEAT, 0.0, 37.0, 33.0).is_err());
        assert!(Body::new(70.0, 0.0, 1.8, 37.0, 33.0).is_err());
    }

    #[test]
    fn body_accepts_below_zero_temperatures() {
        // A frigid skin temperature is physically meaningful as input.
        let b = Body::new(70.0, BODY_SPECIFIC_HEAT, 1.8, 37.0, -5.0).unwrap();
        assert!(approx(b.skin_temp_c, -5.0, EPS));
    }

    #[test]
    fn heat_capacity_is_mass_times_specific_heat() {
        let b = Body::new(70.0, 3492.0, 1.8, 37.0, 33.0).unwrap();
        // 70 * 3492 = 244_440 J/K
        assert!(approx(b.heat_capacity(), 244_440.0, 1e-6));
    }

    #[test]
    fn dubois_matches_known_reference_value() {
        // DuBois for the canonical 70 kg, 170 cm adult ≈ 1.81 m^2.
        let a = Body::dubois_surface_area(70.0, 170.0).unwrap();
        assert!(approx(a, 1.81, 0.02), "got {a}");
    }

    #[test]
    fn dubois_standard_man_71point7_kg_172cm() {
        // The DuBois "standard man" (71.7 kg, 172 cm) is defined to be
        // ~1.85 m^2 — that is what the formula was fit to.
        let a = Body::dubois_surface_area(71.7, 172.0).unwrap();
        assert!(approx(a, 1.85, 0.02), "got {a}");
    }

    #[test]
    fn standard_adult_uses_dubois_and_textbook_specific_heat() {
        let b = Body::standard_adult(70.0, 170.0, 37.0, 33.0).unwrap();
        assert!(approx(b.specific_heat, BODY_SPECIFIC_HEAT, EPS));
        let a = Body::dubois_surface_area(70.0, 170.0).unwrap();
        assert!(approx(b.surface_area_m2, a, EPS));
    }

    #[test]
    fn with_core_shifted_moves_only_core() {
        let b = Body::new(70.0, 3492.0, 1.8, 37.0, 33.0).unwrap();
        let b2 = b.with_core_shifted(0.5);
        assert!(approx(b2.core_temp_c, 37.5, EPS));
        assert!(approx(b2.skin_temp_c, 33.0, EPS));
        assert!(approx(b2.mass_kg, 70.0, EPS));
    }

    #[test]
    fn skin_minus_air_is_the_sensible_gradient() {
        let b = Body::new(70.0, 3492.0, 1.8, 37.0, 33.0).unwrap();
        assert!(approx(b.skin_minus_air(20.0), 13.0, EPS));
        // Hotter air than skin reverses the sign (heat flows in).
        assert!(b.skin_minus_air(40.0) < 0.0);
    }

    #[test]
    fn sweat_rejects_negative_rate_but_allows_zero() {
        assert!(Sweat::new(-0.1, LATENT_HEAT_SWEAT).is_err());
        assert!(Sweat::new(0.0, LATENT_HEAT_SWEAT).is_ok());
        assert!(Sweat::new(1e-4, 0.0).is_err());
    }

    #[test]
    fn evaporative_power_is_rate_times_latent_heat() {
        // 1 g/min evaporated = 1e-3/60 kg/s; * 2_426_000 J/kg ≈ 40.4 W.
        let s = Sweat::from_rate(1.0e-3 / 60.0).unwrap();
        assert!(
            approx(s.evaporative_power(), 40.433, 1e-3),
            "got {}",
            s.evaporative_power()
        );
    }

    #[test]
    fn higher_sweat_rate_gives_more_cooling() {
        let lo = Sweat::from_rate(1.0e-5).unwrap();
        let hi = Sweat::from_rate(2.0e-5).unwrap();
        assert!(hi.evaporative_power() > lo.evaporative_power());
        // Exactly proportional: double the rate, double the power.
        assert!(approx(
            hi.evaporative_power(),
            2.0 * lo.evaporative_power(),
            1e-6
        ));
    }

    #[test]
    fn no_sweat_removes_no_heat() {
        assert!(approx(Sweat::none().evaporative_power(), 0.0, EPS));
    }
}
