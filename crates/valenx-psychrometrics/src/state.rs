//! A fully resolved moist-air state.
//!
//! ## What
//!
//! [`MoistAirState`] bundles a complete psychrometric point — dry-bulb
//! temperature, total pressure, partial vapour pressure, relative
//! humidity, humidity ratio, dew point, and specific enthalpy — derived
//! consistently from a single specification. The two constructors,
//! [`MoistAirState::from_relative_humidity`] and
//! [`MoistAirState::from_vapour_pressure`], are the usual entry points.
//!
//! ## Model
//!
//! The state simply composes the relations in
//! [`crate::saturation`] and [`crate::moist_air`]: it evaluates the
//! saturation pressure once, derives the vapour pressure from the given
//! relative humidity (or vice versa), and fills in the humidity ratio,
//! dew point and enthalpy. Because the same Magnus curve underlies both
//! the saturation pressure and the dew-point inverse, the invariants
//! `dew_point_c <= dry_bulb_c` and `RH = 100%` exactly at the dew point
//! hold by construction.
//!
//! ## Honest scope
//!
//! Same as the rest of the crate: textbook constant-property
//! psychrometrics, suitable for research and teaching, not a certified
//! HVAC or metrology tool.

use serde::{Deserialize, Serialize};

use crate::error::Result;
use crate::moist_air::{
    humidity_ratio, moist_air_enthalpy, relative_humidity, vapour_pressure_from_rh,
    STANDARD_PRESSURE_PA,
};
use crate::saturation::{dew_point, saturation_pressure};

/// A consistent moist-air psychrometric state.
///
/// All fields are derived together so they never disagree; construct
/// one with [`MoistAirState::from_relative_humidity`] or
/// [`MoistAirState::from_vapour_pressure`] rather than assembling the
/// fields by hand.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MoistAirState {
    /// Dry-bulb temperature, degrees Celsius.
    pub dry_bulb_c: f64,
    /// Total (barometric) pressure, pascals.
    pub pressure_pa: f64,
    /// Saturation vapour pressure at the dry-bulb temperature, pascals.
    pub saturation_pressure_pa: f64,
    /// Actual partial vapour pressure, pascals.
    pub vapour_pressure_pa: f64,
    /// Relative humidity, fraction in `0..=1`.
    pub relative_humidity: f64,
    /// Humidity ratio, kilograms of water vapour per kilogram of dry
    /// air.
    pub humidity_ratio: f64,
    /// Dew-point temperature, degrees Celsius.
    pub dew_point_c: f64,
    /// Specific enthalpy, kilojoules per kilogram of dry air.
    pub enthalpy_kj_per_kg: f64,
}

impl MoistAirState {
    /// Resolve a state from dry-bulb temperature, relative humidity and
    /// total pressure.
    ///
    /// # Errors
    ///
    /// Propagates any [`crate::error::PsychroError`] from the underlying
    /// relations — for example a relative humidity outside `[0, 1]`, a
    /// non-positive pressure, or a temperature below absolute zero.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_psychrometrics::state::MoistAirState;
    ///
    /// let s = MoistAirState::from_relative_humidity(25.0, 0.5, 101_325.0).unwrap();
    /// assert!(s.dew_point_c < s.dry_bulb_c);
    /// assert!(s.humidity_ratio > 0.0);
    /// ```
    pub fn from_relative_humidity(dry_bulb_c: f64, rh: f64, pressure_pa: f64) -> Result<Self> {
        let saturation_pressure_pa = saturation_pressure(dry_bulb_c)?;
        let vapour_pressure_pa = vapour_pressure_from_rh(rh, dry_bulb_c)?;
        Self::assemble(
            dry_bulb_c,
            pressure_pa,
            saturation_pressure_pa,
            vapour_pressure_pa,
        )
    }

    /// Resolve a state from dry-bulb temperature, an actual partial
    /// vapour pressure, and total pressure.
    ///
    /// # Errors
    ///
    /// Propagates any [`crate::error::PsychroError`] from the underlying
    /// relations — including [`crate::error::PsychroError::Unphysical`]
    /// when the vapour pressure meets or exceeds the total pressure.
    pub fn from_vapour_pressure(
        dry_bulb_c: f64,
        vapour_pressure_pa: f64,
        pressure_pa: f64,
    ) -> Result<Self> {
        let saturation_pressure_pa = saturation_pressure(dry_bulb_c)?;
        Self::assemble(
            dry_bulb_c,
            pressure_pa,
            saturation_pressure_pa,
            vapour_pressure_pa,
        )
    }

    /// Resolve a state at one standard atmosphere from dry-bulb
    /// temperature and relative humidity.
    ///
    /// # Errors
    ///
    /// As [`MoistAirState::from_relative_humidity`].
    pub fn at_sea_level(dry_bulb_c: f64, rh: f64) -> Result<Self> {
        Self::from_relative_humidity(dry_bulb_c, rh, STANDARD_PRESSURE_PA)
    }

    /// Internal: fill every derived field from the four primaries,
    /// validating each through the topic functions.
    ///
    /// Perfectly dry air (`vapour_pressure_pa == 0`) is a valid state:
    /// its humidity ratio is zero and it has no condensation
    /// temperature, so the dew point is reported as
    /// [`f64::NEG_INFINITY`] (the limit of the Magnus inverse as the
    /// vapour pressure tends to zero), which still satisfies the
    /// `dew_point_c < dry_bulb_c` invariant.
    fn assemble(
        dry_bulb_c: f64,
        pressure_pa: f64,
        saturation_pressure_pa: f64,
        vapour_pressure_pa: f64,
    ) -> Result<Self> {
        let relative_humidity = relative_humidity(vapour_pressure_pa, dry_bulb_c)?;
        let humidity_ratio = humidity_ratio(vapour_pressure_pa, pressure_pa)?;
        let dew_point_c = if vapour_pressure_pa == 0.0 {
            f64::NEG_INFINITY
        } else {
            dew_point(vapour_pressure_pa)?
        };
        let enthalpy_kj_per_kg = moist_air_enthalpy(dry_bulb_c, humidity_ratio)?;
        Ok(Self {
            dry_bulb_c,
            pressure_pa,
            saturation_pressure_pa,
            vapour_pressure_pa,
            relative_humidity,
            humidity_ratio,
            dew_point_c,
            enthalpy_kj_per_kg,
        })
    }

    /// Whether the air is saturated (relative humidity at or above one,
    /// within `tol`).
    pub fn is_saturated(&self, tol: f64) -> bool {
        self.relative_humidity >= 1.0 - tol
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The dew point of any unsaturated state is strictly below the
    /// dry-bulb temperature; for a saturated state they coincide.
    #[test]
    fn dew_point_never_exceeds_dry_bulb() {
        for &t in &[0.0, 10.0, 20.0, 30.0, 40.0] {
            for &rh in &[0.1, 0.5, 0.9] {
                let s = MoistAirState::from_relative_humidity(t, rh, 101_325.0).unwrap();
                assert!(
                    s.dew_point_c < s.dry_bulb_c,
                    "dew {} >= dry {} at t={t}, rh={rh}",
                    s.dew_point_c,
                    s.dry_bulb_c
                );
            }
            // Saturated: dew point equals dry-bulb.
            let s = MoistAirState::from_relative_humidity(t, 1.0, 101_325.0).unwrap();
            assert!(
                (s.dew_point_c - t).abs() < 1e-9,
                "dew {} != {t}",
                s.dew_point_c
            );
        }
    }

    /// Cooling air to its dew point yields a 100%-RH saturated state.
    #[test]
    fn rh_is_unity_at_dew_point() {
        // Start from a humid but unsaturated state.
        let start = MoistAirState::from_relative_humidity(30.0, 0.6, 101_325.0).unwrap();
        // Cool to the dew point at the same vapour pressure.
        let cooled = MoistAirState::from_vapour_pressure(
            start.dew_point_c,
            start.vapour_pressure_pa,
            101_325.0,
        )
        .unwrap();
        assert!(
            (cooled.relative_humidity - 1.0).abs() < 1e-9,
            "RH at dew point = {}",
            cooled.relative_humidity
        );
        assert!(cooled.is_saturated(1e-6));
    }

    /// At saturation the vapour pressure equals the saturation pressure.
    #[test]
    fn vapour_pressure_equals_saturation_at_full_rh() {
        let s = MoistAirState::from_relative_humidity(25.0, 1.0, 101_325.0).unwrap();
        assert!(
            (s.vapour_pressure_pa - s.saturation_pressure_pa).abs() < 1e-9,
            "pv {} != psat {}",
            s.vapour_pressure_pa,
            s.saturation_pressure_pa
        );
    }

    /// Enthalpy rises with relative humidity at fixed temperature and
    /// pressure (more moisture carries more latent heat).
    #[test]
    fn enthalpy_rises_with_relative_humidity() {
        let mut prev = MoistAirState::at_sea_level(25.0, 0.0)
            .unwrap()
            .enthalpy_kj_per_kg;
        for &rh in &[0.2, 0.4, 0.6, 0.8, 1.0] {
            let cur = MoistAirState::at_sea_level(25.0, rh)
                .unwrap()
                .enthalpy_kj_per_kg;
            assert!(cur > prev, "h not increasing at rh={rh}: {cur} <= {prev}");
            prev = cur;
        }
    }

    /// Perfectly dry air is a valid state: zero humidity ratio, zero
    /// relative humidity, and a dew point at negative infinity (still
    /// below the dry-bulb temperature).
    #[test]
    fn dry_air_is_a_valid_state() {
        let s = MoistAirState::at_sea_level(25.0, 0.0).unwrap();
        assert!(
            (s.relative_humidity).abs() < 1e-12,
            "rh {}",
            s.relative_humidity
        );
        assert!((s.humidity_ratio).abs() < 1e-12, "w {}", s.humidity_ratio);
        assert!(s.vapour_pressure_pa.abs() < 1e-12);
        assert_eq!(s.dew_point_c, f64::NEG_INFINITY);
        assert!(s.dew_point_c < s.dry_bulb_c);
        // Enthalpy reduces to the dry-air term 1.006 * T.
        assert!((s.enthalpy_kj_per_kg - 1.006 * 25.0).abs() < 1e-12);
    }

    /// The state round-trips through serde JSON unchanged.
    #[test]
    fn serde_round_trip() {
        let s = MoistAirState::at_sea_level(20.0, 0.55).unwrap();
        let json = serde_json::to_string(&s).unwrap();
        let back: MoistAirState = serde_json::from_str(&json).unwrap();
        assert_eq!(s, back);
    }

    /// A vapour pressure exceeding the total pressure is rejected.
    #[test]
    fn rejects_supersaturated_vapour_pressure() {
        assert!(MoistAirState::from_vapour_pressure(20.0, 200_000.0, 101_325.0).is_err());
    }
}
