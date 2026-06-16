//! # valenx-psychrometrics — moist-air psychrometrics
//!
//! Closed-form properties of moist air: saturation vapour pressure,
//! relative humidity, humidity ratio, dew point and specific enthalpy,
//! plus a fully resolved [`MoistAirState`] that ties them together.
//!
//! ## What
//!
//! Three small topic modules and a convenience state struct:
//!
//! 1. [`saturation`] — the Magnus / Tetens [`saturation_pressure`] of
//!    water and its exact closed-form inverse, [`dew_point`].
//! 2. [`moist_air`] — [`relative_humidity`] and its inverse, the
//!    [`humidity_ratio`] `w = 0.622 pv / (p - pv)`, the
//!    [`saturation_humidity_ratio`] `ws` and the
//!    [`degree_of_saturation`] `mu = w / ws`, and the
//!    [`moist_air_enthalpy`] `h = 1.006 T + w (2501 + 1.86 T)`.
//! 3. [`state`] — [`MoistAirState`], a serde-serialisable point that
//!    derives every quantity above consistently from one specification.
//!
//! Every fallible function returns
//! [`Result<_, PsychroError>`](error::PsychroError), whose
//! [`code`](error::PsychroError::code) and
//! [`category`](error::PsychroError::category) accessors are stable for
//! telemetry.
//!
//! ## Model
//!
//! The saturation vapour pressure follows the Magnus–Tetens fit
//! `psat(T) = 610.78 exp(17.27 T / (T + 237.3))` Pa with `T` in degrees
//! Celsius, whose analytic inverse gives the dew point. Relative
//! humidity is `RH = pv / psat(T)`. The humidity ratio comes from the
//! ideal-gas partial-pressure balance with the water-to-air molar-mass
//! ratio `0.622`; evaluating it at the saturation pressure gives the
//! saturation humidity ratio `ws`, and the degree of saturation
//! `mu = w / ws = RH (p - psat) / (p - pv)` is the humidity-ratio
//! analogue of relative humidity (always `<= RH`). The specific enthalpy
//! per kilogram of dry air uses
//! constant specific heats and the `0 degC` latent heat of
//! vaporisation. Each module documents its own equation and references
//! in detail.
//!
//! ```
//! use valenx_psychrometrics::MoistAirState;
//!
//! // Typical indoor air: 22 degC, 50% RH, sea level.
//! let air = MoistAirState::at_sea_level(22.0, 0.50).expect("valid state");
//! assert!(air.dew_point_c < air.dry_bulb_c);
//! println!(
//!     "dew point {:.1} degC, w {:.4} kg/kg, h {:.1} kJ/kg",
//!     air.dew_point_c, air.humidity_ratio, air.enthalpy_kj_per_kg,
//! );
//! ```
//!
//! ## Honest scope
//!
//! Research / educational grade. The formulas here are textbook
//! closed-form and well-established numerical models — the Magnus fit
//! reproduces the standard saturation table to a fraction of a percent,
//! and the humidity-ratio and enthalpy relations are the canonical
//! constant-property approximations used across HVAC teaching. The
//! deliberate simplifications are a single liquid-water saturation
//! curve (no separate ice / sublimation branch below freezing), a
//! fixed `0.622` molar-mass ratio, and temperature-independent specific
//! heats. None of these makes a result meaningless within ordinary
//! ambient conditions, but this crate is NOT a clinical / medical tool
//! and NOT a certified production engineering or metrology reference;
//! for those, defer to a validated property library such as CoolProp or
//! REFPROP.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod moist_air;
pub mod saturation;
pub mod state;

// --- Convenience re-exports of the most-used items --------------------

pub use error::{ErrorCategory, PsychroError, Result};
pub use moist_air::{
    degree_of_saturation, humidity_ratio, moist_air_enthalpy, relative_humidity,
    saturation_humidity_ratio, vapour_pressure_from_rh, STANDARD_PRESSURE_PA,
};
pub use saturation::{dew_point, saturation_pressure};
pub use state::MoistAirState;

#[cfg(test)]
mod tests {
    use super::*;

    /// End-to-end consistency across the whole crate at a hand-checked
    /// point: 25 degC, 60% RH, one standard atmosphere.
    #[test]
    fn end_to_end_consistency() {
        let t = 25.0;
        let rh = 0.60;
        let air = MoistAirState::at_sea_level(t, rh).unwrap();

        // Relative humidity and vapour pressure agree with the topic
        // functions.
        let psat = saturation_pressure(t).unwrap();
        let pv = vapour_pressure_from_rh(rh, t).unwrap();
        assert!((air.vapour_pressure_pa - pv).abs() < 1e-9);
        assert!((air.saturation_pressure_pa - psat).abs() < 1e-9);
        assert!((relative_humidity(pv, t).unwrap() - rh).abs() < 1e-12);

        // Humidity ratio matches the closed form.
        let w = humidity_ratio(pv, STANDARD_PRESSURE_PA).unwrap();
        assert!((air.humidity_ratio - w).abs() < 1e-12);

        // Enthalpy matches the closed form.
        let h = moist_air_enthalpy(t, w).unwrap();
        assert!((air.enthalpy_kj_per_kg - h).abs() < 1e-12);

        // Dew point is below dry-bulb and inverts the saturation curve.
        assert!(air.dew_point_c < t);
        assert!((dew_point(pv).unwrap() - air.dew_point_c).abs() < 1e-12);
    }
}
