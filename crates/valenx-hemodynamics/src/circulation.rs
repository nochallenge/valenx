//! Whole-circulation lumped relations.
//!
//! These treat the systemic circulation as a single Ohm-analogue
//! resistor driven by the heart as a flow source. They are the standard
//! bedside-physiology identities relating heart rate, stroke volume,
//! cardiac output, systemic vascular resistance and mean arterial
//! pressure.
//!
//! # Models
//!
//! - Cardiac output ([`cardiac_output`]):
//!   `CO = HR * SV` — the heart rate times the per-beat stroke volume.
//! - Mean arterial pressure ([`mean_arterial_pressure`]):
//!   `MAP = CO * SVR` (strictly, `MAP - CVP = CO * SVR`; with the
//!   central venous pressure taken as the zero reference this reduces to
//!   `MAP = CO * SVR`).
//! - Systemic vascular resistance ([`systemic_vascular_resistance`]):
//!   the inverse identity `SVR = MAP / CO`.
//!
//! The functions are unit-agnostic — any self-consistent system works.
//! Two common conventions:
//!
//! - SI: `HR` in 1/s (beats/s), `SV` in m^3, `CO` in m^3/s, `SVR` in
//!   Pa·s/m^3, `MAP` in Pa.
//! - Clinical: `HR` in beats/min, `SV` in mL, `CO` in mL/min and (with
//!   pressure in mmHg) `SVR` in mmHg·min/mL — the so-called "Wood
//!   units" up to the conventional 80 dyn·s·cm^-5 scaling.

use crate::error::{require_non_negative, require_positive};
use crate::HemodynamicsError;

/// Cardiac output from heart rate and stroke volume: `CO = HR * SV`.
///
/// The volume of blood the heart ejects per unit time is the per-beat
/// stroke volume multiplied by the beat frequency.
///
/// # Units
///
/// The result carries the product of the input units: e.g. `HR` in
/// beats/min and `SV` in mL gives `CO` in mL/min; `HR` in 1/s and `SV`
/// in m^3 gives `CO` in m^3/s.
///
/// # Errors
///
/// Returns [`HemodynamicsError::Negative`] / [`HemodynamicsError::NotFinite`]
/// if `heart_rate` or `stroke_volume` is negative or non-finite. Zero
/// is permitted (a stopped heart has zero output).
pub fn cardiac_output(heart_rate: f64, stroke_volume: f64) -> Result<f64, HemodynamicsError> {
    let hr = require_non_negative("heart_rate", heart_rate)?;
    let sv = require_non_negative("stroke_volume", stroke_volume)?;
    Ok(hr * sv)
}

/// Mean arterial pressure from cardiac output and systemic vascular
/// resistance: `MAP = CO * SVR`.
///
/// This is the Ohm analogue `pressure = flow * resistance` applied to
/// the whole systemic circuit, with the venous return pressure taken as
/// the zero reference.
///
/// # Units
///
/// The result carries the product of the input units (e.g. `CO` in
/// m^3/s times `SVR` in Pa·s/m^3 gives `MAP` in Pa).
///
/// # Errors
///
/// Returns [`HemodynamicsError::Negative`] / [`HemodynamicsError::NotFinite`]
/// if `cardiac_output` or `systemic_vascular_resistance` is negative or
/// non-finite.
pub fn mean_arterial_pressure(
    cardiac_output: f64,
    systemic_vascular_resistance: f64,
) -> Result<f64, HemodynamicsError> {
    let co = require_non_negative("cardiac_output", cardiac_output)?;
    let svr = require_non_negative("systemic_vascular_resistance", systemic_vascular_resistance)?;
    Ok(co * svr)
}

/// Systemic vascular resistance from mean arterial pressure and cardiac
/// output: `SVR = MAP / CO`.
///
/// The inverse of [`mean_arterial_pressure`]; recovers the resistance
/// that would produce a given pressure at a given output.
///
/// # Units
///
/// The result carries the ratio of the input units (e.g. `MAP` in Pa
/// over `CO` in m^3/s gives `SVR` in Pa·s/m^3).
///
/// # Errors
///
/// Returns [`HemodynamicsError::NonPositive`] / [`HemodynamicsError::NotFinite`]
/// if `cardiac_output` is not strictly positive (division by zero), or
/// [`HemodynamicsError::Negative`] / [`HemodynamicsError::NotFinite`] if
/// `mean_arterial_pressure` is negative or non-finite.
pub fn systemic_vascular_resistance(
    mean_arterial_pressure: f64,
    cardiac_output: f64,
) -> Result<f64, HemodynamicsError> {
    let map = require_non_negative("mean_arterial_pressure", mean_arterial_pressure)?;
    let co = require_positive("cardiac_output", cardiac_output)?;
    Ok(map / co)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cardiac_output_is_hr_times_sv() {
        // 70 beats/min * 70 mL = 4900 mL/min ~ 4.9 L/min (a normal CO).
        let co = cardiac_output(70.0, 70.0).expect("valid");
        assert!((co - 4900.0).abs() < 1e-9);
    }

    #[test]
    fn cardiac_output_scales_linearly() {
        let co1 = cardiac_output(60.0, 80.0).expect("valid");
        let co2 = cardiac_output(120.0, 80.0).expect("valid");
        assert!((co2 / co1 - 2.0).abs() < 1e-12);
    }

    #[test]
    fn map_equals_co_times_svr() {
        // VALIDATE: MAP = CO * SVR.
        // SI-ish: CO = 9.0e-5 m^3/s (~5.4 L/min), SVR = 1.2e8 Pa·s/m^3.
        let co = 9.0e-5;
        let svr = 1.2e8;
        let map = mean_arterial_pressure(co, svr).expect("valid");
        assert!((map - co * svr).abs() < 1e-9 * (co * svr));
        // ~10800 Pa ~ 81 mmHg, a physiological MAP.
        assert!(map > 0.0);
    }

    #[test]
    fn svr_inverts_map() {
        // VALIDATE: the inverse identity round-trips MAP = CO * SVR.
        let co = 9.0e-5;
        let svr = 1.2e8;
        let map = mean_arterial_pressure(co, svr).expect("valid");
        let svr_back = systemic_vascular_resistance(map, co).expect("valid");
        assert!(
            (svr_back - svr).abs() < 1e-3 * svr,
            "svr={svr}, back={svr_back}"
        );
    }

    #[test]
    fn map_chain_from_hr_sv_svr() {
        // Full chain HR, SV -> CO -> MAP.
        let hr = 1.2; // beats/s (72 bpm)
        let sv = 7.0e-5; // m^3 (70 mL)
        let svr = 1.2e8; // Pa·s/m^3
        let co = cardiac_output(hr, sv).expect("valid");
        let map = mean_arterial_pressure(co, svr).expect("valid");
        let expected = hr * sv * svr;
        assert!((map - expected).abs() < 1e-9 * expected);
    }

    #[test]
    fn doubling_svr_doubles_map_at_fixed_co() {
        let co = 9.0e-5;
        let svr = 1.2e8;
        let map1 = mean_arterial_pressure(co, svr).expect("valid");
        let map2 = mean_arterial_pressure(co, 2.0 * svr).expect("valid");
        assert!((map2 / map1 - 2.0).abs() < 1e-12);
    }

    #[test]
    fn zero_output_gives_zero_pressure() {
        let map = mean_arterial_pressure(0.0, 1.2e8).expect("valid");
        assert!(map.abs() < 1e-12);
    }

    #[test]
    fn invalid_inputs_are_rejected() {
        assert!(cardiac_output(-1.0, 70.0).is_err());
        assert!(cardiac_output(70.0, f64::NAN).is_err());
        assert!(mean_arterial_pressure(-1.0, 1.0).is_err());
        // SVR requires a strictly-positive CO (no divide-by-zero).
        assert!(systemic_vascular_resistance(100.0, 0.0).is_err());
        assert!(systemic_vascular_resistance(-1.0, 5.0).is_err());
    }
}
