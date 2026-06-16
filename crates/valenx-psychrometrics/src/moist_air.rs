//! Moist-air properties: relative humidity, humidity ratio, enthalpy.
//!
//! ## Model
//!
//! Moist air is treated as an ideal mixture of dry air and water
//! vapour. The defining relations are textbook psychrometrics:
//!
//! Relative humidity is the ratio of the actual partial vapour
//! pressure to the saturation pressure at the same dry-bulb
//! temperature,
//!
//! ```text
//! RH = pv / psat(T)              (dimensionless, 0..1)
//! ```
//!
//! The **humidity ratio** (a.k.a. mixing ratio) is the mass of water
//! vapour per unit mass of dry air. From the ideal-gas partial
//! pressures and the molar-mass ratio of water to dry air
//! (`18.015 / 28.966 ~= 0.621945`, conventionally rounded to `0.622`),
//!
//! ```text
//! w = 0.622 * pv / (p - pv)      [kg water / kg dry air]
//! ```
//!
//! where `p` is the total barometric pressure. As `pv` approaches `p`
//! the denominator vanishes and `w` diverges, which is rejected as an
//! [unphysical](crate::error::PsychroError::Unphysical) state.
//!
//! The **specific enthalpy** of moist air, per kilogram of dry air and
//! referenced to `0 degC` dry air and liquid water, is
//!
//! ```text
//! h = 1.006 * T + w * (2501 + 1.86 * T)      [kJ / kg dry air]
//! ```
//!
//! Here `1.006` kJ/(kg·K) is the specific heat of dry air, `2501`
//! kJ/kg is the latent heat of vaporisation of water at `0 degC`, and
//! `1.86` kJ/(kg·K) is the specific heat of water vapour.
//!
//! ## Honest scope
//!
//! Constant specific heats and a fixed `0.622` molar-mass ratio are the
//! standard engineering approximations and are accurate to well within
//! a percent over normal HVAC conditions. This is a research /
//! educational model; it is not a replacement for a validated property
//! library (CoolProp, REFPROP) or a certified HVAC design tool.

use crate::error::{PsychroError, Result};
use crate::saturation::saturation_pressure;

/// Molar-mass ratio of water vapour to dry air, conventionally rounded
/// to `0.622` in psychrometric practice.
pub const RATIO_MW: f64 = 0.622;

/// Specific heat of dry air, in kJ/(kg·K).
pub const CP_DRY_AIR: f64 = 1.006;

/// Latent heat of vaporisation of water at `0 degC`, in kJ/kg.
pub const LATENT_HEAT_0C: f64 = 2501.0;

/// Specific heat of water vapour, in kJ/(kg·K).
pub const CP_VAPOUR: f64 = 1.86;

/// One standard atmosphere, in pascals. A convenient default total
/// pressure for sea-level calculations.
pub const STANDARD_PRESSURE_PA: f64 = 101_325.0;

/// Validate a total (barometric) pressure argument, in pascals.
fn check_pressure(p_pa: f64) -> Result<()> {
    if !p_pa.is_finite() || p_pa <= 0.0 {
        return Err(PsychroError::bad_parameter(
            "pressure_pa",
            "total pressure must be a positive, finite value",
        ));
    }
    Ok(())
}

/// Validate a partial vapour pressure against the total pressure.
fn check_vapour(pv_pa: f64, p_pa: f64) -> Result<()> {
    if !pv_pa.is_finite() || pv_pa < 0.0 {
        return Err(PsychroError::bad_parameter(
            "pv_pa",
            "vapour pressure must be a non-negative, finite value",
        ));
    }
    if pv_pa >= p_pa {
        return Err(PsychroError::unphysical(format!(
            "vapour pressure {pv_pa} Pa meets or exceeds total pressure {p_pa} Pa"
        )));
    }
    Ok(())
}

/// Relative humidity (a dimensionless fraction in `0..=1`) for an
/// actual partial vapour pressure `pv_pa` at dry-bulb temperature
/// `t_c`.
///
/// `RH = pv / psat(T)`.
///
/// # Errors
///
/// Returns [`PsychroError::BadParameter`] if `pv_pa` is negative or
/// non-finite, or if `t_c` is rejected by [`saturation_pressure`].
///
/// # Examples
///
/// ```
/// use valenx_psychrometrics::moist_air::relative_humidity;
/// use valenx_psychrometrics::saturation::saturation_pressure;
///
/// let t = 25.0;
/// // Half of the saturation vapour pressure gives RH = 0.5.
/// let pv = 0.5 * saturation_pressure(t).unwrap();
/// let rh = relative_humidity(pv, t).unwrap();
/// assert!((rh - 0.5).abs() < 1e-12, "got {rh}");
/// ```
pub fn relative_humidity(pv_pa: f64, t_c: f64) -> Result<f64> {
    if !pv_pa.is_finite() || pv_pa < 0.0 {
        return Err(PsychroError::bad_parameter(
            "pv_pa",
            "vapour pressure must be a non-negative, finite value",
        ));
    }
    let psat = saturation_pressure(t_c)?;
    Ok(pv_pa / psat)
}

/// Actual partial vapour pressure (pascals) from a relative humidity
/// `rh` (a fraction in `0..=1`) at dry-bulb temperature `t_c`.
///
/// The inverse of [`relative_humidity`]: `pv = RH * psat(T)`.
///
/// # Errors
///
/// Returns [`PsychroError::BadParameter`] if `rh` is outside `[0, 1]`
/// or non-finite, or if `t_c` is rejected by [`saturation_pressure`].
pub fn vapour_pressure_from_rh(rh: f64, t_c: f64) -> Result<f64> {
    if !rh.is_finite() || !(0.0..=1.0).contains(&rh) {
        return Err(PsychroError::bad_parameter(
            "rh",
            "relative humidity must be a fraction in [0, 1]",
        ));
    }
    let psat = saturation_pressure(t_c)?;
    Ok(rh * psat)
}

/// Humidity ratio `w` (kilograms of water vapour per kilogram of dry
/// air) for a partial vapour pressure `pv_pa` and total pressure
/// `p_pa`.
///
/// `w = 0.622 * pv / (p - pv)`.
///
/// # Errors
///
/// Returns [`PsychroError::BadParameter`] for a non-positive total
/// pressure or a negative / non-finite vapour pressure, and
/// [`PsychroError::Unphysical`] if `pv_pa >= p_pa` (the ratio would
/// diverge).
///
/// # Examples
///
/// ```
/// use valenx_psychrometrics::moist_air::{humidity_ratio, STANDARD_PRESSURE_PA};
///
/// let w = humidity_ratio(1500.0, STANDARD_PRESSURE_PA).unwrap();
/// assert!(w > 0.0 && w < 0.05, "got {w}");
/// ```
pub fn humidity_ratio(pv_pa: f64, p_pa: f64) -> Result<f64> {
    check_pressure(p_pa)?;
    check_vapour(pv_pa, p_pa)?;
    Ok(RATIO_MW * pv_pa / (p_pa - pv_pa))
}

/// Saturation humidity ratio `ws` (kg water / kg dry air): the humidity
/// ratio of *saturated* air at dry-bulb temperature `t_c` and total
/// pressure `p_pa`,
///
/// `ws = 0.622 * psat(T) / (p - psat(T))`.
///
/// This is the ceiling the actual [`humidity_ratio`] approaches as the
/// air saturates (`RH -> 1`); it is [`humidity_ratio`] evaluated at the
/// [`saturation_pressure`].
///
/// # Errors
///
/// Returns [`PsychroError::BadParameter`] for a non-positive total
/// pressure or a `t_c` rejected by [`saturation_pressure`], and
/// [`PsychroError::Unphysical`] if the saturation pressure meets or
/// exceeds the total pressure (saturated air cannot exist — the air is
/// at or above its boiling point for that pressure).
pub fn saturation_humidity_ratio(t_c: f64, p_pa: f64) -> Result<f64> {
    let psat = saturation_pressure(t_c)?;
    humidity_ratio(psat, p_pa)
}

/// Degree of saturation `mu = w / ws` (dimensionless): the actual
/// humidity ratio `w` as a fraction of the [`saturation_humidity_ratio`]
/// at the same dry-bulb temperature and total pressure.
///
/// Also called the *percentage humidity*, it is the humidity-ratio analogue
/// of relative humidity and is always slightly *below* the relative
/// humidity for unsaturated air: with `RH = pv / psat`,
///
/// ```text
/// mu = RH * (p - psat) / (p - RH * psat) <= RH,
/// ```
/// with equality only at `RH = 0` and `RH = 1`.
///
/// For physical sub-saturated air `mu` lies in `[0, 1]`; a value above
/// one signals a supersaturated input `w > ws`.
///
/// # Errors
///
/// Returns [`PsychroError::BadParameter`] if `w` is negative or
/// non-finite, and propagates the errors of
/// [`saturation_humidity_ratio`].
pub fn degree_of_saturation(w: f64, t_c: f64, p_pa: f64) -> Result<f64> {
    if !w.is_finite() || w < 0.0 {
        return Err(PsychroError::bad_parameter(
            "w",
            "humidity ratio must be a non-negative, finite value",
        ));
    }
    let ws = saturation_humidity_ratio(t_c, p_pa)?;
    Ok(w / ws)
}

/// Specific enthalpy of moist air, in kilojoules per kilogram of dry
/// air, for a dry-bulb temperature `t_c` (degrees Celsius) and a
/// humidity ratio `w` (kg/kg).
///
/// `h = 1.006 * T + w * (2501 + 1.86 * T)`.
///
/// # Errors
///
/// Returns [`PsychroError::BadParameter`] if `w` is negative or
/// non-finite, or if `t_c` is non-finite.
///
/// # Examples
///
/// ```
/// use valenx_psychrometrics::moist_air::moist_air_enthalpy;
///
/// // Dry air at 0 degC has (by reference choice) zero enthalpy.
/// let h = moist_air_enthalpy(0.0, 0.0).unwrap();
/// assert!(h.abs() < 1e-12, "got {h}");
/// ```
pub fn moist_air_enthalpy(t_c: f64, w: f64) -> Result<f64> {
    if !t_c.is_finite() {
        return Err(PsychroError::bad_parameter(
            "t_c",
            "temperature must be finite",
        ));
    }
    if !w.is_finite() || w < 0.0 {
        return Err(PsychroError::bad_parameter(
            "w",
            "humidity ratio must be a non-negative, finite value",
        ));
    }
    Ok(CP_DRY_AIR * t_c + w * (LATENT_HEAT_0C + CP_VAPOUR * t_c))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::saturation::saturation_pressure;

    /// At the saturation vapour pressure the relative humidity is 100%.
    #[test]
    fn rh_is_unity_at_saturation() {
        for t in [0.0, 10.0, 20.0, 30.0, 40.0] {
            let psat = saturation_pressure(t).unwrap();
            let rh = relative_humidity(psat, t).unwrap();
            assert!(
                (rh - 1.0).abs() < 1e-12,
                "RH at saturation wrong at {t}: {rh}"
            );
        }
    }

    /// Relative humidity and vapour pressure are exact inverses.
    #[test]
    fn rh_and_vapour_pressure_round_trip() {
        let t = 22.0;
        for rh in [0.0, 0.1, 0.35, 0.6, 0.85, 1.0] {
            let pv = vapour_pressure_from_rh(rh, t).unwrap();
            let back = relative_humidity(pv, t).unwrap();
            assert!(
                (back - rh).abs() < 1e-12,
                "round trip failed at rh={rh}: {back}"
            );
        }
    }

    /// Humidity ratio increases as vapour pressure rises at fixed total
    /// pressure.
    #[test]
    fn humidity_ratio_rises_with_vapour_pressure() {
        let p = STANDARD_PRESSURE_PA;
        let mut prev = humidity_ratio(100.0, p).unwrap();
        for pv in [500.0, 1000.0, 2000.0, 3000.0, 5000.0] {
            let cur = humidity_ratio(pv, p).unwrap();
            assert!(cur > prev, "w not increasing at pv={pv}: {cur} <= {prev}");
            prev = cur;
        }
    }

    /// Known value: w = 0.622 * pv / (p - pv) reproduces the formula
    /// exactly for a hand-checked case.
    #[test]
    fn humidity_ratio_known_value() {
        // pv = 2000 Pa, p = 100000 Pa -> w = 0.622 * 2000 / 98000.
        let w = humidity_ratio(2000.0, 100_000.0).unwrap();
        let expect = 0.622 * 2000.0 / 98_000.0;
        assert!((w - expect).abs() < 1e-12, "got {w}, expected {expect}");
    }

    /// A vapour pressure at or above the total pressure is unphysical.
    #[test]
    fn humidity_ratio_rejects_vapour_at_or_above_total() {
        let err = humidity_ratio(101_325.0, 101_325.0).unwrap_err();
        assert_eq!(err.code(), "psychro.unphysical");
        assert!(humidity_ratio(120_000.0, 101_325.0).is_err());
    }

    /// Enthalpy rises with temperature at a fixed humidity ratio.
    #[test]
    fn enthalpy_rises_with_temperature() {
        let w = 0.010;
        let mut prev = moist_air_enthalpy(-5.0, w).unwrap();
        for t in [0.0, 10.0, 20.0, 30.0, 40.0] {
            let cur = moist_air_enthalpy(t, w).unwrap();
            assert!(cur > prev, "h not increasing at t={t}: {cur} <= {prev}");
            prev = cur;
        }
    }

    /// Enthalpy rises with humidity ratio at a fixed temperature.
    #[test]
    fn enthalpy_rises_with_humidity_ratio() {
        let t = 25.0;
        let mut prev = moist_air_enthalpy(t, 0.0).unwrap();
        for w in [0.002, 0.005, 0.010, 0.015, 0.020] {
            let cur = moist_air_enthalpy(t, w).unwrap();
            assert!(cur > prev, "h not increasing at w={w}: {cur} <= {prev}");
            prev = cur;
        }
    }

    /// Known value: standard textbook state 25 degC, w = 0.010 gives
    /// roughly 50.6 kJ/kg.
    #[test]
    fn enthalpy_known_value() {
        let h = moist_air_enthalpy(25.0, 0.010).unwrap();
        let expect = 1.006 * 25.0 + 0.010 * (2501.0 + 1.86 * 25.0);
        assert!((h - expect).abs() < 1e-12, "got {h}, expected {expect}");
        // Sanity vs the published ~50.6 kJ/kg figure.
        assert!((h - 50.6).abs() < 0.5, "got {h}");
    }

    #[test]
    fn rejects_bad_rh() {
        assert!(vapour_pressure_from_rh(-0.1, 20.0).is_err());
        assert!(vapour_pressure_from_rh(1.5, 20.0).is_err());
    }

    #[test]
    fn rejects_negative_humidity_ratio_in_enthalpy() {
        assert!(moist_air_enthalpy(20.0, -0.001).is_err());
    }

    #[test]
    fn rejects_nonpositive_total_pressure() {
        assert!(humidity_ratio(1000.0, 0.0).is_err());
        assert!(humidity_ratio(1000.0, -1.0).is_err());
    }

    #[test]
    fn saturation_humidity_ratio_known_value() {
        // 25 degC, one atm -> ws ~= 0.0201 kg/kg (textbook).
        let p = STANDARD_PRESSURE_PA;
        let ws = saturation_humidity_ratio(25.0, p).unwrap();
        let psat = saturation_pressure(25.0).unwrap();
        assert!((ws - humidity_ratio(psat, p).unwrap()).abs() < 1e-15);
        assert!((ws - 0.0201).abs() < 5e-4, "ws = {ws}");
    }

    #[test]
    fn actual_w_never_exceeds_saturation_and_equals_it_at_rh_one() {
        // For RH <= 1 the actual humidity ratio stays at or below ws, and
        // reaches it exactly at saturation.
        let (t, p) = (30.0, STANDARD_PRESSURE_PA);
        let ws = saturation_humidity_ratio(t, p).unwrap();
        for rh in [0.0, 0.25, 0.5, 0.8, 1.0] {
            let pv = vapour_pressure_from_rh(rh, t).unwrap();
            let w = humidity_ratio(pv, p).unwrap();
            assert!(w <= ws + 1e-12, "w {w} > ws {ws} at rh={rh}");
            if (rh - 1.0).abs() < 1e-12 {
                assert!((w - ws).abs() < 1e-12, "w != ws at saturation");
            }
        }
    }

    #[test]
    fn degree_of_saturation_matches_closed_form_and_stays_below_rh() {
        // GOLD identity: mu = RH (p - psat) / (p - pv), and mu <= RH with
        // equality only at the endpoints.
        let (t, p) = (35.0, STANDARD_PRESSURE_PA);
        let psat = saturation_pressure(t).unwrap();
        for rh in [0.1, 0.4, 0.7, 0.95] {
            let pv = vapour_pressure_from_rh(rh, t).unwrap();
            let w = humidity_ratio(pv, p).unwrap();
            let mu = degree_of_saturation(w, t, p).unwrap();
            let expected = rh * (p - psat) / (p - pv);
            assert!(
                (mu - expected).abs() < 1e-12,
                "mu {mu} vs {expected} at rh={rh}"
            );
            assert!(mu < rh, "mu {mu} should be below RH {rh}");
            assert!(mu > 0.0 && mu < 1.0);
        }
    }

    #[test]
    fn degree_of_saturation_endpoints() {
        let (t, p) = (20.0, STANDARD_PRESSURE_PA);
        // RH = 0 -> mu = 0.
        assert!(degree_of_saturation(0.0, t, p).unwrap().abs() < 1e-12);
        // RH = 1 -> w = ws -> mu = 1.
        let ws = saturation_humidity_ratio(t, p).unwrap();
        assert!((degree_of_saturation(ws, t, p).unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn degree_of_saturation_rejects_negative_w() {
        assert!(degree_of_saturation(-0.001, 20.0, STANDARD_PRESSURE_PA).is_err());
        // Non-positive total pressure propagates from ws.
        assert!(degree_of_saturation(0.01, 20.0, 0.0).is_err());
    }
}
