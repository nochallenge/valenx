//! Saturation vapour pressure and dew point (Magnus / Tetens form).
//!
//! ## Model
//!
//! The saturation vapour pressure of water over a flat liquid surface
//! is given by the **Magnus–Tetens** empirical formula
//!
//! ```text
//! psat(T) = C * exp( A * T / (T + B) )      [Pa],  T in degrees Celsius
//! ```
//!
//! with the widely used coefficients
//!
//! ```text
//! A = 17.27,   B = 237.3 degC,   C = 610.78 Pa
//! ```
//!
//! These reproduce the standard saturation table to within a fraction
//! of a percent across the meteorological range (roughly 0–60 degC).
//! At the triple point `T = 0.01 degC` the formula returns ~611 Pa,
//! the textbook value.
//!
//! ## Dew point
//!
//! Because the Magnus form is analytically invertible, the **dew
//! point** — the temperature to which moist air must be cooled (at
//! constant pressure and humidity ratio) for the water vapour to begin
//! to condense — has a closed form. Given an actual partial vapour
//! pressure `pv`, set `gamma = ln(pv / C)`; then
//!
//! ```text
//! Td = B * gamma / (A - gamma)      [degC]
//! ```
//!
//! This inverse is *exact* with respect to [`saturation_pressure`]:
//! feeding `psat(T)` back through [`dew_point`] recovers `T`, so the
//! dew point of fully saturated air equals its temperature.
//!
//! ## Honest scope
//!
//! This is the single-curve Magnus fit over *liquid* water. It does not
//! switch to an ice (sublimation) curve below freezing, and it carries
//! the few-tenths-of-a-percent error inherent in the empirical
//! coefficients versus a full formulation such as IAPWS. It is a
//! research / educational model, not a metrology reference.

use crate::error::{PsychroError, Result};

/// Magnus dimensionless coefficient `A`.
pub const MAGNUS_A: f64 = 17.27;

/// Magnus temperature offset `B`, in degrees Celsius.
pub const MAGNUS_B: f64 = 237.3;

/// Magnus reference pressure `C`, in pascals (saturation pressure at
/// `0 degC`).
pub const MAGNUS_C: f64 = 610.78;

/// Absolute-zero temperature in degrees Celsius. Used to reject
/// physically impossible inputs.
pub const ABSOLUTE_ZERO_C: f64 = -273.15;

/// Saturation vapour pressure of water over a liquid surface, in
/// pascals, for a dry-bulb temperature `t_c` in degrees Celsius.
///
/// Implements the Magnus–Tetens formula described in the
/// [module documentation](self).
///
/// # Errors
///
/// Returns [`PsychroError::BadParameter`] if `t_c` is below absolute
/// zero or is not finite.
///
/// # Examples
///
/// ```
/// use valenx_psychrometrics::saturation::saturation_pressure;
///
/// // ~611 Pa at the triple point.
/// let p = saturation_pressure(0.01).unwrap();
/// assert!((p - 611.0).abs() < 2.0, "got {p}");
/// ```
pub fn saturation_pressure(t_c: f64) -> Result<f64> {
    if !t_c.is_finite() {
        return Err(PsychroError::bad_parameter(
            "t_c",
            "temperature must be finite",
        ));
    }
    if t_c <= ABSOLUTE_ZERO_C {
        return Err(PsychroError::bad_parameter(
            "t_c",
            format!("temperature {t_c} degC is at or below absolute zero"),
        ));
    }
    Ok(MAGNUS_C * (MAGNUS_A * t_c / (t_c + MAGNUS_B)).exp())
}

/// Dew-point temperature in degrees Celsius for an actual partial
/// vapour pressure `pv_pa` (pascals).
///
/// Implements the closed-form inverse of [`saturation_pressure`]
/// described in the [module documentation](self).
///
/// # Errors
///
/// Returns [`PsychroError::BadParameter`] if `pv_pa` is not strictly
/// positive or is not finite. (A zero or negative vapour pressure has
/// no dew point: the air is perfectly dry.)
///
/// # Examples
///
/// ```
/// use valenx_psychrometrics::saturation::{dew_point, saturation_pressure};
///
/// // Dew point of saturated air equals its temperature.
/// let t = 25.0_f64;
/// let pv = saturation_pressure(t).unwrap();
/// let td = dew_point(pv).unwrap();
/// assert!((td - t).abs() < 1e-9, "got {td}");
/// ```
pub fn dew_point(pv_pa: f64) -> Result<f64> {
    if !pv_pa.is_finite() {
        return Err(PsychroError::bad_parameter(
            "pv_pa",
            "vapour pressure must be finite",
        ));
    }
    if pv_pa <= 0.0 {
        return Err(PsychroError::bad_parameter(
            "pv_pa",
            "vapour pressure must be strictly positive to have a dew point",
        ));
    }
    let gamma = (pv_pa / MAGNUS_C).ln();
    Ok(MAGNUS_B * gamma / (MAGNUS_A - gamma))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Saturation pressure at the triple point is ~611 Pa.
    #[test]
    fn psat_at_zero_celsius() {
        let p = saturation_pressure(0.0).unwrap();
        // Magnus reference value at exactly 0 degC.
        assert!((p - MAGNUS_C).abs() < 1e-9, "got {p}");
    }

    /// Known reference: psat(20 degC) is close to the textbook 2.34 kPa.
    #[test]
    fn psat_at_twenty_celsius_matches_table() {
        let p = saturation_pressure(20.0).unwrap();
        // ASHRAE table value ~2339 Pa; Magnus fit lands within ~10 Pa.
        assert!((p - 2339.0).abs() < 15.0, "got {p}");
    }

    /// Known reference: psat(100 degC) is close to one standard
    /// atmosphere (boiling).
    #[test]
    fn psat_at_hundred_celsius_near_atmospheric() {
        let p = saturation_pressure(100.0).unwrap();
        // Magnus over-/under-shoots at the high end; loose tolerance.
        assert!((p - 101_325.0).abs() < 4000.0, "got {p}");
    }

    /// Saturation pressure increases monotonically with temperature.
    #[test]
    fn psat_monotonic_in_temperature() {
        let mut prev = saturation_pressure(-10.0).unwrap();
        for t in [-5.0, 0.0, 10.0, 20.0, 30.0, 40.0, 50.0] {
            let cur = saturation_pressure(t).unwrap();
            assert!(cur > prev, "psat not increasing at {t}: {cur} <= {prev}");
            prev = cur;
        }
    }

    /// The dew point is the exact inverse of the saturation curve.
    #[test]
    fn dew_point_inverts_saturation() {
        for t in [-5.0, 0.0, 5.0, 15.0, 25.0, 35.0, 45.0] {
            let pv = saturation_pressure(t).unwrap();
            let td = dew_point(pv).unwrap();
            assert!((td - t).abs() < 1e-9, "round trip failed at {t}: {td}");
        }
    }

    /// Dew point rises with vapour pressure (more moisture, warmer dew).
    #[test]
    fn dew_point_monotonic_in_vapour_pressure() {
        let mut prev = dew_point(200.0).unwrap();
        for pv in [400.0, 800.0, 1200.0, 2000.0, 3000.0] {
            let cur = dew_point(pv).unwrap();
            assert!(
                cur > prev,
                "dew point not increasing at {pv}: {cur} <= {prev}"
            );
            prev = cur;
        }
    }

    #[test]
    fn rejects_below_absolute_zero() {
        assert!(saturation_pressure(-300.0).is_err());
    }

    #[test]
    fn rejects_nonpositive_vapour_pressure() {
        assert!(dew_point(0.0).is_err());
        assert!(dew_point(-5.0).is_err());
    }
}
