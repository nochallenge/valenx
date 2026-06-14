//! Larson-Miller time-temperature parameter for stress-rupture life.
//!
//! ## Model
//!
//! The Larson-Miller parameter (LMP) collapses the family of
//! stress-rupture curves measured at many temperatures and many
//! exposure times onto a single master curve of stress versus LMP. For
//! a temperature `T` (kelvin) and a time to rupture `t_r` (hours) it is
//!
//! ```text
//!   LMP = T * (C + log10(t_r))
//! ```
//!
//! where `C` is the Larson-Miller constant — a material-dependent fit
//! parameter that is `~20` for many ferrous alloys (the classic default
//! Larson and Miller proposed in 1952). The grouping says that the same
//! amount of creep damage (and therefore rupture) is reached either
//! quickly at high temperature or slowly at low temperature, along a
//! constant-LMP locus.
//!
//! Because a material's master curve fixes LMP for a given applied
//! stress, the relation inverts to predict the **time to rupture** at
//! any service temperature:
//!
//! ```text
//!   log10(t_r) = LMP / T - C
//!   t_r        = 10^(LMP / T - C)
//! ```
//!
//! This is the standard engineering use: read LMP off the master curve
//! at the design stress, then solve for the life at the operating
//! temperature.
//!
//! ## Honest scope
//!
//! These are the textbook closed-form Larson-Miller relations. They do
//! not themselves contain any material data — you must supply LMP (or a
//! measured rupture point) and the constant `C` from a qualified data
//! source. The parameter is an empirical extrapolation tool; real
//! life-assessment must account for scatter, multiaxiality, oxidation,
//! microstructural evolution and the validated bounds of the master
//! curve. Research / educational grade only.

use crate::error::{require_finite, require_positive, CreepError};
use serde::{Deserialize, Serialize};

/// The classic default Larson-Miller constant for many ferrous /
/// low-alloy steels (Larson & Miller, 1952). Always prefer a value
/// fitted to the specific material when one is available.
pub const DEFAULT_C: f64 = 20.0;

/// A single stress-rupture data point: a temperature and the observed
/// (or predicted) time to rupture at that temperature.
///
/// Construct one with [`RupturePoint::new`], which validates that both
/// quantities are finite and strictly positive.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RupturePoint {
    /// Absolute temperature (kelvin).
    pub temperature_k: f64,
    /// Time to rupture (hours).
    pub time_hours: f64,
}

impl RupturePoint {
    /// Build a validated rupture point.
    ///
    /// # Errors
    ///
    /// Returns [`CreepError`] if either argument is non-finite or not
    /// strictly positive. A temperature in kelvin and a rupture time in
    /// hours are both physically positive.
    pub fn new(temperature_k: f64, time_hours: f64) -> Result<Self, CreepError> {
        Ok(Self {
            temperature_k: require_positive("temperature_k", temperature_k)?,
            time_hours: require_positive("time_hours", time_hours)?,
        })
    }

    /// The Larson-Miller parameter of this point for the constant `c`.
    ///
    /// Convenience wrapper around [`larson_miller_parameter`] using this
    /// point's own temperature and time.
    ///
    /// # Errors
    ///
    /// Returns [`CreepError`] if `c` is non-finite.
    pub fn parameter(&self, c: f64) -> Result<f64, CreepError> {
        larson_miller_parameter(self.temperature_k, self.time_hours, c)
    }
}

/// Evaluate the Larson-Miller parameter `LMP = T * (C + log10(t_r))`.
///
/// `temperature_k` is `T` in kelvin, `time_hours` is the rupture time
/// `t_r` in hours, and `c` is the Larson-Miller constant.
///
/// # Errors
///
/// Returns [`CreepError`] if `temperature_k` or `time_hours` is
/// non-finite or not strictly positive, or if `c` is non-finite.
///
/// # Examples
///
/// ```
/// use valenx_creep::larson_miller::larson_miller_parameter;
///
/// // 1 hour at 1000 K with C = 20: log10(1) = 0, so LMP = 1000 * 20.
/// let lmp = larson_miller_parameter(1000.0, 1.0, 20.0).unwrap();
/// assert!((lmp - 20_000.0).abs() < 1e-9);
/// ```
pub fn larson_miller_parameter(
    temperature_k: f64,
    time_hours: f64,
    c: f64,
) -> Result<f64, CreepError> {
    let temperature_k = require_positive("temperature_k", temperature_k)?;
    let time_hours = require_positive("time_hours", time_hours)?;
    let c = require_finite("c", c)?;
    Ok(temperature_k * (c + time_hours.log10()))
}

/// Solve the Larson-Miller relation for the time to rupture.
///
/// Inverts `LMP = T * (C + log10(t_r))` to give
/// `t_r = 10^(LMP / T - C)` hours at temperature `temperature_k`
/// (kelvin) for the constant `c`.
///
/// This is the headline engineering query: given the LMP read off a
/// material's master curve at the design stress, how long will the part
/// last at the operating temperature.
///
/// # Errors
///
/// Returns [`CreepError`] if `temperature_k` is non-finite or not
/// strictly positive, or if `lmp` or `c` is non-finite.
///
/// # Examples
///
/// ```
/// use valenx_creep::larson_miller::rupture_time_hours;
///
/// // LMP = 20000 at 1000 K with C = 20 inverts to log10(t) = 0 → 1 h.
/// let t = rupture_time_hours(20_000.0, 1000.0, 20.0).unwrap();
/// assert!((t - 1.0).abs() < 1e-9);
/// ```
pub fn rupture_time_hours(lmp: f64, temperature_k: f64, c: f64) -> Result<f64, CreepError> {
    let lmp = require_finite("lmp", lmp)?;
    let temperature_k = require_positive("temperature_k", temperature_k)?;
    let c = require_finite("c", c)?;
    let log10_t = lmp / temperature_k - c;
    Ok(10f64.powf(log10_t))
}

/// Solve the Larson-Miller relation for the service temperature at
/// which a given life is reached.
///
/// Inverts the relation for `T`: from `LMP = T * (C + log10(t_r))`,
///
/// ```text
///   T = LMP / (C + log10(t_r))
/// ```
///
/// kelvin, for a target rupture time `time_hours` and constant `c`.
///
/// # Errors
///
/// Returns [`CreepError`] if `time_hours` is non-finite or not strictly
/// positive, if `lmp` or `c` is non-finite, or if the denominator
/// `C + log10(t_r)` is zero (the temperature would be undefined).
pub fn rupture_temperature_k(lmp: f64, time_hours: f64, c: f64) -> Result<f64, CreepError> {
    let lmp = require_finite("lmp", lmp)?;
    let time_hours = require_positive("time_hours", time_hours)?;
    let c = require_finite("c", c)?;
    let denom = c + time_hours.log10();
    if denom == 0.0 {
        return Err(CreepError::InvalidValue {
            name: "c + log10(time_hours)",
            value: denom,
            reason: "denominator is zero; temperature is undefined",
        });
    }
    Ok(lmp / denom)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for floating-point comparisons. The closed forms are
    /// exact up to round-off, so a tight epsilon is appropriate.
    const EPS: f64 = 1e-9;

    #[test]
    fn lmp_matches_defining_formula() {
        // Ground truth: directly evaluate T * (C + log10(t)).
        let t = 1100.0;
        let time = 1000.0; // log10(1000) = 3 exactly.
        let c = 20.0;
        let expected = t * (c + 3.0);
        let got = larson_miller_parameter(t, time, c).unwrap();
        assert!(
            (got - expected).abs() < EPS,
            "LMP mismatch: got {got}, expected {expected}"
        );
        // And 1100 * 23 = 25300.
        assert!((got - 25_300.0).abs() < 1e-6);
    }

    #[test]
    fn lmp_with_unit_time_is_t_times_c() {
        // log10(1) = 0, so LMP collapses to T * C.
        let lmp = larson_miller_parameter(950.0, 1.0, DEFAULT_C).unwrap();
        assert!((lmp - 950.0 * DEFAULT_C).abs() < EPS, "got {lmp}");
    }

    #[test]
    fn higher_temperature_shortens_rupture_for_same_lmp() {
        // Same LMP (same damage / same applied stress on the master
        // curve), two temperatures: the hotter part ruptures sooner.
        let lmp = 25_000.0;
        let c = DEFAULT_C;
        let t_cool = rupture_time_hours(lmp, 1000.0, c).unwrap();
        let t_hot = rupture_time_hours(lmp, 1100.0, c).unwrap();
        assert!(
            t_hot < t_cool,
            "hotter should rupture sooner: t_hot {t_hot}, t_cool {t_cool}"
        );
    }

    #[test]
    fn rupture_time_round_trips_through_lmp() {
        // Compute LMP from a known point, then invert it back to the
        // time at the same temperature: must recover the original time.
        let t = 1023.0;
        let time = 5_000.0;
        let c = 19.0;
        let lmp = larson_miller_parameter(t, time, c).unwrap();
        let back = rupture_time_hours(lmp, t, c).unwrap();
        let rel = (back - time).abs() / time;
        assert!(rel < 1e-9, "round-trip failed: back {back}, time {time}");
    }

    #[test]
    fn rupture_point_parameter_round_trips() {
        let p = RupturePoint::new(1000.0, 100.0).unwrap();
        let c = DEFAULT_C;
        let lmp = p.parameter(c).unwrap();
        let back = rupture_time_hours(lmp, p.temperature_k, c).unwrap();
        assert!((back - p.time_hours).abs() < 1e-6, "got {back}");
    }

    #[test]
    fn temperature_inversion_round_trips() {
        // Invert for T, then recompute LMP at that T and time: recover.
        let lmp = 27_000.0;
        let time = 2_000.0;
        let c = DEFAULT_C;
        let t = rupture_temperature_k(lmp, time, c).unwrap();
        let back = larson_miller_parameter(t, time, c).unwrap();
        assert!((back - lmp).abs() < 1e-6, "got {back}, want {lmp}");
    }

    #[test]
    fn larger_lmp_gives_longer_life_at_fixed_temperature() {
        // On a master curve, lower stress → larger LMP → longer life.
        let c = DEFAULT_C;
        let t = 1000.0;
        let life_low = rupture_time_hours(24_000.0, t, c).unwrap();
        let life_high = rupture_time_hours(26_000.0, t, c).unwrap();
        assert!(
            life_high > life_low,
            "larger LMP should mean longer life: {life_high} vs {life_low}"
        );
    }

    #[test]
    fn ten_kelvin_step_changes_log_life_predictably() {
        // From t = 10^(LMP/T - C), holding LMP and C fixed, the log10
        // life is exactly LMP/T - C; check two temperatures agree with
        // the analytic value.
        let lmp = 25_000.0;
        let c = DEFAULT_C;
        for &t in &[900.0, 1000.0, 1200.0_f64] {
            let life = rupture_time_hours(lmp, t, c).unwrap();
            let expected_log = lmp / t - c;
            assert!(
                (life.log10() - expected_log).abs() < 1e-9,
                "log life mismatch at {t} K"
            );
        }
    }

    #[test]
    fn rejects_non_positive_temperature() {
        assert!(larson_miller_parameter(0.0, 1.0, 20.0).is_err());
        assert!(larson_miller_parameter(-5.0, 1.0, 20.0).is_err());
        assert!(rupture_time_hours(25_000.0, 0.0, 20.0).is_err());
        assert!(RupturePoint::new(-1.0, 10.0).is_err());
    }

    #[test]
    fn rejects_non_positive_time() {
        assert!(larson_miller_parameter(1000.0, 0.0, 20.0).is_err());
        assert!(RupturePoint::new(1000.0, 0.0).is_err());
        assert!(rupture_temperature_k(25_000.0, -3.0, 20.0).is_err());
    }

    #[test]
    fn rejects_non_finite_inputs() {
        assert!(larson_miller_parameter(1000.0, 10.0, f64::NAN).is_err());
        assert!(rupture_time_hours(f64::INFINITY, 1000.0, 20.0).is_err());
        assert!(rupture_temperature_k(25_000.0, 10.0, f64::NAN).is_err());
    }

    #[test]
    fn temperature_inversion_rejects_zero_denominator() {
        // C + log10(t) = 0 when t = 10^(-C); choose C = 3, t = 1e-3.
        let err = rupture_temperature_k(1000.0, 1e-3, 3.0);
        assert!(err.is_err());
    }
}
