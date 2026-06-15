//! Isentropic choked-flow detection through a nozzle / orifice.
//!
//! ## Model
//!
//! When compressed air discharges through a restriction (a valve port, an
//! orifice, a nozzle), the *velocity* at the throat cannot exceed the
//! local speed of sound. Once the downstream-to-upstream **absolute**
//! pressure ratio falls to or below a critical value, the flow at the
//! throat becomes sonic ("choked") and the mass-flow rate stops increasing
//! as the downstream pressure is lowered further.
//!
//! For an isentropic ideal gas the critical pressure ratio depends only on
//! the ratio of specific heats `k = c_p / c_v`:
//!
//! ```text
//! (p_down / p_up)_crit = (2 / (k + 1)) ^ ( k / (k - 1) )
//! ```
//!
//! For diatomic air, `k = 1.4`, this evaluates to about `0.528`. So:
//!
//! ```text
//! flow is choked  <=>  p_down / p_up <= 0.528   (air)
//! ```
//!
//! The condition is on **absolute** pressures (the gas laws know nothing of
//! gauges); convert with [`crate::compression::absolute_pressure`] first if
//! your readings are gauge.
//!
//! ## Honest scope
//!
//! This module decides *whether* flow is choked and at *what* ratio the
//! transition happens — it is the textbook isentropic threshold. It does
//! not compute the actual mass-flow rate (which needs the throat area, the
//! upstream temperature and a discharge coefficient), nor does it model
//! friction, heat transfer, real-gas behaviour, or condensation of the
//! water vapour in real shop air. It assumes `k` is constant.

use crate::error::{PneumaticsError, Result};

/// Ratio of specific heats for dry air (diatomic ideal gas), `k = 1.4`.
pub const GAMMA_AIR: f64 = 1.4;

/// The critical (choking) downstream-to-upstream absolute pressure ratio
/// for air, `(2/(k+1))^(k/(k-1))` with `k = 1.4`, which is approximately
/// `0.5283`.
pub const CRITICAL_PRESSURE_RATIO_AIR: f64 = 0.528_281_787_717_174_2;

/// The critical pressure ratio `(2/(k+1))^(k/(k-1))` for a gas with ratio
/// of specific heats `k` (dimensionless). Below this downstream/upstream
/// absolute-pressure ratio the flow is choked.
///
/// # Errors
///
/// [`PneumaticsError::NonPositive`] if `k` is not finite and `> 1`
/// (physically `k > 1` for any gas; `k <= 1` would make the exponent
/// singular or negative).
///
/// # Examples
///
/// ```
/// use valenx_pneumatics::flow::{critical_pressure_ratio, GAMMA_AIR, CRITICAL_PRESSURE_RATIO_AIR};
/// let r = critical_pressure_ratio(GAMMA_AIR).unwrap();
/// assert!((r - CRITICAL_PRESSURE_RATIO_AIR).abs() < 1e-12);
/// // Monatomic gas (k = 5/3) chokes at a lower ratio (~0.487).
/// let mono = critical_pressure_ratio(5.0 / 3.0).unwrap();
/// assert!((mono - 0.487_139_).abs() < 1e-3);
/// ```
pub fn critical_pressure_ratio(k: f64) -> Result<f64> {
    // `k` must be finite and strictly greater than 1: the exponent
    // k/(k-1) is singular at k = 1 and the model is meaningless for k <= 1.
    let k = PneumaticsError::positive("k", k)?;
    if k <= 1.0 {
        return Err(PneumaticsError::NonPositive {
            what: "k (must exceed 1)",
            value: k,
        });
    }
    Ok((2.0 / (k + 1.0)).powf(k / (k - 1.0)))
}

/// Whether flow from `upstream_abs` to `downstream_abs` (both **absolute**
/// pressures, pascals) is choked for a gas with ratio of specific heats
/// `k`. Choked when `p_down / p_up <= critical_pressure_ratio(k)`.
///
/// # Errors
///
/// - [`PneumaticsError::NonPositive`] if `upstream_abs` is not finite and
///   `> 0`, or if `k` is invalid (see [`critical_pressure_ratio`]).
/// - [`PneumaticsError::Negative`] if `downstream_abs` is negative or
///   non-finite (a perfect vacuum, zero, is allowed and is always choked).
/// - [`PneumaticsError::Geometry`] if `downstream_abs > upstream_abs`
///   (back-flow: the pressure ratio model does not apply).
///
/// # Examples
///
/// ```
/// use valenx_pneumatics::flow::{is_choked, GAMMA_AIR};
/// // 7 bar abs venting to 1 bar abs: ratio 0.143 < 0.528 -> choked.
/// assert!(is_choked(700_000.0, 100_000.0, GAMMA_AIR).unwrap());
/// // 7 bar abs to 6 bar abs: ratio 0.857 > 0.528 -> not choked.
/// assert!(!is_choked(700_000.0, 600_000.0, GAMMA_AIR).unwrap());
/// ```
pub fn is_choked(upstream_abs: f64, downstream_abs: f64, k: f64) -> Result<bool> {
    let ratio = pressure_ratio(upstream_abs, downstream_abs)?;
    let crit = critical_pressure_ratio(k)?;
    Ok(ratio <= crit)
}

/// Convenience wrapper for [`is_choked`] specialised to air (`k = 1.4`).
///
/// # Errors
///
/// Same as [`is_choked`].
///
/// # Examples
///
/// ```
/// use valenx_pneumatics::flow::is_choked_air;
/// assert!(is_choked_air(700_000.0, 100_000.0).unwrap());
/// assert!(!is_choked_air(700_000.0, 600_000.0).unwrap());
/// ```
pub fn is_choked_air(upstream_abs: f64, downstream_abs: f64) -> Result<bool> {
    is_choked(upstream_abs, downstream_abs, GAMMA_AIR)
}

/// The downstream-to-upstream absolute-pressure ratio `p_down / p_up`
/// (dimensionless, in `[0, 1]`).
///
/// # Errors
///
/// - [`PneumaticsError::NonPositive`] if `upstream_abs` is not finite and
///   `> 0`.
/// - [`PneumaticsError::Negative`] if `downstream_abs` is negative or
///   non-finite.
/// - [`PneumaticsError::Geometry`] if `downstream_abs > upstream_abs`.
///
/// # Examples
///
/// ```
/// use valenx_pneumatics::flow::pressure_ratio;
/// assert!((pressure_ratio(800_000.0, 200_000.0).unwrap() - 0.25).abs() < 1e-12);
/// ```
pub fn pressure_ratio(upstream_abs: f64, downstream_abs: f64) -> Result<f64> {
    let up = PneumaticsError::positive("upstream_abs", upstream_abs)?;
    let down = PneumaticsError::non_negative("downstream_abs", downstream_abs)?;
    if down > up {
        return Err(PneumaticsError::Geometry(
            "downstream pressure exceeds upstream (back-flow); ratio undefined",
        ));
    }
    Ok(down / up)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute-difference tolerance for floating comparisons.
    const EPS: f64 = 1e-12;

    #[test]
    fn critical_ratio_for_air_is_about_half() {
        // Closed-form check against (2/2.4)^(1.4/0.4).
        let r = critical_pressure_ratio(GAMMA_AIR).unwrap();
        let expected = (2.0_f64 / 2.4).powf(1.4 / 0.4);
        assert!((r - expected).abs() < EPS);
        assert!((r - CRITICAL_PRESSURE_RATIO_AIR).abs() < EPS);
        // Sanity: physically it sits near 0.528.
        assert!((r - 0.528).abs() < 1e-3);
    }

    #[test]
    fn published_constant_matches_function() {
        // The exported constant must equal the computed value for air.
        assert!(
            (CRITICAL_PRESSURE_RATIO_AIR - critical_pressure_ratio(GAMMA_AIR).unwrap()).abs() < EPS
        );
    }

    #[test]
    fn monatomic_gas_chokes_lower_than_air() {
        // k = 5/3 -> (2/(8/3))^((5/3)/(2/3)) = 0.75^2.5 ~= 0.487.
        let mono = critical_pressure_ratio(5.0 / 3.0).unwrap();
        let expected = (0.75_f64).powf(2.5);
        assert!((mono - expected).abs() < EPS);
        assert!(mono < CRITICAL_PRESSURE_RATIO_AIR);
    }

    #[test]
    fn choked_below_critical_ratio() {
        // 7 bar abs to 1 bar abs -> ratio 0.1429 < 0.528 -> choked.
        assert!(is_choked_air(700_000.0, 100_000.0).unwrap());
    }

    #[test]
    fn not_choked_above_critical_ratio() {
        // 7 bar abs to 6 bar abs -> ratio 0.857 > 0.528 -> subsonic.
        assert!(!is_choked_air(700_000.0, 600_000.0).unwrap());
    }

    #[test]
    fn boundary_is_choked_inclusive() {
        // At exactly the critical ratio the flow is choked (<=).
        let up = 700_000.0;
        let down = up * CRITICAL_PRESSURE_RATIO_AIR;
        assert!(is_choked_air(up, down).unwrap());

        // A hair above the critical ratio is NOT choked.
        let down_hi = up * (CRITICAL_PRESSURE_RATIO_AIR + 1e-6);
        assert!(!is_choked_air(up, down_hi).unwrap());

        // A hair below the critical ratio IS choked.
        let down_lo = up * (CRITICAL_PRESSURE_RATIO_AIR - 1e-6);
        assert!(is_choked_air(up, down_lo).unwrap());
    }

    #[test]
    fn venting_to_vacuum_is_choked() {
        // Zero downstream pressure -> ratio 0 -> always choked.
        assert!(is_choked_air(700_000.0, 0.0).unwrap());
    }

    #[test]
    fn pressure_ratio_is_quotient() {
        let r = pressure_ratio(800_000.0, 200_000.0).unwrap();
        assert!((r - 0.25).abs() < EPS);
    }

    #[test]
    fn equal_pressures_give_ratio_one_and_not_choked() {
        let r = pressure_ratio(500_000.0, 500_000.0).unwrap();
        assert!((r - 1.0).abs() < EPS);
        assert!(!is_choked_air(500_000.0, 500_000.0).unwrap());
    }

    #[test]
    fn rejects_backflow() {
        // Downstream above upstream is rejected, not silently flipped.
        assert!(pressure_ratio(100_000.0, 200_000.0).is_err());
        assert!(is_choked_air(100_000.0, 200_000.0).is_err());
    }

    #[test]
    fn rejects_nonpositive_upstream() {
        assert!(pressure_ratio(0.0, 0.0).is_err());
        assert!(is_choked_air(-1.0, 0.0).is_err());
    }

    #[test]
    fn rejects_negative_downstream() {
        let err = pressure_ratio(700_000.0, -1.0).unwrap_err();
        assert_eq!(err.code(), "pneumatics.negative");
    }

    #[test]
    fn rejects_invalid_k() {
        // k must exceed 1.
        assert!(critical_pressure_ratio(1.0).is_err());
        assert!(critical_pressure_ratio(0.5).is_err());
        assert!(critical_pressure_ratio(f64::NAN).is_err());
    }
}
