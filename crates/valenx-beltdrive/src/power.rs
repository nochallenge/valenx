//! Transmitted power, centrifugal tension, and drive capacity.
//!
//! ## Transmitted power
//!
//! The net driving force at the belt is the difference between the
//! tight- and slack-side tensions, and power is force times belt speed:
//!
//! ```text
//! P = (T1 - T2) * v
//! ```
//!
//! With `T1`, `T2` in newtons and `v` in m/s, `P` is in watts.
//!
//! ## Centrifugal tension
//!
//! A belt of linear mass density `m` (kg/m) running at speed `v`
//! develops a tension that acts equally on both sides as it is flung
//! outward around the pulleys:
//!
//! ```text
//! Tc = m * v^2
//! ```
//!
//! Centrifugal tension does no useful work but consumes part of the
//! allowable belt tension, so it reduces the power a drive can transmit
//! at high speed. The effective tensions available for driving are
//! `T1 - Tc` (tight) and `T2 - Tc` (slack), and these are what enter the
//! capstan limit `(T1 - Tc)/(T2 - Tc) = exp(mu*theta)`.

use crate::error::BeltError;

/// Transmitted power `P = (T1 - T2) * v`.
///
/// `t1` is the tight-side and `t2` the slack-side tension (consistent
/// force unit, e.g. newtons); `belt_speed` is the belt linear speed in
/// the matching length-per-time unit (e.g. m/s). With newtons and m/s
/// the result is watts.
///
/// # Errors
///
/// Returns [`BeltError::BadParameter`] if any input is non-finite, if
/// `belt_speed` is negative, or if `t1 < t2` (the tight side must carry
/// at least the slack-side tension for a driving belt).
pub fn transmitted_power(t1: f64, t2: f64, belt_speed: f64) -> Result<f64, BeltError> {
    require_finite("t1", t1)?;
    require_finite("t2", t2)?;
    require_non_negative("belt_speed", belt_speed)?;
    if t1 < t2 {
        return Err(BeltError::bad_parameter(
            "t1",
            format!("tight-side tension {t1} must be >= slack-side tension {t2}"),
        ));
    }
    Ok((t1 - t2) * belt_speed)
}

/// Centrifugal tension `Tc = m * v^2`.
///
/// `linear_density` is the belt mass per unit length (kg/m) and
/// `belt_speed` the belt linear speed (m/s); the result is in newtons.
///
/// # Errors
///
/// Returns [`BeltError::BadParameter`] if `linear_density` is negative
/// or `belt_speed` is negative (or either is non-finite).
pub fn centrifugal_tension(linear_density: f64, belt_speed: f64) -> Result<f64, BeltError> {
    require_non_negative("linear_density", linear_density)?;
    require_non_negative("belt_speed", belt_speed)?;
    Ok(linear_density * belt_speed * belt_speed)
}

/// Maximum power a drive can transmit at the point of slipping for a
/// given maximum allowable tight-side tension.
///
/// At impending slip the effective tensions obey the capstan relation
/// `(T1 - Tc)/(T2 - Tc) = exp(mu*theta)`, so the driving force is
///
/// ```text
/// T1 - T2 = (T1 - Tc) * (1 - 1/k),   k = exp(mu*theta)
/// ```
///
/// and `P_max = (T1 - T2) * v`. The centrifugal term `Tc = m*v^2` is
/// computed from `linear_density` and `belt_speed`.
///
/// # Errors
///
/// Returns [`BeltError::BadParameter`] for out-of-domain inputs (see
/// [`centrifugal_tension`] and [`crate::friction::tension_ratio`]), and
/// when the centrifugal tension equals or exceeds `t1_max`, in which
/// case no effective tension is left to drive the load.
pub fn max_power(
    t1_max: f64,
    linear_density: f64,
    belt_speed: f64,
    mu: f64,
    wrap_angle: f64,
) -> Result<f64, BeltError> {
    require_non_negative("t1_max", t1_max)?;
    let tc = centrifugal_tension(linear_density, belt_speed)?;
    if tc >= t1_max {
        return Err(BeltError::bad_parameter(
            "belt_speed",
            format!("centrifugal tension {tc} leaves no driving tension below t1_max {t1_max}"),
        ));
    }
    let k = crate::friction::tension_ratio(mu, wrap_angle)?;
    // Effective driving force at impending slip.
    let driving_force = (t1_max - tc) * (1.0 - 1.0 / k);
    Ok(driving_force * belt_speed)
}

/// Belt linear speed that maximises transmitted power for a belt with a
/// fixed maximum allowable tension, `v* = sqrt(T1_max / (3*m))`.
///
/// Differentiating `P_max(v)` (with `Tc = m*v^2`) and setting the
/// derivative to zero gives the classic result that maximum power is
/// reached when the centrifugal tension equals one third of the maximum
/// allowable tension, i.e. `Tc = T1_max / 3`.
///
/// # Errors
///
/// Returns [`BeltError::BadParameter`] if `t1_max` is not strictly
/// positive or `linear_density` is not strictly positive.
pub fn speed_for_max_power(t1_max: f64, linear_density: f64) -> Result<f64, BeltError> {
    if !t1_max.is_finite() || t1_max <= 0.0 {
        return Err(BeltError::bad_parameter(
            "t1_max",
            format!("must be finite and > 0, got {t1_max}"),
        ));
    }
    if !linear_density.is_finite() || linear_density <= 0.0 {
        return Err(BeltError::bad_parameter(
            "linear_density",
            format!("must be finite and > 0, got {linear_density}"),
        ));
    }
    Ok((t1_max / (3.0 * linear_density)).sqrt())
}

/// Validate that `value` is finite.
fn require_finite(name: &'static str, value: f64) -> Result<(), BeltError> {
    if !value.is_finite() {
        return Err(BeltError::bad_parameter(
            name,
            format!("must be finite, got {value}"),
        ));
    }
    Ok(())
}

/// Validate that `value` is finite and not negative.
fn require_non_negative(name: &'static str, value: f64) -> Result<(), BeltError> {
    if !value.is_finite() || value < 0.0 {
        return Err(BeltError::bad_parameter(
            name,
            format!("must be finite and >= 0, got {value}"),
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::friction::tension_ratio;
    use std::f64::consts::PI;

    const EPS: f64 = 1e-9;

    #[test]
    fn transmitted_power_matches_force_times_speed() {
        // (1000 - 400) N at 10 m/s -> 6000 W.
        let p = transmitted_power(1000.0, 400.0, 10.0).unwrap();
        assert!((p - 6000.0).abs() < EPS, "power was {p}");
    }

    #[test]
    fn transmitted_power_zero_when_tensions_equal() {
        let p = transmitted_power(500.0, 500.0, 12.0).unwrap();
        assert!(p.abs() < EPS, "power was {p}");
    }

    #[test]
    fn transmitted_power_scales_with_speed() {
        let p1 = transmitted_power(900.0, 300.0, 5.0).unwrap();
        let p2 = transmitted_power(900.0, 300.0, 10.0).unwrap();
        assert!((p2 - 2.0 * p1).abs() < EPS);
    }

    #[test]
    fn transmitted_power_rejects_tight_below_slack() {
        let err = transmitted_power(300.0, 900.0, 5.0).unwrap_err();
        assert_eq!(err.code(), "beltdrive.bad_parameter");
    }

    #[test]
    fn transmitted_power_rejects_negative_speed() {
        assert!(transmitted_power(900.0, 300.0, -1.0).is_err());
    }

    #[test]
    fn centrifugal_tension_matches_m_v_squared() {
        // m = 0.5 kg/m, v = 20 m/s -> Tc = 0.5 * 400 = 200 N.
        let tc = centrifugal_tension(0.5, 20.0).unwrap();
        assert!((tc - 200.0).abs() < EPS, "Tc was {tc}");
    }

    #[test]
    fn centrifugal_tension_grows_quadratically() {
        let tc1 = centrifugal_tension(0.5, 10.0).unwrap();
        let tc2 = centrifugal_tension(0.5, 20.0).unwrap();
        // Doubling speed quadruples Tc.
        assert!((tc2 - 4.0 * tc1).abs() < EPS);
    }

    #[test]
    fn centrifugal_tension_zero_at_rest() {
        let tc = centrifugal_tension(0.5, 0.0).unwrap();
        assert!(tc.abs() < EPS);
    }

    #[test]
    fn centrifugal_tension_rejects_negative_density() {
        assert!(centrifugal_tension(-0.1, 10.0).is_err());
    }

    #[test]
    fn max_power_recovers_force_times_speed_no_centrifugal() {
        // With m = 0, Tc = 0; effective force = T1_max*(1 - 1/k).
        let t1 = 1000.0;
        let mu = 0.3;
        let theta = PI;
        let v = 10.0;
        let p = max_power(t1, 0.0, v, mu, theta).unwrap();
        let k = tension_ratio(mu, theta).unwrap();
        let expected = t1 * (1.0 - 1.0 / k) * v;
        assert!((p - expected).abs() < 1e-6, "power was {p}");
    }

    #[test]
    fn max_power_consistent_with_capstan_tensions() {
        // Build T1, T2 from the capstan limit and centrifugal tension,
        // then check max_power == (T1 - T2)*v.
        let t1 = 1200.0;
        let m = 0.4;
        let v = 15.0;
        let mu = 0.25;
        let theta = 2.8;
        let tc = centrifugal_tension(m, v).unwrap();
        let k = tension_ratio(mu, theta).unwrap();
        // (T1 - Tc)/(T2 - Tc) = k  ->  T2 = Tc + (T1 - Tc)/k.
        let t2 = tc + (t1 - tc) / k;
        let direct = transmitted_power(t1, t2, v).unwrap();
        let viamax = max_power(t1, m, v, mu, theta).unwrap();
        assert!(
            (direct - viamax).abs() < 1e-6,
            "direct={direct} max={viamax}"
        );
    }

    #[test]
    fn higher_wrap_raises_max_power() {
        // More wrap -> more capacity at the same allowable tension.
        let low = max_power(1000.0, 0.3, 12.0, 0.3, PI - 0.4).unwrap();
        let high = max_power(1000.0, 0.3, 12.0, 0.3, PI + 0.4).unwrap();
        assert!(high > low, "high={high} low={low}");
    }

    #[test]
    fn higher_friction_raises_max_power() {
        let low = max_power(1000.0, 0.3, 12.0, 0.2, PI).unwrap();
        let high = max_power(1000.0, 0.3, 12.0, 0.4, PI).unwrap();
        assert!(high > low, "high={high} low={low}");
    }

    #[test]
    fn max_power_rejects_centrifugal_exceeding_t1max() {
        // m*v^2 >= T1_max -> no driving tension left.
        let err = max_power(100.0, 1.0, 20.0, 0.3, PI).unwrap_err();
        assert_eq!(err.code(), "beltdrive.bad_parameter");
    }

    #[test]
    fn optimum_speed_puts_tc_at_one_third_of_t1max() {
        // v* = sqrt(T1_max/(3m)) -> Tc(v*) = T1_max/3.
        let t1 = 1500.0;
        let m = 0.5;
        let v_star = speed_for_max_power(t1, m).unwrap();
        let tc = centrifugal_tension(m, v_star).unwrap();
        assert!((tc - t1 / 3.0).abs() < 1e-6, "Tc at optimum was {tc}");
    }

    #[test]
    fn optimum_speed_is_a_local_maximum_of_power() {
        // Power at v* should beat power slightly either side.
        let t1 = 1500.0;
        let m = 0.5;
        let mu = 0.3;
        let theta = PI;
        let v_star = speed_for_max_power(t1, m).unwrap();
        let p_star = max_power(t1, m, v_star, mu, theta).unwrap();
        let p_lo = max_power(t1, m, v_star * 0.9, mu, theta).unwrap();
        let p_hi = max_power(t1, m, v_star * 1.1, mu, theta).unwrap();
        assert!(p_star >= p_lo, "p*={p_star} p_lo={p_lo}");
        assert!(p_star >= p_hi, "p*={p_star} p_hi={p_hi}");
    }

    #[test]
    fn optimum_speed_rejects_nonpositive_inputs() {
        assert!(speed_for_max_power(0.0, 0.5).is_err());
        assert!(speed_for_max_power(1000.0, 0.0).is_err());
    }
}
