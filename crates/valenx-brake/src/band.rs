//! Band / drum (capstan) brake tension and torque.
//!
//! A flexible band wraps a rotating drum over a contact (wrap) angle
//! `theta`. Integrating the differential Coulomb friction along the
//! contact gives the classic capstan (Euler–Eytelwein / belt-friction)
//! equation relating the tight-side tension `T1` to the slack-side
//! tension `T2`:
//!
//! ```text
//! T1 / T2 = exp(mu * theta)
//! ```
//!
//! with `theta` in **radians**. The net braking force the band exerts
//! at the drum surface is `T1 - T2`, and the braking torque about the
//! drum axis is `(T1 - T2) * r`.
//!
//! The capstan relation inverts cleanly for the wrap angle:
//! [`wrap_angle_for_ratio`] solves `theta = ln(ratio) / mu`, the
//! design question of how much wrap (how many turns `theta / 2*pi`) a
//! given friction needs to reach a target tension ratio.

use crate::error::{check_finite, check_friction, check_non_negative, check_positive, BrakeError};

/// Tension ratio `T1 / T2 = exp(mu * theta)` of a band/drum brake.
///
/// This is the dimensionless mechanical-advantage factor the band's
/// friction multiplies across the wrap.
///
/// # Parameters
/// - `mu` — friction coefficient, in `(0, MU_MAX]`.
/// - `wrap_angle_rad` — contact (wrap) angle `theta`, in **radians** (>= 0).
///
/// # Errors
/// [`BrakeError`] if `mu` is out of range or `wrap_angle_rad` is
/// non-finite or negative.
///
/// # Examples
/// ```
/// use valenx_brake::band::tension_ratio;
/// // A full 360° wrap (2π rad) with mu = 0.25.
/// let r = tension_ratio(0.25, std::f64::consts::TAU).unwrap();
/// assert!((r - (0.25 * std::f64::consts::TAU).exp()).abs() < 1e-12);
/// ```
pub fn tension_ratio(mu: f64, wrap_angle_rad: f64) -> Result<f64, BrakeError> {
    let mu = check_friction(mu)?;
    let theta = check_non_negative("wrap_angle_rad", wrap_angle_rad)?;
    Ok((mu * theta).exp())
}

/// The wrap angle `theta` (radians) needed to reach a target tension
/// ratio, inverting the capstan equation `T1/T2 = exp(mu*theta)`:
///
/// ```text
/// theta = ln(ratio) / mu
/// ```
///
/// This is the classic capstan-design question — how much wrap (e.g. how
/// many turns `theta / 2*pi` of rope around a bollard) a given friction
/// coefficient needs to multiply tension by `ratio`. Because the ratio
/// grows *exponentially* with wrap, the required wrap grows only
/// *logarithmically* with the ratio. A `target_ratio` of `1` needs no
/// wrap (`theta = 0`), and every ratio `>= 1` is reachable.
///
/// # Parameters
/// - `mu` — friction coefficient, in `(0, MU_MAX]`.
/// - `target_ratio` — desired `T1/T2`, dimensionless and `>= 1`.
///
/// # Errors
/// [`BrakeError`] if `mu` is out of range, `target_ratio` is non-finite,
/// or `target_ratio < 1` (a ratio below 1 would require a negative wrap
/// angle, which is unphysical).
///
/// # Examples
/// ```
/// use valenx_brake::band::wrap_angle_for_ratio;
/// // With mu = ln 2, doubling the tension (ratio = 2) needs theta = 1 rad.
/// let theta = wrap_angle_for_ratio(std::f64::consts::LN_2, 2.0).unwrap();
/// assert!((theta - 1.0).abs() < 1e-12);
/// ```
pub fn wrap_angle_for_ratio(mu: f64, target_ratio: f64) -> Result<f64, BrakeError> {
    let mu = check_friction(mu)?;
    let ratio = check_finite("target_ratio", target_ratio)?;
    if ratio < 1.0 {
        // ratio = exp(mu*theta) >= 1 for any non-negative wrap; a smaller
        // target is unreachable. Report it as `ratio - 1` being negative.
        return Err(BrakeError::Negative {
            name: "target_ratio - 1",
            value: ratio - 1.0,
        });
    }
    Ok(ratio.ln() / mu)
}

/// Tight-side tension `T1` given the slack-side tension `T2`.
///
/// `T1 = T2 * exp(mu * theta)`.
///
/// # Errors
/// [`BrakeError`] if `mu` is out of range, `slack_tension_n` is
/// non-finite or negative, or `wrap_angle_rad` is non-finite or
/// negative.
///
/// # Examples
/// ```
/// use valenx_brake::band::tight_side_tension;
/// let t1 = tight_side_tension(100.0, 0.3, 90.0_f64.to_radians()).unwrap();
/// assert!(t1 > 100.0);
/// ```
pub fn tight_side_tension(
    slack_tension_n: f64,
    mu: f64,
    wrap_angle_rad: f64,
) -> Result<f64, BrakeError> {
    let t2 = check_non_negative("slack_tension_n", slack_tension_n)?;
    let ratio = tension_ratio(mu, wrap_angle_rad)?;
    Ok(t2 * ratio)
}

/// Slack-side tension `T2` given the tight-side tension `T1`.
///
/// `T2 = T1 / exp(mu * theta) = T1 * exp(-mu * theta)`.
///
/// # Errors
/// [`BrakeError`] if `mu` is out of range, `tight_tension_n` is
/// non-finite or negative, or `wrap_angle_rad` is non-finite or
/// negative.
///
/// # Examples
/// ```
/// use valenx_brake::band::slack_side_tension;
/// let t2 = slack_side_tension(271.83, 1.0, 1.0).unwrap();
/// assert!((t2 - 100.0).abs() < 0.01);
/// ```
pub fn slack_side_tension(
    tight_tension_n: f64,
    mu: f64,
    wrap_angle_rad: f64,
) -> Result<f64, BrakeError> {
    let t1 = check_non_negative("tight_tension_n", tight_tension_n)?;
    let ratio = tension_ratio(mu, wrap_angle_rad)?;
    Ok(t1 / ratio)
}

/// Net braking force at the drum surface, `T1 - T2`, in newtons.
///
/// Given the slack-side tension and the friction/wrap, this is the
/// effective tangential drag the band applies to the drum.
///
/// # Errors
/// [`BrakeError`] if `mu` is out of range, `slack_tension_n` is
/// non-finite or negative, or `wrap_angle_rad` is non-finite or
/// negative.
///
/// # Examples
/// ```
/// use valenx_brake::band::braking_force;
/// // With T2 = 100 N, mu*theta = ln(2): T1 = 200, so force = 100.
/// let f = braking_force(100.0, std::f64::consts::LN_2, 1.0).unwrap();
/// assert!((f - 100.0).abs() < 1e-9);
/// ```
pub fn braking_force(
    slack_tension_n: f64,
    mu: f64,
    wrap_angle_rad: f64,
) -> Result<f64, BrakeError> {
    let t2 = check_non_negative("slack_tension_n", slack_tension_n)?;
    let t1 = tight_side_tension(t2, mu, wrap_angle_rad)?;
    Ok(t1 - t2)
}

/// Braking torque about the drum axis, `(T1 - T2) * r`, in newton-metres.
///
/// # Errors
/// [`BrakeError`] if `mu` is out of range, `slack_tension_n` is
/// non-finite or negative, `wrap_angle_rad` is non-finite or negative,
/// or `drum_radius_m` is non-finite or non-positive.
///
/// # Examples
/// ```
/// use valenx_brake::band::braking_torque;
/// // force 100 N at r = 0.2 m -> 20 N·m.
/// let t = braking_torque(100.0, std::f64::consts::LN_2, 1.0, 0.2).unwrap();
/// assert!((t - 20.0).abs() < 1e-9);
/// ```
pub fn braking_torque(
    slack_tension_n: f64,
    mu: f64,
    wrap_angle_rad: f64,
    drum_radius_m: f64,
) -> Result<f64, BrakeError> {
    let force = braking_force(slack_tension_n, mu, wrap_angle_rad)?;
    let r = check_positive("drum_radius_m", drum_radius_m)?;
    Ok(force * r)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{LN_2, PI, TAU};

    const EPS: f64 = 1e-12;

    #[test]
    fn ratio_matches_exponential() {
        // T1/T2 = exp(mu*theta) for a 180° wrap.
        let r = tension_ratio(0.3, PI).unwrap();
        let expected = (0.3 * PI).exp();
        assert!((r - expected).abs() < EPS, "got {r} expected {expected}");
    }

    #[test]
    fn zero_wrap_gives_unit_ratio() {
        // exp(mu*0) = 1: no wrap, no friction multiplication.
        let r = tension_ratio(0.5, 0.0).unwrap();
        assert!((r - 1.0).abs() < EPS, "got {r}");
    }

    #[test]
    fn ln2_special_case() {
        // Choosing mu*theta = ln 2 makes the ratio exactly 2.
        let r = tension_ratio(LN_2, 1.0).unwrap();
        assert!((r - 2.0).abs() < EPS, "got {r}");
        // And with mu = 1, theta = ln 2 gives the same.
        let r2 = tension_ratio(1.0, LN_2).unwrap();
        assert!((r2 - 2.0).abs() < EPS, "got {r2}");
    }

    #[test]
    fn full_turn_ratio() {
        // A full 2π wrap with mu = 0.2.
        let r = tension_ratio(0.2, TAU).unwrap();
        let expected = (0.2 * TAU).exp();
        assert!((r - expected).abs() < EPS, "got {r}");
    }

    #[test]
    fn ratio_increases_with_wrap_angle() {
        // More wrap -> strictly larger ratio (monotonic in theta).
        let a = tension_ratio(0.3, PI * 0.5).unwrap();
        let b = tension_ratio(0.3, PI).unwrap();
        let c = tension_ratio(0.3, PI * 1.5).unwrap();
        assert!(a < b && b < c, "a {a} b {b} c {c}");
    }

    #[test]
    fn tight_and_slack_round_trip() {
        // T1 from T2, then T2 back from T1, recovers the original.
        let t2 = 100.0;
        let t1 = tight_side_tension(t2, 0.3, PI).unwrap();
        let back = slack_side_tension(t1, 0.3, PI).unwrap();
        assert!((back - t2).abs() < 1e-9, "t1 {t1} back {back}");
        // And T1 should be T2 * exp(mu*theta).
        assert!((t1 - t2 * (0.3 * PI).exp()).abs() < 1e-9, "t1 {t1}");
    }

    #[test]
    fn braking_force_is_difference() {
        // mu*theta = ln 2 -> T1 = 2*T2 -> force = T2.
        let f = braking_force(100.0, LN_2, 1.0).unwrap();
        assert!((f - 100.0).abs() < 1e-9, "got {f}");
    }

    #[test]
    fn braking_torque_is_force_times_radius() {
        // force 100 N (ln2 case) at 0.2 m radius -> 20 N·m.
        let t = braking_torque(100.0, LN_2, 1.0, 0.2).unwrap();
        assert!((t - 20.0).abs() < 1e-9, "got {t}");
        // Cross-check against the explicit force*radius product.
        let f = braking_force(100.0, LN_2, 1.0).unwrap();
        assert!((t - f * 0.2).abs() < 1e-9, "t {t} f {f}");
    }

    #[test]
    fn rejects_bad_inputs() {
        assert_eq!(
            tension_ratio(0.0, PI).unwrap_err().code(),
            "brake.friction_out_of_range"
        );
        assert_eq!(
            tension_ratio(0.3, -0.1).unwrap_err().code(),
            "brake.negative"
        );
        assert_eq!(
            braking_torque(100.0, 0.3, PI, 0.0).unwrap_err().code(),
            "brake.non_positive"
        );
        assert_eq!(
            tension_ratio(0.3, f64::INFINITY).unwrap_err().code(),
            "brake.not_finite"
        );
    }

    #[test]
    fn wrap_angle_inverts_tension_ratio_both_ways() {
        let mu = 0.3;
        // ratio -> theta -> ratio.
        let r0 = 2.5;
        let theta = wrap_angle_for_ratio(mu, r0).unwrap();
        let r_back = tension_ratio(mu, theta).unwrap();
        assert!((r_back - r0).abs() < 1e-12, "ratio round-trip {r_back}");
        // theta -> ratio -> theta.
        let theta0 = PI;
        let r = tension_ratio(mu, theta0).unwrap();
        let theta_back = wrap_angle_for_ratio(mu, r).unwrap();
        assert!(
            (theta_back - theta0).abs() < 1e-12,
            "theta round-trip {theta_back}"
        );
    }

    #[test]
    fn wrap_angle_hand_value() {
        // mu*theta = ln(ratio). mu = ln 2, ratio = 2 -> theta = 1.
        let theta = wrap_angle_for_ratio(LN_2, 2.0).unwrap();
        assert!((theta - 1.0).abs() < EPS, "got {theta}");
        // mu = 0.3, ratio = 2 -> theta = ln(2) / 0.3.
        let theta2 = wrap_angle_for_ratio(0.3, 2.0).unwrap();
        assert!((theta2 - LN_2 / 0.3).abs() < EPS, "got {theta2}");
    }

    #[test]
    fn unit_ratio_needs_no_wrap() {
        // ratio = 1 -> theta = 0 (ln 1 = 0), for any mu.
        for mu in [0.1_f64, 0.5, 1.0] {
            let theta = wrap_angle_for_ratio(mu, 1.0).unwrap();
            assert!(theta.abs() < EPS, "mu {mu} theta {theta}");
        }
    }

    #[test]
    fn wrap_angle_monotonic() {
        // Larger target ratio -> larger wrap; larger mu -> smaller wrap.
        let a = wrap_angle_for_ratio(0.3, 2.0).unwrap();
        let b = wrap_angle_for_ratio(0.3, 4.0).unwrap();
        assert!(b > a, "more ratio should need more wrap: {a} {b}");
        let lo_mu = wrap_angle_for_ratio(0.2, 3.0).unwrap();
        let hi_mu = wrap_angle_for_ratio(0.6, 3.0).unwrap();
        assert!(
            hi_mu < lo_mu,
            "more friction should need less wrap: {lo_mu} {hi_mu}"
        );
    }

    #[test]
    fn capstan_wrap_grows_logarithmically() {
        // Holding ratio grows exponentially with turns, so the wrap grows
        // only logarithmically with ratio: ratio 100 needs exactly twice
        // the wrap of ratio 10 (since 100 = 10^2).
        let mu = 0.25;
        let w10 = wrap_angle_for_ratio(mu, 10.0).unwrap();
        let w100 = wrap_angle_for_ratio(mu, 100.0).unwrap();
        assert!((w100 - 2.0 * w10).abs() < 1e-9, "w10 {w10} w100 {w100}");
    }

    #[test]
    fn wrap_angle_rejects_bad_inputs() {
        // Friction out of range.
        assert_eq!(
            wrap_angle_for_ratio(0.0, 2.0).unwrap_err().code(),
            "brake.friction_out_of_range"
        );
        // Ratio below 1 would require a negative wrap angle.
        assert_eq!(
            wrap_angle_for_ratio(0.3, 0.5).unwrap_err().code(),
            "brake.negative"
        );
        // Non-finite ratio.
        assert_eq!(
            wrap_angle_for_ratio(0.3, f64::NAN).unwrap_err().code(),
            "brake.not_finite"
        );
        assert_eq!(
            wrap_angle_for_ratio(0.3, f64::INFINITY).unwrap_err().code(),
            "brake.not_finite"
        );
    }
}
