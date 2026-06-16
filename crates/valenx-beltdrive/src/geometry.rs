//! Belt-drive kinematics and pulley geometry.
//!
//! Closed-form relations for an idealised flat-belt drive with rigid
//! pulleys and an inextensible, non-slipping belt:
//!
//! ## Speed ratio
//!
//! With no slip the belt linear speed is common to both pulleys, so the
//! surface speeds match: `pi * D_driver * N_driver = pi * D_driven *
//! N_driven`. The transmission (velocity) ratio is therefore
//!
//! ```text
//! i = D_driven / D_driver = N_driver / N_driven
//! ```
//!
//! ## Belt linear speed
//!
//! The belt travels at the pulley rim speed:
//!
//! ```text
//! v = pi * D * N
//! ```
//!
//! where `D` is the pitch diameter and `N` the rotational speed. With
//! `D` in metres and `N` in revolutions per second, `v` is in m/s.
//!
//! ## Open-belt wrap geometry
//!
//! For an open belt (pulleys turning the same way) on parallel shafts a
//! distance `C` apart, the angle of wrap on the small pulley is less
//! than 180 degrees and on the large pulley more than 180 degrees:
//!
//! ```text
//! alpha     = asin((R_large - R_small) / C)
//! theta_sml = pi - 2 * alpha
//! theta_lrg = pi + 2 * alpha
//! ```
//!
//! and the total belt length is
//!
//! ```text
//! L = 2*C*cos(alpha) + R_sml*theta_sml + R_lrg*theta_lrg.
//! ```
//!
//! ## Crossed-belt wrap geometry
//!
//! A crossed belt figure-eights so the shafts counter-rotate; both pulleys
//! then wrap by the *same* angle, set by the **sum** of the radii:
//!
//! ```text
//! gamma = asin((R_small + R_large) / C)
//! theta = pi + 2 * gamma        (on both pulleys)
//! L     = 2*C*cos(gamma) + (R_small + R_large) * theta
//! ```
//!
//! ([`wrap_angle_crossed`], [`belt_length_crossed`]). A crossed belt is
//! always a little longer than the open belt of the same geometry.

use crate::error::BeltError;
use std::f64::consts::PI;

/// Transmission (velocity) ratio of a belt drive,
/// `i = D_driven / D_driver`.
///
/// A ratio greater than one is a speed *reduction* (the driven pulley
/// is larger and turns slower); less than one is a speed *increase*.
///
/// # Errors
///
/// Returns [`BeltError::BadParameter`] if either diameter is not
/// strictly positive.
pub fn speed_ratio(driver_diameter: f64, driven_diameter: f64) -> Result<f64, BeltError> {
    require_positive("driver_diameter", driver_diameter)?;
    require_positive("driven_diameter", driven_diameter)?;
    Ok(driven_diameter / driver_diameter)
}

/// Belt linear (rim) speed `v = pi * D * N`.
///
/// `diameter` is the pulley pitch diameter and `rev_per_sec` its
/// rotational speed in revolutions per second; the result carries the
/// same length unit as `diameter`, per unit time (e.g. m given metres
/// and rev/s). Use [`rpm_to_rev_per_sec`] to convert from RPM.
///
/// # Errors
///
/// Returns [`BeltError::BadParameter`] if `diameter` is not strictly
/// positive or `rev_per_sec` is negative.
pub fn belt_speed(diameter: f64, rev_per_sec: f64) -> Result<f64, BeltError> {
    require_positive("diameter", diameter)?;
    require_non_negative("rev_per_sec", rev_per_sec)?;
    Ok(PI * diameter * rev_per_sec)
}

/// Convert a rotational speed from revolutions per minute to
/// revolutions per second (divide by 60).
///
/// # Errors
///
/// Returns [`BeltError::BadParameter`] if `rpm` is negative.
pub fn rpm_to_rev_per_sec(rpm: f64) -> Result<f64, BeltError> {
    require_non_negative("rpm", rpm)?;
    Ok(rpm / 60.0)
}

/// Driven-pulley rotational speed implied by no-slip kinematics,
/// `N_driven = N_driver * D_driver / D_driven`.
///
/// Returned in the same time unit as `driver_speed`.
///
/// # Errors
///
/// Returns [`BeltError::BadParameter`] if either diameter is not
/// strictly positive or `driver_speed` is negative.
pub fn driven_speed(
    driver_diameter: f64,
    driven_diameter: f64,
    driver_speed: f64,
) -> Result<f64, BeltError> {
    require_non_negative("driver_speed", driver_speed)?;
    let i = speed_ratio(driver_diameter, driven_diameter)?;
    Ok(driver_speed / i)
}

/// Angle of wrap on each pulley of an **open**-belt drive, returned as
/// `(theta_small, theta_large)` in radians.
///
/// The small pulley sees `pi - 2*alpha` and the large pulley
/// `pi + 2*alpha`, where `alpha = asin((R_large - R_small) / C)`. The
/// two angles always sum to `2*pi`.
///
/// `r_small` and `r_large` are pulley radii (any consistent length
/// unit) and `center_distance` is the shaft-to-shaft spacing in the
/// same unit. The order of the two radii does not matter; the smaller
/// is treated as the small pulley.
///
/// # Errors
///
/// Returns [`BeltError::BadParameter`] for non-positive radii or a
/// non-positive centre distance, and [`BeltError::DegenerateGeometry`]
/// if the centre distance is too small for the radius difference
/// (`|R_large - R_small| >= C`), which has no real wrap geometry.
pub fn wrap_angles_open(
    r_small: f64,
    r_large: f64,
    center_distance: f64,
) -> Result<(f64, f64), BeltError> {
    require_positive("r_small", r_small)?;
    require_positive("r_large", r_large)?;
    require_positive("center_distance", center_distance)?;

    let (rs, rl) = if r_small <= r_large {
        (r_small, r_large)
    } else {
        (r_large, r_small)
    };
    let dr = rl - rs;
    if dr >= center_distance {
        return Err(BeltError::DegenerateGeometry(format!(
            "centre distance {center_distance} too small for radius difference {dr}"
        )));
    }
    let alpha = (dr / center_distance).asin();
    let theta_small = PI - 2.0 * alpha;
    let theta_large = PI + 2.0 * alpha;
    Ok((theta_small, theta_large))
}

/// Total belt length of an **open**-belt drive,
/// `L = 2*C*cos(alpha) + R_small*theta_small + R_large*theta_large`.
///
/// This is the exact geometric length (two straight tangents plus the
/// two wrapped arcs), not the common small-angle approximation. Units
/// follow the inputs.
///
/// # Errors
///
/// Same conditions as [`wrap_angles_open`].
pub fn belt_length_open(
    r_small: f64,
    r_large: f64,
    center_distance: f64,
) -> Result<f64, BeltError> {
    let (rs, rl) = if r_small <= r_large {
        (r_small, r_large)
    } else {
        (r_large, r_small)
    };
    let (theta_small, theta_large) = wrap_angles_open(r_small, r_large, center_distance)?;
    // alpha recovered from the small-pulley wrap: theta_small = pi - 2*alpha.
    let alpha = (PI - theta_small) / 2.0;
    let tangent = 2.0 * center_distance * alpha.cos();
    Ok(tangent + rs * theta_small + rl * theta_large)
}

/// Angle of wrap on **each** pulley of a **crossed**-belt drive (the belt
/// figure-eights so the shafts counter-rotate), in radians.
///
/// Unlike the open belt, a crossed belt wraps both pulleys by the *same*
/// angle, which is always greater than `pi`:
///
/// ```text
/// gamma = asin((r_small + r_large) / C)
/// theta = pi + 2*gamma   (on both pulleys)
/// ```
///
/// Note the wrap is set by the **sum** of the radii (vs the *difference*
/// for the open belt), so a crossed belt needs `C > r_small + r_large`.
///
/// # Errors
///
/// Returns [`BeltError::BadParameter`] for non-positive radii or centre
/// distance, and [`BeltError::DegenerateGeometry`] if
/// `r_small + r_large >= C` (the belt would have to pass through the
/// pulleys).
pub fn wrap_angle_crossed(
    r_small: f64,
    r_large: f64,
    center_distance: f64,
) -> Result<f64, BeltError> {
    require_positive("r_small", r_small)?;
    require_positive("r_large", r_large)?;
    require_positive("center_distance", center_distance)?;

    let sum = r_small + r_large;
    if sum >= center_distance {
        return Err(BeltError::DegenerateGeometry(format!(
            "centre distance {center_distance} too small for radius sum {sum}"
        )));
    }
    let gamma = (sum / center_distance).asin();
    Ok(PI + 2.0 * gamma)
}

/// Total belt length of a **crossed**-belt drive,
/// `L = 2*C*cos(gamma) + (r_small + r_large) * theta`, with
/// `theta = pi + 2*gamma` the common wrap angle.
///
/// The exact geometric length (two crossing tangents plus the two equal
/// wrapped arcs). For the same pulleys and centre distance a crossed belt
/// is always longer than the [open](belt_length_open) one (by `~4 r1 r2 /
/// C`), which is why reversing rotation costs a little extra belt.
///
/// # Errors
///
/// Same conditions as [`wrap_angle_crossed`].
pub fn belt_length_crossed(
    r_small: f64,
    r_large: f64,
    center_distance: f64,
) -> Result<f64, BeltError> {
    let theta = wrap_angle_crossed(r_small, r_large, center_distance)?;
    // gamma recovered from the wrap: theta = pi + 2*gamma.
    let gamma = (theta - PI) / 2.0;
    let tangent = 2.0 * center_distance * gamma.cos();
    Ok(tangent + (r_small + r_large) * theta)
}

/// Validate that `value` is finite and strictly greater than zero.
fn require_positive(name: &'static str, value: f64) -> Result<(), BeltError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(BeltError::bad_parameter(
            name,
            format!("must be finite and > 0, got {value}"),
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

    const EPS: f64 = 1e-9;

    #[test]
    fn speed_ratio_is_driven_over_driver() {
        // 200 mm driven on a 100 mm driver -> 2:1 reduction.
        let i = speed_ratio(100.0, 200.0).unwrap();
        assert!((i - 2.0).abs() < EPS, "ratio was {i}");
    }

    #[test]
    fn speed_ratio_below_one_is_an_overdrive() {
        let i = speed_ratio(200.0, 100.0).unwrap();
        assert!((i - 0.5).abs() < EPS, "ratio was {i}");
    }

    #[test]
    fn speed_ratio_rejects_nonpositive_diameter() {
        assert_eq!(
            speed_ratio(0.0, 100.0).unwrap_err().code(),
            "beltdrive.bad_parameter"
        );
        assert!(speed_ratio(100.0, -5.0).is_err());
    }

    #[test]
    fn belt_speed_matches_pi_d_n() {
        // D = 0.1 m, N = 10 rev/s -> v = pi * 0.1 * 10 = pi.
        let v = belt_speed(0.1, 10.0).unwrap();
        assert!((v - PI).abs() < EPS, "speed was {v}");
    }

    #[test]
    fn belt_speed_scales_linearly_with_speed() {
        let v1 = belt_speed(0.2, 5.0).unwrap();
        let v2 = belt_speed(0.2, 10.0).unwrap();
        assert!((v2 - 2.0 * v1).abs() < EPS);
    }

    #[test]
    fn belt_speed_zero_speed_is_zero() {
        let v = belt_speed(0.25, 0.0).unwrap();
        assert!(v.abs() < EPS);
    }

    #[test]
    fn rpm_conversion_divides_by_sixty() {
        let n = rpm_to_rev_per_sec(1800.0).unwrap();
        assert!((n - 30.0).abs() < EPS, "rev/s was {n}");
    }

    #[test]
    fn rim_speeds_match_at_both_pulleys() {
        // No-slip: pi*D1*N1 must equal pi*D2*N2.
        let d1 = 0.1;
        let d2 = 0.25;
        let n1 = 20.0; // rev/s
        let n2 = driven_speed(d1, d2, n1).unwrap();
        let v1 = belt_speed(d1, n1).unwrap();
        let v2 = belt_speed(d2, n2).unwrap();
        assert!((v1 - v2).abs() < EPS, "v1={v1} v2={v2}");
    }

    #[test]
    fn driven_speed_obeys_inverse_ratio() {
        // 2:1 reduction -> driven turns at half the driver speed.
        let n2 = driven_speed(100.0, 200.0, 1000.0).unwrap();
        assert!((n2 - 500.0).abs() < EPS, "driven speed was {n2}");
    }

    #[test]
    fn equal_pulleys_wrap_half_turn_each() {
        // R_small == R_large -> alpha = 0 -> both wraps = pi.
        let (ts, tl) = wrap_angles_open(0.05, 0.05, 0.5).unwrap();
        assert!((ts - PI).abs() < EPS, "small wrap was {ts}");
        assert!((tl - PI).abs() < EPS, "large wrap was {tl}");
    }

    #[test]
    fn open_belt_wrap_angles_sum_to_two_pi() {
        let (ts, tl) = wrap_angles_open(0.05, 0.15, 0.6).unwrap();
        assert!((ts + tl - 2.0 * PI).abs() < EPS, "sum was {}", ts + tl);
        // Small pulley wraps less than half a turn, large more.
        assert!(ts < PI);
        assert!(tl > PI);
    }

    #[test]
    fn open_belt_wrap_matches_closed_form_alpha() {
        // R_large - R_small = 0.1, C = 0.2 -> alpha = asin(0.5) = 30 deg.
        let (ts, tl) = wrap_angles_open(0.1, 0.2, 0.2).unwrap();
        let alpha = (0.5f64).asin();
        assert!((ts - (PI - 2.0 * alpha)).abs() < EPS, "small wrap {ts}");
        assert!((tl - (PI + 2.0 * alpha)).abs() < EPS, "large wrap {tl}");
    }

    #[test]
    fn wrap_radius_order_is_symmetric() {
        // Swapping the two radii must not change the (small, large) pair.
        let a = wrap_angles_open(0.05, 0.15, 0.6).unwrap();
        let b = wrap_angles_open(0.15, 0.05, 0.6).unwrap();
        assert!((a.0 - b.0).abs() < EPS);
        assert!((a.1 - b.1).abs() < EPS);
    }

    #[test]
    fn wrap_rejects_too_small_center_distance() {
        // |R_large - R_small| = 0.1 but C = 0.05 < 0.1 -> degenerate
        // (the centre distance cannot clear the radius difference).
        let err = wrap_angles_open(0.05, 0.15, 0.05).unwrap_err();
        assert_eq!(err.code(), "beltdrive.degenerate_geometry");
    }

    #[test]
    fn equal_pulley_belt_length_is_two_arcs_plus_two_spans() {
        // R = 0.05, C = 0.5: L = 2*C + 2*pi*R (alpha = 0).
        let l = belt_length_open(0.05, 0.05, 0.5).unwrap();
        let expected = 2.0 * 0.5 + 2.0 * PI * 0.05;
        assert!((l - expected).abs() < EPS, "length was {l}");
    }

    #[test]
    fn belt_length_exceeds_straight_span() {
        // The belt must be longer than two straight centre spans.
        let l = belt_length_open(0.05, 0.15, 0.6).unwrap();
        assert!(l > 2.0 * 0.6, "length {l} not above 2*C");
    }

    #[test]
    fn crossed_wrap_exceeds_pi_and_matches_closed_form() {
        // r_small + r_large = 0.3, C = 0.6 -> gamma = asin(0.5) = 30 deg,
        // theta = pi + 2*(pi/6) = 4*pi/3.
        let theta = wrap_angle_crossed(0.1, 0.2, 0.6).unwrap();
        let gamma = (0.5f64).asin();
        assert!((theta - (PI + 2.0 * gamma)).abs() < EPS, "theta {theta}");
        assert!(theta > PI, "crossed wrap must exceed pi");
    }

    #[test]
    fn crossed_belt_is_longer_than_open_for_same_geometry() {
        // Reversing rotation (crossed) costs ~4 r1 r2 / C extra belt.
        let (rs, rl, c) = (0.05, 0.15, 0.6);
        let open = belt_length_open(rs, rl, c).unwrap();
        let crossed = belt_length_crossed(rs, rl, c).unwrap();
        assert!(
            crossed > open,
            "crossed {crossed} should exceed open {open}"
        );
        // The exact difference matches the textbook ~4 r1 r2 / C to leading
        // order (loose bound covers the higher-order terms).
        assert!((crossed - open - 4.0 * rs * rl / c).abs() < 1e-3);
    }

    #[test]
    fn crossed_equal_pulleys_length_matches_construction() {
        // R = 0.05, C = 0.5: gamma = asin(0.2), L = 2C cos(gamma) + 2R*theta.
        let l = belt_length_crossed(0.05, 0.05, 0.5).unwrap();
        let gamma = (0.2f64).asin();
        let expected = 2.0 * 0.5 * gamma.cos() + 0.1 * (PI + 2.0 * gamma);
        assert!((l - expected).abs() < EPS, "length {l} vs {expected}");
        assert!(l > 2.0 * 0.5, "crossed belt below 2*C");
    }

    #[test]
    fn crossed_rejects_center_distance_below_radius_sum() {
        // r_small + r_large = 0.2 but C = 0.15 < 0.2 -> degenerate.
        let err = wrap_angle_crossed(0.05, 0.15, 0.15).unwrap_err();
        assert_eq!(err.code(), "beltdrive.degenerate_geometry");
        assert!(belt_length_crossed(0.05, 0.15, 0.15).is_err());
    }
}
