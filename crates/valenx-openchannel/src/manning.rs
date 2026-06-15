//! Manning's equation for steady uniform open-channel flow.
//!
//! ## Model
//!
//! Under the steady-uniform-flow (normal-flow) assumption the mean
//! velocity in SI units is given by Manning's equation
//!
//! `v = (1/n) R^(2/3) S^(1/2)`
//!
//! and the volumetric discharge by
//!
//! `Q = v A = (1/n) A R^(2/3) S^(1/2)`
//!
//! where
//!
//! - `n` is Manning's roughness coefficient (dimensionless in the SI
//!   form, s·m^(-1/3) when carrying units),
//! - `A` is the flow area in m²,
//! - `R = A / P` is the hydraulic radius in m,
//! - `S` is the channel-bed (energy-line) slope, dimensionless (m/m).
//!
//! The `1/n` form here is the SI form; the US-customary form carries an
//! extra `1.49` factor, which this crate deliberately does not use — all
//! inputs and outputs are SI.
//!
//! ## Honest scope
//!
//! This is the textbook empirical normal-flow relation. It assumes a
//! prismatic channel, fully turbulent rough flow, a constant roughness
//! around the wetted perimeter, and a mild bed slope so the depth-normal
//! and vertical directions coincide. It is a learning / first-estimate
//! tool, not a calibrated hydraulic-design or flood model.

use crate::error::OpenChannelError;
use crate::geometry::Channel;

/// Manning mean velocity `v = (1/n) R^(2/3) S^(1/2)` in m/s.
///
/// Takes the hydraulic radius and slope directly (no geometry), so it
/// can be reused for any cross-section whose `R` is already known.
///
/// # Errors
///
/// Returns [`OpenChannelError`] if `manning_n` is not finite-positive,
/// if `hydraulic_radius_m` is not finite-positive, or if `slope` is
/// negative / non-finite (a zero slope yields a zero velocity).
pub fn velocity_from_radius(
    manning_n: f64,
    hydraulic_radius_m: f64,
    slope: f64,
) -> Result<f64, OpenChannelError> {
    let n = OpenChannelError::non_positive("manning_n", manning_n)?;
    let r = OpenChannelError::non_positive("hydraulic_radius_m", hydraulic_radius_m)?;
    let s = OpenChannelError::negative("slope", slope)?;
    Ok((1.0 / n) * r.powf(2.0 / 3.0) * s.sqrt())
}

/// Manning mean velocity in m/s for a [`Channel`] flowing at `depth_m`.
///
/// Equivalent to [`velocity_from_radius`] with `R = A/P` taken from the
/// channel geometry at the given depth.
///
/// # Errors
///
/// Returns [`OpenChannelError`] if any of `depth_m`, `manning_n` or
/// `slope` violates its domain (see [`velocity_from_radius`] and the
/// geometry accessors).
pub fn velocity(
    channel: &Channel,
    depth_m: f64,
    manning_n: f64,
    slope: f64,
) -> Result<f64, OpenChannelError> {
    let r = channel.hydraulic_radius_m(depth_m)?;
    velocity_from_radius(manning_n, r, slope)
}

/// Manning discharge `Q = (1/n) A R^(2/3) S^(1/2)` in m³/s for a
/// [`Channel`] flowing at `depth_m`.
///
/// # Errors
///
/// Returns [`OpenChannelError`] if any of `depth_m`, `manning_n` or
/// `slope` violates its domain.
pub fn discharge(
    channel: &Channel,
    depth_m: f64,
    manning_n: f64,
    slope: f64,
) -> Result<f64, OpenChannelError> {
    let a = channel.area_m2(depth_m)?;
    let v = velocity(channel, depth_m, manning_n, slope)?;
    Ok(v * a)
}

/// Solve for the **normal depth** `y_n`: the flow depth at which
/// Manning's equation delivers the target discharge `target_q_m3s` in
/// the given channel for the given roughness and slope.
///
/// Discharge increases monotonically with depth for these prismatic
/// shapes, so the root is unique and found by bisection.
///
/// # Errors
///
/// Returns [`OpenChannelError::NonPositive`] / [`OpenChannelError::Negative`]
/// for an out-of-domain input, or [`OpenChannelError::Convergence`] if
/// the target discharge cannot be bracketed below `max_depth_m`.
pub fn normal_depth(
    channel: &Channel,
    target_q_m3s: f64,
    manning_n: f64,
    slope: f64,
    max_depth_m: f64,
) -> Result<f64, OpenChannelError> {
    let q_target = OpenChannelError::non_positive("target_q_m3s", target_q_m3s)?;
    // Validate the slope / roughness once up front via a probe evaluation
    // at a small depth (also surfaces a zero slope, which makes Q == 0
    // and the target unreachable).
    let _ = OpenChannelError::non_positive("manning_n", manning_n)?;
    let s = OpenChannelError::negative("slope", slope)?;
    let max_depth = OpenChannelError::non_positive("max_depth_m", max_depth_m)?;
    if s == 0.0 {
        return Err(OpenChannelError::Convergence(format!(
            "slope is zero so discharge is zero everywhere; cannot reach Q = {q_target}"
        )));
    }

    let f = |y: f64| -> Result<f64, OpenChannelError> {
        Ok(discharge(channel, y, manning_n, slope)? - q_target)
    };

    let mut lo = 1e-9_f64;
    let mut hi = max_depth;
    let f_hi = f(hi)?;
    if f_hi < 0.0 {
        return Err(OpenChannelError::Convergence(format!(
            "target discharge {q_target} not reached at max_depth_m = {max_depth}"
        )));
    }
    // f(lo) is negative for any positive target (Q -> 0 as y -> 0), so the
    // root is bracketed in (lo, hi].
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        let f_mid = f(mid)?;
        if f_mid.abs() < 1e-12 || (hi - lo) < 1e-12 {
            return Ok(mid);
        }
        if f_mid > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Ok(0.5 * (lo + hi))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn discharge_matches_hand_computed_value() {
        // Rectangular b = 3, y = 1  ->  A = 3, P = 5, R = 0.6.
        // n = 0.013, S = 0.001.
        // v = (1/0.013) * 0.6^(2/3) * 0.001^0.5
        // Q = v * A.
        let ch = Channel::rectangular(3.0).unwrap();
        let n: f64 = 0.013;
        let s: f64 = 0.001;
        let r: f64 = 0.6;
        let v_expected = (1.0 / n) * r.powf(2.0 / 3.0) * s.sqrt();
        let q_expected = v_expected * 3.0;
        let q = discharge(&ch, 1.0, n, s).unwrap();
        assert!((q - q_expected).abs() < 1e-9);
    }

    #[test]
    fn velocity_from_radius_matches_formula() {
        // Direct evaluation of v = (1/n) R^(2/3) S^(1/2).
        let n = 0.025;
        let r: f64 = 1.2;
        let s = 0.0009;
        let v = velocity_from_radius(n, r, s).unwrap();
        let expected = (1.0 / n) * r.powf(2.0 / 3.0) * s.sqrt();
        assert!((v - expected).abs() < 1e-12);
    }

    #[test]
    fn discharge_grows_with_sqrt_of_slope() {
        // Quadrupling the slope multiplies Q by exactly 2 (Q ∝ sqrt(S)).
        let ch = Channel::rectangular(3.0).unwrap();
        let q1 = discharge(&ch, 1.0, 0.013, 0.001).unwrap();
        let q4 = discharge(&ch, 1.0, 0.013, 0.004).unwrap();
        assert!((q4 / q1 - 2.0).abs() < 1e-9);
    }

    #[test]
    fn discharge_scales_inversely_with_roughness() {
        // Doubling n halves Q (Q ∝ 1/n).
        let ch = Channel::rectangular(3.0).unwrap();
        let q_smooth = discharge(&ch, 1.0, 0.013, 0.001).unwrap();
        let q_rough = discharge(&ch, 1.0, 0.026, 0.001).unwrap();
        assert!((q_smooth / q_rough - 2.0).abs() < 1e-9);
    }

    #[test]
    fn discharge_is_velocity_times_area() {
        let ch = Channel::trapezoidal(2.0, 1.0).unwrap();
        let v = velocity(&ch, 1.5, 0.02, 0.0016).unwrap();
        let a = ch.area_m2(1.5).unwrap();
        let q = discharge(&ch, 1.5, 0.02, 0.0016).unwrap();
        assert!((q - v * a).abs() < 1e-12);
    }

    #[test]
    fn zero_slope_gives_zero_velocity() {
        let ch = Channel::rectangular(3.0).unwrap();
        let v = velocity(&ch, 1.0, 0.013, 0.0).unwrap();
        assert!(v.abs() < 1e-15);
    }

    #[test]
    fn normal_depth_inverts_discharge() {
        // Pick a depth, compute Q, then recover the depth from Q.
        let ch = Channel::trapezoidal(4.0, 2.0).unwrap();
        let (n, s, y) = (0.015, 0.0005, 1.234);
        let q = discharge(&ch, y, n, s).unwrap();
        let y_back = normal_depth(&ch, q, n, s, 10.0).unwrap();
        assert!((y_back - y).abs() < 1e-6);
    }

    #[test]
    fn normal_depth_unreachable_target_errors() {
        let ch = Channel::rectangular(1.0).unwrap();
        // Absurdly large target that no depth below 0.01 m can deliver.
        let err = normal_depth(&ch, 1.0e6, 0.013, 0.001, 0.01).unwrap_err();
        assert!(matches!(err, OpenChannelError::Convergence(_)));
    }

    #[test]
    fn rejects_bad_inputs() {
        let ch = Channel::rectangular(3.0).unwrap();
        assert!(velocity(&ch, 1.0, 0.0, 0.001).is_err());
        assert!(velocity(&ch, 1.0, 0.013, -0.001).is_err());
        assert!(velocity_from_radius(0.013, -1.0, 0.001).is_err());
    }
}
