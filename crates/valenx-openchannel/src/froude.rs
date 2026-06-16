//! Froude number, flow regime, specific energy and critical depth.
//!
//! ## Model
//!
//! The Froude number compares flow inertia to the speed of a small
//! gravity (surface) wave. For an open channel it is
//!
//! `Fr = v / sqrt(g D)`
//!
//! where `v` is the mean velocity, `g` the gravitational acceleration,
//! and `D = A / T` the hydraulic (mean) depth (which equals the flow
//! depth `y` for a rectangular channel). The flow is
//!
//! - **subcritical** when `Fr < 1` (tranquil, downstream-controlled),
//! - **critical** when `Fr = 1`,
//! - **supercritical** when `Fr > 1` (rapid, upstream-controlled).
//!
//! The specific energy (energy head referenced to the channel bed) is
//!
//! `E = y + v^2 / (2 g)`
//!
//! and, for a fixed discharge `Q`, is minimised at the **critical
//! depth** `y_c`, where `Fr = 1`. Critical flow satisfies the general
//! relation
//!
//! `Q^2 T(y_c) = g A(y_c)^3`
//!
//! which this module solves by bisection. For a rectangular channel of
//! width `b` it reduces to the closed form
//! `y_c = (q^2 / g)^(1/3)` with unit discharge `q = Q / b`.

use crate::error::OpenChannelError;
use crate::geometry::Channel;

/// Standard gravitational acceleration in m/s² used throughout this
/// crate (`9.81`, the conventional textbook value).
pub const GRAVITY_M_S2: f64 = 9.81;

/// Open-channel flow regime classified by the Froude number.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FlowRegime {
    /// `Fr < 1` — tranquil, gravity waves can travel upstream.
    Subcritical,
    /// `Fr` within tolerance of `1` — critical flow.
    Critical,
    /// `Fr > 1` — rapid, disturbances are swept downstream only.
    Supercritical,
}

/// Froude number `Fr = v / sqrt(g * hydraulic_depth_m)` from a velocity
/// and a length scale directly.
///
/// Pass the hydraulic depth `D = A/T` as the length scale (it is the
/// flow depth `y` for a rectangular channel).
///
/// # Errors
///
/// Returns [`OpenChannelError`] if `velocity_m_s` is not finite, if
/// `hydraulic_depth_m` is not finite-positive, or if `gravity_m_s2` is
/// not finite-positive.
pub fn froude_number(
    velocity_m_s: f64,
    hydraulic_depth_m: f64,
    gravity_m_s2: f64,
) -> Result<f64, OpenChannelError> {
    let v = OpenChannelError::not_finite("velocity_m_s", velocity_m_s)?;
    let d = OpenChannelError::non_positive("hydraulic_depth_m", hydraulic_depth_m)?;
    let g = OpenChannelError::non_positive("gravity_m_s2", gravity_m_s2)?;
    Ok(v / (g * d).sqrt())
}

/// Classify a Froude number into a [`FlowRegime`], treating values
/// within `tol` of `1.0` as [`FlowRegime::Critical`].
///
/// `tol` must be non-negative; a typical choice is `1e-6`.
///
/// # Errors
///
/// Returns [`OpenChannelError`] if `froude` is not finite or if `tol` is
/// negative / non-finite.
pub fn classify_regime(froude: f64, tol: f64) -> Result<FlowRegime, OpenChannelError> {
    let fr = OpenChannelError::not_finite("froude", froude)?;
    let tol = OpenChannelError::negative("tol", tol)?;
    if (fr - 1.0).abs() <= tol {
        Ok(FlowRegime::Critical)
    } else if fr < 1.0 {
        Ok(FlowRegime::Subcritical)
    } else {
        Ok(FlowRegime::Supercritical)
    }
}

/// Specific energy `E = y + v^2 / (2 g)` in metres for a known velocity
/// at flow depth `depth_m`.
///
/// # Errors
///
/// Returns [`OpenChannelError`] if `depth_m` is not finite-positive,
/// `velocity_m_s` is not finite, or `gravity_m_s2` is not
/// finite-positive.
pub fn specific_energy_m(
    depth_m: f64,
    velocity_m_s: f64,
    gravity_m_s2: f64,
) -> Result<f64, OpenChannelError> {
    let y = OpenChannelError::non_positive("depth_m", depth_m)?;
    let v = OpenChannelError::not_finite("velocity_m_s", velocity_m_s)?;
    let g = OpenChannelError::non_positive("gravity_m_s2", gravity_m_s2)?;
    Ok(y + v * v / (2.0 * g))
}

/// Specific energy in metres for a [`Channel`] carrying discharge
/// `discharge_m3s` at flow depth `depth_m`, using `v = Q / A`.
///
/// # Errors
///
/// Returns [`OpenChannelError`] if `depth_m` is not finite-positive,
/// `discharge_m3s` is not finite-positive, or `gravity_m_s2` is not
/// finite-positive.
pub fn specific_energy_for_discharge_m(
    channel: &Channel,
    discharge_m3s: f64,
    depth_m: f64,
    gravity_m_s2: f64,
) -> Result<f64, OpenChannelError> {
    let q = OpenChannelError::non_positive("discharge_m3s", discharge_m3s)?;
    let a = channel.area_m2(depth_m)?;
    let v = q / a;
    specific_energy_m(depth_m, v, gravity_m_s2)
}

/// Froude number for a [`Channel`] carrying discharge `discharge_m3s` at
/// flow depth `depth_m`.
///
/// Uses `v = Q / A` and the hydraulic depth `D = A / T` as the wave
/// length scale.
///
/// # Errors
///
/// Returns [`OpenChannelError`] if `depth_m` is not finite-positive,
/// `discharge_m3s` is not finite-positive, or `gravity_m_s2` is not
/// finite-positive.
pub fn froude_for_discharge(
    channel: &Channel,
    discharge_m3s: f64,
    depth_m: f64,
    gravity_m_s2: f64,
) -> Result<f64, OpenChannelError> {
    let q = OpenChannelError::non_positive("discharge_m3s", discharge_m3s)?;
    let a = channel.area_m2(depth_m)?;
    let d = channel.hydraulic_depth_m(depth_m)?;
    let v = q / a;
    froude_number(v, d, gravity_m_s2)
}

/// Critical depth `y_c` for a [`Channel`] carrying discharge
/// `discharge_m3s`: the depth at which `Fr = 1`, found by solving
/// `Q^2 T = g A^3` via bisection.
///
/// The function `phi(y) = g A(y)^3 - Q^2 T(y)` is negative at small `y`
/// and positive at large `y` for these prismatic shapes, so the root is
/// unique and bracketed in `(0, max_depth_m]`.
///
/// # Errors
///
/// Returns [`OpenChannelError`] for an out-of-domain input, or
/// [`OpenChannelError::Convergence`] if critical flow cannot be
/// bracketed below `max_depth_m`.
pub fn critical_depth(
    channel: &Channel,
    discharge_m3s: f64,
    gravity_m_s2: f64,
    max_depth_m: f64,
) -> Result<f64, OpenChannelError> {
    let q = OpenChannelError::non_positive("discharge_m3s", discharge_m3s)?;
    let g = OpenChannelError::non_positive("gravity_m_s2", gravity_m_s2)?;
    let max_depth = OpenChannelError::non_positive("max_depth_m", max_depth_m)?;
    let q2 = q * q;

    // phi(y) = g A^3 - Q^2 T. Root at Fr = 1.
    let phi = |y: f64| -> Result<f64, OpenChannelError> {
        let a = channel.area_m2(y)?;
        let t = channel.top_width_m(y)?;
        Ok(g * a * a * a - q2 * t)
    };

    let mut lo = 1e-9_f64;
    let mut hi = max_depth;
    let phi_hi = phi(hi)?;
    if phi_hi < 0.0 {
        return Err(OpenChannelError::Convergence(format!(
            "critical flow not reached at max_depth_m = {max_depth} for Q = {q}"
        )));
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        let phi_mid = phi(mid)?;
        if phi_mid.abs() < 1e-12 || (hi - lo) < 1e-12 {
            return Ok(mid);
        }
        if phi_mid > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Ok(0.5 * (lo + hi))
}

/// Sequent (conjugate) depth of a hydraulic jump in a **rectangular**
/// channel, from the Bélanger momentum equation:
///
/// `y2 = (y1 / 2) * (sqrt(1 + 8 * Fr1^2) - 1)`
///
/// where `upstream_depth_m` is the upstream flow depth `y1` and
/// `upstream_froude` the upstream Froude number `Fr1`. A physical jump
/// runs from a supercritical state (`Fr1 > 1`, giving `y2 > y1` — the
/// flow deepens and slows) toward subcritical; at `Fr1 = 1` the two
/// depths coincide (`y2 = y1`, no jump). The two depths are *conjugate*:
/// feeding the downstream state back through the same relation recovers
/// `y1`.
///
/// This is the rectangular (constant-width) closed form, which follows
/// from equating the specific force on either side of the jump;
/// trapezoidal or other sections need the full momentum-function balance
/// and are not covered here.
///
/// # Errors
///
/// Returns [`OpenChannelError`] if `upstream_depth_m` or `upstream_froude`
/// is not finite and strictly positive.
pub fn sequent_depth(upstream_depth_m: f64, upstream_froude: f64) -> Result<f64, OpenChannelError> {
    let y1 = OpenChannelError::non_positive("upstream_depth_m", upstream_depth_m)?;
    let fr1 = OpenChannelError::non_positive("upstream_froude", upstream_froude)?;
    Ok(0.5 * y1 * ((1.0 + 8.0 * fr1 * fr1).sqrt() - 1.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn froude_number_matches_formula() {
        // v = 2, D = 1, g = 9.81  ->  Fr = 2 / sqrt(9.81).
        let fr = froude_number(2.0, 1.0, GRAVITY_M_S2).unwrap();
        let expected = 2.0 / GRAVITY_M_S2.sqrt();
        assert!((fr - expected).abs() < EPS);
    }

    #[test]
    fn froude_unity_when_velocity_equals_wave_speed() {
        // v = sqrt(g D) makes Fr exactly 1.
        let g = GRAVITY_M_S2;
        let d = 0.75;
        let v = (g * d).sqrt();
        let fr = froude_number(v, d, g).unwrap();
        assert!((fr - 1.0).abs() < EPS);
    }

    #[test]
    fn regime_classification() {
        assert_eq!(classify_regime(0.5, 1e-6).unwrap(), FlowRegime::Subcritical);
        assert_eq!(classify_regime(1.0, 1e-6).unwrap(), FlowRegime::Critical);
        assert_eq!(
            classify_regime(2.0, 1e-6).unwrap(),
            FlowRegime::Supercritical
        );
        // Within tolerance counts as critical.
        assert_eq!(
            classify_regime(1.0 + 5e-7, 1e-6).unwrap(),
            FlowRegime::Critical
        );
    }

    #[test]
    fn specific_energy_matches_formula() {
        // E = y + v^2 / (2 g).
        let e = specific_energy_m(1.5, 3.0, GRAVITY_M_S2).unwrap();
        let expected = 1.5 + 9.0 / (2.0 * GRAVITY_M_S2);
        assert!((e - expected).abs() < EPS);
    }

    #[test]
    fn rectangular_critical_depth_matches_closed_form() {
        // Rectangular b = 5, Q = 10  ->  q = 2 m²/s (unit discharge).
        // y_c = (q^2 / g)^(1/3) = (4 / 9.81)^(1/3).
        let ch = Channel::rectangular(5.0).unwrap();
        let q = 10.0;
        let yc = critical_depth(&ch, q, GRAVITY_M_S2, 10.0).unwrap();
        let unit_q = q / 5.0;
        let yc_closed = (unit_q * unit_q / GRAVITY_M_S2).powf(1.0 / 3.0);
        assert!((yc - yc_closed).abs() < 1e-7);
    }

    #[test]
    fn froude_is_unity_at_critical_depth() {
        // The depth returned by `critical_depth` must give Fr ≈ 1.
        let ch = Channel::trapezoidal(3.0, 1.5).unwrap();
        let q = 12.0;
        let yc = critical_depth(&ch, q, GRAVITY_M_S2, 10.0).unwrap();
        let fr = froude_for_discharge(&ch, q, yc, GRAVITY_M_S2).unwrap();
        assert!((fr - 1.0).abs() < 1e-6);
    }

    #[test]
    fn specific_energy_is_minimised_at_critical_depth() {
        // E(y_c) should be a local minimum: lower than E at y_c ± δ.
        let ch = Channel::rectangular(4.0).unwrap();
        let q = 8.0;
        let yc = critical_depth(&ch, q, GRAVITY_M_S2, 10.0).unwrap();
        let e_c = specific_energy_for_discharge_m(&ch, q, yc, GRAVITY_M_S2).unwrap();
        let e_lo = specific_energy_for_discharge_m(&ch, q, yc - 0.05, GRAVITY_M_S2).unwrap();
        let e_hi = specific_energy_for_discharge_m(&ch, q, yc + 0.05, GRAVITY_M_S2).unwrap();
        assert!(e_c < e_lo);
        assert!(e_c < e_hi);
    }

    #[test]
    fn rectangular_critical_specific_energy_is_three_halves_yc() {
        // For a rectangular channel, E_min = (3/2) y_c exactly.
        let ch = Channel::rectangular(6.0).unwrap();
        let q = 9.0;
        let yc = critical_depth(&ch, q, GRAVITY_M_S2, 10.0).unwrap();
        let e_c = specific_energy_for_discharge_m(&ch, q, yc, GRAVITY_M_S2).unwrap();
        assert!((e_c - 1.5 * yc).abs() < 1e-6);
    }

    #[test]
    fn critical_depth_unreachable_errors() {
        let ch = Channel::rectangular(1.0).unwrap();
        let err = critical_depth(&ch, 1.0e6, GRAVITY_M_S2, 0.01).unwrap_err();
        assert!(matches!(err, OpenChannelError::Convergence(_)));
    }

    #[test]
    fn rejects_bad_inputs() {
        let ch = Channel::rectangular(3.0).unwrap();
        assert!(froude_number(2.0, 0.0, GRAVITY_M_S2).is_err());
        assert!(froude_number(f64::NAN, 1.0, GRAVITY_M_S2).is_err());
        assert!(specific_energy_m(1.0, 2.0, -9.81).is_err());
        assert!(critical_depth(&ch, -1.0, GRAVITY_M_S2, 10.0).is_err());
        assert!(classify_regime(1.0, -1e-6).is_err());
    }

    #[test]
    fn sequent_depth_matches_belanger_closed_form() {
        // y1 = 0.5 m, Fr1 = 3 -> y2 = 0.5 * 0.5 * (sqrt(73) - 1).
        let y2 = sequent_depth(0.5, 3.0).unwrap();
        let expected = 0.5 * 0.5 * ((1.0 + 8.0 * 9.0_f64).sqrt() - 1.0);
        assert!((y2 - expected).abs() < EPS, "y2 = {y2}, want {expected}");
        // Supercritical upstream deepens the flow.
        assert!(y2 > 0.5, "jump should deepen: {y2}");
    }

    #[test]
    fn no_jump_at_critical_froude() {
        // GOLD: Fr1 = 1 is the conjugate of itself, so y2 == y1.
        let y2 = sequent_depth(1.3, 1.0).unwrap();
        assert!((y2 - 1.3).abs() < EPS, "y2 = {y2} should equal y1 at Fr=1");
        // Subcritical upstream gives a shallower (supercritical) conjugate.
        let y2_sub = sequent_depth(1.3, 0.5).unwrap();
        assert!(
            y2_sub < 1.3,
            "subcritical conjugate should be shallower: {y2_sub}"
        );
    }

    #[test]
    fn sequent_depth_is_conjugate_reciprocal() {
        // GOLD reciprocity: from continuity (q = V*y const) the downstream
        // Froude is Fr2 = Fr1 * (y1/y2)^(3/2), and the Bélanger relation
        // applied to the downstream state returns y1.
        for &fr1 in &[1.5, 2.0, 4.0, 6.0] {
            let y1 = 0.8;
            let y2 = sequent_depth(y1, fr1).unwrap();
            let fr2 = fr1 * (y1 / y2).powf(1.5);
            assert!(fr2 < 1.0, "downstream of a jump must be subcritical: {fr2}");
            let y1_back = sequent_depth(y2, fr2).unwrap();
            assert!(
                (y1_back - y1).abs() < 1e-9,
                "conjugate round-trip {y1_back} vs {y1} at Fr1={fr1}"
            );
        }
    }

    #[test]
    fn sequent_depth_increases_with_upstream_froude() {
        let mut prev = 0.0;
        for k in 1..20 {
            let fr1 = 1.0 + k as f64 * 0.25;
            let y2 = sequent_depth(0.6, fr1).unwrap();
            assert!(y2 > prev, "y2 not increasing at Fr1={fr1}: {y2} <= {prev}");
            prev = y2;
        }
    }

    #[test]
    fn sequent_depth_rejects_bad_inputs() {
        assert!(sequent_depth(0.0, 3.0).is_err());
        assert!(sequent_depth(-1.0, 3.0).is_err());
        assert!(sequent_depth(0.5, 0.0).is_err());
        assert!(sequent_depth(0.5, f64::NAN).is_err());
    }
}
