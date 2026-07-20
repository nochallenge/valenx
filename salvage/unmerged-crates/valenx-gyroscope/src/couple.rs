//! Core gyroscopic relations for a symmetric rigid rotor.
//!
//! A rotor with polar mass moment of inertia `I` spinning at angular
//! velocity `w` carries spin angular momentum `L = I*w` directed along
//! its axis. If that axis is swung (precessed) at angular velocity `wp`
//! about an axis perpendicular to the spin axis, the angular-momentum
//! vector must change at the rate
//!
//! `C = dL/dt = I * w * wp`
//!
//! `C` is the magnitude of the **active** gyroscopic couple — the
//! external couple that must be applied to force the precession. By
//! Newton's third law the rotor exerts an equal and opposite **reactive**
//! couple on its bearings / frame.

use crate::error::{non_negative, positive, Result};

/// Convert a rotational speed in revolutions per minute to angular
/// velocity in radians per second: `w = 2*pi*N / 60`.
///
/// # Errors
///
/// Returns [`crate::GyroError::NonFinite`] if `rpm` is `NaN`/`Inf`.
/// Negative rpm (reverse rotation) is accepted and simply flips the sign
/// of the returned angular velocity.
#[inline]
pub fn rpm_to_rad_s(rpm: f64) -> Result<f64> {
    let rpm = crate::error::finite("rpm", rpm)?;
    Ok(rpm * (2.0 * core::f64::consts::PI) / 60.0)
}

/// Convert angular velocity in radians per second back to revolutions per
/// minute: `N = 60*w / (2*pi)`.
///
/// # Errors
///
/// Returns [`crate::GyroError::NonFinite`] if `rad_s` is `NaN`/`Inf`.
#[inline]
pub fn rad_s_to_rpm(rad_s: f64) -> Result<f64> {
    let rad_s = crate::error::finite("rad_s", rad_s)?;
    Ok(rad_s * 60.0 / (2.0 * core::f64::consts::PI))
}

/// Spin angular momentum `L = I * w` (kg*m^2/s) along the rotor axis.
///
/// # Errors
///
/// Returns an error if `inertia` is not strictly positive, or if
/// `spin` is non-finite.
#[inline]
pub fn angular_momentum(inertia: f64, spin: f64) -> Result<f64> {
    let inertia = positive("inertia", inertia)?;
    let spin = crate::error::finite("spin", spin)?;
    Ok(inertia * spin)
}

/// Gyroscopic couple magnitude `C = I * w * wp` (N*m).
///
/// All three arguments are treated as magnitudes: `inertia` and `spin`
/// must be strictly positive, and `precession` must be non-negative. In
/// the limit `wp -> 0` the couple correctly collapses to `C -> 0`.
///
/// # Errors
///
/// Returns [`crate::GyroError::NonPositive`] if `inertia` or `spin` is
/// not `> 0`, [`crate::GyroError::Negative`] if `precession` is `< 0`,
/// or [`crate::GyroError::NonFinite`] for any `NaN`/`Inf` input.
#[inline]
pub fn gyroscopic_couple(inertia: f64, spin: f64, precession: f64) -> Result<f64> {
    let inertia = positive("inertia", inertia)?;
    let spin = positive("spin", spin)?;
    let precession = non_negative("precession", precession)?;
    Ok(inertia * spin * precession)
}

/// Precession rate `wp = C / (I * w)` (rad/s) implied by an applied
/// couple `C` on a rotor with momentum `I*w`. This is the inverse of
/// [`gyroscopic_couple`].
///
/// # Errors
///
/// Returns [`crate::GyroError::NonPositive`] if `inertia` or `spin`
/// is not `> 0` (they sit in the denominator), [`crate::GyroError::Negative`]
/// if `couple` is `< 0`, or [`crate::GyroError::NonFinite`] for any
/// `NaN`/`Inf` input.
#[inline]
pub fn precession_rate(couple: f64, inertia: f64, spin: f64) -> Result<f64> {
    let couple = non_negative("couple", couple)?;
    let inertia = positive("inertia", inertia)?;
    let spin = positive("spin", spin)?;
    Ok(couple / (inertia * spin))
}

/// Spin angular velocity `w = C / (I * wp)` (rad/s) required to react a
/// given couple `C` while precessing at `wp`. Inverse of
/// [`gyroscopic_couple`] solved for the spin rate.
///
/// # Errors
///
/// Returns [`crate::GyroError::NonPositive`] if `inertia` or
/// `precession` is not `> 0` (they sit in the denominator),
/// [`crate::GyroError::Negative`] if `couple` is `< 0`, or
/// [`crate::GyroError::NonFinite`] for any `NaN`/`Inf` input.
#[inline]
pub fn spin_for_couple(couple: f64, inertia: f64, precession: f64) -> Result<f64> {
    let couple = non_negative("couple", couple)?;
    let inertia = positive("inertia", inertia)?;
    let precession = positive("precession", precession)?;
    Ok(couple / (inertia * precession))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::GyroError;
    use core::f64::consts::PI;

    const EPS: f64 = 1e-9;

    #[test]
    fn rpm_round_trips_through_rad_s() {
        // 1000 rpm = 2*pi*1000/60 rad/s, hand-worked.
        let w = rpm_to_rad_s(1000.0).unwrap();
        assert!((w - 104.719_755_119_659_77).abs() < 1e-9);
        // Inverse must return the original rpm exactly (to eps).
        let n = rad_s_to_rpm(w).unwrap();
        assert!((n - 1000.0).abs() < 1e-9);
    }

    #[test]
    fn one_rev_per_second_is_two_pi_rad_s() {
        // 60 rpm = 1 rev/s = 2*pi rad/s.
        let w = rpm_to_rad_s(60.0).unwrap();
        assert!((w - 2.0 * PI).abs() < EPS);
    }

    #[test]
    fn negative_rpm_flips_sign() {
        let w = rpm_to_rad_s(-60.0).unwrap();
        assert!((w + 2.0 * PI).abs() < EPS);
    }

    #[test]
    fn angular_momentum_is_inertia_times_spin() {
        // I = 20, w = 104.71975511965977 => L = 2094.3951... hand-worked.
        let w = rpm_to_rad_s(1000.0).unwrap();
        let l = angular_momentum(20.0, w).unwrap();
        assert!((l - 2094.3951023931954).abs() < 1e-6);
    }

    #[test]
    fn couple_matches_textbook_aeroplane_number() {
        // I = 20 kg*m^2, N = 1000 rpm, wp = v/R = 50/100 = 0.5 rad/s.
        // C = I*w*wp = 20 * 104.71975511965977 * 0.5 = 1047.1975512 N*m.
        let w = rpm_to_rad_s(1000.0).unwrap();
        let c = gyroscopic_couple(20.0, w, 0.5).unwrap();
        assert!((c - 1047.1975511965977).abs() < 1e-4);
    }

    #[test]
    fn couple_vanishes_as_precession_goes_to_zero() {
        // Limiting case wp -> 0 gives C -> 0.
        let w = rpm_to_rad_s(1000.0).unwrap();
        let c = gyroscopic_couple(20.0, w, 0.0).unwrap();
        assert!(c.abs() < EPS);
        // And a tiny precession gives a correspondingly tiny couple.
        let c_small = gyroscopic_couple(20.0, w, 1e-12).unwrap();
        assert!(c_small.abs() < 1e-6);
        assert!(c_small > 0.0);
    }

    #[test]
    fn precession_rate_is_exact_inverse_of_couple() {
        let w = rpm_to_rad_s(1000.0).unwrap();
        let c = gyroscopic_couple(20.0, w, 0.5).unwrap();
        let wp = precession_rate(c, 20.0, w).unwrap();
        assert!((wp - 0.5).abs() < EPS);
    }

    #[test]
    fn spin_for_couple_is_exact_inverse() {
        // Recover the spin rate from couple, inertia, and precession.
        let w = rpm_to_rad_s(1000.0).unwrap();
        let c = gyroscopic_couple(20.0, w, 0.5).unwrap();
        let w_back = spin_for_couple(c, 20.0, 0.5).unwrap();
        assert!((w_back - w).abs() < 1e-9);
    }

    #[test]
    fn zero_couple_gives_zero_precession() {
        let w = rpm_to_rad_s(1000.0).unwrap();
        let wp = precession_rate(0.0, 20.0, w).unwrap();
        assert!(wp.abs() < EPS);
    }

    #[test]
    fn couple_rejects_bad_inputs() {
        // Non-positive inertia / spin.
        assert_eq!(
            gyroscopic_couple(0.0, 100.0, 0.5).unwrap_err().code(),
            "gyroscope.non_positive"
        );
        assert_eq!(
            gyroscopic_couple(20.0, -1.0, 0.5).unwrap_err().code(),
            "gyroscope.non_positive"
        );
        // Negative precession (magnitude must be >= 0).
        assert_eq!(
            gyroscopic_couple(20.0, 100.0, -0.5).unwrap_err().code(),
            "gyroscope.negative"
        );
        // Non-finite.
        assert!(matches!(
            gyroscopic_couple(20.0, f64::NAN, 0.5),
            Err(GyroError::NonFinite { .. })
        ));
        assert!(matches!(
            gyroscopic_couple(f64::INFINITY, 100.0, 0.5),
            Err(GyroError::NonFinite { .. })
        ));
    }

    #[test]
    fn inverses_reject_zero_denominators() {
        // precession_rate divides by I*w.
        assert_eq!(
            precession_rate(100.0, 0.0, 100.0).unwrap_err().code(),
            "gyroscope.non_positive"
        );
        assert_eq!(
            precession_rate(100.0, 20.0, 0.0).unwrap_err().code(),
            "gyroscope.non_positive"
        );
        // spin_for_couple divides by I*wp.
        assert_eq!(
            spin_for_couple(100.0, 20.0, 0.0).unwrap_err().code(),
            "gyroscope.non_positive"
        );
        // Negative couple rejected.
        assert_eq!(
            precession_rate(-1.0, 20.0, 100.0).unwrap_err().code(),
            "gyroscope.negative"
        );
    }

    #[test]
    fn rpm_conversions_reject_non_finite() {
        assert!(matches!(
            rpm_to_rad_s(f64::NAN),
            Err(GyroError::NonFinite { .. })
        ));
        assert!(matches!(
            rad_s_to_rpm(f64::INFINITY),
            Err(GyroError::NonFinite { .. })
        ));
    }
}
