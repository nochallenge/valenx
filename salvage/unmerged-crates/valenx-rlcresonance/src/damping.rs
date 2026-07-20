//! Transient damping of the series RLC step/impulse response.
//!
//! The source-free series RLC loop obeys the second-order ODE
//!
//! ```text
//! L d2i/dt2 + R di/dt + i/C = 0
//! ```
//!
//! whose characteristic roots are set by the neper (damping) frequency
//! `alpha = R / (2 L)` and the resonant frequency `omega0 = 1/sqrt(L C)`.
//! The dimensionless **damping ratio**
//!
//! ```text
//! zeta = alpha / omega0 = (R / 2) sqrt(C / L) = 1 / (2 Q)
//! ```
//!
//! classifies the regime, and the boundary `zeta = 1` happens at the
//! **critical resistance** `R_crit = 2 sqrt(L / C)`. These are the
//! standard textbook results and are exactly the ground truth the tests
//! pin to.

use crate::error::{require_non_negative, require_positive, Result};

/// Which transient regime a series RLC loop is in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DampingRegime {
    /// `zeta < 1`: complex-conjugate roots, the response rings (decaying
    /// sinusoid) at the damped frequency before settling.
    Underdamped,
    /// `zeta == 1`: a repeated real root, the fastest non-oscillatory
    /// settling. Occurs exactly at `R = 2 sqrt(L/C)`.
    CriticallyDamped,
    /// `zeta > 1`: two distinct real roots, a sluggish non-oscillatory
    /// decay.
    Overdamped,
}

/// Neper (damping) frequency `alpha = R / (2 L)`, in radians per second.
///
/// This is the exponential decay rate of the transient envelope.
///
/// # Errors
///
/// Returns [`crate::RlcError::Negative`] for a negative or non-finite
/// `resistance_ohm`, and [`crate::RlcError::NonPositive`] for a
/// non-positive `inductance_h`.
///
/// # Example
///
/// ```
/// use valenx_rlcresonance::neper_frequency_rad_s;
/// // R = 10, L = 1 mH -> alpha = 10/(2 * 1e-3) = 5000 rad/s.
/// let a = neper_frequency_rad_s(10.0, 1e-3).unwrap();
/// assert!((a - 5000.0).abs() < 1e-9);
/// ```
pub fn neper_frequency_rad_s(resistance_ohm: f64, inductance_h: f64) -> Result<f64> {
    let r = require_non_negative("resistance", resistance_ohm)?;
    let l = require_positive("inductance", inductance_h)?;
    Ok(r / (2.0 * l))
}

/// Damping ratio `zeta = (R / 2) sqrt(C / L)` of a series RLC loop.
///
/// Equivalently `zeta = alpha / omega0 = 1 / (2 Q)`. Dimensionless; the
/// regime boundary is `zeta = 1`.
///
/// # Errors
///
/// Returns [`crate::RlcError::Negative`] for a negative or non-finite
/// `resistance_ohm`, and [`crate::RlcError::NonPositive`] for a
/// non-positive `inductance_h` or `capacitance_f`.
///
/// # Example
///
/// ```
/// use valenx_rlcresonance::damping_ratio;
/// // R = 10, L = 1 mH, C = 1 uF -> zeta = 5 * sqrt(1e-6/1e-3) = 0.158114.
/// let z = damping_ratio(10.0, 1e-3, 1e-6).unwrap();
/// assert!((z - 0.15811388300841897).abs() < 1e-12);
/// ```
pub fn damping_ratio(resistance_ohm: f64, inductance_h: f64, capacitance_f: f64) -> Result<f64> {
    let r = require_non_negative("resistance", resistance_ohm)?;
    let l = require_positive("inductance", inductance_h)?;
    let c = require_positive("capacitance", capacitance_f)?;
    Ok((r / 2.0) * (c / l).sqrt())
}

/// Critical resistance `R_crit = 2 sqrt(L / C)` that yields `zeta = 1`.
///
/// At this resistance the series loop is critically damped: the fastest
/// settling with no overshoot.
///
/// # Errors
///
/// Returns [`crate::RlcError::NonPositive`] if `inductance_h` or
/// `capacitance_f` is not finite and `> 0`.
///
/// # Example
///
/// ```
/// use valenx_rlcresonance::critical_resistance_ohm;
/// // L = 1 mH, C = 1 uF -> R_crit = 2 sqrt(1e-3/1e-6) = 63.2456 ohm.
/// let rc = critical_resistance_ohm(1e-3, 1e-6).unwrap();
/// assert!((rc - 63.245553203367585).abs() < 1e-9);
/// ```
pub fn critical_resistance_ohm(inductance_h: f64, capacitance_f: f64) -> Result<f64> {
    let l = require_positive("inductance", inductance_h)?;
    let c = require_positive("capacitance", capacitance_f)?;
    Ok(2.0 * (l / c).sqrt())
}

/// Classify the transient [`DampingRegime`] of a series RLC loop.
///
/// Compares the damping ratio against `1` with an absolute tolerance so
/// that a resistance numerically equal to [`critical_resistance_ohm`]
/// reports [`DampingRegime::CriticallyDamped`] rather than tipping into
/// one side by a rounding bit.
///
/// # Errors
///
/// Same domain checks as [`damping_ratio`].
///
/// # Example
///
/// ```
/// use valenx_rlcresonance::{classify_regime, DampingRegime};
/// // Low R -> underdamped (rings).
/// assert_eq!(classify_regime(10.0, 1e-3, 1e-6).unwrap(), DampingRegime::Underdamped);
/// // At the critical resistance -> critically damped.
/// assert_eq!(
///     classify_regime(63.245553203367585, 1e-3, 1e-6).unwrap(),
///     DampingRegime::CriticallyDamped,
/// );
/// ```
pub fn classify_regime(
    resistance_ohm: f64,
    inductance_h: f64,
    capacitance_f: f64,
) -> Result<DampingRegime> {
    let zeta = damping_ratio(resistance_ohm, inductance_h, capacitance_f)?;
    // Relative tolerance on the unit-magnitude boundary: 1 ulp-ish band.
    const TOL: f64 = 1e-9;
    Ok(if (zeta - 1.0).abs() <= TOL {
        DampingRegime::CriticallyDamped
    } else if zeta < 1.0 {
        DampingRegime::Underdamped
    } else {
        DampingRegime::Overdamped
    })
}

/// Damped ringing frequency of an underdamped series loop, in rad/s.
///
/// `omega_d = omega0 sqrt(1 - zeta^2)`. This is the angular frequency at
/// which an underdamped circuit actually oscillates, always at or below
/// `omega0`. Defined only for `zeta < 1`; at or beyond critical damping
/// the response does not oscillate, which is reported as
/// [`crate::RlcError::Singular`].
///
/// # Errors
///
/// Same domain checks as [`damping_ratio`], plus
/// [`crate::RlcError::Singular`] when the circuit is not underdamped
/// (`zeta >= 1`), where no real damped frequency exists.
///
/// # Example
///
/// ```
/// use valenx_rlcresonance::{damped_frequency_rad_s, resonant_angular_frequency_rad_s};
/// let l = 1e-3;
/// let c = 1e-6;
/// let wd = damped_frequency_rad_s(10.0, l, c).unwrap();
/// let w0 = resonant_angular_frequency_rad_s(l, c).unwrap();
/// // omega_d is just below omega0 for light damping (zeta ~= 0.158).
/// assert!(wd < w0 && wd > 0.98 * w0);
/// ```
pub fn damped_frequency_rad_s(
    resistance_ohm: f64,
    inductance_h: f64,
    capacitance_f: f64,
) -> Result<f64> {
    let zeta = damping_ratio(resistance_ohm, inductance_h, capacitance_f)?;
    if zeta >= 1.0 {
        return Err(crate::RlcError::Singular(
            "not underdamped (zeta >= 1): no real damped frequency",
        ));
    }
    let w0 = crate::resonant_angular_frequency_rad_s(inductance_h, capacitance_f)?;
    Ok(w0 * (1.0 - zeta * zeta).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{resonant_angular_frequency_rad_s, series_quality_factor, RlcError};

    const EPS: f64 = 1e-9;

    #[test]
    fn neper_frequency_matches_textbook() {
        // R = 10, L = 1 mH -> alpha = 10/(2 * 1e-3) = 5000.
        let a = neper_frequency_rad_s(10.0, 1e-3).unwrap();
        assert!((a - 5000.0).abs() < EPS, "alpha = {a}");
    }

    #[test]
    fn damping_ratio_matches_textbook() {
        // R = 10, L = 1 mH, C = 1 uF.
        // zeta = (10/2) sqrt(1e-6/1e-3) = 5 * sqrt(1e-3) = 5*0.0316227766.
        let z = damping_ratio(10.0, 1e-3, 1e-6).unwrap();
        assert!((z - 0.15811388300841897).abs() < 1e-12, "zeta = {z}");
    }

    #[test]
    fn damping_ratio_is_alpha_over_omega0() {
        // zeta = alpha / omega0 — cross-check against the two frequencies.
        let (r, l, c) = (33.0, 4.7e-3, 1e-8);
        let z = damping_ratio(r, l, c).unwrap();
        let a = neper_frequency_rad_s(r, l).unwrap();
        let w0 = resonant_angular_frequency_rad_s(l, c).unwrap();
        assert!((z - a / w0).abs() < EPS, "zeta = {z}");
    }

    #[test]
    fn damping_ratio_is_reciprocal_of_twice_q() {
        // zeta = 1/(2 Q) for the series circuit.
        let (r, l, c) = (10.0, 1e-3, 1e-6);
        let z = damping_ratio(r, l, c).unwrap();
        let q = series_quality_factor(r, l, c).unwrap();
        assert!((z - 1.0 / (2.0 * q)).abs() < EPS, "zeta = {z}");
    }

    #[test]
    fn critical_resistance_matches_textbook() {
        // L = 1 mH, C = 1 uF -> R_crit = 2 sqrt(1e-3/1e-6) = 2*31.6227766.
        let rc = critical_resistance_ohm(1e-3, 1e-6).unwrap();
        assert!((rc - 63.245553203367585).abs() < EPS, "R_crit = {rc}");
    }

    #[test]
    fn damping_ratio_is_unity_at_critical_resistance() {
        // By construction zeta(R_crit) = 1 exactly.
        let (l, c) = (2.2e-3, 4.7e-7);
        let rc = critical_resistance_ohm(l, c).unwrap();
        let z = damping_ratio(rc, l, c).unwrap();
        assert!((z - 1.0).abs() < EPS, "zeta = {z}");
    }

    #[test]
    fn classify_underdamped_critical_overdamped() {
        let (l, c) = (1e-3, 1e-6);
        let rc = critical_resistance_ohm(l, c).unwrap();
        // Well below critical -> rings.
        assert_eq!(
            classify_regime(rc * 0.1, l, c).unwrap(),
            DampingRegime::Underdamped
        );
        // Exactly critical.
        assert_eq!(
            classify_regime(rc, l, c).unwrap(),
            DampingRegime::CriticallyDamped
        );
        // Well above critical -> sluggish.
        assert_eq!(
            classify_regime(rc * 10.0, l, c).unwrap(),
            DampingRegime::Overdamped
        );
    }

    #[test]
    fn zero_resistance_is_undamped_underdamped() {
        // R = 0 -> zeta = 0 -> lossless oscillator, classed underdamped.
        let z = damping_ratio(0.0, 1e-3, 1e-6).unwrap();
        assert!(z.abs() < EPS, "zeta = {z}");
        assert_eq!(
            classify_regime(0.0, 1e-3, 1e-6).unwrap(),
            DampingRegime::Underdamped
        );
    }

    #[test]
    fn damped_frequency_below_omega0_and_correct_value() {
        // R = 10, L = 1 mH, C = 1 uF. zeta = 0.158113883.
        // omega_d = omega0 sqrt(1 - zeta^2).
        let (r, l, c) = (10.0, 1e-3, 1e-6);
        let w0 = resonant_angular_frequency_rad_s(l, c).unwrap();
        let z = damping_ratio(r, l, c).unwrap();
        let expected = w0 * (1.0 - z * z).sqrt();
        let wd = damped_frequency_rad_s(r, l, c).unwrap();
        assert!((wd - expected).abs() < 1e-6, "omega_d = {wd}");
        assert!(wd < w0, "omega_d should be below omega0");
    }

    #[test]
    fn damped_frequency_undamped_limit_equals_omega0() {
        // R = 0 -> zeta = 0 -> omega_d = omega0.
        let (l, c) = (1e-3, 1e-6);
        let w0 = resonant_angular_frequency_rad_s(l, c).unwrap();
        let wd = damped_frequency_rad_s(0.0, l, c).unwrap();
        assert!((wd - w0).abs() < 1e-6, "omega_d = {wd}");
    }

    #[test]
    fn damped_frequency_rejects_critical_and_overdamped() {
        let (l, c) = (1e-3, 1e-6);
        let rc = critical_resistance_ohm(l, c).unwrap();
        // At critical damping there is no real damped frequency.
        let err = damped_frequency_rad_s(rc, l, c).unwrap_err();
        assert_eq!(err.code(), "rlc.singular");
        // And overdamped likewise.
        assert!(matches!(
            damped_frequency_rad_s(rc * 5.0, l, c),
            Err(RlcError::Singular(_))
        ));
    }

    #[test]
    fn damping_constructors_reject_bad_input() {
        assert_eq!(
            neper_frequency_rad_s(-1.0, 1e-3).unwrap_err().code(),
            "rlc.negative"
        );
        assert_eq!(
            damping_ratio(10.0, 0.0, 1e-6).unwrap_err().code(),
            "rlc.non_positive"
        );
        assert_eq!(
            critical_resistance_ohm(1e-3, f64::NAN).unwrap_err().code(),
            "rlc.non_positive"
        );
        assert!(matches!(
            classify_regime(10.0, 1e-3, -1.0),
            Err(RlcError::NonPositive {
                name: "capacitance",
                ..
            })
        ));
    }
}
