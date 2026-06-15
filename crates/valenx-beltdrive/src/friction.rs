//! Capstan / Euler belt-friction relations.
//!
//! For a belt about to slip on a pulley, the tight-side tension `T1`
//! and slack-side tension `T2` satisfy the Euler (capstan) equation
//!
//! ```text
//! T1 / T2 = exp(mu * theta)
//! ```
//!
//! where `mu` is the coefficient of friction between belt and pulley
//! and `theta` is the angle of wrap in radians. The ratio grows
//! exponentially with both the friction coefficient and the wrap angle,
//! so a larger wrap or a grippier belt raises the load a drive can
//! carry before slipping.
//!
//! ## V-belts
//!
//! A V-belt seats in a groove of total included angle `2*beta`, which
//! wedges the belt and multiplies the effective friction. The capstan
//! relation then uses an *effective* coefficient
//! `mu_eff = mu / sin(beta)`:
//!
//! ```text
//! T1 / T2 = exp((mu / sin(beta)) * theta).
//! ```

use crate::error::BeltError;

/// Capstan tension ratio `T1 / T2 = exp(mu * theta)` for a flat belt at
/// the point of slipping.
///
/// `mu` is the belt/pulley coefficient of friction (dimensionless) and
/// `wrap_angle` the angle of wrap in radians. With `mu = 0` or
/// `wrap_angle = 0` the ratio is exactly one (no friction advantage).
///
/// # Errors
///
/// Returns [`BeltError::BadParameter`] if `mu` is negative or
/// `wrap_angle` is negative (or either is non-finite).
pub fn tension_ratio(mu: f64, wrap_angle: f64) -> Result<f64, BeltError> {
    require_non_negative("mu", mu)?;
    require_non_negative("wrap_angle", wrap_angle)?;
    Ok((mu * wrap_angle).exp())
}

/// Effective coefficient of friction for a V-belt in a groove of
/// half-angle `beta` (radians), `mu_eff = mu / sin(beta)`.
///
/// The groove wedging makes `mu_eff` larger than the flat-belt `mu`,
/// which is why V-belts transmit more power for the same wrap and
/// nominal friction.
///
/// # Errors
///
/// Returns [`BeltError::BadParameter`] if `mu` is negative or `beta` is
/// not in the open interval `(0, pi/2]` (a groove half-angle must be
/// strictly positive so `sin(beta) > 0`).
pub fn v_belt_effective_mu(mu: f64, beta: f64) -> Result<f64, BeltError> {
    require_non_negative("mu", mu)?;
    if !beta.is_finite() || beta <= 0.0 || beta > std::f64::consts::FRAC_PI_2 {
        return Err(BeltError::bad_parameter(
            "beta",
            format!("groove half-angle must be in (0, pi/2], got {beta}"),
        ));
    }
    Ok(mu / beta.sin())
}

/// Capstan tension ratio for a V-belt,
/// `T1 / T2 = exp((mu / sin(beta)) * theta)`.
///
/// Convenience composition of [`v_belt_effective_mu`] and
/// [`tension_ratio`]. For the same `mu`, `theta`, and a typical groove
/// (`beta` ~ 17–19 degrees) this exceeds the flat-belt ratio.
///
/// # Errors
///
/// Returns [`BeltError::BadParameter`] under the union of the
/// conditions in [`v_belt_effective_mu`] and [`tension_ratio`].
pub fn v_belt_tension_ratio(mu: f64, beta: f64, wrap_angle: f64) -> Result<f64, BeltError> {
    let mu_eff = v_belt_effective_mu(mu, beta)?;
    tension_ratio(mu_eff, wrap_angle)
}

/// Slack-side tension `T2` that holds a given tight-side tension `T1` at
/// the point of slipping, `T2 = T1 / exp(mu * theta)`.
///
/// # Errors
///
/// Returns [`BeltError::BadParameter`] if `t1` is negative or the
/// friction inputs are out of domain (see [`tension_ratio`]).
pub fn slack_tension(t1: f64, mu: f64, wrap_angle: f64) -> Result<f64, BeltError> {
    require_non_negative("t1", t1)?;
    let ratio = tension_ratio(mu, wrap_angle)?;
    Ok(t1 / ratio)
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
    use std::f64::consts::PI;

    const EPS: f64 = 1e-9;

    #[test]
    fn tension_ratio_matches_exp_mu_theta() {
        // mu = 0.3, theta = pi (half wrap) -> exp(0.3*pi).
        let r = tension_ratio(0.3, PI).unwrap();
        assert!((r - (0.3 * PI).exp()).abs() < EPS, "ratio was {r}");
    }

    #[test]
    fn zero_friction_gives_unit_ratio() {
        let r = tension_ratio(0.0, PI).unwrap();
        assert!((r - 1.0).abs() < EPS, "ratio was {r}");
    }

    #[test]
    fn zero_wrap_gives_unit_ratio() {
        let r = tension_ratio(0.4, 0.0).unwrap();
        assert!((r - 1.0).abs() < EPS, "ratio was {r}");
    }

    #[test]
    fn larger_friction_raises_the_ratio() {
        // Higher mu -> higher capacity (larger T1/T2).
        let low = tension_ratio(0.2, PI).unwrap();
        let high = tension_ratio(0.4, PI).unwrap();
        assert!(high > low, "high={high} low={low}");
    }

    #[test]
    fn larger_wrap_raises_the_ratio() {
        // Higher wrap angle -> higher capacity.
        let small = tension_ratio(0.3, PI - 0.5).unwrap();
        let big = tension_ratio(0.3, PI + 0.5).unwrap();
        assert!(big > small, "big={big} small={small}");
    }

    #[test]
    fn ratio_is_monotonic_and_above_one() {
        // For mu>0, theta>0 the ratio strictly exceeds 1.
        let r = tension_ratio(0.25, 2.5).unwrap();
        assert!(r > 1.0, "ratio was {r}");
    }

    #[test]
    fn tension_ratio_rejects_negative_mu() {
        let err = tension_ratio(-0.1, PI).unwrap_err();
        assert_eq!(err.code(), "beltdrive.bad_parameter");
    }

    #[test]
    fn tension_ratio_rejects_negative_wrap() {
        assert!(tension_ratio(0.3, -1.0).is_err());
    }

    #[test]
    fn v_belt_effective_mu_is_amplified() {
        // beta = 30 deg -> sin = 0.5 -> mu_eff = 2*mu.
        let mu_eff = v_belt_effective_mu(0.2, PI / 6.0).unwrap();
        assert!((mu_eff - 0.4).abs() < EPS, "mu_eff was {mu_eff}");
    }

    #[test]
    fn v_belt_ratio_beats_flat_belt() {
        // Same mu, wrap; V-groove wedging gives a higher ratio.
        let flat = tension_ratio(0.2, PI).unwrap();
        let vee = v_belt_tension_ratio(0.2, PI / 9.0, PI).unwrap(); // beta = 20 deg
        assert!(vee > flat, "vee={vee} flat={flat}");
    }

    #[test]
    fn v_belt_effective_mu_rejects_zero_beta() {
        let err = v_belt_effective_mu(0.2, 0.0).unwrap_err();
        assert_eq!(err.code(), "beltdrive.bad_parameter");
    }

    #[test]
    fn v_belt_effective_mu_rejects_beta_above_ninety() {
        assert!(v_belt_effective_mu(0.2, PI).is_err());
    }

    #[test]
    fn slack_tension_recovers_the_ratio() {
        // T2 = T1 / exp(mu*theta); then T1/T2 must equal the ratio.
        let t1 = 1000.0;
        let t2 = slack_tension(t1, 0.3, PI).unwrap();
        let ratio = tension_ratio(0.3, PI).unwrap();
        assert!((t1 / t2 - ratio).abs() < 1e-6, "t1/t2={}", t1 / t2);
    }

    #[test]
    fn slack_tension_is_below_tight_tension() {
        let t2 = slack_tension(800.0, 0.35, PI).unwrap();
        assert!(t2 < 800.0, "t2 was {t2}");
    }
}
