//! Disc / caliper brake friction torque.
//!
//! A caliper clamps the rotor with a normal (clamp) force `F`. Each of
//! the `n_pads` friction faces develops a Coulomb friction force
//! `mu * F` acting at the effective (mean) radius `r_eff`, so the
//! retarding torque is
//!
//! ```text
//! T = mu * F * n_pads * r_eff
//! ```
//!
//! The relation is linear in `F`, `n_pads` and `r_eff`, so it inverts
//! cleanly: [`clamp_force_for_torque`] and [`effective_radius_for_torque`]
//! solve for the clamp force or radius that yields a target torque.

use crate::error::{check_count, check_friction, check_positive, BrakeError};

/// Retarding torque of a disc / caliper brake, in newton-metres.
///
/// `T = mu * F * n_pads * r_eff`.
///
/// # Parameters
/// - `mu` — friction coefficient, in `(0, MU_MAX]`.
/// - `clamp_force_n` — normal clamp force per pad face `F`, in newtons (> 0).
/// - `n_pads` — number of friction faces (>= 1; a floating single-piston
///   caliper has 2).
/// - `r_eff_m` — effective (mean) friction radius `r_eff`, in metres (> 0).
///
/// # Errors
/// Returns a [`BrakeError`] if `mu` is out of range, any continuous
/// parameter is non-finite or non-positive, or `n_pads == 0`.
///
/// # Examples
/// ```
/// use valenx_brake::disc::disc_torque;
/// let t = disc_torque(0.4, 8_000.0, 2, 0.12).unwrap();
/// assert!((t - 768.0).abs() < 1e-9);
/// ```
pub fn disc_torque(
    mu: f64,
    clamp_force_n: f64,
    n_pads: u32,
    r_eff_m: f64,
) -> Result<f64, BrakeError> {
    let mu = check_friction(mu)?;
    let f = check_positive("clamp_force_n", clamp_force_n)?;
    let n = check_count("n_pads", n_pads)?;
    let r = check_positive("r_eff_m", r_eff_m)?;
    Ok(mu * f * f64::from(n) * r)
}

/// Per-face friction force `mu * F` of a disc brake, in newtons.
///
/// This is the tangential drag at one pad face, before it is multiplied
/// by the face count and radius to give torque.
///
/// # Errors
/// [`BrakeError`] if `mu` is out of range or `clamp_force_n` is
/// non-finite or non-positive.
///
/// # Examples
/// ```
/// use valenx_brake::disc::friction_force;
/// let f = friction_force(0.4, 8_000.0).unwrap();
/// assert!((f - 3_200.0).abs() < 1e-9);
/// ```
pub fn friction_force(mu: f64, clamp_force_n: f64) -> Result<f64, BrakeError> {
    let mu = check_friction(mu)?;
    let f = check_positive("clamp_force_n", clamp_force_n)?;
    Ok(mu * f)
}

/// Clamp force `F` required to develop a target disc-brake torque.
///
/// Inverts `T = mu * F * n_pads * r_eff` for `F`:
/// `F = T / (mu * n_pads * r_eff)`.
///
/// # Errors
/// [`BrakeError`] if `mu` is out of range, `target_torque_nm` or
/// `r_eff_m` is non-finite or non-positive, or `n_pads == 0`.
///
/// # Examples
/// ```
/// use valenx_brake::disc::clamp_force_for_torque;
/// let f = clamp_force_for_torque(768.0, 0.4, 2, 0.12).unwrap();
/// assert!((f - 8_000.0).abs() < 1e-6);
/// ```
pub fn clamp_force_for_torque(
    target_torque_nm: f64,
    mu: f64,
    n_pads: u32,
    r_eff_m: f64,
) -> Result<f64, BrakeError> {
    let t = check_positive("target_torque_nm", target_torque_nm)?;
    let mu = check_friction(mu)?;
    let n = check_count("n_pads", n_pads)?;
    let r = check_positive("r_eff_m", r_eff_m)?;
    Ok(t / (mu * f64::from(n) * r))
}

/// Effective radius `r_eff` required to develop a target disc-brake torque.
///
/// Inverts `T = mu * F * n_pads * r_eff` for `r_eff`:
/// `r_eff = T / (mu * F * n_pads)`.
///
/// # Errors
/// [`BrakeError`] if `mu` is out of range, `target_torque_nm` or
/// `clamp_force_n` is non-finite or non-positive, or `n_pads == 0`.
///
/// # Examples
/// ```
/// use valenx_brake::disc::effective_radius_for_torque;
/// let r = effective_radius_for_torque(768.0, 0.4, 8_000.0, 2).unwrap();
/// assert!((r - 0.12).abs() < 1e-9);
/// ```
pub fn effective_radius_for_torque(
    target_torque_nm: f64,
    mu: f64,
    clamp_force_n: f64,
    n_pads: u32,
) -> Result<f64, BrakeError> {
    let t = check_positive("target_torque_nm", target_torque_nm)?;
    let mu = check_friction(mu)?;
    let f = check_positive("clamp_force_n", clamp_force_n)?;
    let n = check_count("n_pads", n_pads)?;
    Ok(t / (mu * f * f64::from(n)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::ErrorCategory;

    const EPS: f64 = 1e-9;

    #[test]
    fn matches_closed_form() {
        // T = mu*F*N*r = 0.4 * 8000 * 2 * 0.12 = 768.
        let t = disc_torque(0.4, 8_000.0, 2, 0.12).unwrap();
        assert!((t - 768.0).abs() < EPS, "got {t}");
    }

    #[test]
    fn single_pad_known_value() {
        // mu*F*N*r = 0.35 * 5000 * 1 * 0.15 = 262.5.
        let t = disc_torque(0.35, 5_000.0, 1, 0.15).unwrap();
        assert!((t - 262.5).abs() < EPS, "got {t}");
    }

    #[test]
    fn torque_scales_linearly_with_clamp_force() {
        // Doubling F must double T (linearity in clamp force).
        let base = disc_torque(0.4, 4_000.0, 2, 0.12).unwrap();
        let doubled = disc_torque(0.4, 8_000.0, 2, 0.12).unwrap();
        assert!(
            (doubled - 2.0 * base).abs() < EPS,
            "base {base} doubled {doubled}"
        );
    }

    #[test]
    fn torque_scales_linearly_with_radius() {
        // Tripling r_eff must triple T (linearity in radius).
        let base = disc_torque(0.4, 8_000.0, 2, 0.05).unwrap();
        let tripled = disc_torque(0.4, 8_000.0, 2, 0.15).unwrap();
        assert!(
            (tripled - 3.0 * base).abs() < EPS,
            "base {base} tripled {tripled}"
        );
    }

    #[test]
    fn torque_scales_with_pad_count() {
        // Two faces give exactly twice a single face, all else equal.
        let one = disc_torque(0.42, 6_000.0, 1, 0.1).unwrap();
        let two = disc_torque(0.42, 6_000.0, 2, 0.1).unwrap();
        assert!((two - 2.0 * one).abs() < EPS, "one {one} two {two}");
    }

    #[test]
    fn friction_force_is_mu_times_clamp() {
        let f = friction_force(0.4, 8_000.0).unwrap();
        assert!((f - 3_200.0).abs() < EPS, "got {f}");
    }

    #[test]
    fn clamp_force_inverts_torque() {
        // Round-trip: torque -> required clamp force -> torque.
        let t = disc_torque(0.4, 8_000.0, 2, 0.12).unwrap();
        let f = clamp_force_for_torque(t, 0.4, 2, 0.12).unwrap();
        assert!((f - 8_000.0).abs() < 1e-6, "got {f}");
        let back = disc_torque(0.4, f, 2, 0.12).unwrap();
        assert!((back - t).abs() < 1e-6, "got {back}");
    }

    #[test]
    fn radius_inverts_torque() {
        let t = disc_torque(0.4, 8_000.0, 2, 0.12).unwrap();
        let r = effective_radius_for_torque(t, 0.4, 8_000.0, 2).unwrap();
        assert!((r - 0.12).abs() < EPS, "got {r}");
    }

    #[test]
    fn rejects_bad_inputs() {
        assert_eq!(
            disc_torque(0.0, 8_000.0, 2, 0.12).unwrap_err().code(),
            "brake.friction_out_of_range"
        );
        assert_eq!(
            disc_torque(0.4, 0.0, 2, 0.12).unwrap_err().code(),
            "brake.non_positive"
        );
        assert_eq!(
            disc_torque(0.4, 8_000.0, 0, 0.12).unwrap_err().code(),
            "brake.zero_count"
        );
        assert_eq!(
            disc_torque(0.4, 8_000.0, 2, -0.1).unwrap_err().code(),
            "brake.non_positive"
        );
        assert_eq!(
            disc_torque(0.4, f64::NAN, 2, 0.12).unwrap_err().category(),
            ErrorCategory::Input
        );
    }
}
