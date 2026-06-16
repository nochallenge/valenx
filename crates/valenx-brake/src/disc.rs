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
//!
//! The effective radius itself can be derived from the annular pad
//! geometry under the two standard contact assumptions:
//! [`effective_radius_uniform_wear`] (a worn, bedded-in brake,
//! `(r_i + r_o)/2`) and [`effective_radius_uniform_pressure`] (a fresh
//! pad, `(2/3)(r_o^3 - r_i^3)/(r_o^2 - r_i^2)`).

use crate::error::{check_count, check_friction, check_non_negative, check_positive, BrakeError};

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

/// Effective friction radius of an annular pad under the **uniform-wear**
/// assumption, in metres.
///
/// As a brake beds in, the contact pressure redistributes so that the
/// product `p * r` is constant (the pad wears fastest where it slides
/// fastest). Integrating the friction torque under that assumption puts
/// the effective radius at the arithmetic mean of the pad's inner and
/// outer radii:
///
/// ```text
/// r_eff = (r_inner + r_outer) / 2
/// ```
///
/// This is the smaller of the two standard estimates and the one usually
/// used for a worn, in-service brake. Feed the result to [`disc_torque`]
/// as `r_eff_m`.
///
/// # Parameters
/// - `r_inner_m` — inner pad radius, in metres (>= 0).
/// - `r_outer_m` — outer pad radius, in metres (> 0, and > `r_inner_m`).
///
/// # Errors
/// [`BrakeError`] if either radius is non-finite, `r_inner_m` is
/// negative, `r_outer_m` is non-positive, or `r_outer_m <= r_inner_m`
/// (a degenerate annulus).
///
/// # Examples
/// ```
/// use valenx_brake::disc::effective_radius_uniform_wear;
/// // Pad from 100 mm to 140 mm -> mean radius 120 mm.
/// let r = effective_radius_uniform_wear(0.10, 0.14).unwrap();
/// assert!((r - 0.12).abs() < 1e-12);
/// ```
pub fn effective_radius_uniform_wear(r_inner_m: f64, r_outer_m: f64) -> Result<f64, BrakeError> {
    let r_i = check_non_negative("r_inner_m", r_inner_m)?;
    let r_o = check_positive("r_outer_m", r_outer_m)?;
    // A non-degenerate annulus needs the outer radius to exceed the inner.
    check_positive("r_outer_m - r_inner_m", r_o - r_i)?;
    Ok(0.5 * (r_i + r_o))
}

/// Effective friction radius of an annular pad under the
/// **uniform-pressure** assumption, in metres.
///
/// For a fresh pad pressing with uniform contact pressure, integrating
/// the friction torque `r * (p * 2*pi*r dr)` over the annulus and
/// dividing by the total clamp force gives
///
/// ```text
/// r_eff = (2/3) * (r_outer^3 - r_inner^3) / (r_outer^2 - r_inner^2)
/// ```
///
/// This is the larger of the two standard estimates and applies to a new
/// brake before it has worn in; it always exceeds the uniform-wear value
/// [`effective_radius_uniform_wear`] by exactly
/// `(r_outer - r_inner)^2 / (6*(r_outer + r_inner))`. Feed the result to
/// [`disc_torque`] as `r_eff_m`.
///
/// # Parameters
/// - `r_inner_m` — inner pad radius, in metres (>= 0).
/// - `r_outer_m` — outer pad radius, in metres (> 0, and > `r_inner_m`).
///
/// # Errors
/// [`BrakeError`] if either radius is non-finite, `r_inner_m` is
/// negative, `r_outer_m` is non-positive, or `r_outer_m <= r_inner_m`
/// (a degenerate annulus).
///
/// # Examples
/// ```
/// use valenx_brake::disc::effective_radius_uniform_pressure;
/// // Solid disc (r_inner = 0): r_eff = (2/3)*r_outer.
/// let r = effective_radius_uniform_pressure(0.0, 0.15).unwrap();
/// assert!((r - 0.1).abs() < 1e-12);
/// ```
pub fn effective_radius_uniform_pressure(
    r_inner_m: f64,
    r_outer_m: f64,
) -> Result<f64, BrakeError> {
    let r_i = check_non_negative("r_inner_m", r_inner_m)?;
    let r_o = check_positive("r_outer_m", r_outer_m)?;
    // A non-degenerate annulus needs the outer radius to exceed the inner.
    check_positive("r_outer_m - r_inner_m", r_o - r_i)?;
    let num = r_o.powi(3) - r_i.powi(3);
    let den = r_o.powi(2) - r_i.powi(2);
    Ok((2.0 / 3.0) * num / den)
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

    #[test]
    fn uniform_wear_is_mean_radius() {
        // (0.10 + 0.14) / 2 = 0.12 m.
        let r = effective_radius_uniform_wear(0.10, 0.14).unwrap();
        assert!((r - 0.12).abs() < EPS, "got {r}");
    }

    #[test]
    fn uniform_pressure_closed_form() {
        // (2/3)·(0.14³ − 0.10³)/(0.14² − 0.10²)
        //   = (2/3)·0.001744/0.0096 = 0.121111… m.
        let r = effective_radius_uniform_pressure(0.10, 0.14).unwrap();
        let num = 0.14_f64.powi(3) - 0.10_f64.powi(3);
        let den = 0.14_f64.powi(2) - 0.10_f64.powi(2);
        let expected = (2.0 / 3.0) * num / den;
        assert!((r - expected).abs() < EPS, "got {r} vs {expected}");
        assert!((r - 0.121_111_111).abs() < 1e-6, "got {r}");
    }

    #[test]
    fn pressure_exceeds_wear_by_exact_gap() {
        // Derivable identity: r_p − r_w = (r_o − r_i)² / (6·(r_o + r_i)),
        // strictly positive for any proper annulus. A fresh (uniform
        // pressure) brake has a larger effective radius than a worn one.
        let (r_i, r_o) = (0.10, 0.14);
        let r_w = effective_radius_uniform_wear(r_i, r_o).unwrap();
        let r_p = effective_radius_uniform_pressure(r_i, r_o).unwrap();
        assert!(r_p > r_w, "pressure {r_p} should exceed wear {r_w}");
        let gap = (r_o - r_i).powi(2) / (6.0 * (r_o + r_i));
        assert!((r_p - r_w - gap).abs() < EPS, "gap mismatch");
    }

    #[test]
    fn solid_disc_limits() {
        // r_inner = 0: wear → r_o/2, pressure → (2/3)·r_o.
        let r_o = 0.15;
        let r_w = effective_radius_uniform_wear(0.0, r_o).unwrap();
        let r_p = effective_radius_uniform_pressure(0.0, r_o).unwrap();
        assert!((r_w - r_o / 2.0).abs() < EPS, "wear {r_w}");
        assert!((r_p - (2.0 / 3.0) * r_o).abs() < EPS, "pressure {r_p}");
    }

    #[test]
    fn thin_annulus_both_approach_mean() {
        // As the band narrows (r_i → r_o) both estimates collapse onto the
        // common radius (the (r_o−r_i)² gap vanishes).
        let (r_i, r_o) = (0.1399, 0.1400);
        let r_w = effective_radius_uniform_wear(r_i, r_o).unwrap();
        let r_p = effective_radius_uniform_pressure(r_i, r_o).unwrap();
        assert!((r_p - r_w).abs() < 1e-7, "w {r_w} p {r_p}");
    }

    #[test]
    fn effective_radius_feeds_disc_torque() {
        // The computed r_eff is a valid input to disc_torque.
        let r = effective_radius_uniform_wear(0.10, 0.14).unwrap();
        let t = disc_torque(0.4, 8_000.0, 2, r).unwrap();
        assert!((t - 0.4 * 8_000.0 * 2.0 * 0.12).abs() < 1e-6, "got {t}");
    }

    #[test]
    fn effective_radius_rejects_bad_annulus() {
        // Outer must strictly exceed inner.
        assert_eq!(
            effective_radius_uniform_wear(0.14, 0.10)
                .unwrap_err()
                .code(),
            "brake.non_positive"
        );
        assert_eq!(
            effective_radius_uniform_pressure(0.12, 0.12)
                .unwrap_err()
                .code(),
            "brake.non_positive"
        );
        // Negative inner radius.
        assert_eq!(
            effective_radius_uniform_wear(-0.01, 0.14)
                .unwrap_err()
                .code(),
            "brake.negative"
        );
        // Non-positive outer radius.
        assert_eq!(
            effective_radius_uniform_pressure(0.0, 0.0)
                .unwrap_err()
                .code(),
            "brake.non_positive"
        );
        // Non-finite input.
        assert_eq!(
            effective_radius_uniform_wear(f64::NAN, 0.14)
                .unwrap_err()
                .code(),
            "brake.not_finite"
        );
    }
}
