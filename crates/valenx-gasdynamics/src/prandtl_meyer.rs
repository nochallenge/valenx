//! Prandtl–Meyer expansion: the supersonic expansion-fan relation.
//!
//! The Prandtl–Meyer function `ν(M)` is the angle through which a uniform sonic
//! (`M = 1`) flow must turn to expand isentropically to Mach `M ≥ 1`:
//!
//! ```text
//! ν(M) = sqrt((γ+1)/(γ-1)) · atan( sqrt( (γ-1)/(γ+1) · (M²-1) ) ) − atan( sqrt(M²-1) )
//! ```
//!
//! For a centred expansion that turns the flow through a deflection angle `θ`,
//! the downstream Mach follows from `ν(M₂) = ν(M₁) + θ`
//! ([`mach_after_expansion`]). As `M → ∞`, `ν` approaches the finite maximum
//! `ν_max = (π/2)·(sqrt((γ+1)/(γ-1)) − 1)` ([`nu_max`]) — for air (`γ = 1.4`),
//! `ν_max ≈ 130.45°`.
//!
//! Reference: Anderson, *Modern Compressible Flow*; NACA Report 1135. Same
//! perfect-gas scope and honest limitations as the rest of the crate.

use crate::error::{check_gamma, GasError, Result};

/// The Prandtl–Meyer function `ν(M)` in **radians** for Mach `mach` (`≥ 1`) and
/// specific-heat ratio `gamma`.
///
/// # Errors
///
/// [`GasError::BadGamma`] if `gamma` is not finite and `> 1`;
/// [`GasError::BadMach`] if `mach` is not finite or `< 1` (the function is only
/// defined for supersonic flow).
pub fn prandtl_meyer_angle(mach: f64, gamma: f64) -> Result<f64> {
    let gamma = check_gamma(gamma)?;
    if !mach.is_finite() || mach < 1.0 {
        return Err(GasError::bad_mach(
            mach,
            "prandtl_meyer",
            "must be finite and >= 1 (supersonic)",
        ));
    }
    let gp = gamma + 1.0;
    let gm = gamma - 1.0;
    let k = (gp / gm).sqrt();
    let m2_minus_1 = mach * mach - 1.0;
    Ok(k * ((gm / gp) * m2_minus_1).sqrt().atan() - m2_minus_1.sqrt().atan())
}

/// The asymptotic maximum turning angle `ν_max = (π/2)(sqrt((γ+1)/(γ-1)) − 1)`
/// (radians), approached as `M → ∞`.
///
/// # Errors
///
/// [`GasError::BadGamma`] if `gamma` is not finite and `> 1`.
pub fn nu_max(gamma: f64) -> Result<f64> {
    let gamma = check_gamma(gamma)?;
    let k = ((gamma + 1.0) / (gamma - 1.0)).sqrt();
    Ok(std::f64::consts::FRAC_PI_2 * (k - 1.0))
}

/// Invert the Prandtl–Meyer function: the Mach number whose expansion angle is
/// `nu` radians (`0 ≤ nu < ν_max`). Found by bisection, since `ν(M)` is strictly
/// increasing.
///
/// # Errors
///
/// [`GasError::BadGamma`] if `gamma` is invalid; [`GasError::BadMach`] if `nu`
/// is not finite, negative, or at/beyond [`nu_max`] (unreachable at finite `M`).
pub fn mach_from_prandtl_meyer(nu: f64, gamma: f64) -> Result<f64> {
    let gamma = check_gamma(gamma)?;
    let nu_limit = nu_max(gamma)?;
    if !nu.is_finite() || nu < 0.0 || nu >= nu_limit {
        return Err(GasError::bad_mach(
            nu,
            "prandtl_meyer_inverse",
            "nu must be finite and in [0, nu_max)",
        ));
    }
    // Bracket: ν(1) = 0; grow the upper Mach until ν(hi) ≥ nu, then bisect.
    let mut lo = 1.0_f64;
    let mut hi = 2.0_f64;
    while prandtl_meyer_angle(hi, gamma)? < nu {
        hi *= 2.0;
        if hi > 1e9 {
            break;
        }
    }
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if prandtl_meyer_angle(mid, gamma)? < nu {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(0.5 * (lo + hi))
}

/// The downstream Mach number after a centred expansion that turns supersonic
/// flow at Mach `m1` through a deflection `deflection` (radians, `≥ 0`), via
/// `ν(M₂) = ν(M₁) + deflection`.
///
/// # Errors
///
/// As [`prandtl_meyer_angle`] / [`mach_from_prandtl_meyer`]; also
/// [`GasError::BadMach`] if `deflection` is not finite or negative, or if the
/// requested total turning reaches [`nu_max`] (the flow would expand to vacuum).
pub fn mach_after_expansion(m1: f64, deflection: f64, gamma: f64) -> Result<f64> {
    let nu1 = prandtl_meyer_angle(m1, gamma)?;
    if !deflection.is_finite() || deflection < 0.0 {
        return Err(GasError::bad_mach(
            deflection,
            "expansion_deflection",
            "must be finite and >= 0",
        ));
    }
    mach_from_prandtl_meyer(nu1 + deflection, gamma)
}

#[cfg(test)]
mod tests {
    use super::*;

    // NACA Report 1135 reference values for air (γ = 1.4), ν in degrees.
    #[test]
    fn matches_naca_1135_table_for_air() {
        let g = 1.4;
        let cases = [
            (1.0, 0.0),
            (2.0, 26.380),
            (3.0, 49.757),
            (4.0, 65.785),
            (5.0, 76.920),
        ];
        for (mach, nu_deg) in cases {
            let got = prandtl_meyer_angle(mach, g).unwrap().to_degrees();
            assert!(
                (got - nu_deg).abs() < 0.01,
                "ν({mach}) = {got:.3}° vs NACA {nu_deg}°"
            );
        }
    }

    #[test]
    fn nu_max_for_air_is_about_130_45_degrees() {
        let got = nu_max(1.4).unwrap().to_degrees();
        assert!((got - 130.454).abs() < 0.01, "ν_max = {got:.3}°");
    }

    #[test]
    fn monotonic_increasing_in_mach() {
        let g = 1.4;
        let mut prev = -1.0;
        for i in 0..=40 {
            let m = 1.0 + f64::from(i) * 0.25;
            let nu = prandtl_meyer_angle(m, g).unwrap();
            assert!(nu > prev, "ν not increasing at M={m}");
            prev = nu;
        }
    }

    #[test]
    fn inverse_round_trips() {
        let g = 1.4;
        for &m in &[1.5, 2.0, 2.5, 4.0, 6.0] {
            let nu = prandtl_meyer_angle(m, g).unwrap();
            let recovered = mach_from_prandtl_meyer(nu, g).unwrap();
            assert!(
                (recovered - m).abs() < 1e-6,
                "round-trip M {m} -> {recovered}"
            );
        }
    }

    #[test]
    fn expansion_through_a_deflection() {
        // Air at M=2 (ν≈26.38°) expanded through 10° -> ν≈36.38° -> M≈2.385
        // (NACA 1135: ν=36.38° corresponds to M ≈ 2.38).
        let g = 1.4;
        let m2 = mach_after_expansion(2.0, 10.0_f64.to_radians(), g).unwrap();
        assert!(
            (2.35..2.42).contains(&m2),
            "M after 10° expansion = {m2:.4}"
        );
        // The total angle is exactly ν(M1) + θ.
        let total = prandtl_meyer_angle(m2, g).unwrap().to_degrees();
        let expected = prandtl_meyer_angle(2.0, g).unwrap().to_degrees() + 10.0;
        assert!((total - expected).abs() < 1e-6);
    }

    #[test]
    fn rejects_subsonic_and_bad_inputs() {
        assert!(prandtl_meyer_angle(0.8, 1.4).is_err());
        assert!(prandtl_meyer_angle(2.0, 1.0).is_err()); // gamma must be > 1
        assert!(mach_from_prandtl_meyer(-0.1, 1.4).is_err());
        // nu at/above the asymptote is unreachable.
        let nu_limit = nu_max(1.4).unwrap();
        assert!(mach_from_prandtl_meyer(nu_limit, 1.4).is_err());
    }
}
