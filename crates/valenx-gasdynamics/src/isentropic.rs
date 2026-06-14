//! Isentropic (adiabatic, reversible) compressible-flow relations for a
//! calorically-perfect ideal gas.
//!
//! For a steady one-dimensional flow that is brought to rest
//! isentropically, the local static state `(T, p, rho)` and the
//! stagnation state `(T0, p0, rho0)` are linked to the local Mach number
//! `M` through (Anderson, *Modern Compressible Flow*, ch. 3; NACA Report
//! 1135):
//!
//! ```text
//! T0 / T    = 1 + (gamma - 1)/2 * M^2
//! p0 / p    = (T0 / T) ^ ( gamma / (gamma - 1) )
//! rho0/rho  = (T0 / T) ^ (   1   / (gamma - 1) )
//! ```
//!
//! and the area a stream-tube must have, relative to the sonic
//! (`M = 1`) throat area `A*`, follows the **area-Mach relation**
//!
//! ```text
//! A / A* = (1 / M) * [ (2 / (gamma + 1)) * (1 + (gamma - 1)/2 * M^2) ]
//!                      ^ ( (gamma + 1) / (2 (gamma - 1)) )
//! ```
//!
//! All four are exact algebraic functions of `M` and `gamma`; no
//! iteration is involved.
//!
//! ## Honest scope
//!
//! These are the textbook closed-form relations for a *single* perfect
//! gas with constant specific heats. They ignore real-gas effects
//! (variable `cp`, dissociation, vibrational excitation), viscosity, and
//! heat addition. They are intended for study and first-order
//! engineering estimates, not certified design.

use serde::{Deserialize, Serialize};

use crate::error::{check_gamma, check_mach_nonneg, check_mach_pos, Result};

/// The three isentropic stagnation-to-static ratios at a given Mach
/// number.
///
/// Each field is the *stagnation over static* ratio (so every value is
/// `>= 1`, equal to one only at `M = 0`). Construct it with
/// [`stagnation_ratios`].
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StagnationRatios {
    /// Temperature ratio `T0 / T = 1 + (gamma - 1)/2 * M^2`.
    pub t0_over_t: f64,
    /// Pressure ratio `p0 / p = (T0/T)^(gamma/(gamma-1))`.
    pub p0_over_p: f64,
    /// Density ratio `rho0 / rho = (T0/T)^(1/(gamma-1))`.
    pub rho0_over_rho: f64,
}

/// Isentropic stagnation temperature ratio
/// `T0 / T = 1 + (gamma - 1)/2 * M^2`.
///
/// Equals one at `M = 0` and increases monotonically with `M`. Accepts
/// any finite `M >= 0`.
///
/// # Errors
///
/// Returns [`GasError::BadMach`](crate::error::GasError::BadMach) if `M`
/// is non-finite or negative, and
/// [`GasError::BadGamma`](crate::error::GasError::BadGamma) if `gamma`
/// is non-finite or `<= 1`.
pub fn temperature_ratio(mach: f64, gamma: f64) -> Result<f64> {
    let g = check_gamma(gamma)?;
    let m = check_mach_nonneg(mach, "temperature_ratio")?;
    Ok(1.0 + (g - 1.0) / 2.0 * m * m)
}

/// Isentropic stagnation pressure ratio
/// `p0 / p = (1 + (gamma - 1)/2 * M^2)^(gamma/(gamma-1))`.
///
/// Equals one at `M = 0` and increases monotonically with `M`.
///
/// # Errors
///
/// As [`temperature_ratio`].
pub fn pressure_ratio(mach: f64, gamma: f64) -> Result<f64> {
    let g = check_gamma(gamma)?;
    let t = temperature_ratio(mach, g)?;
    Ok(t.powf(g / (g - 1.0)))
}

/// Isentropic stagnation density ratio
/// `rho0 / rho = (1 + (gamma - 1)/2 * M^2)^(1/(gamma-1))`.
///
/// Equals one at `M = 0` and increases monotonically with `M`.
///
/// # Errors
///
/// As [`temperature_ratio`].
pub fn density_ratio(mach: f64, gamma: f64) -> Result<f64> {
    let g = check_gamma(gamma)?;
    let t = temperature_ratio(mach, g)?;
    Ok(t.powf(1.0 / (g - 1.0)))
}

/// All three isentropic stagnation ratios at once.
///
/// Cheaper and more convenient than three separate calls when a caller
/// needs the full set; the temperature ratio is computed once and reused.
///
/// # Errors
///
/// As [`temperature_ratio`].
pub fn stagnation_ratios(mach: f64, gamma: f64) -> Result<StagnationRatios> {
    let g = check_gamma(gamma)?;
    let t = temperature_ratio(mach, g)?;
    Ok(StagnationRatios {
        t0_over_t: t,
        p0_over_p: t.powf(g / (g - 1.0)),
        rho0_over_rho: t.powf(1.0 / (g - 1.0)),
    })
}

/// Isentropic area-Mach ratio `A / A*`, the local stream-tube area
/// relative to the sonic-throat area.
///
/// ```text
/// A / A* = (1 / M) * [ (2 / (gamma + 1)) * (1 + (gamma - 1)/2 * M^2) ]
///                      ^ ( (gamma + 1) / (2 (gamma - 1)) )
/// ```
///
/// Equal to one exactly at `M = 1`; greater than one for every other
/// Mach number (the ratio is double-valued in `M` — one subsonic and one
/// supersonic branch share each `A/A* > 1`). The flow must be moving, so
/// `M = 0` is rejected (the expression divides by `M`).
///
/// # Errors
///
/// Returns [`GasError::BadMach`](crate::error::GasError::BadMach) if `M`
/// is non-finite or `<= 0`, and
/// [`GasError::BadGamma`](crate::error::GasError::BadGamma) if `gamma`
/// is non-finite or `<= 1`.
pub fn area_mach_ratio(mach: f64, gamma: f64) -> Result<f64> {
    let g = check_gamma(gamma)?;
    let m = check_mach_pos(mach, "area_mach_ratio")?;
    let t = 1.0 + (g - 1.0) / 2.0 * m * m;
    let base = 2.0 / (g + 1.0) * t;
    let exp = (g + 1.0) / (2.0 * (g - 1.0));
    Ok(base.powf(exp) / m)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute tolerance for analytic comparisons. The relations are
    /// closed-form, so agreement is limited only by f64 round-off; this
    /// bound is comfortably tighter than the published-table precision.
    const EPS: f64 = 1e-9;

    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn stagnation_ratios_are_unity_at_zero_mach() {
        let r = stagnation_ratios(0.0, 1.4).unwrap();
        assert!(close(r.t0_over_t, 1.0, EPS), "T0/T = {}", r.t0_over_t);
        assert!(close(r.p0_over_p, 1.0, EPS), "p0/p = {}", r.p0_over_p);
        assert!(
            close(r.rho0_over_rho, 1.0, EPS),
            "rho0/rho = {}",
            r.rho0_over_rho
        );
    }

    #[test]
    fn stagnation_ratios_rise_with_mach() {
        // Strictly increasing in M for every ratio.
        let lo = stagnation_ratios(0.5, 1.4).unwrap();
        let hi = stagnation_ratios(2.5, 1.4).unwrap();
        assert!(hi.t0_over_t > lo.t0_over_t);
        assert!(hi.p0_over_p > lo.p0_over_p);
        assert!(hi.rho0_over_rho > lo.rho0_over_rho);
        // And every ratio exceeds one once the flow is moving.
        assert!(lo.t0_over_t > 1.0 && lo.p0_over_p > 1.0 && lo.rho0_over_rho > 1.0);
    }

    #[test]
    fn stagnation_ratios_match_known_values_m2_gamma14() {
        // Ground truth (NACA Report 1135 / Anderson App. A, gamma = 1.4):
        //   M = 2.0 -> T0/T = 1 + 0.2*4 = 1.8 (exact)
        //              p0/p = 1.8^3.5 = 7.82445...
        //              rho0/rho = 1.8^2.5 = 4.34692...
        // The stagnation-temperature ratio is exact; the pressure and
        // density ratios are `1.8` raised to `gamma/(gamma-1) = 3.5` and
        // `1/(gamma-1) = 2.5`. We compare against those closed forms (to
        // avoid pinning a platform-specific last-ULP decimal) and also
        // tabulate the rounded table values for the human reader.
        let r = stagnation_ratios(2.0, 1.4).unwrap();
        assert!(close(r.t0_over_t, 1.8, EPS), "T0/T = {}", r.t0_over_t);
        assert!(
            close(r.p0_over_p, 1.8_f64.powf(3.5), EPS),
            "p0/p = {}",
            r.p0_over_p
        );
        assert!(
            close(r.rho0_over_rho, 1.8_f64.powf(2.5), EPS),
            "rho0/rho = {}",
            r.rho0_over_rho
        );
        // Independent four-significant-figure cross-check against the
        // printed NACA-1135 table entries.
        assert!(close(r.p0_over_p, 7.8244, 1e-3), "p0/p = {}", r.p0_over_p);
        assert!(
            close(r.rho0_over_rho, 4.3469, 1e-3),
            "rho0/rho = {}",
            r.rho0_over_rho
        );
        // The separate accessors agree with the bundled struct.
        assert!(close(
            temperature_ratio(2.0, 1.4).unwrap(),
            r.t0_over_t,
            EPS
        ));
        assert!(close(pressure_ratio(2.0, 1.4).unwrap(), r.p0_over_p, EPS));
        assert!(close(
            density_ratio(2.0, 1.4).unwrap(),
            r.rho0_over_rho,
            EPS
        ));
    }

    #[test]
    fn pressure_ratio_equals_temp_ratio_to_the_gamma_power() {
        // p0/p = (T0/T)^(gamma/(gamma-1)) and rho0/rho = (T0/T)^(1/(gamma-1))
        // must hold identically at an arbitrary off-table state.
        let g = 1.3;
        let m = 1.37;
        let t = temperature_ratio(m, g).unwrap();
        let p = pressure_ratio(m, g).unwrap();
        let rho = density_ratio(m, g).unwrap();
        assert!(close(p, t.powf(g / (g - 1.0)), EPS), "p={p}");
        assert!(close(rho, t.powf(1.0 / (g - 1.0)), EPS), "rho={rho}");
        // Consistency of the equation of state: p0/p = (T0/T)*(rho0/rho).
        assert!(close(p, t * rho, EPS), "p={p}, t*rho={}", t * rho);
    }

    #[test]
    fn area_ratio_is_one_at_mach_one() {
        for &g in &[1.2_f64, 1.33, 1.4, 5.0 / 3.0] {
            let ar = area_mach_ratio(1.0, g).unwrap();
            assert!(close(ar, 1.0, EPS), "gamma={g}: A/A* = {ar}");
        }
    }

    #[test]
    fn area_ratio_matches_known_value_m2_gamma14() {
        // Ground truth: gamma = 1.4, M = 2.0 -> A/A* = 1.6875 exactly
        //   = (1/2) * (1.5)^3  (see Anderson App. A; matches the value
        //   baked into valenx-astro's nozzle solver test).
        let ar = area_mach_ratio(2.0, 1.4).unwrap();
        assert!(close(ar, 1.6875, EPS), "A/A* = {ar}");
    }

    #[test]
    fn area_ratio_exceeds_one_away_from_sonic_and_is_double_valued() {
        // Both a subsonic and a supersonic Mach number share each
        // A/A* > 1; pick the classic gamma = 1.4 pair M = 0.5 and the
        // matching supersonic root, and confirm both branches are > 1.
        let sub = area_mach_ratio(0.5, 1.4).unwrap();
        let sup = area_mach_ratio(2.0, 1.4).unwrap();
        assert!(sub > 1.0, "subsonic A/A* = {sub}");
        assert!(sup > 1.0, "supersonic A/A* = {sup}");
        // M = 0.5, gamma = 1.4 -> A/A* = 1.33984375 (NACA 1135 table),
        // which is exactly representable in binary f64.
        assert!(close(sub, 1.339_843_75, EPS), "A/A* = {sub}");
    }

    #[test]
    fn domain_errors() {
        // gamma must be > 1.
        assert!(temperature_ratio(1.0, 1.0).is_err());
        assert!(area_mach_ratio(1.0, 0.9).is_err());
        // Mach must be finite and non-negative for stagnation ratios.
        assert!(temperature_ratio(-0.1, 1.4).is_err());
        assert!(temperature_ratio(f64::NAN, 1.4).is_err());
        // Area-Mach divides by M, so M = 0 is rejected.
        assert!(area_mach_ratio(0.0, 1.4).is_err());
    }
}
