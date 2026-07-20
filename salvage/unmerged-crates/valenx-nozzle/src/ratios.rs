//! Stagnation-to-static property ratios for steady isentropic flow of a
//! calorically perfect gas.
//!
//! For a gas with constant ratio of specific heats `gamma`, every quantity
//! below is a pure function of the local Mach number `M`. The temperature
//! ratio comes straight from conservation of stagnation enthalpy; the
//! pressure and density ratios follow by the isentropic relation
//! `p / rho^gamma = const`.
//!
//! Governing equations (`g` denotes `gamma`):
//!
//! ```text
//! T0/T   = 1 + (g-1)/2 * M^2
//! p0/p   = (T0/T)^( g / (g-1) )
//! rho0/rho = (T0/T)^( 1 / (g-1) )
//! ```
//!
//! All three reduce to `1` at `M = 0` (the stagnation point itself) and
//! increase monotonically with `M`.

use crate::error::{check_gamma, check_mach, NozzleError, Result};

/// Stagnation-to-static **temperature** ratio `T0/T`.
///
/// `T0/T = 1 + (gamma-1)/2 * M^2`.
///
/// # Errors
///
/// Returns a [`NozzleError`] if `gamma` is outside `(1, 5/3]` or `mach` is
/// negative / non-finite (see [`check_gamma`] / [`check_mach`]).
///
/// # Example
///
/// ```
/// // Air (gamma = 1.4) at M = 1: T0/T = 1 + 0.2 = 1.2.
/// let r = valenx_nozzle::temperature_ratio(1.0, 1.4).unwrap();
/// assert!((r - 1.2).abs() < 1e-12);
/// ```
pub fn temperature_ratio(mach: f64, gamma: f64) -> Result<f64> {
    check_gamma(gamma)?;
    check_mach(mach)?;
    Ok(1.0 + 0.5 * (gamma - 1.0) * mach * mach)
}

/// Stagnation-to-static **pressure** ratio `p0/p`.
///
/// `p0/p = (T0/T)^( gamma / (gamma-1) )`.
///
/// # Errors
///
/// Returns a [`NozzleError`] for an out-of-range `gamma` or a negative /
/// non-finite `mach`.
///
/// # Example
///
/// ```
/// // Air at M = 1: p0/p = 1.2^3.5 ~= 1.892929.
/// let r = valenx_nozzle::pressure_ratio(1.0, 1.4).unwrap();
/// assert!((r - 1.892929_158_98).abs() < 1e-8);
/// ```
pub fn pressure_ratio(mach: f64, gamma: f64) -> Result<f64> {
    let t = temperature_ratio(mach, gamma)?;
    Ok(t.powf(gamma / (gamma - 1.0)))
}

/// Stagnation-to-static **density** ratio `rho0/rho`.
///
/// `rho0/rho = (T0/T)^( 1 / (gamma-1) )`.
///
/// # Errors
///
/// Returns a [`NozzleError`] for an out-of-range `gamma` or a negative /
/// non-finite `mach`.
///
/// # Example
///
/// ```
/// // Air at M = 1: rho0/rho = 1.2^2.5 ~= 1.577439.
/// let r = valenx_nozzle::density_ratio(1.0, 1.4).unwrap();
/// assert!((r - 1.577_440).abs() < 1e-5);
/// ```
pub fn density_ratio(mach: f64, gamma: f64) -> Result<f64> {
    let t = temperature_ratio(mach, gamma)?;
    Ok(t.powf(1.0 / (gamma - 1.0)))
}

/// The **critical pressure ratio** `p*/p0` — the static-to-stagnation
/// pressure at the sonic throat (`M = 1`).
///
/// `p*/p0 = (2 / (gamma+1))^( gamma / (gamma-1) )`.
///
/// This is the reciprocal of [`pressure_ratio`] evaluated at `M = 1`. For
/// air (`gamma = 1.4`) it is the familiar `0.5283`: when the back-pressure
/// ratio falls to this value, the throat chokes.
///
/// # Errors
///
/// Returns [`NozzleError`] if `gamma` is outside `(1, 5/3]`.
///
/// # Example
///
/// ```
/// let r = valenx_nozzle::critical_pressure_ratio(1.4).unwrap();
/// assert!((r - 0.528_281_8).abs() < 1e-6);
/// ```
pub fn critical_pressure_ratio(gamma: f64) -> Result<f64> {
    check_gamma(gamma)?;
    Ok((2.0 / (gamma + 1.0)).powf(gamma / (gamma - 1.0)))
}

/// Recover the Mach number from a static-to-stagnation **pressure** ratio
/// `p/p0`, the inverse of [`pressure_ratio`].
///
/// Inverting `p0/p = (1 + (g-1)/2 M^2)^(g/(g-1))` gives the closed form
///
/// ```text
/// M = sqrt( 2/(g-1) * ( (p0/p)^((g-1)/g) - 1 ) ).
/// ```
///
/// The argument here is `p/p0` (the reciprocal), so `1.0` maps to `M = 0`
/// and smaller values to faster flow.
///
/// # Errors
///
/// Returns [`NozzleError::PressureRatioOutOfRange`] if `p_over_p0` is not in
/// `(0, 1]`, or a validation error for a bad `gamma`.
///
/// # Example
///
/// ```
/// // Round-trip: p/p0 at M = 2 then back.
/// let p_over_p0 = 1.0 / valenx_nozzle::pressure_ratio(2.0, 1.4).unwrap();
/// let m = valenx_nozzle::mach_from_pressure_ratio(p_over_p0, 1.4).unwrap();
/// assert!((m - 2.0).abs() < 1e-9);
/// ```
pub fn mach_from_pressure_ratio(p_over_p0: f64, gamma: f64) -> Result<f64> {
    check_gamma(gamma)?;
    if !p_over_p0.is_finite() {
        return Err(NozzleError::NotFinite {
            name: "p_over_p0",
            value: p_over_p0,
        });
    }
    if p_over_p0 <= 0.0 || p_over_p0 > 1.0 {
        return Err(NozzleError::PressureRatioOutOfRange { value: p_over_p0 });
    }
    let p0_over_p = 1.0 / p_over_p0;
    let inner = p0_over_p.powf((gamma - 1.0) / gamma) - 1.0;
    // `inner` is >= 0 across the valid domain; clamp away a tiny negative
    // round-off at p/p0 == 1 so the sqrt stays real.
    let inner = inner.max(0.0);
    Ok((2.0 / (gamma - 1.0) * inner).sqrt())
}
