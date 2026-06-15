//! Hydrostatic pressure and the gauge ↔ absolute relationship.
//!
//! In a fluid at rest of constant density `rho` under gravity `g`, the
//! pressure increases linearly with depth `h` below the free surface:
//!
//! ```text
//! P_gauge(h) = rho * g * h
//! P_abs(h)   = P_surface + rho * g * h
//! ```
//!
//! where `P_surface` is the absolute pressure acting on the free
//! surface (one atmosphere for an open tank). *Gauge* pressure is
//! measured relative to the local ambient (surface) pressure; *absolute*
//! pressure is measured relative to a perfect vacuum. The two differ by
//! exactly the surface pressure:
//!
//! ```text
//! P_abs = P_gauge + P_surface
//! ```
//!
//! This module is the foundational building block — buoyancy,
//! submerged-plate forces and manometers all reduce to the linear
//! `rho * g * h` law applied at the right depths.

use crate::error::{require_finite, require_non_negative, require_positive, Result};
use crate::fluid::{Fluid, STANDARD_ATMOSPHERE_PA, STANDARD_GRAVITY};

/// Gauge pressure at depth `depth_m` below the free surface of `fluid`,
/// in pascals: `P = rho * g * h`.
///
/// The result is measured relative to the pressure at the free surface,
/// so it is zero at the surface and grows linearly with depth.
///
/// # Errors
///
/// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
/// if `gravity` is not strictly positive or `depth_m` is negative /
/// non-finite.
pub fn gauge_pressure(fluid: &Fluid, gravity: f64, depth_m: f64) -> Result<f64> {
    let gravity = require_positive("gravity", gravity)?;
    let depth_m = require_non_negative("depth_m", depth_m)?;
    Ok(fluid.density() * gravity * depth_m)
}

/// Absolute pressure at depth `depth_m` below a free surface that is
/// itself at absolute pressure `surface_pressure_pa`, in pascals:
/// `P_abs = P_surface + rho * g * h`.
///
/// For an open tank pass [`STANDARD_ATMOSPHERE_PA`] as the surface
/// pressure; see [`absolute_pressure_open`] for that common case.
///
/// # Errors
///
/// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
/// if `surface_pressure_pa` is negative / non-finite, `gravity` is not
/// strictly positive, or `depth_m` is negative / non-finite.
pub fn absolute_pressure(
    fluid: &Fluid,
    gravity: f64,
    depth_m: f64,
    surface_pressure_pa: f64,
) -> Result<f64> {
    let surface = require_non_negative("surface_pressure_pa", surface_pressure_pa)?;
    let gauge = gauge_pressure(fluid, gravity, depth_m)?;
    Ok(surface + gauge)
}

/// Absolute pressure at depth `depth_m` below an open free surface
/// exposed to one standard atmosphere, in pascals:
/// `P_abs = 1 atm + rho * g * h`.
///
/// # Errors
///
/// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
/// if `gravity` is not strictly positive or `depth_m` is negative /
/// non-finite.
pub fn absolute_pressure_open(fluid: &Fluid, gravity: f64, depth_m: f64) -> Result<f64> {
    absolute_pressure(fluid, gravity, depth_m, STANDARD_ATMOSPHERE_PA)
}

/// Convert an absolute pressure to a gauge pressure given the ambient
/// (reference) pressure: `P_gauge = P_abs - P_ambient`.
///
/// The result may be negative — a partial vacuum has a gauge pressure
/// below zero — so this returns any finite value.
///
/// # Errors
///
/// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
/// if either argument is non-finite or negative (an absolute pressure
/// cannot be below a perfect vacuum).
pub fn gauge_from_absolute(absolute_pa: f64, ambient_pa: f64) -> Result<f64> {
    let absolute = require_non_negative("absolute_pa", absolute_pa)?;
    let ambient = require_non_negative("ambient_pa", ambient_pa)?;
    Ok(absolute - ambient)
}

/// Convert a gauge pressure to an absolute pressure given the ambient
/// (reference) pressure: `P_abs = P_gauge + P_ambient`.
///
/// # Errors
///
/// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
/// if `gauge_pa` is non-finite, `ambient_pa` is negative / non-finite,
/// or the implied absolute pressure would fall below a perfect vacuum
/// (`P_gauge + P_ambient < 0`).
pub fn absolute_from_gauge(gauge_pa: f64, ambient_pa: f64) -> Result<f64> {
    let gauge = require_finite("gauge_pa", gauge_pa)?;
    let ambient = require_non_negative("ambient_pa", ambient_pa)?;
    let absolute = gauge + ambient;
    if absolute < 0.0 {
        return Err(crate::error::FluidStaticsError::invalid(
            "gauge_pa",
            format!(
                "implied absolute pressure {absolute} Pa is below vacuum \
                 (gauge {gauge} + ambient {ambient})"
            ),
        ));
    }
    Ok(absolute)
}

/// The vertical depth, in metres, at which `fluid` under gravity
/// `gravity` produces gauge pressure `target_pa`: `h = P / (rho * g)`.
///
/// This is the inverse of [`gauge_pressure`] and the basis of reading a
/// pressure as a "head" of fluid (e.g. millimetres of mercury).
///
/// # Errors
///
/// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
/// if `gravity` is not strictly positive or `target_pa` is negative /
/// non-finite.
pub fn depth_for_gauge_pressure(fluid: &Fluid, gravity: f64, target_pa: f64) -> Result<f64> {
    let gravity = require_positive("gravity", gravity)?;
    let target_pa = require_non_negative("target_pa", target_pa)?;
    Ok(target_pa / (fluid.density() * gravity))
}

/// Pressure head, in metres of the given `fluid`, equivalent to a gauge
/// pressure `gauge_pa`: `h = P / (rho * g)`.
///
/// A convenience alias for [`depth_for_gauge_pressure`] phrased in the
/// language of fluid head.
///
/// # Errors
///
/// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
/// if `gravity` is not strictly positive or `gauge_pa` is negative /
/// non-finite.
pub fn pressure_head(fluid: &Fluid, gravity: f64, gauge_pa: f64) -> Result<f64> {
    depth_for_gauge_pressure(fluid, gravity, gauge_pa)
}

/// Gauge pressure at depth `depth_m` of `fluid` under
/// [`STANDARD_GRAVITY`], in pascals — a convenience wrapper around
/// [`gauge_pressure`].
///
/// # Errors
///
/// Returns [`Invalid`](crate::FluidStaticsError::Invalid)
/// if `depth_m` is negative / non-finite.
pub fn gauge_pressure_standard(fluid: &Fluid, depth_m: f64) -> Result<f64> {
    gauge_pressure(fluid, STANDARD_GRAVITY, depth_m)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    #[test]
    fn surface_gauge_pressure_is_zero() {
        let p = gauge_pressure(&Fluid::water(), STANDARD_GRAVITY, 0.0).unwrap();
        assert!(p.abs() < EPS, "got {p}");
    }

    #[test]
    fn ten_metres_of_water_is_about_one_atmosphere() {
        // rho*g*h = 1000 * 9.80665 * 10 = 98066.5 Pa, ~0.97 atm — the
        // textbook "10 m of water ≈ 1 atm" rule of thumb.
        let p = gauge_pressure(&Fluid::water(), STANDARD_GRAVITY, 10.0).unwrap();
        assert!((p - 98_066.5).abs() < 1e-3, "got {p}");
        assert!((p / STANDARD_ATMOSPHERE_PA - 0.9678).abs() < 1e-3);
    }

    #[test]
    fn pressure_is_linear_in_depth() {
        // P(2h) = 2 * P(h) and P(h1+h2) = P(h1)+P(h2): linearity.
        let f = Fluid::seawater();
        let p1 = gauge_pressure(&f, STANDARD_GRAVITY, 3.0).unwrap();
        let p2 = gauge_pressure(&f, STANDARD_GRAVITY, 6.0).unwrap();
        assert!((p2 - 2.0 * p1).abs() < EPS, "p1={p1} p2={p2}");

        let pa = gauge_pressure(&f, STANDARD_GRAVITY, 4.0).unwrap();
        let pb = gauge_pressure(&f, STANDARD_GRAVITY, 7.0).unwrap();
        let psum = gauge_pressure(&f, STANDARD_GRAVITY, 11.0).unwrap();
        assert!((psum - (pa + pb)).abs() < EPS, "pa={pa} pb={pb} sum={psum}");
    }

    #[test]
    fn gauge_pressure_slope_equals_specific_weight() {
        // dP/dh = rho*g exactly: a finite difference over any interval
        // recovers the specific weight.
        let f = Fluid::water();
        let gamma = f.specific_weight(STANDARD_GRAVITY).unwrap();
        let p_lo = gauge_pressure(&f, STANDARD_GRAVITY, 2.0).unwrap();
        let p_hi = gauge_pressure(&f, STANDARD_GRAVITY, 5.0).unwrap();
        let slope = (p_hi - p_lo) / (5.0 - 2.0);
        assert!((slope - gamma).abs() < EPS, "slope={slope} gamma={gamma}");
    }

    #[test]
    fn absolute_is_gauge_plus_surface() {
        let f = Fluid::water();
        let gauge = gauge_pressure(&f, STANDARD_GRAVITY, 5.0).unwrap();
        let abs_open = absolute_pressure_open(&f, STANDARD_GRAVITY, 5.0).unwrap();
        assert!(
            (abs_open - (gauge + STANDARD_ATMOSPHERE_PA)).abs() < EPS,
            "gauge={gauge} abs={abs_open}"
        );
    }

    #[test]
    fn gauge_absolute_round_trip() {
        let abs_pa = 150_000.0;
        let ambient = STANDARD_ATMOSPHERE_PA;
        let gauge = gauge_from_absolute(abs_pa, ambient).unwrap();
        let back = absolute_from_gauge(gauge, ambient).unwrap();
        assert!((back - abs_pa).abs() < EPS, "back={back}");
        // The gauge value is the expected difference.
        assert!((gauge - (abs_pa - ambient)).abs() < EPS, "gauge={gauge}");
    }

    #[test]
    fn partial_vacuum_has_negative_gauge() {
        // An absolute pressure below ambient reads negative on a gauge.
        let gauge = gauge_from_absolute(80_000.0, STANDARD_ATMOSPHERE_PA).unwrap();
        assert!(gauge < 0.0, "got {gauge}");
    }

    #[test]
    fn absolute_from_gauge_rejects_below_vacuum() {
        // -2 atm gauge at 1 atm ambient would be -1 atm absolute: impossible.
        let err = absolute_from_gauge(-2.0 * STANDARD_ATMOSPHERE_PA, STANDARD_ATMOSPHERE_PA);
        assert!(err.is_err());
    }

    #[test]
    fn depth_for_pressure_inverts_pressure_for_depth() {
        let f = Fluid::seawater();
        let h = 7.3;
        let p = gauge_pressure(&f, STANDARD_GRAVITY, h).unwrap();
        let h_back = depth_for_gauge_pressure(&f, STANDARD_GRAVITY, p).unwrap();
        assert!((h_back - h).abs() < 1e-9, "h_back={h_back}");
    }

    #[test]
    fn one_atmosphere_of_mercury_head_is_about_760_mm() {
        // The classic barometer: 1 atm ≈ 760 mmHg.
        let head =
            pressure_head(&Fluid::mercury(), STANDARD_GRAVITY, STANDARD_ATMOSPHERE_PA).unwrap();
        assert!((head - 0.7637).abs() < 1e-3, "got {head} m");
    }

    #[test]
    fn rejects_bad_arguments() {
        assert!(gauge_pressure(&Fluid::water(), 0.0, 1.0).is_err());
        assert!(gauge_pressure(&Fluid::water(), STANDARD_GRAVITY, -1.0).is_err());
        assert!(absolute_pressure(&Fluid::water(), STANDARD_GRAVITY, 1.0, -1.0).is_err());
    }

    #[test]
    fn standard_gravity_wrapper_matches_explicit() {
        let f = Fluid::water();
        let a = gauge_pressure_standard(&f, 3.0).unwrap();
        let b = gauge_pressure(&f, STANDARD_GRAVITY, 3.0).unwrap();
        assert!((a - b).abs() < EPS);
    }
}
