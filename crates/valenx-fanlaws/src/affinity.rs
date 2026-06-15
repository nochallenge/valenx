//! Fan affinity laws: scaling a single fan between two operating points.
//!
//! For one geometrically fixed fan, changing only its rotational speed
//! from `n1` to `n2` (with the same gas density and an unchanged system
//! curve) scales the operating point by powers of the speed ratio
//! `r = n2 / n1`:
//!
//! ```text
//! Q2 / Q1 = r          (volumetric flow is linear in speed)
//! P2 / P1 = r^2        (pressure rise grows with the square of speed)
//! W2 / W1 = r^3        (shaft power grows with the cube of speed)
//! ```
//!
//! When instead only the gas density changes from `rho1` to `rho2` at
//! fixed speed (e.g. the same fan moving hotter, thinner air), the
//! volumetric flow is unchanged but pressure and power scale linearly
//! with density:
//!
//! ```text
//! Q2 / Q1 = 1
//! P2 / P1 = rho2 / rho1
//! W2 / W1 = rho2 / rho1
//! ```
//!
//! ## Honest scope
//!
//! These are idealised similarity relations. They assume the fan
//! efficiency is unchanged across the speed step (true only over modest
//! ratios where Reynolds-number and tip-clearance effects are small),
//! incompressible flow (so density is treated as a free multiplier),
//! and that the system the fan works against follows a fixed
//! quadratic resistance curve so the operating point stays on the same
//! similarity ray. Real fans deviate; treat the output as a first-order
//! estimate, not a guarantee.

use crate::error::{require_positive, FanLawError};

/// A single fan operating point: volumetric flow, pressure rise, shaft
/// power, rotational speed, and the gas density at which it was
/// measured.
///
/// Units are not fixed by this type — it is dimension-agnostic — but
/// they must be **self-consistent** across the two points you scale
/// between (e.g. both flows in m^3/s, both speeds in rev/min). A common
/// SI choice is flow in m^3/s, pressure in pascals, power in watts,
/// speed in rad/s or rev/min, and density in kg/m^3.
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct OperatingPoint {
    /// Volumetric flow rate `Q` (non-negative).
    pub flow: f64,
    /// Static (or total) pressure rise `dP` across the fan
    /// (non-negative).
    pub pressure: f64,
    /// Shaft (mechanical) power `W` delivered to the fan (non-negative).
    pub power: f64,
    /// Rotational speed `N` of the impeller (strictly positive).
    pub speed: f64,
    /// Gas density `rho` at this point (strictly positive).
    pub density: f64,
}

impl OperatingPoint {
    /// Construct a validated operating point.
    ///
    /// `flow`, `pressure`, and `power` may be zero (shut-off / no-flow)
    /// but not negative; `speed` and `density` must be strictly
    /// positive because they appear in ratios.
    ///
    /// # Errors
    ///
    /// Returns [`FanLawError`] if any value is non-finite, if a
    /// non-negative quantity is negative, or if a strictly-positive
    /// quantity is zero or negative.
    pub fn new(
        flow: f64,
        pressure: f64,
        power: f64,
        speed: f64,
        density: f64,
    ) -> Result<Self, FanLawError> {
        use crate::error::require_non_negative;
        Ok(Self {
            flow: require_non_negative("flow", flow)?,
            pressure: require_non_negative("pressure", pressure)?,
            power: require_non_negative("power", power)?,
            speed: require_positive("speed", speed)?,
            density: require_positive("density", density)?,
        })
    }
}

/// Scale only the **flow** of a fan from speed `n1` to speed `n2`.
///
/// Flow is linear in speed, so `q2 = q1 * (n2 / n1)`.
///
/// # Errors
///
/// Returns [`FanLawError`] if `q1` is negative or non-finite, or if
/// either speed is not strictly positive.
pub fn scale_flow(q1: f64, n1: f64, n2: f64) -> Result<f64, FanLawError> {
    use crate::error::require_non_negative;
    let q1 = require_non_negative("q1", q1)?;
    let n1 = require_positive("n1", n1)?;
    let n2 = require_positive("n2", n2)?;
    Ok(q1 * (n2 / n1))
}

/// Scale only the **pressure rise** of a fan from speed `n1` to speed
/// `n2` at fixed density.
///
/// Pressure grows with the square of speed: `p2 = p1 * (n2 / n1)^2`.
///
/// # Errors
///
/// Returns [`FanLawError`] if `p1` is negative or non-finite, or if
/// either speed is not strictly positive.
pub fn scale_pressure(p1: f64, n1: f64, n2: f64) -> Result<f64, FanLawError> {
    use crate::error::require_non_negative;
    let p1 = require_non_negative("p1", p1)?;
    let n1 = require_positive("n1", n1)?;
    let n2 = require_positive("n2", n2)?;
    let r = n2 / n1;
    Ok(p1 * r * r)
}

/// Scale only the **shaft power** of a fan from speed `n1` to speed
/// `n2` at fixed density.
///
/// Power grows with the cube of speed: `w2 = w1 * (n2 / n1)^3`.
///
/// # Errors
///
/// Returns [`FanLawError`] if `w1` is negative or non-finite, or if
/// either speed is not strictly positive.
pub fn scale_power(w1: f64, n1: f64, n2: f64) -> Result<f64, FanLawError> {
    use crate::error::require_non_negative;
    let w1 = require_non_negative("w1", w1)?;
    let n1 = require_positive("n1", n1)?;
    let n2 = require_positive("n2", n2)?;
    let r = n2 / n1;
    Ok(w1 * r * r * r)
}

/// Density correction for **pressure** at fixed speed: a fan develops
/// pressure in proportion to the gas density, so moving the same fan in
/// denser air raises the pressure rise linearly.
///
/// `p2 = p1 * (rho2 / rho1)`.
///
/// # Errors
///
/// Returns [`FanLawError`] if `p1` is negative or non-finite, or if
/// either density is not strictly positive.
pub fn correct_pressure_for_density(p1: f64, rho1: f64, rho2: f64) -> Result<f64, FanLawError> {
    use crate::error::require_non_negative;
    let p1 = require_non_negative("p1", p1)?;
    let rho1 = require_positive("rho1", rho1)?;
    let rho2 = require_positive("rho2", rho2)?;
    Ok(p1 * (rho2 / rho1))
}

/// Density correction for **power** at fixed speed: shaft power also
/// scales linearly with density (the work rate to move a denser fluid
/// rises in the same proportion as the pressure it develops).
///
/// `w2 = w1 * (rho2 / rho1)`.
///
/// # Errors
///
/// Returns [`FanLawError`] if `w1` is negative or non-finite, or if
/// either density is not strictly positive.
pub fn correct_power_for_density(w1: f64, rho1: f64, rho2: f64) -> Result<f64, FanLawError> {
    use crate::error::require_non_negative;
    let w1 = require_non_negative("w1", w1)?;
    let rho1 = require_positive("rho1", rho1)?;
    let rho2 = require_positive("rho2", rho2)?;
    Ok(w1 * (rho2 / rho1))
}

/// Apply the full affinity transform to an [`OperatingPoint`]: change
/// the impeller speed to `n2` **and** the gas density to `rho2` in one
/// step, combining the speed powers with the linear density factor.
///
/// The returned point carries the scaled flow, pressure, and power, the
/// new `speed = n2`, and the new `density = rho2`:
///
/// ```text
/// r       = n2 / point.speed
/// d       = rho2 / point.density
/// flow    = Q * r            (density does not affect volumetric flow)
/// pressure= P * r^2 * d
/// power   = W * r^3 * d
/// ```
///
/// # Errors
///
/// Returns [`FanLawError`] if `n2` or `rho2` is not strictly positive
/// or non-finite. (The `point` was already validated at construction.)
pub fn scale_operating_point(
    point: &OperatingPoint,
    n2: f64,
    rho2: f64,
) -> Result<OperatingPoint, FanLawError> {
    let n2 = require_positive("n2", n2)?;
    let rho2 = require_positive("rho2", rho2)?;
    let r = n2 / point.speed;
    let d = rho2 / point.density;
    Ok(OperatingPoint {
        flow: point.flow * r,
        pressure: point.pressure * r * r * d,
        power: point.power * r * r * r * d,
        speed: n2,
        density: rho2,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for analytic float comparisons.
    const EPS: f64 = 1e-9;

    #[test]
    fn doubling_speed_doubles_flow_quadruples_pressure_octuples_power() {
        // Ground truth: at r = 2, Q->2Q, dP->4dP, W->8W.
        let n1 = 1000.0;
        let n2 = 2000.0;
        assert!((scale_flow(10.0, n1, n2).unwrap() - 20.0).abs() < EPS);
        assert!((scale_pressure(50.0, n1, n2).unwrap() - 200.0).abs() < EPS);
        assert!((scale_power(3.0, n1, n2).unwrap() - 24.0).abs() < EPS);
    }

    #[test]
    fn halving_speed_halves_flow_quarters_pressure_eighths_power() {
        // r = 1/2: Q->Q/2, dP->dP/4, W->W/8.
        let n1 = 1800.0;
        let n2 = 900.0;
        assert!((scale_flow(8.0, n1, n2).unwrap() - 4.0).abs() < EPS);
        assert!((scale_pressure(80.0, n1, n2).unwrap() - 20.0).abs() < EPS);
        assert!((scale_power(64.0, n1, n2).unwrap() - 8.0).abs() < EPS);
    }

    #[test]
    fn flow_is_exactly_linear_in_speed() {
        // Q(k*n)/Q(n) == k for a sweep of ratios.
        let q1 = 12.5;
        let n1 = 1450.0;
        for k in [0.25_f64, 0.5, 1.0, 1.7, 3.0, 10.0] {
            let n2 = k * n1;
            let q2 = scale_flow(q1, n1, n2).unwrap();
            assert!(
                (q2 - k * q1).abs() < EPS,
                "k={k}: expected {expected}, got {q2}",
                expected = k * q1
            );
        }
    }

    #[test]
    fn pressure_scales_with_square_and_power_with_cube_over_a_sweep() {
        let p1 = 250.0;
        let w1 = 5.0;
        let n1 = 1200.0;
        for k in [0.5_f64, 0.8, 1.0, 1.25, 2.0, 3.5] {
            let n2 = k * n1;
            let p2 = scale_pressure(p1, n1, n2).unwrap();
            let w2 = scale_power(w1, n1, n2).unwrap();
            assert!(
                (p2 - p1 * k * k).abs() < 1e-7,
                "k={k}: pressure {p2} != {expected}",
                expected = p1 * k * k
            );
            assert!(
                (w2 - w1 * k * k * k).abs() < 1e-7,
                "k={k}: power {w2} != {expected}",
                expected = w1 * k * k * k
            );
        }
    }

    #[test]
    fn pressure_and_power_scale_linearly_with_density() {
        // rho doubled -> pressure and power doubled, at fixed speed.
        assert!((correct_pressure_for_density(100.0, 1.0, 2.0).unwrap() - 200.0).abs() < EPS);
        assert!((correct_power_for_density(4.0, 1.0, 2.0).unwrap() - 8.0).abs() < EPS);
        // Standard-air -> high-altitude thinner air halves both.
        assert!((correct_pressure_for_density(300.0, 1.2, 0.6).unwrap() - 150.0).abs() < EPS);
        assert!((correct_power_for_density(10.0, 1.2, 0.6).unwrap() - 5.0).abs() < EPS);
    }

    #[test]
    fn identity_transforms_are_no_ops() {
        // Same speed -> unchanged everywhere.
        let n = 1500.0;
        assert!((scale_flow(7.0, n, n).unwrap() - 7.0).abs() < EPS);
        assert!((scale_pressure(7.0, n, n).unwrap() - 7.0).abs() < EPS);
        assert!((scale_power(7.0, n, n).unwrap() - 7.0).abs() < EPS);
        // Same density -> unchanged.
        assert!((correct_pressure_for_density(7.0, 1.2, 1.2).unwrap() - 7.0).abs() < EPS);
    }

    #[test]
    fn full_operating_point_transform_combines_speed_and_density() {
        // Q=10, P=100, W=2 at N=1000, rho=1.2.
        // Go to N=2000 (r=2), rho=0.6 (d=0.5).
        // flow  = 10 * 2          = 20
        // press = 100 * 4 * 0.5   = 200
        // power = 2  * 8 * 0.5    = 8
        let p0 = OperatingPoint::new(10.0, 100.0, 2.0, 1000.0, 1.2).unwrap();
        let p1 = scale_operating_point(&p0, 2000.0, 0.6).unwrap();
        assert!((p1.flow - 20.0).abs() < EPS);
        assert!((p1.pressure - 200.0).abs() < EPS);
        assert!((p1.power - 8.0).abs() < EPS);
        assert!((p1.speed - 2000.0).abs() < EPS);
        assert!((p1.density - 0.6).abs() < EPS);
    }

    #[test]
    fn full_transform_matches_individual_laws_at_fixed_density() {
        // With rho unchanged the combined transform must reproduce the
        // three single-quantity scalings exactly.
        let p0 = OperatingPoint::new(15.0, 220.0, 6.5, 980.0, 1.18).unwrap();
        let n2 = 1450.0;
        let scaled = scale_operating_point(&p0, n2, p0.density).unwrap();
        assert!((scaled.flow - scale_flow(p0.flow, p0.speed, n2).unwrap()).abs() < EPS);
        assert!((scaled.pressure - scale_pressure(p0.pressure, p0.speed, n2).unwrap()).abs() < EPS);
        assert!((scaled.power - scale_power(p0.power, p0.speed, n2).unwrap()).abs() < EPS);
    }

    #[test]
    fn round_trip_speed_change_recovers_original() {
        // Scaling up then back down by the inverse ratio returns the
        // start (within float tolerance).
        let q = 9.0;
        let n1 = 1100.0;
        let n2 = 2750.0;
        let up = scale_power(q, n1, n2).unwrap();
        let back = scale_power(up, n2, n1).unwrap();
        assert!((back - q).abs() < 1e-7, "round trip gave {back}");
    }

    #[test]
    fn rejects_non_positive_speed_and_density() {
        assert!(scale_flow(1.0, 0.0, 100.0).is_err());
        assert!(scale_pressure(1.0, 100.0, -5.0).is_err());
        assert!(correct_pressure_for_density(1.0, 0.0, 1.0).is_err());
        assert!(OperatingPoint::new(1.0, 1.0, 1.0, 0.0, 1.2).is_err());
        assert!(OperatingPoint::new(-1.0, 1.0, 1.0, 1.0, 1.2).is_err());
    }
}
