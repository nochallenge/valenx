//! The centrifugal-pump **affinity laws**.
//!
//! For one fixed-geometry pump running the same fluid, changing the shaft
//! speed `N` scales the duty point by simple power laws:
//!
//! ```text
//!   Q2 / Q1 = (N2 / N1)
//!   H2 / H1 = (N2 / N1)^2
//!   P2 / P1 = (N2 / N1)^3
//! ```
//!
//! so doubling the speed doubles the flow, quadruples the head, and
//! octuples the shaft power. Only the speed *ratio* matters, so `N` may
//! be given in any consistent unit (rpm, rad/s, …).
//!
//! The laws assume ideal scaling — constant efficiency, fully turbulent
//! flow, and unchanged geometry. They are the standard first-order tool
//! for re-rating a measured curve to a new speed; see the crate-level
//! "Honest scope" for the caveats.

use serde::{Deserialize, Serialize};

use crate::error::{require_non_negative, require_positive, PumpError};

/// A single measured (or rated) pump operating point at a known speed.
///
/// All four quantities are tied together by the affinity laws when the
/// speed changes; [`scale_to_speed`] produces the corresponding point at
/// any other speed.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct DutyPoint {
    /// Shaft speed `N` (rpm, rad/s, … — any unit; only ratios are used).
    pub speed: f64,
    /// Volumetric flow rate `Q`, in m³/s.
    pub flow_m3s: f64,
    /// Total head `H`, in metres of pumped fluid.
    pub head_m: f64,
    /// Shaft power `P`, in watts.
    pub power_w: f64,
}

impl DutyPoint {
    /// Build a validated duty point.
    ///
    /// `speed` must be finite and strictly positive (a ratio is taken
    /// against it). `flow_m3s`, `head_m` and `power_w` must be finite and
    /// non-negative.
    ///
    /// # Errors
    ///
    /// Returns [`PumpError::BadParameter`] if any value is out of domain.
    pub fn new(speed: f64, flow_m3s: f64, head_m: f64, power_w: f64) -> Result<Self, PumpError> {
        Ok(Self {
            speed: require_positive("speed", speed)?,
            flow_m3s: require_non_negative("flow_m3s", flow_m3s)?,
            head_m: require_non_negative("head_m", head_m)?,
            power_w: require_non_negative("power_w", power_w)?,
        })
    }
}

/// Scale a measured [`DutyPoint`] to a new shaft speed using the affinity
/// laws.
///
/// With `r = new_speed / base.speed`, the returned point has
/// `Q = base.Q · r`, `H = base.H · r²` and `P = base.P · r³`.
///
/// # Errors
///
/// Returns [`PumpError::BadParameter`] if `new_speed` is not finite and
/// strictly positive.
///
/// # Examples
///
/// ```
/// use valenx_pump::affinity::{scale_to_speed, DutyPoint};
///
/// let base = DutyPoint::new(1000.0, 0.02, 12.0, 5000.0).unwrap();
/// let doubled = scale_to_speed(&base, 2000.0).unwrap();
/// assert!((doubled.flow_m3s - 0.04).abs() < 1e-12);
/// assert!((doubled.head_m - 48.0).abs() < 1e-12);
/// assert!((doubled.power_w - 40_000.0).abs() < 1e-9);
/// ```
pub fn scale_to_speed(base: &DutyPoint, new_speed: f64) -> Result<DutyPoint, PumpError> {
    let new_speed = require_positive("new_speed", new_speed)?;
    let r = new_speed / base.speed;
    Ok(DutyPoint {
        speed: new_speed,
        flow_m3s: base.flow_m3s * r,
        head_m: base.head_m * r * r,
        power_w: base.power_w * r * r * r,
    })
}

/// Scale a flow rate by the affinity law `Q ∝ N`.
///
/// Returns `flow_m3s · (new_speed / base_speed)`.
///
/// # Errors
///
/// Returns [`PumpError::BadParameter`] if either speed is not finite and
/// strictly positive, or if `flow_m3s` is not finite and non-negative.
pub fn scale_flow(flow_m3s: f64, base_speed: f64, new_speed: f64) -> Result<f64, PumpError> {
    let flow_m3s = require_non_negative("flow_m3s", flow_m3s)?;
    let base_speed = require_positive("base_speed", base_speed)?;
    let new_speed = require_positive("new_speed", new_speed)?;
    Ok(flow_m3s * (new_speed / base_speed))
}

/// Scale a head by the affinity law `H ∝ N²`.
///
/// Returns `head_m · (new_speed / base_speed)²`.
///
/// # Errors
///
/// Returns [`PumpError::BadParameter`] if either speed is not finite and
/// strictly positive, or if `head_m` is not finite and non-negative.
pub fn scale_head(head_m: f64, base_speed: f64, new_speed: f64) -> Result<f64, PumpError> {
    let head_m = require_non_negative("head_m", head_m)?;
    let base_speed = require_positive("base_speed", base_speed)?;
    let new_speed = require_positive("new_speed", new_speed)?;
    let r = new_speed / base_speed;
    Ok(head_m * r * r)
}

/// Scale a power by the affinity law `P ∝ N³`.
///
/// Returns `power_w · (new_speed / base_speed)³`.
///
/// # Errors
///
/// Returns [`PumpError::BadParameter`] if either speed is not finite and
/// strictly positive, or if `power_w` is not finite and non-negative.
pub fn scale_power(power_w: f64, base_speed: f64, new_speed: f64) -> Result<f64, PumpError> {
    let power_w = require_non_negative("power_w", power_w)?;
    let base_speed = require_positive("base_speed", base_speed)?;
    let new_speed = require_positive("new_speed", new_speed)?;
    let r = new_speed / base_speed;
    Ok(power_w * r * r * r)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn double_speed_doubles_flow_quadruples_head_octuples_power() {
        let base = DutyPoint::new(1450.0, 0.05, 30.0, 18_000.0).unwrap();
        let fast = scale_to_speed(&base, 2900.0).unwrap();
        // Q ×2, H ×4, P ×8 — the canonical affinity ground truth.
        assert!((fast.flow_m3s - 0.10).abs() < EPS);
        assert!((fast.head_m - 120.0).abs() < EPS);
        assert!((fast.power_w - 144_000.0).abs() < 1e-6);
        assert!((fast.speed - 2900.0).abs() < EPS);
    }

    #[test]
    fn half_speed_halves_flow_quarters_head_eighths_power() {
        let base = DutyPoint::new(3000.0, 0.08, 40.0, 24_000.0).unwrap();
        let slow = scale_to_speed(&base, 1500.0).unwrap();
        assert!((slow.flow_m3s - 0.04).abs() < EPS);
        assert!((slow.head_m - 10.0).abs() < EPS);
        assert!((slow.power_w - 3_000.0).abs() < EPS);
    }

    #[test]
    fn identity_speed_is_a_fixed_point() {
        let base = DutyPoint::new(1750.0, 0.033, 22.5, 9_100.0).unwrap();
        let same = scale_to_speed(&base, 1750.0).unwrap();
        assert!((same.flow_m3s - base.flow_m3s).abs() < EPS);
        assert!((same.head_m - base.head_m).abs() < EPS);
        assert!((same.power_w - base.power_w).abs() < EPS);
    }

    #[test]
    fn scalar_helpers_match_arbitrary_ratio() {
        // 1450 -> 1740 rpm is a ratio of exactly 1.2.
        let r: f64 = 1740.0 / 1450.0;
        let q = scale_flow(0.05, 1450.0, 1740.0).unwrap();
        let h = scale_head(30.0, 1450.0, 1740.0).unwrap();
        let p = scale_power(18_000.0, 1450.0, 1740.0).unwrap();
        assert!((q - 0.05 * r).abs() < EPS);
        assert!((h - 30.0 * r * r).abs() < EPS);
        assert!((p - 18_000.0 * r * r * r).abs() < 1e-6);
        // And r is what we think it is.
        assert!((r - 1.2).abs() < EPS);
    }

    #[test]
    fn power_grows_as_cube_of_speed_ratio() {
        // Tripling the speed must multiply power by 27.
        let base = DutyPoint::new(500.0, 0.01, 4.0, 200.0).unwrap();
        let fast = scale_to_speed(&base, 1500.0).unwrap();
        assert!((fast.power_w - 200.0 * 27.0).abs() < EPS);
        assert!((fast.head_m - 4.0 * 9.0).abs() < EPS);
        assert!((fast.flow_m3s - 0.01 * 3.0).abs() < EPS);
    }

    #[test]
    fn rejects_non_positive_speed() {
        assert!(DutyPoint::new(0.0, 0.01, 1.0, 1.0).is_err());
        let base = DutyPoint::new(1000.0, 0.01, 1.0, 1.0).unwrap();
        let err = scale_to_speed(&base, -10.0).unwrap_err();
        assert_eq!(err.code(), "pump.bad_parameter");
    }

    #[test]
    fn rejects_non_finite_inputs() {
        assert!(DutyPoint::new(1000.0, f64::NAN, 1.0, 1.0).is_err());
        assert!(DutyPoint::new(1000.0, 0.01, f64::INFINITY, 1.0).is_err());
    }
}
