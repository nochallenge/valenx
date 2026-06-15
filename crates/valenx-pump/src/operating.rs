//! The pump/system **operating point**.
//!
//! A pump runs where its head curve crosses the system curve — the one
//! flow at which the head the pump *produces* equals the head the system
//! *demands*. Modelling the pump curve as the idealised downward parabola
//!
//! ```text
//!   H_pump(Q) = H0 - a * Q^2        (a > 0, shut-off head H0 at Q=0)
//! ```
//!
//! and the system curve as `H_sys(Q) = H_static + K·Q²`
//! ([`crate::system::SystemCurve`]), equating the two heads gives a
//! quadratic with a single non-negative root:
//!
//! ```text
//!   H0 - a*Q^2 = H_static + K*Q^2
//!   Q* = sqrt( (H0 - H_static) / (a + K) )
//!   H* = H_static + K * Q*^2
//! ```
//!
//! A solution exists only when the pump's shut-off head `H0` exceeds the
//! system's static head `H_static`; otherwise the pump cannot even start
//! moving liquid and there is no operating point.

use serde::{Deserialize, Serialize};

use crate::error::{require_finite, require_positive, PumpError};
use crate::system::SystemCurve;

/// An idealised centrifugal-pump head curve `H = H₀ − a·Q²`.
///
/// `H₀` is the shut-off (zero-flow) head and `a > 0` sets how fast the
/// head droops with flow. This is the parabolic teaching model, not a
/// fitted manufacturer curve (see the crate-level "Honest scope").
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct PumpCurve {
    /// Shut-off head `H₀`, in metres of fluid (the head at zero flow).
    pub shutoff_head_m: f64,
    /// Droop coefficient `a > 0`, in m·s²/m⁶ (metres of head lost per
    /// (m³/s)²).
    pub droop_a: f64,
}

impl PumpCurve {
    /// Build a validated pump curve.
    ///
    /// `shutoff_head_m` must be finite. `droop_a` must be finite and
    /// strictly positive (a flat pump curve `a = 0` has no unique
    /// intersection with a flat system).
    ///
    /// # Errors
    ///
    /// Returns [`PumpError::BadParameter`] for an out-of-domain value.
    pub fn new(shutoff_head_m: f64, droop_a: f64) -> Result<Self, PumpError> {
        Ok(Self {
            shutoff_head_m: require_finite("shutoff_head_m", shutoff_head_m)?,
            droop_a: require_positive("droop_a", droop_a)?,
        })
    }

    /// The head the pump produces at volumetric flow `flow_m3s`, in
    /// metres: `H₀ − a·Q²`.
    pub fn head_m(&self, flow_m3s: f64) -> f64 {
        self.shutoff_head_m - self.droop_a * flow_m3s * flow_m3s
    }
}

/// The flow/head pair where a pump and a system run together.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct OperatingPoint {
    /// Operating flow `Q*`, in m³/s.
    pub flow_m3s: f64,
    /// Operating head `H*`, in metres of fluid.
    pub head_m: f64,
}

/// Solve for the operating point where `pump` and `system` heads are
/// equal.
///
/// Returns the unique non-negative-flow intersection
/// `Q* = sqrt((H₀ − H_static)/(a + K))`, `H* = H_static + K·Q*²`.
///
/// # Errors
///
/// Returns [`PumpError::Inconsistent`] if the pump's shut-off head does
/// not exceed the system's static head (`H₀ ≤ H_static`), in which case
/// the pump cannot deliver any flow and no operating point exists.
///
/// # Examples
///
/// ```
/// use valenx_pump::{
///     operating::{operating_point, PumpCurve},
///     system::SystemCurve,
/// };
///
/// // Pump: 50 m shut-off, droop a = 1000. System: 10 m static, K = 4000.
/// let pump = PumpCurve::new(50.0, 1000.0).unwrap();
/// let sys = SystemCurve::new(10.0, 4000.0).unwrap();
/// let op = operating_point(&pump, &sys).unwrap();
/// // Q* = sqrt((50-10)/(1000+4000)) = sqrt(0.008) = 0.08944...
/// assert!((op.flow_m3s - (40.0_f64 / 5000.0).sqrt()).abs() < 1e-12);
/// // Both curves give the same head there.
/// assert!((pump.head_m(op.flow_m3s) - sys.head_m(op.flow_m3s)).abs() < 1e-9);
/// ```
pub fn operating_point(
    pump: &PumpCurve,
    system: &SystemCurve,
) -> Result<OperatingPoint, PumpError> {
    let head_excess = pump.shutoff_head_m - system.static_head_m;
    if head_excess <= 0.0 {
        return Err(PumpError::Inconsistent(format!(
            "pump shut-off head {h0} does not exceed system static head {hs}; no flow",
            h0 = pump.shutoff_head_m,
            hs = system.static_head_m
        )));
    }
    // a + K is strictly positive: a > 0 by construction, K >= 0.
    let denom = pump.droop_a + system.resistance_k;
    let flow_m3s = (head_excess / denom).sqrt();
    let head_m = system.head_m(flow_m3s);
    Ok(OperatingPoint { flow_m3s, head_m })
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn operating_point_balances_both_heads() {
        let pump = PumpCurve::new(50.0, 1000.0).unwrap();
        let sys = SystemCurve::new(10.0, 4000.0).unwrap();
        let op = operating_point(&pump, &sys).unwrap();
        // The defining property: pump head == system head at Q*.
        assert!((pump.head_m(op.flow_m3s) - sys.head_m(op.flow_m3s)).abs() < EPS);
        assert!((op.head_m - sys.head_m(op.flow_m3s)).abs() < EPS);
    }

    #[test]
    fn flow_matches_closed_form_root() {
        let pump = PumpCurve::new(50.0, 1000.0).unwrap();
        let sys = SystemCurve::new(10.0, 4000.0).unwrap();
        let op = operating_point(&pump, &sys).unwrap();
        // Q* = sqrt((50-10)/(1000+4000)) = sqrt(0.008).
        let expected_q = (40.0_f64 / 5000.0).sqrt();
        assert!((op.flow_m3s - expected_q).abs() < EPS);
        // H* = 10 + 4000 * 0.008 = 42.
        assert!((op.head_m - 42.0).abs() < EPS);
    }

    #[test]
    fn clean_integer_case() {
        // Choose numbers so Q* is exactly 0.1: need (H0-Hs)/(a+K) = 0.01.
        // H0=20, Hs=4 -> excess 16; a+K = 1600 -> a=600, K=1000.
        let pump = PumpCurve::new(20.0, 600.0).unwrap();
        let sys = SystemCurve::new(4.0, 1000.0).unwrap();
        let op = operating_point(&pump, &sys).unwrap();
        assert!((op.flow_m3s - 0.1).abs() < EPS);
        // H* = 4 + 1000*0.01 = 14; check pump side too: 20 - 600*0.01 = 14.
        assert!((op.head_m - 14.0).abs() < EPS);
        assert!((pump.head_m(0.1) - 14.0).abs() < EPS);
    }

    #[test]
    fn stiffer_system_throttles_flow_and_raises_head() {
        let pump = PumpCurve::new(40.0, 800.0).unwrap();
        let soft = SystemCurve::new(5.0, 1000.0).unwrap();
        let stiff = SystemCurve::new(5.0, 4000.0).unwrap();
        let op_soft = operating_point(&pump, &soft).unwrap();
        let op_stiff = operating_point(&pump, &stiff).unwrap();
        // A more restrictive system moves the operating point back up the
        // pump curve: less flow, more head.
        assert!(op_stiff.flow_m3s < op_soft.flow_m3s);
        assert!(op_stiff.head_m > op_soft.head_m);
    }

    #[test]
    fn zero_resistance_system_meets_at_static_head() {
        // K = 0: system is flat at H_static, pump meets it where
        // H0 - a Q^2 = H_static -> Q = sqrt((H0-Hs)/a), H* = H_static.
        let pump = PumpCurve::new(30.0, 1200.0).unwrap();
        let sys = SystemCurve::new(6.0, 0.0).unwrap();
        let op = operating_point(&pump, &sys).unwrap();
        assert!((op.head_m - 6.0).abs() < EPS);
        let expected_q = (24.0_f64 / 1200.0).sqrt();
        assert!((op.flow_m3s - expected_q).abs() < EPS);
    }

    #[test]
    fn no_solution_when_static_head_exceeds_shutoff() {
        let pump = PumpCurve::new(10.0, 1000.0).unwrap();
        let sys = SystemCurve::new(15.0, 2000.0).unwrap();
        let err = operating_point(&pump, &sys).unwrap_err();
        assert_eq!(err.code(), "pump.inconsistent");
    }

    #[test]
    fn equal_heads_give_zero_flow_solution_is_rejected() {
        // H0 == H_static is the degenerate boundary: excess is zero.
        let pump = PumpCurve::new(12.0, 1000.0).unwrap();
        let sys = SystemCurve::new(12.0, 2000.0).unwrap();
        assert!(operating_point(&pump, &sys).is_err());
    }

    #[test]
    fn rejects_non_positive_droop() {
        assert!(PumpCurve::new(20.0, 0.0).is_err());
        assert!(PumpCurve::new(20.0, -5.0).is_err());
    }
}
