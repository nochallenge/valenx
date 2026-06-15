//! The piping **system resistance curve**.
//!
//! The head a piping system demands to pass a flow `Q` is a fixed static
//! lift plus a friction loss that grows with the square of flow:
//!
//! ```text
//!   H_sys(Q) = H_static + K * Q^2
//! ```
//!
//! `H_static` is the elevation (and any pressure) difference the pump
//! must overcome at zero flow; it can be **negative** for a flooded /
//! downhill arrangement where the static head helps the flow. `K ≥ 0` is
//! a lumped resistance coefficient (units m·s²/m⁶) that rolls together
//! pipe friction and fittings. Because the loss term is quadratic, the
//! curve steepens as flow rises — the defining shape of a system curve.

use serde::{Deserialize, Serialize};

use crate::error::{require_finite, require_non_negative, PumpError};

/// A quadratic system head curve `H = H_static + K·Q²`.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct SystemCurve {
    /// Static head `H_static`, in metres. May be negative (downhill /
    /// flooded suction where gravity assists the flow).
    pub static_head_m: f64,
    /// Resistance coefficient `K ≥ 0`, in m·s²/m⁶ (i.e. metres of head
    /// per (m³/s)²). Larger `K` is a more restrictive system.
    pub resistance_k: f64,
}

impl SystemCurve {
    /// Build a validated system curve.
    ///
    /// `static_head_m` must be finite (any sign). `resistance_k` must be
    /// finite and non-negative.
    ///
    /// # Errors
    ///
    /// Returns [`PumpError::BadParameter`] if `resistance_k` is negative
    /// or either value is non-finite.
    pub fn new(static_head_m: f64, resistance_k: f64) -> Result<Self, PumpError> {
        Ok(Self {
            static_head_m: require_finite("static_head_m", static_head_m)?,
            resistance_k: require_non_negative("resistance_k", resistance_k)?,
        })
    }

    /// The head the system requires to pass volumetric flow `flow_m3s`,
    /// in metres: `H_static + K·Q²`.
    ///
    /// The flow is squared, so the sign of `flow_m3s` does not matter;
    /// passing a negative flow returns the same head as its magnitude.
    pub fn head_m(&self, flow_m3s: f64) -> f64 {
        self.static_head_m + self.resistance_k * flow_m3s * flow_m3s
    }

    /// The flow at which this system demands `target_head_m` of head, the
    /// inverse of [`SystemCurve::head_m`].
    ///
    /// Solving `H_static + K·Q² = H_target` for the non-negative root
    /// gives `Q = sqrt((H_target − H_static) / K)`.
    ///
    /// # Errors
    ///
    /// Returns [`PumpError::BadParameter`] if `target_head_m` is
    /// non-finite. Returns [`PumpError::Inconsistent`] if `K` is zero (a
    /// flat curve has no single flow for a head) or if `target_head_m`
    /// lies below `static_head_m` (the system can never demand less than
    /// its static head).
    pub fn flow_at_head_m3s(&self, target_head_m: f64) -> Result<f64, PumpError> {
        let target_head_m = require_finite("target_head_m", target_head_m)?;
        if self.resistance_k == 0.0 {
            return Err(PumpError::Inconsistent(
                "system curve has zero resistance, so head is independent of flow".to_string(),
            ));
        }
        let numerator = target_head_m - self.static_head_m;
        if numerator < 0.0 {
            return Err(PumpError::Inconsistent(format!(
                "target head {target_head_m} is below the static head {static}",
                static = self.static_head_m
            )));
        }
        Ok((numerator / self.resistance_k).sqrt())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    #[test]
    fn head_is_static_plus_k_q_squared() {
        let sys = SystemCurve::new(10.0, 2000.0).unwrap();
        // At Q = 0 the head is purely static.
        assert!((sys.head_m(0.0) - 10.0).abs() < EPS);
        // At Q = 0.1: 10 + 2000*0.01 = 30.
        assert!((sys.head_m(0.1) - 30.0).abs() < EPS);
    }

    #[test]
    fn loss_rises_with_square_of_flow() {
        let sys = SystemCurve::new(0.0, 500.0).unwrap();
        // Doubling the flow must quadruple the (purely frictional) head.
        let h1 = sys.head_m(0.2);
        let h2 = sys.head_m(0.4);
        assert!((h2 / h1 - 4.0).abs() < 1e-9);
        // Tripling -> 9x.
        let h3 = sys.head_m(0.6);
        assert!((h3 / h1 - 9.0).abs() < 1e-9);
    }

    #[test]
    fn negative_static_head_is_allowed() {
        // Flooded suction: 5 m of assisting head, plus friction.
        let sys = SystemCurve::new(-5.0, 100.0).unwrap();
        // At Q=0.1: -5 + 100*0.01 = -4.
        assert!((sys.head_m(0.1) - (-4.0)).abs() < EPS);
    }

    #[test]
    fn flow_at_head_inverts_head_at_flow() {
        let sys = SystemCurve::new(8.0, 1500.0).unwrap();
        let q = 0.12;
        let h = sys.head_m(q);
        let recovered = sys.flow_at_head_m3s(h).unwrap();
        assert!((recovered - q).abs() < 1e-9);
    }

    #[test]
    fn flow_sign_does_not_matter() {
        let sys = SystemCurve::new(3.0, 250.0).unwrap();
        assert!((sys.head_m(0.3) - sys.head_m(-0.3)).abs() < EPS);
    }

    #[test]
    fn flat_curve_has_no_flow_for_head() {
        let sys = SystemCurve::new(5.0, 0.0).unwrap();
        let err = sys.flow_at_head_m3s(5.0).unwrap_err();
        assert_eq!(err.code(), "pump.inconsistent");
    }

    #[test]
    fn head_below_static_is_unreachable() {
        let sys = SystemCurve::new(12.0, 800.0).unwrap();
        let err = sys.flow_at_head_m3s(10.0).unwrap_err();
        assert_eq!(err.code(), "pump.inconsistent");
    }

    #[test]
    fn rejects_negative_resistance() {
        let err = SystemCurve::new(1.0, -1.0).unwrap_err();
        assert_eq!(err.code(), "pump.bad_parameter");
    }
}
