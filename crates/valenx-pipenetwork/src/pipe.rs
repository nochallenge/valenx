//! Pipes: the edges of a hydraulic network.
//!
//! A [`Pipe`] carries a signed assumed flow `q` and a resistance
//! coefficient `k` such that the *signed* head loss along the pipe (in the
//! direction of positive `q`) is
//!
//! ```text
//! h = k * q * |q|
//! ```
//!
//! and its magnitude is `k * q^2`. Writing the loss as `q * |q|` rather
//! than `q^2` keeps the sign of the head loss tied to the sign of the
//! flow, which is what makes the Hardy-Cross loop sum vanish at
//! convergence regardless of flow direction.
//!
//! This module also provides two textbook ways to *compute* `k` from pipe
//! geometry and a friction model, so callers do not have to derive it by
//! hand:
//!
//! - [`darcy_weisbach_k`] for the Darcy-Weisbach loss law, and
//! - [`hazen_williams_k`] for the Hazen-Williams empirical law.

use crate::error::NetworkError;
use serde::{Deserialize, Serialize};

/// Acceleration due to gravity used by the Darcy-Weisbach resistance
/// helper, in metres per second squared (standard gravity).
pub const GRAVITY_M_S2: f64 = 9.806_65;

/// A single pipe (edge) in a hydraulic network.
///
/// The pipe stores an *assumed* volumetric flow `q` and a resistance
/// coefficient `k`. The Hardy-Cross solver mutates `q` in place as it
/// iterates; `k` is treated as a constant of the pipe (it does not depend
/// on the flow under the quadratic loss model used here).
///
/// Units are deliberately left to the caller, but must be *self-consistent*
/// across the whole network: if `q` is in m^3/s and `k` is chosen so that
/// `k * q^2` is in metres of head, then every head loss in the network is
/// in metres. The [`darcy_weisbach_k`] and [`hazen_williams_k`] helpers
/// document the unit system they assume.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Pipe {
    /// Resistance coefficient `k >= 0` in the loss law `h = k * q * |q|`.
    pub k: f64,
    /// Current (assumed, then corrected) signed volumetric flow `q`.
    ///
    /// Positive `q` flows in the pipe's own reference direction; negative
    /// `q` flows the other way. The Hardy-Cross iteration overwrites this
    /// field with successively better estimates.
    pub q: f64,
}

impl Pipe {
    /// Build a pipe from a resistance coefficient and an initial assumed
    /// flow.
    ///
    /// # Errors
    ///
    /// Returns [`NetworkError::BadParameter`] if `k` is not finite or is
    /// negative, or if `q` is not finite. A zero `k` is permitted (an
    /// idealised loss-free pipe) but note that a loop made entirely of
    /// loss-free pipes has an undefined Hardy-Cross correction and the
    /// solver will reject it at solve time.
    pub fn new(k: f64, q: f64) -> Result<Self, NetworkError> {
        if !k.is_finite() || k < 0.0 {
            return Err(NetworkError::bad(
                "k",
                format!("resistance coefficient must be finite and >= 0, got {k}"),
            ));
        }
        if !q.is_finite() {
            return Err(NetworkError::bad(
                "q",
                format!("assumed flow must be finite, got {q}"),
            ));
        }
        Ok(Self { k, q })
    }

    /// Signed head loss along this pipe in its reference direction:
    /// `h = k * q * |q|`.
    ///
    /// The sign follows the sign of `q`, so reversing the flow reverses
    /// the head loss.
    pub fn head_loss(&self) -> f64 {
        self.k * self.q * self.q.abs()
    }

    /// Derivative of the head loss with respect to flow, `dh/dq = 2 k |q|`.
    ///
    /// This is the per-pipe term that goes into the denominator of the
    /// Hardy-Cross loop correction.
    pub fn head_loss_slope(&self) -> f64 {
        2.0 * self.k * self.q.abs()
    }
}

/// Resistance coefficient `k` for the **Darcy-Weisbach** loss law.
///
/// The Darcy-Weisbach head loss for a full circular pipe is
///
/// ```text
/// h_f = f * (L / D) * V^2 / (2 g)
/// ```
///
/// where `V = q / A` and `A = pi D^2 / 4`. Substituting and collecting
/// constants gives `h_f = k * q^2` with
///
/// ```text
/// k = 8 f L / (pi^2 g D^5).
/// ```
///
/// # Parameters
///
/// - `friction_factor` — the dimensionless Darcy friction factor `f`
///   (`>= 0`; e.g. ~0.02 for turbulent flow in commercial steel).
/// - `length_m` — pipe length `L` in metres (`> 0`).
/// - `diameter_m` — internal diameter `D` in metres (`> 0`).
///
/// With these SI inputs the resulting `k` makes `k * q^2` a head in metres
/// when `q` is in cubic metres per second.
///
/// # Errors
///
/// Returns [`NetworkError::BadParameter`] if any argument is non-finite,
/// if `friction_factor` is negative, or if `length_m` or `diameter_m` is
/// not strictly positive.
pub fn darcy_weisbach_k(
    friction_factor: f64,
    length_m: f64,
    diameter_m: f64,
) -> Result<f64, NetworkError> {
    if !friction_factor.is_finite() || friction_factor < 0.0 {
        return Err(NetworkError::bad(
            "friction_factor",
            format!("must be finite and >= 0, got {friction_factor}"),
        ));
    }
    if !length_m.is_finite() || length_m <= 0.0 {
        return Err(NetworkError::bad(
            "length_m",
            format!("must be finite and > 0, got {length_m}"),
        ));
    }
    if !diameter_m.is_finite() || diameter_m <= 0.0 {
        return Err(NetworkError::bad(
            "diameter_m",
            format!("must be finite and > 0, got {diameter_m}"),
        ));
    }
    let pi = std::f64::consts::PI;
    let k = 8.0 * friction_factor * length_m / (pi * pi * GRAVITY_M_S2 * diameter_m.powi(5));
    Ok(k)
}

/// Resistance coefficient `k` for the **Hazen-Williams** empirical law,
/// linearised to the quadratic `h = k * q^2` form used by this solver.
///
/// The true Hazen-Williams head loss varies as `q^1.852`, not `q^2`. The
/// classic Hardy-Cross formulation in this crate uses an exponent of 2, so
/// this helper returns the `k` that reproduces the Hazen-Williams loss
/// *exactly at a chosen reference flow* `q_ref` and approximates it nearby:
///
/// ```text
/// k = h_HW(q_ref) / q_ref^2
/// ```
///
/// where, in SI units,
///
/// ```text
/// h_HW = 10.67 * L * q^1.852 / (C^1.852 * D^4.87).
/// ```
///
/// This is a deliberate, documented simplification — see the crate-level
/// "Honest scope" note. For a network whose flows stay near `q_ref` it is
/// a reasonable engineering approximation; far from `q_ref` the quadratic
/// model diverges from Hazen-Williams.
///
/// # Parameters
///
/// - `c` — Hazen-Williams roughness coefficient `C` (`> 0`; ~100-150 for
///   common pipe materials).
/// - `length_m` — pipe length in metres (`> 0`).
/// - `diameter_m` — internal diameter in metres (`> 0`).
/// - `q_ref` — reference flow in m^3/s at which the quadratic fit is exact
///   (`> 0`).
///
/// # Errors
///
/// Returns [`NetworkError::BadParameter`] if any argument is non-finite,
/// or if any of `c`, `length_m`, `diameter_m`, or `q_ref` is not strictly
/// positive.
pub fn hazen_williams_k(
    c: f64,
    length_m: f64,
    diameter_m: f64,
    q_ref: f64,
) -> Result<f64, NetworkError> {
    for (name, value) in [
        ("c", c),
        ("length_m", length_m),
        ("diameter_m", diameter_m),
        ("q_ref", q_ref),
    ] {
        if !value.is_finite() || value <= 0.0 {
            return Err(NetworkError::bad(
                name,
                format!("must be finite and > 0, got {value}"),
            ));
        }
    }
    let h_hw = 10.67 * length_m * q_ref.powf(1.852) / (c.powf(1.852) * diameter_m.powf(4.87));
    Ok(h_hw / (q_ref * q_ref))
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    #[test]
    fn head_loss_is_k_q_squared_in_magnitude() {
        let p = Pipe::new(3.0, 2.0).unwrap();
        // h = k q |q| = 3 * 2 * 2 = 12; magnitude equals k q^2 = 3 * 4 = 12.
        assert!((p.head_loss() - 12.0).abs() < EPS);
        assert!((p.head_loss().abs() - p.k * p.q * p.q).abs() < EPS);
    }

    #[test]
    fn head_loss_sign_follows_flow_sign() {
        let forward = Pipe::new(3.0, 2.0).unwrap();
        let reverse = Pipe::new(3.0, -2.0).unwrap();
        // Same magnitude, opposite sign — this is why we use q|q|, not q^2.
        assert!((forward.head_loss() + reverse.head_loss()).abs() < EPS);
        assert!(forward.head_loss() > 0.0);
        assert!(reverse.head_loss() < 0.0);
    }

    #[test]
    fn head_loss_slope_is_two_k_abs_q() {
        let p = Pipe::new(3.0, -2.0).unwrap();
        // dh/dq = 2 k |q| = 2 * 3 * 2 = 12, always non-negative.
        assert!((p.head_loss_slope() - 12.0).abs() < EPS);
        assert!(p.head_loss_slope() >= 0.0);
    }

    #[test]
    fn slope_matches_finite_difference_of_loss() {
        // Verify head_loss_slope() really is d/dq [k q|q|] at q != 0.
        let k = 1.7;
        let q = 2.3;
        let p = Pipe::new(k, q).unwrap();
        let dq = 1e-6;
        let hp = k * (q + dq) * (q + dq).abs();
        let hm = k * (q - dq) * (q - dq).abs();
        let numeric = (hp - hm) / (2.0 * dq);
        assert!((p.head_loss_slope() - numeric).abs() < 1e-6);
    }

    #[test]
    fn zero_flow_has_zero_loss_and_zero_slope() {
        let p = Pipe::new(5.0, 0.0).unwrap();
        assert!(p.head_loss().abs() < EPS);
        assert!(p.head_loss_slope().abs() < EPS);
    }

    #[test]
    fn negative_k_is_rejected() {
        let err = Pipe::new(-1.0, 1.0).unwrap_err();
        assert!(matches!(err, NetworkError::BadParameter { name: "k", .. }));
    }

    #[test]
    fn non_finite_inputs_are_rejected() {
        assert!(Pipe::new(f64::NAN, 1.0).is_err());
        assert!(Pipe::new(f64::INFINITY, 1.0).is_err());
        assert!(Pipe::new(1.0, f64::NAN).is_err());
        assert!(Pipe::new(1.0, f64::INFINITY).is_err());
    }

    #[test]
    fn darcy_weisbach_k_matches_closed_form() {
        // k = 8 f L / (pi^2 g D^5).  Choose round numbers and hand-check.
        let f = 0.02;
        let length = 100.0;
        let diameter = 0.3;
        let k = darcy_weisbach_k(f, length, diameter).unwrap();
        let pi = std::f64::consts::PI;
        let expected = 8.0 * f * length / (pi * pi * GRAVITY_M_S2 * diameter.powi(5));
        assert!((k - expected).abs() < EPS * expected.max(1.0));
        // Sanity: a real 300 mm, 100 m, f=0.02 pipe has a sizeable k.
        // h at 0.1 m^3/s should be a few metres of head.
        let h = k * 0.1 * 0.1;
        assert!(h > 0.5 && h < 50.0, "head {h} out of physical band");
    }

    #[test]
    fn darcy_weisbach_k_scales_with_length_and_diameter() {
        let base = darcy_weisbach_k(0.02, 100.0, 0.3).unwrap();
        // Doubling length doubles k.
        let twice_len = darcy_weisbach_k(0.02, 200.0, 0.3).unwrap();
        assert!((twice_len - 2.0 * base).abs() < EPS * base);
        // Halving diameter multiplies k by 2^5 = 32 (k ~ 1/D^5).
        let half_dia = darcy_weisbach_k(0.02, 100.0, 0.15).unwrap();
        assert!((half_dia - 32.0 * base).abs() < 1e-6 * half_dia);
    }

    #[test]
    fn darcy_weisbach_k_rejects_bad_geometry() {
        assert!(darcy_weisbach_k(-0.01, 1.0, 0.1).is_err());
        assert!(darcy_weisbach_k(0.02, 0.0, 0.1).is_err());
        assert!(darcy_weisbach_k(0.02, 1.0, 0.0).is_err());
        assert!(darcy_weisbach_k(f64::NAN, 1.0, 0.1).is_err());
    }

    #[test]
    fn hazen_williams_k_reproduces_loss_at_reference_flow() {
        // By construction k = h_HW(q_ref) / q_ref^2, so k * q_ref^2 must
        // equal the Hazen-Williams loss at q_ref exactly.
        let c = 130.0;
        let length = 500.0;
        let diameter = 0.25;
        let q_ref = 0.05;
        let k = hazen_williams_k(c, length, diameter, q_ref).unwrap();
        let h_hw = 10.67 * length * q_ref.powf(1.852) / (c.powf(1.852) * diameter.powf(4.87));
        assert!((k * q_ref * q_ref - h_hw).abs() < 1e-12 * h_hw.max(1.0));
    }

    #[test]
    fn hazen_williams_k_rejects_nonpositive() {
        assert!(hazen_williams_k(0.0, 1.0, 0.1, 0.01).is_err());
        assert!(hazen_williams_k(130.0, 0.0, 0.1, 0.01).is_err());
        assert!(hazen_williams_k(130.0, 1.0, 0.0, 0.01).is_err());
        assert!(hazen_williams_k(130.0, 1.0, 0.1, 0.0).is_err());
    }
}
