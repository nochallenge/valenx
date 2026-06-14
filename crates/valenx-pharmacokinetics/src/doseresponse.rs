//! Emax / Hill sigmoidal dose-response (pharmacodynamics).
//!
//! Where the one- and two-compartment models describe *how much* drug is
//! present over time, the dose-response model describes *what it does*.
//! The sigmoidal Emax (Hill) equation relates the effect `E` to the dose
//! (or concentration) `d`:
//!
//! `E(d) = Emax · d^n / (EC50^n + d^n)`,
//!
//! where `Emax` is the maximal attainable effect (the asymptote as
//! `d → ∞`), `EC50` is the dose producing half of `Emax`, and `n` is the
//! Hill coefficient (steepness of the curve; `n = 1` is the classic
//! hyperbolic / Michaelis-Menten shape, `n > 1` is sigmoidal).
//!
//! Two exact properties anchor the model and are exercised by the tests:
//! at the half-maximal dose the effect is exactly half the maximum,
//! `E(EC50) = Emax/2`, and for `n > 0` the effect is a strictly
//! increasing function of dose (a monotone agonist response).

use crate::error::{require_non_negative, require_positive, Result};
use serde::{Deserialize, Serialize};

/// An Emax / Hill sigmoidal dose-response model.
///
/// Construct with [`DoseResponse::new`], which validates that `emax`,
/// `ec50`, and the Hill coefficient `n` are strictly positive.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct DoseResponse {
    emax: f64,
    ec50: f64,
    n: f64,
}

impl DoseResponse {
    /// Build the model `E(d) = Emax·d^n / (EC50^n + d^n)`.
    ///
    /// `emax` is the maximal effect, `ec50` the half-maximal dose, and
    /// `n` the Hill coefficient — all strictly positive.
    ///
    /// # Errors
    ///
    /// Returns [`PkError::NotPositive`](crate::PkError::NotPositive) if
    /// any of `emax`, `ec50`, or `n` is not strictly positive (or
    /// non-finite).
    pub fn new(emax: f64, ec50: f64, n: f64) -> Result<Self> {
        let emax = require_positive("emax", emax)?;
        let ec50 = require_positive("ec50", ec50)?;
        let n = require_positive("n", n)?;
        Ok(Self { emax, ec50, n })
    }

    /// The maximal attainable effect `Emax` (the asymptote as `d → ∞`).
    #[inline]
    pub fn emax(&self) -> f64 {
        self.emax
    }

    /// The half-maximal dose `EC50`.
    #[inline]
    pub fn ec50(&self) -> f64 {
        self.ec50
    }

    /// The Hill coefficient `n` (curve steepness).
    #[inline]
    pub fn hill(&self) -> f64 {
        self.n
    }

    /// Effect `E(d) = Emax·d^n / (EC50^n + d^n)` at dose `d >= 0`.
    ///
    /// At `d = 0` the effect is `0`; as `d → ∞` it approaches `Emax`.
    ///
    /// # Errors
    ///
    /// Returns [`PkError::Negative`](crate::PkError::Negative) if `d` is
    /// negative or non-finite.
    pub fn effect(&self, d: f64) -> Result<f64> {
        let d = require_non_negative("d", d)?;
        let dn = d.powf(self.n);
        let ec50n = self.ec50.powf(self.n);
        Ok(self.emax * dn / (ec50n + dn))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-12;

    #[test]
    fn rejects_bad_parameters() {
        assert!(DoseResponse::new(0.0, 5.0, 1.0).is_err());
        assert!(DoseResponse::new(100.0, 0.0, 1.0).is_err());
        assert!(DoseResponse::new(100.0, 5.0, 0.0).is_err());
        assert!(DoseResponse::new(-100.0, 5.0, 1.0).is_err());
        assert!(DoseResponse::new(100.0, 5.0, f64::NAN).is_err());
    }

    #[test]
    fn effect_at_zero_is_zero() {
        let m = DoseResponse::new(100.0, 5.0, 1.0).unwrap();
        let e = m.effect(0.0).unwrap();
        assert!(e.abs() < TOL, "E(0) = {e}");
    }

    #[test]
    fn effect_at_ec50_is_half_emax() {
        // ANALYTIC: E(EC50) = Emax/2, for any Hill coefficient.
        for &n in &[0.5_f64, 1.0, 2.0, 4.0] {
            let m = DoseResponse::new(100.0, 5.0, n).unwrap();
            let e = m.effect(m.ec50()).unwrap();
            assert!(
                (e - 50.0).abs() < 1e-10,
                "E(EC50) = {e} for n = {n}, expected 50"
            );
        }
    }

    #[test]
    fn effect_is_strictly_increasing_in_dose() {
        // ANALYTIC: for n > 0 the agonist response is monotone increasing.
        let m = DoseResponse::new(100.0, 5.0, 2.0).unwrap();
        let mut prev = -1.0;
        for i in 0..=50 {
            let d = i as f64 * 0.5;
            let e = m.effect(d).unwrap();
            assert!(e > prev, "E should increase: E({d}) = {e}, prev = {prev}");
            prev = e;
        }
    }

    #[test]
    fn effect_approaches_emax_for_large_dose() {
        let m = DoseResponse::new(100.0, 5.0, 1.0).unwrap();
        let e = m.effect(1.0e6).unwrap();
        assert!((e - 100.0).abs() < 1e-3, "E(large) = {e}");
        // Never exceeds Emax.
        assert!(e < 100.0, "E must stay below Emax, got {e}");
    }

    #[test]
    fn hyperbolic_case_matches_closed_form() {
        // n = 1 reduces to E = Emax·d / (EC50 + d); check d = EC50/3.
        let m = DoseResponse::new(80.0, 6.0, 1.0).unwrap();
        let d = 2.0; // EC50/3.
        let expected = 80.0 * d / (6.0 + d);
        let e = m.effect(d).unwrap();
        assert!((e - expected).abs() < TOL, "E = {e}, expected {expected}");
    }

    #[test]
    fn sigmoidicity_steepens_with_hill_coefficient() {
        // Below EC50 a steeper Hill coefficient gives a smaller fractional
        // effect; above EC50 a larger one. Compare n = 1 vs n = 4 at
        // d = EC50/2 (below) and d = 2·EC50 (above).
        let shallow = DoseResponse::new(100.0, 5.0, 1.0).unwrap();
        let steep = DoseResponse::new(100.0, 5.0, 4.0).unwrap();
        let below = 2.5;
        let above = 10.0;
        assert!(
            steep.effect(below).unwrap() < shallow.effect(below).unwrap(),
            "steep curve should be lower below EC50"
        );
        assert!(
            steep.effect(above).unwrap() > shallow.effect(above).unwrap(),
            "steep curve should be higher above EC50"
        );
    }

    #[test]
    fn rejects_negative_dose() {
        let m = DoseResponse::new(100.0, 5.0, 1.0).unwrap();
        assert!(m.effect(-1.0).is_err());
        assert!(m.effect(f64::INFINITY).is_err());
    }
}
