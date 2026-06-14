//! The two-parameter Weibull lifetime model.
//!
//! The Weibull distribution is the workhorse of reliability engineering
//! because its single *shape* parameter `beta` lets one family describe
//! all three regions of the bathtub curve:
//!
//! - `beta < 1` — a *decreasing* hazard (infant mortality / burn-in),
//! - `beta = 1` — a *constant* hazard; the model degenerates exactly to
//!   the [exponential](crate::exponential) with `lambda = 1 / eta`,
//! - `beta > 1` — an *increasing* hazard (wear-out); `beta = 2` is the
//!   Rayleigh distribution (linearly rising hazard).
//!
//! With shape `beta > 0` and scale (characteristic life) `eta > 0`:
//!
//! - reliability   `R(t) = exp(-(t / eta)^beta)`,
//! - unreliability `F(t) = 1 - exp(-(t / eta)^beta)`,
//! - density       `f(t) = (beta / eta) (t / eta)^(beta - 1) exp(-(t / eta)^beta)`,
//! - hazard        `h(t) = (beta / eta) (t / eta)^(beta - 1)`,
//! - mean life     `MTTF = eta * Gamma(1 + 1 / beta)`.
//!
//! The scale `eta` is the *characteristic life*: `R(eta) = exp(-1) = 1/e`
//! for every shape, so about 63.2% of items have failed by `t = eta`
//! regardless of `beta`.

use serde::{Deserialize, Serialize};

use crate::error::{require_finite, require_time, ReliabilityError};

/// A two-parameter Weibull lifetime model.
///
/// Construct one with [`Weibull::new`]. Both parameters are validated to
/// be strictly positive at construction.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Weibull {
    /// Shape parameter `beta` (dimensionless), `> 0`.
    shape: f64,
    /// Scale / characteristic life `eta` (time units), `> 0`.
    scale: f64,
}

impl Weibull {
    /// Build a Weibull model from a shape `beta` and scale `eta`.
    ///
    /// # Errors
    ///
    /// Returns [`ReliabilityError::NonPositiveShape`] if `beta <= 0`,
    /// [`ReliabilityError::NonPositiveScale`] if `eta <= 0`, or
    /// [`ReliabilityError::NotFinite`] if either argument is not finite.
    pub fn new(shape: f64, scale: f64) -> Result<Self, ReliabilityError> {
        let shape = require_finite("beta", shape)?;
        let scale = require_finite("eta", scale)?;
        if shape <= 0.0 {
            return Err(ReliabilityError::NonPositiveShape { value: shape });
        }
        if scale <= 0.0 {
            return Err(ReliabilityError::NonPositiveScale { value: scale });
        }
        Ok(Self { shape, scale })
    }

    /// The shape parameter `beta`.
    #[must_use]
    pub fn shape(&self) -> f64 {
        self.shape
    }

    /// The scale parameter (characteristic life) `eta`.
    #[must_use]
    pub fn scale(&self) -> f64 {
        self.scale
    }

    /// The reliability (survival probability)
    /// `R(t) = exp(-(t / eta)^beta)`.
    ///
    /// `R(0) = 1` and `R(t) -> 0` as `t -> infinity`.
    ///
    /// # Errors
    ///
    /// Returns [`ReliabilityError::NegativeTime`] if `t < 0`.
    pub fn reliability(&self, t: f64) -> Result<f64, ReliabilityError> {
        let t = require_time(t)?;
        Ok((-(t / self.scale).powf(self.shape)).exp())
    }

    /// The unreliability / cumulative failure probability
    /// `F(t) = 1 - R(t)`.
    ///
    /// # Errors
    ///
    /// Returns [`ReliabilityError::NegativeTime`] if `t < 0`.
    pub fn unreliability(&self, t: f64) -> Result<f64, ReliabilityError> {
        Ok(1.0 - self.reliability(t)?)
    }

    /// The instantaneous hazard rate
    /// `h(t) = (beta / eta) (t / eta)^(beta - 1)`.
    ///
    /// At `t = 0` this is `0` for `beta > 1`, the constant `1 / eta` for
    /// `beta = 1`, and diverges to `+infinity` for `beta < 1`.
    ///
    /// # Errors
    ///
    /// Returns [`ReliabilityError::NegativeTime`] if `t < 0`.
    pub fn hazard(&self, t: f64) -> Result<f64, ReliabilityError> {
        let t = require_time(t)?;
        Ok((self.shape / self.scale) * (t / self.scale).powf(self.shape - 1.0))
    }

    /// The probability-density function
    /// `f(t) = h(t) R(t)`.
    ///
    /// # Errors
    ///
    /// Returns [`ReliabilityError::NegativeTime`] if `t < 0`.
    pub fn pdf(&self, t: f64) -> Result<f64, ReliabilityError> {
        let t = require_time(t)?;
        let h = (self.shape / self.scale) * (t / self.scale).powf(self.shape - 1.0);
        let r = (-(t / self.scale).powf(self.shape)).exp();
        Ok(h * r)
    }

    /// The mean time to failure `MTTF = eta * Gamma(1 + 1 / beta)`.
    ///
    /// For `beta = 1` this reduces to `eta` (matching the exponential
    /// `MTBF = 1 / lambda = eta`). The gamma factor is evaluated with a
    /// Lanczos approximation (see [`gamma`]).
    #[must_use]
    pub fn mttf(&self) -> f64 {
        self.scale * gamma(1.0 + 1.0 / self.shape)
    }

    /// The time `t` at which the reliability drops to `target`, i.e. the
    /// quantile solving `R(t) = target`.
    ///
    /// Inverting `R` gives `t = eta * (-ln(target))^(1 / beta)`. The
    /// special case `target = exp(-1)` returns `eta` exactly (the
    /// characteristic life).
    ///
    /// # Errors
    ///
    /// Returns [`ReliabilityError::ProbabilityOutOfRange`] if `target`
    /// is not in `[0, 1]`. A `target` of `0.0` yields `+infinity`.
    pub fn time_for_reliability(&self, target: f64) -> Result<f64, ReliabilityError> {
        let target = crate::error::require_probability(target)?;
        Ok(self.scale * (-target.ln()).powf(1.0 / self.shape))
    }
}

/// The gamma function `Gamma(x)` for real `x > 0`, via the Lanczos
/// approximation (g = 7, nine coefficients).
///
/// Reliability uses this only to evaluate `Gamma(1 + 1 / beta)` for the
/// Weibull mean, where the argument lies in `(1, 2]` for `beta >= 1`;
/// the implementation is nonetheless valid for any positive real `x`
/// (the reflection formula extends it below `0.5`). Relative error is
/// better than `1e-13` across the range that matters here.
///
/// Reference values it must reproduce: `Gamma(1) = 1`, `Gamma(2) = 1`,
/// `Gamma(n) = (n - 1)!`, and `Gamma(1/2) = sqrt(pi)`.
#[must_use]
pub fn gamma(x: f64) -> f64 {
    // Lanczos coefficients for g = 7.
    const G: f64 = 7.0;
    const COEFFS: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_1,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];

    if x < 0.5 {
        // Reflection: Gamma(x) Gamma(1 - x) = pi / sin(pi x).
        std::f64::consts::PI / ((std::f64::consts::PI * x).sin() * gamma(1.0 - x))
    } else {
        let x = x - 1.0;
        let mut a = COEFFS[0];
        let t = x + G + 0.5;
        for (i, &c) in COEFFS.iter().enumerate().skip(1) {
            a += c / (x + i as f64);
        }
        (2.0 * std::f64::consts::PI).sqrt() * t.powf(x + 0.5) * (-t).exp() * a
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::exponential::Exponential;

    /// Tolerance for floating-point ground-truth comparisons.
    const EPS: f64 = 1e-10;

    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn constructor_validates_both_parameters() {
        assert!(matches!(
            Weibull::new(0.0, 1.0),
            Err(ReliabilityError::NonPositiveShape { .. })
        ));
        assert!(matches!(
            Weibull::new(-1.0, 1.0),
            Err(ReliabilityError::NonPositiveShape { .. })
        ));
        assert!(matches!(
            Weibull::new(1.0, 0.0),
            Err(ReliabilityError::NonPositiveScale { .. })
        ));
        assert!(matches!(
            Weibull::new(f64::NAN, 1.0),
            Err(ReliabilityError::NotFinite { .. })
        ));
        assert!(Weibull::new(1.5, 2.0).is_ok());
    }

    #[test]
    fn reliability_at_zero_is_one() {
        // GROUND TRUTH: R(0) = exp(0) = 1 (for beta >= 1; for beta < 1
        // the limit at 0 is still 1).
        for &(b, e) in &[(0.5, 1.0), (1.0, 3.0), (2.0, 5.0), (3.5, 0.5)] {
            let w = Weibull::new(b, e).unwrap();
            let r0 = w.reliability(0.0).unwrap();
            assert!(close(r0, 1.0, EPS), "R(0) = {r0} for beta = {b}");
        }
    }

    #[test]
    fn reliability_tends_to_zero_at_large_time() {
        // GROUND TRUTH: R(t) -> 0 as t -> infinity.
        let w = Weibull::new(2.0, 1.0).unwrap();
        let r_far = w.reliability(50.0).unwrap();
        assert!(r_far < EPS, "R(50) = {r_far} should be ~0");
        // Deep in the tail exp(-(t/eta)^beta) underflows to exactly 0.0:
        // (1000/1)^2 = 1e6, and exp(-1e6) == 0.0 in f64.
        let r_underflow = w.reliability(1000.0).unwrap();
        assert_eq!(r_underflow, 0.0);
    }

    #[test]
    fn reliability_at_characteristic_life_is_one_over_e() {
        // GROUND TRUTH: R(eta) = exp(-1) = 1/e for EVERY shape beta.
        let one_over_e = (-1.0_f64).exp();
        for &b in &[0.5, 1.0, 1.7, 2.0, 4.0] {
            let eta = 12.5;
            let w = Weibull::new(b, eta).unwrap();
            let r = w.reliability(eta).unwrap();
            assert!(
                close(r, one_over_e, EPS),
                "R(eta) = {r} expected 1/e = {one_over_e} (beta = {b})"
            );
        }
    }

    #[test]
    fn reliability_matches_closed_form() {
        // GROUND TRUTH: R(t) = exp(-(t/eta)^beta). beta=2, eta=2, t=2
        //   -> exp(-(1)^2) = exp(-1).
        let w = Weibull::new(2.0, 2.0).unwrap();
        assert!(close(w.reliability(2.0).unwrap(), (-1.0_f64).exp(), EPS));
        // beta=2, eta=2, t=4 -> exp(-(2)^2) = exp(-4).
        assert!(close(w.reliability(4.0).unwrap(), (-4.0_f64).exp(), EPS));
    }

    #[test]
    fn unreliability_complements_reliability() {
        let w = Weibull::new(1.3, 4.0).unwrap();
        for k in 0..=20 {
            let t = k as f64 * 0.5;
            let r = w.reliability(t).unwrap();
            let f = w.unreliability(t).unwrap();
            assert!(close(r + f, 1.0, EPS), "R + F = {} at t = {t}", r + f);
        }
    }

    #[test]
    fn beta_one_reduces_to_exponential() {
        // GROUND TRUTH (validation requirement): a Weibull with beta = 1
        // is exactly an exponential with lambda = 1/eta. R, F, pdf and
        // hazard must all coincide pointwise.
        let eta = 7.0;
        let w = Weibull::new(1.0, eta).unwrap();
        let e = Exponential::new(1.0 / eta).unwrap();
        for k in 0..=30 {
            let t = k as f64 * 0.4;
            assert!(
                close(w.reliability(t).unwrap(), e.reliability(t).unwrap(), EPS),
                "R mismatch at t = {t}"
            );
            assert!(
                close(w.pdf(t).unwrap(), e.pdf(t).unwrap(), EPS),
                "pdf mismatch at t = {t}"
            );
            assert!(
                close(w.hazard(t).unwrap(), e.hazard(t).unwrap(), EPS),
                "hazard mismatch at t = {t}"
            );
        }
        // And the means agree: MTTF = eta * Gamma(2) = eta = MTBF.
        assert!(close(w.mttf(), e.mtbf(), EPS), "mean mismatch");
    }

    #[test]
    fn hazard_is_decreasing_increasing_or_constant_by_shape() {
        // GROUND TRUTH: sign of d h / d t follows beta vs 1.
        let t1 = 1.0;
        let t2 = 3.0;
        // beta < 1: decreasing hazard.
        let burn_in = Weibull::new(0.5, 2.0).unwrap();
        assert!(burn_in.hazard(t1).unwrap() > burn_in.hazard(t2).unwrap());
        // beta = 1: constant hazard = 1/eta.
        let flat = Weibull::new(1.0, 2.0).unwrap();
        assert!(close(flat.hazard(t1).unwrap(), 1.0 / 2.0, EPS));
        assert!(close(flat.hazard(t2).unwrap(), 1.0 / 2.0, EPS));
        // beta > 1: increasing (wear-out) hazard.
        let wear = Weibull::new(2.5, 2.0).unwrap();
        assert!(wear.hazard(t1).unwrap() < wear.hazard(t2).unwrap());
    }

    #[test]
    fn pdf_equals_hazard_times_reliability() {
        // GROUND TRUTH: f(t) = h(t) * R(t).
        let w = Weibull::new(1.8, 3.0).unwrap();
        for k in 1..=20 {
            let t = k as f64 * 0.3;
            let expected = w.hazard(t).unwrap() * w.reliability(t).unwrap();
            assert!(
                close(w.pdf(t).unwrap(), expected, EPS),
                "pdf mismatch at t = {t}"
            );
        }
    }

    #[test]
    fn mttf_matches_known_values() {
        // GROUND TRUTH: beta=2 (Rayleigh) -> MTTF = eta*Gamma(1.5)
        //   = eta * sqrt(pi)/2.
        let eta = 10.0;
        let w = Weibull::new(2.0, eta).unwrap();
        let expected = eta * std::f64::consts::PI.sqrt() / 2.0;
        assert!(
            close(w.mttf(), expected, 1e-9),
            "MTTF = {} expected {expected}",
            w.mttf()
        );
    }

    #[test]
    fn time_for_reliability_inverts_reliability() {
        let w = Weibull::new(2.3, 6.0).unwrap();
        for &target in &[0.99, 0.5, 0.1, 0.02] {
            let t = w.time_for_reliability(target).unwrap();
            let back = w.reliability(t).unwrap();
            assert!(close(back, target, EPS), "R(t({target})) = {back}");
        }
    }

    #[test]
    fn characteristic_life_quantile_returns_eta() {
        // GROUND TRUTH: time for R = 1/e is exactly eta.
        let eta = 9.0;
        let w = Weibull::new(3.0, eta).unwrap();
        let t = w.time_for_reliability((-1.0_f64).exp()).unwrap();
        assert!(close(t, eta, EPS), "quantile(1/e) = {t} expected {eta}");
    }

    #[test]
    fn negative_time_rejected() {
        let w = Weibull::new(2.0, 1.0).unwrap();
        assert!(matches!(
            w.reliability(-1.0),
            Err(ReliabilityError::NegativeTime { .. })
        ));
        assert!(matches!(
            w.hazard(-0.1),
            Err(ReliabilityError::NegativeTime { .. })
        ));
    }

    #[test]
    fn gamma_reproduces_factorials() {
        // GROUND TRUTH: Gamma(n) = (n-1)!.
        assert!(close(gamma(1.0), 1.0, 1e-12), "Gamma(1) = {}", gamma(1.0));
        assert!(close(gamma(2.0), 1.0, 1e-12), "Gamma(2) = {}", gamma(2.0));
        assert!(close(gamma(3.0), 2.0, 1e-12), "Gamma(3) = {}", gamma(3.0));
        assert!(close(gamma(4.0), 6.0, 1e-11), "Gamma(4) = {}", gamma(4.0));
        assert!(close(gamma(5.0), 24.0, 1e-10), "Gamma(5) = {}", gamma(5.0));
        assert!(close(gamma(6.0), 120.0, 1e-9), "Gamma(6) = {}", gamma(6.0));
    }

    #[test]
    fn gamma_half_is_sqrt_pi() {
        // GROUND TRUTH: Gamma(1/2) = sqrt(pi).
        assert!(
            close(gamma(0.5), std::f64::consts::PI.sqrt(), 1e-12),
            "Gamma(0.5) = {}",
            gamma(0.5)
        );
        // Gamma(1.5) = sqrt(pi)/2.
        assert!(
            close(gamma(1.5), std::f64::consts::PI.sqrt() / 2.0, 1e-12),
            "Gamma(1.5) = {}",
            gamma(1.5)
        );
    }

    #[test]
    fn gamma_recurrence_holds() {
        // GROUND TRUTH: Gamma(x+1) = x * Gamma(x).
        for &x in &[0.3, 0.75, 1.2, 2.6, 4.1] {
            assert!(
                close(
                    gamma(x + 1.0),
                    x * gamma(x),
                    1e-9 * gamma(x + 1.0).abs().max(1.0)
                ),
                "recurrence fails at x = {x}"
            );
        }
    }

    #[test]
    fn serde_round_trips() {
        let w = Weibull::new(1.6, 4.2).unwrap();
        let json = serde_json::to_string(&w).unwrap();
        let back: Weibull = serde_json::from_str(&json).unwrap();
        assert_eq!(w, back);
    }
}
