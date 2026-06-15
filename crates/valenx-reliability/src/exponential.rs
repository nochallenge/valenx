//! The exponential (constant-hazard) lifetime model.
//!
//! The exponential distribution is the canonical model for a component
//! whose instantaneous failure rate (hazard) is *constant* in time — the
//! "useful life", flat-bottom region of the bathtub curve, where failures
//! arrive as a memoryless Poisson process. With a constant failure rate
//! `lambda > 0` (failures per unit time):
//!
//! - reliability (survival)   `R(t) = exp(-lambda t)`,
//! - unreliability (CDF)      `F(t) = 1 - exp(-lambda t)`,
//! - probability density      `f(t) = lambda exp(-lambda t)`,
//! - hazard rate              `h(t) = lambda`  (constant),
//! - mean time between/to failure  `MTBF = 1 / lambda`.
//!
//! A defining property is *memorylessness*: a used component that has
//! survived to time `s` is statistically as good as new — its remaining
//! life has the same exponential distribution it started with.

use serde::{Deserialize, Serialize};

use crate::error::{require_finite, require_time, ReliabilityError};

/// An exponential lifetime model parameterised by its constant failure
/// rate `lambda`.
///
/// Construct one with [`Exponential::new`] (from a failure rate) or
/// [`Exponential::from_mtbf`] (from a mean time between failures). The
/// invariant `lambda > 0` is enforced at construction, so every accessor
/// is total (infallible) except those taking a time argument, which
/// validate `t >= 0`.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Exponential {
    /// Constant failure rate `lambda` (failures per unit time), `> 0`.
    lambda: f64,
}

impl Exponential {
    /// Build a model from a constant failure rate `lambda`.
    ///
    /// # Errors
    ///
    /// Returns [`ReliabilityError::NonPositiveRate`] if `lambda <= 0`,
    /// or [`ReliabilityError::NotFinite`] if `lambda` is `NaN` or
    /// infinite.
    pub fn new(lambda: f64) -> Result<Self, ReliabilityError> {
        let lambda = require_finite("lambda", lambda)?;
        if lambda <= 0.0 {
            return Err(ReliabilityError::NonPositiveRate { value: lambda });
        }
        Ok(Self { lambda })
    }

    /// Build a model from a mean time between failures `mtbf = 1 / lambda`.
    ///
    /// # Errors
    ///
    /// Returns [`ReliabilityError::NonPositiveRate`] if `mtbf <= 0`
    /// (the implied rate would be non-positive), or
    /// [`ReliabilityError::NotFinite`] if `mtbf` is not finite.
    pub fn from_mtbf(mtbf: f64) -> Result<Self, ReliabilityError> {
        let mtbf = require_finite("mtbf", mtbf)?;
        if mtbf <= 0.0 {
            // Echo the implied rate so the error talks about lambda.
            return Err(ReliabilityError::NonPositiveRate { value: 1.0 / mtbf });
        }
        Self::new(1.0 / mtbf)
    }

    /// The constant failure rate `lambda` (failures per unit time).
    #[must_use]
    pub fn rate(&self) -> f64 {
        self.lambda
    }

    /// The mean time between failures, `MTBF = 1 / lambda`.
    #[must_use]
    pub fn mtbf(&self) -> f64 {
        1.0 / self.lambda
    }

    /// The instantaneous hazard (failure) rate at time `t`.
    ///
    /// For the exponential model this is the constant `lambda` for all
    /// `t >= 0`.
    ///
    /// # Errors
    ///
    /// Returns [`ReliabilityError::NegativeTime`] if `t < 0` (or
    /// [`ReliabilityError::NotFinite`] if `t` is not finite).
    pub fn hazard(&self, t: f64) -> Result<f64, ReliabilityError> {
        require_time(t)?;
        Ok(self.lambda)
    }

    /// The reliability (survival probability) `R(t) = exp(-lambda t)`.
    ///
    /// `R(0) = 1` and `R(t) -> 0` as `t -> infinity`.
    ///
    /// # Errors
    ///
    /// Returns [`ReliabilityError::NegativeTime`] if `t < 0`.
    pub fn reliability(&self, t: f64) -> Result<f64, ReliabilityError> {
        let t = require_time(t)?;
        Ok((-self.lambda * t).exp())
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

    /// The probability-density function `f(t) = lambda exp(-lambda t)`.
    ///
    /// # Errors
    ///
    /// Returns [`ReliabilityError::NegativeTime`] if `t < 0`.
    pub fn pdf(&self, t: f64) -> Result<f64, ReliabilityError> {
        let t = require_time(t)?;
        Ok(self.lambda * (-self.lambda * t).exp())
    }

    /// The time `t` at which the reliability first drops to `target`,
    /// i.e. the quantile solving `R(t) = target`.
    ///
    /// Because `R` is strictly decreasing, this is the inverse
    /// `t = -ln(target) / lambda`. For example, the *median* life is
    /// `quantile(0.5)`.
    ///
    /// # Errors
    ///
    /// Returns [`ReliabilityError::ProbabilityOutOfRange`] if `target`
    /// is not in `[0, 1]`. A `target` of `0.0` yields `+infinity`
    /// (reliability never truly reaches zero).
    pub fn time_for_reliability(&self, target: f64) -> Result<f64, ReliabilityError> {
        let target = crate::error::require_probability(target)?;
        Ok(-target.ln() / self.lambda)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for floating-point ground-truth comparisons.
    const EPS: f64 = 1e-12;

    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn constructor_rejects_non_positive_and_non_finite_rate() {
        assert!(matches!(
            Exponential::new(0.0),
            Err(ReliabilityError::NonPositiveRate { .. })
        ));
        assert!(matches!(
            Exponential::new(-1.0),
            Err(ReliabilityError::NonPositiveRate { .. })
        ));
        assert!(matches!(
            Exponential::new(f64::NAN),
            Err(ReliabilityError::NotFinite { .. })
        ));
        assert!(Exponential::new(0.001).is_ok());
    }

    #[test]
    fn reliability_at_zero_is_one() {
        // GROUND TRUTH: R(0) = exp(0) = 1 for every lambda.
        for &lambda in &[1e-6, 0.5, 1.0, 42.0] {
            let m = Exponential::new(lambda).unwrap();
            let r0 = m.reliability(0.0).unwrap();
            assert!(close(r0, 1.0, EPS), "R(0) = {r0} for lambda = {lambda}");
        }
    }

    #[test]
    fn reliability_tends_to_zero_at_large_time() {
        // GROUND TRUTH: R(t) -> 0 as t -> infinity.
        let m = Exponential::new(2.0).unwrap();
        let r_far = m.reliability(1e6).unwrap();
        assert!(r_far < EPS, "R(1e6) = {r_far} should be ~0");
        // Deep in the tail, exp(-lambda t) underflows to exactly 0.0:
        // lambda * t = 2 * 1000 = 2000, and exp(-2000) == 0.0 in f64.
        let r_underflow = m.reliability(1000.0).unwrap();
        assert_eq!(r_underflow, 0.0);
    }

    #[test]
    fn reliability_is_strictly_decreasing() {
        let m = Exponential::new(0.3).unwrap();
        let mut prev = m.reliability(0.0).unwrap();
        for k in 1..=50 {
            let t = k as f64 * 0.5;
            let r = m.reliability(t).unwrap();
            assert!(
                r < prev,
                "R must strictly decrease: R({t}) = {r}, prev = {prev}"
            );
            prev = r;
        }
    }

    #[test]
    fn mtbf_is_reciprocal_of_rate() {
        // GROUND TRUTH: MTBF = 1 / lambda.
        for &lambda in &[0.01, 0.25, 1.0, 7.5] {
            let m = Exponential::new(lambda).unwrap();
            assert!(
                close(m.mtbf(), 1.0 / lambda, EPS),
                "MTBF = {} expected {}",
                m.mtbf(),
                1.0 / lambda
            );
        }
    }

    #[test]
    fn from_mtbf_round_trips_rate() {
        let m = Exponential::from_mtbf(500.0).unwrap();
        assert!(close(m.rate(), 1.0 / 500.0, EPS), "rate = {}", m.rate());
        assert!(close(m.mtbf(), 500.0, EPS), "mtbf = {}", m.mtbf());
    }

    #[test]
    fn from_mtbf_rejects_non_positive() {
        assert!(matches!(
            Exponential::from_mtbf(0.0),
            Err(ReliabilityError::NonPositiveRate { .. })
        ));
        assert!(matches!(
            Exponential::from_mtbf(-10.0),
            Err(ReliabilityError::NonPositiveRate { .. })
        ));
    }

    #[test]
    fn reliability_at_mtbf_is_one_over_e() {
        // GROUND TRUTH: R(MTBF) = R(1/lambda) = exp(-1) = 1/e.
        // ~63.2% of items have failed by the MTBF for the exponential.
        let one_over_e = (-1.0_f64).exp();
        for &lambda in &[0.05, 1.0, 13.0] {
            let m = Exponential::new(lambda).unwrap();
            let r = m.reliability(m.mtbf()).unwrap();
            assert!(
                close(r, one_over_e, EPS),
                "R(MTBF) = {r} expected 1/e = {one_over_e} (lambda = {lambda})"
            );
        }
    }

    #[test]
    fn reliability_matches_closed_form() {
        // GROUND TRUTH: spot values of exp(-lambda t).
        let m = Exponential::new(0.1).unwrap();
        // lambda*t = 0.1*10 = 1 -> exp(-1).
        assert!(close(m.reliability(10.0).unwrap(), (-1.0_f64).exp(), EPS));
        // lambda*t = 0.1*20 = 2 -> exp(-2).
        assert!(close(m.reliability(20.0).unwrap(), (-2.0_f64).exp(), EPS));
    }

    #[test]
    fn unreliability_complements_reliability() {
        let m = Exponential::new(0.7).unwrap();
        for k in 0..=20 {
            let t = k as f64 * 0.3;
            let r = m.reliability(t).unwrap();
            let f = m.unreliability(t).unwrap();
            assert!(close(r + f, 1.0, EPS), "R + F = {} at t = {t}", r + f);
        }
    }

    #[test]
    fn hazard_is_constant_and_equals_lambda() {
        // GROUND TRUTH: h(t) = lambda for all t (defining property).
        let lambda = 0.42;
        let m = Exponential::new(lambda).unwrap();
        for k in 0..=10 {
            let t = k as f64;
            assert!(close(m.hazard(t).unwrap(), lambda, EPS), "h({t}) != lambda");
        }
    }

    #[test]
    fn pdf_equals_lambda_times_reliability() {
        // GROUND TRUTH: f(t) = lambda * R(t) = h(t) * R(t).
        let m = Exponential::new(0.9).unwrap();
        for k in 0..=15 {
            let t = k as f64 * 0.4;
            let expected = m.rate() * m.reliability(t).unwrap();
            assert!(close(m.pdf(t).unwrap(), expected, EPS), "f({t}) mismatch");
        }
    }

    #[test]
    fn pdf_at_zero_is_lambda() {
        // GROUND TRUTH: f(0) = lambda * exp(0) = lambda.
        let m = Exponential::new(3.3).unwrap();
        assert!(
            close(m.pdf(0.0).unwrap(), 3.3, EPS),
            "f(0) = {}",
            m.pdf(0.0).unwrap()
        );
    }

    #[test]
    fn memorylessness_holds() {
        // GROUND TRUTH: P(T > s + t | T > s) = R(s + t)/R(s) = R(t).
        let m = Exponential::new(0.6).unwrap();
        let s = 4.0;
        let t = 2.5;
        let conditional = m.reliability(s + t).unwrap() / m.reliability(s).unwrap();
        assert!(
            close(conditional, m.reliability(t).unwrap(), EPS),
            "memoryless property violated: {conditional}"
        );
    }

    #[test]
    fn time_for_reliability_inverts_reliability() {
        // GROUND TRUTH: t = -ln(target)/lambda, and R(t(target)) = target.
        let m = Exponential::new(0.2).unwrap();
        for &target in &[0.9, 0.5, 0.1, 0.01] {
            let t = m.time_for_reliability(target).unwrap();
            let back = m.reliability(t).unwrap();
            assert!(close(back, target, EPS), "R(t({target})) = {back}");
        }
    }

    #[test]
    fn median_life_is_ln2_over_lambda() {
        // GROUND TRUTH: median solves R(t)=0.5 -> t = ln(2)/lambda.
        let lambda = 0.05;
        let m = Exponential::new(lambda).unwrap();
        let median = m.time_for_reliability(0.5).unwrap();
        assert!(
            close(median, 2.0_f64.ln() / lambda, EPS),
            "median = {median} expected ln2/lambda = {}",
            2.0_f64.ln() / lambda
        );
        // Median is shorter than the mean (MTBF) for the exponential.
        assert!(median < m.mtbf());
    }

    #[test]
    fn time_for_reliability_target_zero_is_infinite() {
        let m = Exponential::new(1.0).unwrap();
        assert_eq!(m.time_for_reliability(0.0).unwrap(), f64::INFINITY);
    }

    #[test]
    fn negative_time_is_rejected_everywhere() {
        let m = Exponential::new(1.0).unwrap();
        assert!(matches!(
            m.reliability(-1.0),
            Err(ReliabilityError::NegativeTime { .. })
        ));
        assert!(matches!(
            m.pdf(-0.001),
            Err(ReliabilityError::NegativeTime { .. })
        ));
        assert!(matches!(
            m.hazard(-5.0),
            Err(ReliabilityError::NegativeTime { .. })
        ));
    }

    #[test]
    fn serde_round_trips() {
        let m = Exponential::new(0.123).unwrap();
        let json = serde_json::to_string(&m).unwrap();
        let back: Exponential = serde_json::from_str(&json).unwrap();
        assert_eq!(m, back);
    }
}
