//! Two-compartment biexponential disposition.
//!
//! After an IV bolus the plasma concentration of many drugs is described
//! not by a single exponential but by the sum of two: a fast
//! *distribution* phase (the drug spreading from the central/plasma
//! compartment into a peripheral tissue compartment) and a slow
//! *elimination* phase. In terms of the so-called **macro constants** the
//! curve is
//!
//! `C(t) = A·exp(-α·t) + B·exp(-β·t)`,
//!
//! with `α > β > 0` by convention (`α` the distribution rate, `β` the
//! terminal rate). This module models the curve directly from its four
//! macro constants — `A`, `α`, `B`, `β` — which is exactly the form one
//! fits to plasma data, and does not attempt to re-derive the underlying
//! micro-rate constants (`k10`, `k12`, `k21`) or the inter-compartmental
//! volumes.
//!
//! The intercept `C(0) = A + B` is the (back-extrapolated) concentration
//! at the instant of the bolus, and the terminal half-life is
//! `ln(2)/β`.

use crate::error::{require_non_negative, require_positive, Result};
use serde::{Deserialize, Serialize};

/// A two-compartment biexponential disposition curve, parameterised by
/// its macro constants.
///
/// Construct with [`TwoCompartment::new`], which validates that the two
/// rate constants are strictly positive and the two coefficients are
/// non-negative.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct TwoCompartment {
    a: f64,
    alpha: f64,
    b: f64,
    beta: f64,
}

impl TwoCompartment {
    /// Build the curve `C(t) = A·exp(-α·t) + B·exp(-β·t)` from its macro
    /// constants.
    ///
    /// `a` and `b` are concentration coefficients (non-negative); `alpha`
    /// and `beta` are rate constants (strictly positive). No ordering of
    /// `alpha` vs `beta` is imposed — the caller may pass them in either
    /// order — but the conventional choice is `alpha > beta`.
    ///
    /// # Errors
    ///
    /// Returns [`PkError::NotPositive`](crate::PkError::NotPositive) if
    /// `alpha` or `beta` is not strictly positive (or non-finite), and
    /// [`PkError::Negative`](crate::PkError::Negative) if `a` or `b` is
    /// negative (or non-finite).
    pub fn new(a: f64, alpha: f64, b: f64, beta: f64) -> Result<Self> {
        let a = require_non_negative("a", a)?;
        let alpha = require_positive("alpha", alpha)?;
        let b = require_non_negative("b", b)?;
        let beta = require_positive("beta", beta)?;
        Ok(Self { a, alpha, b, beta })
    }

    /// The distribution-phase coefficient `A`.
    #[inline]
    pub fn a(&self) -> f64 {
        self.a
    }

    /// The distribution-phase rate constant `α`.
    #[inline]
    pub fn alpha(&self) -> f64 {
        self.alpha
    }

    /// The elimination-phase coefficient `B`.
    #[inline]
    pub fn b(&self) -> f64 {
        self.b
    }

    /// The elimination-phase rate constant `β`.
    #[inline]
    pub fn beta(&self) -> f64 {
        self.beta
    }

    /// Concentration `C(t) = A·exp(-α·t) + B·exp(-β·t)` at time `t >= 0`.
    ///
    /// # Errors
    ///
    /// Returns [`PkError::Negative`](crate::PkError::Negative) if `t` is
    /// negative or non-finite.
    pub fn concentration(&self, t: f64) -> Result<f64> {
        let t = require_non_negative("t", t)?;
        Ok(self.a * (-self.alpha * t).exp() + self.b * (-self.beta * t).exp())
    }

    /// Back-extrapolated intercept `C(0) = A + B`.
    #[inline]
    pub fn intercept(&self) -> f64 {
        self.a + self.b
    }

    /// Terminal half-life `ln(2) / β`, governed by the slow elimination
    /// phase.
    #[inline]
    pub fn terminal_half_life(&self) -> f64 {
        std::f64::consts::LN_2 / self.beta
    }

    /// Total exposure `AUC(0→∞) = A/α + B/β`, the integral of
    /// [`concentration`](Self::concentration) over `[0, ∞)`.
    #[inline]
    pub fn auc(&self) -> f64 {
        self.a / self.alpha + self.b / self.beta
    }

    /// Partial area under the concentration-time curve from time 0 to `t`,
    /// `AUC(0→t) = (A/α)(1 - exp(-α·t)) + (B/β)(1 - exp(-β·t))` — the
    /// cumulative exposure accrued by time `t`.
    ///
    /// Each exponential phase contributes its own saturating term, so this
    /// is the integral of [`concentration`](Self::concentration) over
    /// `[0, t]`: it is zero at `t = 0`, rises monotonically, and approaches
    /// the total [`auc`](Self::auc) (`= A/α + B/β`) as `t → ∞`.
    ///
    /// # Errors
    ///
    /// Returns [`PkError::Negative`](crate::PkError::Negative) if `t` is
    /// negative or non-finite.
    pub fn auc_to(&self, t: f64) -> Result<f64> {
        let t = require_non_negative("t", t)?;
        let dist = self.a / self.alpha * (1.0 - (-self.alpha * t).exp());
        let elim = self.b / self.beta * (1.0 - (-self.beta * t).exp());
        Ok(dist + elim)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-12;

    fn model() -> TwoCompartment {
        // A 8, α 1.0; B 2, β 0.1 → C(0) = 10.
        TwoCompartment::new(8.0, 1.0, 2.0, 0.1).unwrap()
    }

    #[test]
    fn rejects_bad_parameters() {
        assert!(TwoCompartment::new(8.0, 0.0, 2.0, 0.1).is_err());
        assert!(TwoCompartment::new(8.0, 1.0, 2.0, -0.1).is_err());
        assert!(TwoCompartment::new(-8.0, 1.0, 2.0, 0.1).is_err());
        assert!(TwoCompartment::new(8.0, 1.0, f64::NAN, 0.1).is_err());
        assert!(TwoCompartment::new(8.0, f64::INFINITY, 2.0, 0.1).is_err());
    }

    #[test]
    fn concentration_at_zero_is_a_plus_b() {
        // ANALYTIC: C(0) = A + B.
        let m = model();
        let c0 = m.concentration(0.0).unwrap();
        assert!((c0 - 10.0).abs() < TOL, "C(0) = {c0}");
        assert!((c0 - m.intercept()).abs() < TOL);
        assert!((c0 - (m.a() + m.b())).abs() < TOL);
    }

    #[test]
    fn concentration_matches_explicit_biexponential() {
        // Check a non-trivial time against a hand-evaluated value.
        // A=8, α=1, B=2, β=0.1 at t=2 → 8·exp(-2) + 2·exp(-0.2).
        let m = model();
        let t = 2.0;
        let expected = 8.0 * (-2.0_f64).exp() + 2.0 * (-0.2_f64).exp();
        let c = m.concentration(t).unwrap();
        assert!(
            (c - expected).abs() < TOL,
            "C(2) = {c}, expected {expected}"
        );
    }

    #[test]
    fn terminal_phase_dominates_at_late_times() {
        // At large t the fast (α) term is negligible, so C(t) ≈ B·exp(-β·t)
        // and the local decay rate approaches β.
        let m = model();
        let t1 = 50.0;
        let t2 = 60.0;
        let c1 = m.concentration(t1).unwrap();
        let c2 = m.concentration(t2).unwrap();
        // -ln(c2/c1)/(t2-t1) should approach β = 0.1.
        let apparent_rate = -(c2 / c1).ln() / (t2 - t1);
        assert!(
            (apparent_rate - 0.1).abs() < 1e-6,
            "apparent terminal rate = {apparent_rate}"
        );
    }

    #[test]
    fn terminal_half_life_is_ln2_over_beta() {
        let m = model();
        let expected = std::f64::consts::LN_2 / 0.1;
        assert!(
            (m.terminal_half_life() - expected).abs() < TOL,
            "t½ = {}",
            m.terminal_half_life()
        );
    }

    #[test]
    fn auc_is_sum_of_term_ratios() {
        // ANALYTIC: AUC = A/α + B/β = 8/1 + 2/0.1 = 8 + 20 = 28.
        let m = model();
        assert!((m.auc() - 28.0).abs() < TOL, "AUC = {}", m.auc());
    }

    #[test]
    fn concentration_is_monotonically_decreasing() {
        let m = model();
        let mut prev = f64::INFINITY;
        for i in 0..30 {
            let t = i as f64;
            let c = m.concentration(t).unwrap();
            assert!(c < prev, "C should decrease: C({t}) = {c}");
            prev = c;
        }
    }

    #[test]
    fn rejects_negative_time() {
        assert!(model().concentration(-0.5).is_err());
    }

    #[test]
    fn partial_auc_is_zero_at_origin_and_total_at_infinity() {
        let m = model();
        assert!(m.auc_to(0.0).unwrap().abs() < TOL, "AUC(0→0) = 0");
        // The terminal half-life is ln2/0.1 ≈ 6.93 h; by t = 1000 both
        // phases are fully accrued.
        let far = m.auc_to(1000.0).unwrap();
        assert!((far - m.auc()).abs() < 1e-9, "AUC(0→∞) → 28, got {far}");
    }

    #[test]
    fn partial_auc_matches_a_trapezoidal_integral_of_concentration() {
        // Numerically integrate the biexponential over [0, T] and compare
        // to the closed form.
        let m = model();
        let big_t = 90.0;
        let n = 180_000;
        let dt = big_t / n as f64;
        let mut area = 0.0;
        let mut prev = m.concentration(0.0).unwrap();
        for i in 1..=n {
            let c = m.concentration(i as f64 * dt).unwrap();
            area += 0.5 * (prev + c) * dt;
            prev = c;
        }
        let closed = m.auc_to(big_t).unwrap();
        assert!(
            (area - closed).abs() < 1e-5,
            "trapz {area} vs closed {closed}"
        );
    }

    #[test]
    fn partial_auc_increases_monotonically_and_rejects_bad_time() {
        let m = model();
        let mut prev = -1.0;
        for i in 0..40 {
            let a = m.auc_to(i as f64).unwrap();
            assert!(a > prev, "AUC(0→t) should increase: {a} after {prev}");
            assert!(a <= m.auc() + TOL, "never exceeds the total");
            prev = a;
        }
        assert!(m.auc_to(-1.0).is_err());
        assert!(m.auc_to(f64::NAN).is_err());
    }
}
