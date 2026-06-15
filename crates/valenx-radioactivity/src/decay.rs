//! Single-nuclide exponential decay.
//!
//! A radioactive nuclide decays by a first-order process: in any short
//! interval each surviving nucleus has a fixed, time-independent
//! probability per unit time of decaying. That probability is the
//! **decay constant** `lambda` (units: inverse time). Everything else in
//! this module is an algebraic consequence of `lambda`.
//!
//! Letting `N(t)` be the number of un-decayed nuclei at time `t`:
//!
//! - Number law: `N(t) = N0 * exp(-lambda * t)`.
//! - Half-life: the time for the population to halve,
//!   `t_half = ln(2) / lambda`.
//! - Mean (average) life: the expected lifetime of one nucleus,
//!   `tau = 1 / lambda = t_half / ln(2)`.
//! - Activity: the decay *rate*, `A(t) = lambda * N(t)`, measured in
//!   becquerel (one decay per second). Because activity is just `lambda`
//!   times the population, it obeys the same exponential law,
//!   `A(t) = A0 * exp(-lambda * t)`.
//!
//! The [`Nuclide`] type stores `lambda` and exposes these as methods,
//! with constructors that accept whichever of `lambda`, `t_half` or `tau`
//! you happen to have.

use serde::{Deserialize, Serialize};

use crate::error::{RadioactivityError, Result};

/// Natural logarithm of 2, `ln(2)`.
///
/// The conversion factor between half-life and mean life:
/// `t_half = LN_2 * tau` and `tau = t_half / LN_2`. Exposed because it is
/// the single physical constant tying the two timescales together and
/// callers frequently want it for their own checks.
pub const LN_2: f64 = std::f64::consts::LN_2;

/// A radioactive nuclide characterised by its decay constant.
///
/// `lambda` is the probability per unit time that any one surviving
/// nucleus decays; it is the only stored quantity because half-life,
/// mean life and activity all derive from it. All time-valued outputs
/// share whatever time unit `lambda` was expressed in (e.g. if `lambda`
/// is per-year, [`Nuclide::half_life`] is in years).
///
/// Construct with [`Nuclide::from_decay_constant`],
/// [`Nuclide::from_half_life`], or [`Nuclide::from_mean_life`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Nuclide {
    /// Decay constant `lambda` (strictly positive, inverse time).
    lambda: f64,
}

impl Nuclide {
    /// Builds a nuclide from its decay constant `lambda` (inverse time).
    ///
    /// # Errors
    ///
    /// [`RadioactivityError::NonPositive`] if `lambda` is not a strictly
    /// positive, finite number — a non-positive decay constant has no
    /// physical meaning (a stable nuclide is the `lambda -> 0` limit,
    /// which this model does not represent).
    pub fn from_decay_constant(lambda: f64) -> Result<Self> {
        let lambda = RadioactivityError::require_positive("lambda", lambda)?;
        Ok(Self { lambda })
    }

    /// Builds a nuclide from its half-life `t_half`.
    ///
    /// Uses `lambda = ln(2) / t_half`.
    ///
    /// # Errors
    ///
    /// [`RadioactivityError::NonPositive`] if `half_life` is not a
    /// strictly positive, finite number.
    pub fn from_half_life(half_life: f64) -> Result<Self> {
        let half_life = RadioactivityError::require_positive("half_life", half_life)?;
        Ok(Self {
            lambda: LN_2 / half_life,
        })
    }

    /// Builds a nuclide from its mean (average) life `tau`.
    ///
    /// Uses `lambda = 1 / tau`.
    ///
    /// # Errors
    ///
    /// [`RadioactivityError::NonPositive`] if `mean_life` is not a
    /// strictly positive, finite number.
    pub fn from_mean_life(mean_life: f64) -> Result<Self> {
        let mean_life = RadioactivityError::require_positive("mean_life", mean_life)?;
        Ok(Self {
            lambda: 1.0 / mean_life,
        })
    }

    /// The decay constant `lambda` (inverse time).
    #[inline]
    pub fn decay_constant(&self) -> f64 {
        self.lambda
    }

    /// The half-life `t_half = ln(2) / lambda`: the time for the
    /// population (or activity) to fall to one half of any starting value.
    #[inline]
    pub fn half_life(&self) -> f64 {
        LN_2 / self.lambda
    }

    /// The mean (average) life `tau = 1 / lambda`: the expected lifetime
    /// of a single nucleus. Equals `half_life / ln(2)`, so it is always
    /// longer than the half-life.
    #[inline]
    pub fn mean_life(&self) -> f64 {
        1.0 / self.lambda
    }

    /// Number of un-decayed nuclei remaining at time `t`,
    /// `N(t) = N0 * exp(-lambda * t)`.
    ///
    /// # Errors
    ///
    /// [`RadioactivityError::NonPositive`] if `n0` is not strictly
    /// positive, or if `t` is negative or non-finite (the law is only
    /// posed forward from the `t = 0` reference instant).
    pub fn remaining(&self, n0: f64, t: f64) -> Result<f64> {
        let n0 = RadioactivityError::require_positive("n0", n0)?;
        let t = RadioactivityError::require_non_negative("t", t)?;
        Ok(n0 * (-self.lambda * t).exp())
    }

    /// Fraction of the original population still un-decayed at time `t`,
    /// `N(t) / N0 = exp(-lambda * t)`. Always in `(0, 1]`.
    ///
    /// # Errors
    ///
    /// [`RadioactivityError::NonPositive`] if `t` is negative or
    /// non-finite.
    pub fn remaining_fraction(&self, t: f64) -> Result<f64> {
        let t = RadioactivityError::require_non_negative("t", t)?;
        Ok((-self.lambda * t).exp())
    }

    /// Activity `A = lambda * N` for a population of `n` nuclei, in the
    /// reciprocal of `lambda`'s time unit (becquerel when `lambda` is
    /// per-second). Activity is the number of decays per unit time.
    ///
    /// # Errors
    ///
    /// [`RadioactivityError::NonPositive`] if `n` is not strictly
    /// positive.
    pub fn activity(&self, n: f64) -> Result<f64> {
        let n = RadioactivityError::require_positive("n", n)?;
        Ok(self.lambda * n)
    }

    /// Activity at time `t` given an initial population `n0`,
    /// `A(t) = lambda * N0 * exp(-lambda * t)`. Equivalently
    /// `A(t) = A0 * exp(-lambda * t)` with `A0 = lambda * N0`: activity
    /// decays with the same constant as the population.
    ///
    /// # Errors
    ///
    /// [`RadioactivityError::NonPositive`] if `n0` is not strictly
    /// positive, or if `t` is negative or non-finite.
    pub fn activity_at(&self, n0: f64, t: f64) -> Result<f64> {
        let n_t = self.remaining(n0, t)?;
        Ok(self.lambda * n_t)
    }

    /// Time required for the population (or activity) to fall to a given
    /// `fraction` of its starting value, `t = -ln(fraction) / lambda`.
    ///
    /// `fraction == 0.5` returns the half-life; `fraction == 1.0` returns
    /// `0`. This is the inverse of [`remaining_fraction`](Self::remaining_fraction).
    ///
    /// # Errors
    ///
    /// [`RadioactivityError::OutOfRange`] if `fraction` is not in the
    /// half-open interval `(0, 1]` (the population can never grow, and
    /// reaching exactly zero takes infinite time).
    pub fn time_to_fraction(&self, fraction: f64) -> Result<f64> {
        if !(fraction.is_finite() && fraction > 0.0 && fraction <= 1.0) {
            return Err(RadioactivityError::OutOfRange {
                what: "fraction",
                value: fraction,
                interval: "(0, 1]",
            });
        }
        Ok(-fraction.ln() / self.lambda)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute-difference assertion helper for float comparisons.
    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn lambda_from_half_life_matches_ln2_over_t_half() {
        // GROUND TRUTH: lambda = ln(2) / t_half.
        let t_half = 8.0;
        let nuc = Nuclide::from_half_life(t_half).unwrap();
        assert!(
            close(nuc.decay_constant(), LN_2 / t_half, 1e-15),
            "lambda = {}",
            nuc.decay_constant()
        );
    }

    #[test]
    fn half_life_round_trips_through_decay_constant() {
        let t_half = 5730.0; // carbon-14, years
        let nuc = Nuclide::from_half_life(t_half).unwrap();
        assert!(
            close(nuc.half_life(), t_half, 1e-9),
            "t_half = {}",
            nuc.half_life()
        );
        // Rebuild from the recovered lambda and check t_half survives.
        let nuc2 = Nuclide::from_decay_constant(nuc.decay_constant()).unwrap();
        assert!(close(nuc2.half_life(), t_half, 1e-9));
    }

    #[test]
    fn population_halves_at_one_half_life() {
        // GROUND TRUTH: N(t_half) = N0 / 2.
        let nuc = Nuclide::from_half_life(10.0).unwrap();
        let n0 = 1.0e6;
        let n = nuc.remaining(n0, nuc.half_life()).unwrap();
        assert!(close(n, n0 / 2.0, 1e-6), "N(t_half) = {n}");
    }

    #[test]
    fn population_quarters_at_two_half_lives() {
        // GROUND TRUTH: N(2 t_half) = N0 / 4, and 1/8 at three.
        let nuc = Nuclide::from_half_life(3.0).unwrap();
        let n0 = 8192.0;
        let two = nuc.remaining(n0, 2.0 * nuc.half_life()).unwrap();
        let three = nuc.remaining(n0, 3.0 * nuc.half_life()).unwrap();
        assert!(close(two, n0 / 4.0, 1e-6), "N(2 t_half) = {two}");
        assert!(close(three, n0 / 8.0, 1e-6), "N(3 t_half) = {three}");
    }

    #[test]
    fn remaining_fraction_is_one_at_zero_and_half_at_half_life() {
        let nuc = Nuclide::from_decay_constant(0.25).unwrap();
        assert!(close(nuc.remaining_fraction(0.0).unwrap(), 1.0, 1e-15));
        assert!(close(
            nuc.remaining_fraction(nuc.half_life()).unwrap(),
            0.5,
            1e-12
        ));
    }

    #[test]
    fn activity_equals_lambda_times_n() {
        // GROUND TRUTH: A = lambda * N.
        let lambda = 0.04;
        let nuc = Nuclide::from_decay_constant(lambda).unwrap();
        let n = 2.5e8;
        assert!(
            close(nuc.activity(n).unwrap(), lambda * n, 1e-3),
            "A = {}",
            nuc.activity(n).unwrap()
        );
    }

    #[test]
    fn activity_halves_at_one_half_life() {
        // GROUND TRUTH: A(t_half) = A0 / 2.
        let nuc = Nuclide::from_half_life(12.0).unwrap();
        let n0 = 1.0e9;
        let a0 = nuc.activity(n0).unwrap();
        let a_half = nuc.activity_at(n0, nuc.half_life()).unwrap();
        assert!(
            close(a_half, a0 / 2.0, 1e-3),
            "A(t_half) = {a_half}, A0 = {a0}"
        );
    }

    #[test]
    fn activity_obeys_same_exponential_as_population() {
        // A(t) / A0 must equal N(t) / N0 for every t.
        let nuc = Nuclide::from_half_life(7.0).unwrap();
        let n0 = 3.3e7;
        let a0 = nuc.activity(n0).unwrap();
        for &t in &[0.0, 1.5, 7.0, 21.0, 50.0] {
            let a_ratio = nuc.activity_at(n0, t).unwrap() / a0;
            let n_ratio = nuc.remaining_fraction(t).unwrap();
            assert!(
                close(a_ratio, n_ratio, 1e-12),
                "t = {t}: {a_ratio} vs {n_ratio}"
            );
        }
    }

    #[test]
    fn mean_life_equals_half_life_over_ln2() {
        // GROUND TRUTH: tau = t_half / ln(2), and tau = 1 / lambda.
        let t_half = 22.3;
        let nuc = Nuclide::from_half_life(t_half).unwrap();
        assert!(
            close(nuc.mean_life(), t_half / LN_2, 1e-12),
            "tau = {}",
            nuc.mean_life()
        );
        assert!(close(nuc.mean_life(), 1.0 / nuc.decay_constant(), 1e-12));
        // Mean life is always longer than the half-life (1/ln2 ~ 1.4427).
        assert!(nuc.mean_life() > nuc.half_life());
    }

    #[test]
    fn fraction_remaining_at_one_mean_life_is_one_over_e() {
        // After one mean life exactly 1/e of the sample survives.
        let nuc = Nuclide::from_mean_life(4.0).unwrap();
        let frac = nuc.remaining_fraction(nuc.mean_life()).unwrap();
        assert!(
            close(frac, 1.0 / std::f64::consts::E, 1e-12),
            "frac = {frac}"
        );
    }

    #[test]
    fn time_to_fraction_inverts_remaining_fraction() {
        let nuc = Nuclide::from_half_life(15.0).unwrap();
        // Half should come back as exactly the half-life.
        assert!(close(
            nuc.time_to_fraction(0.5).unwrap(),
            nuc.half_life(),
            1e-12
        ));
        // Round-trip an arbitrary fraction.
        let f = 0.137;
        let t = nuc.time_to_fraction(f).unwrap();
        assert!(
            close(nuc.remaining_fraction(t).unwrap(), f, 1e-12),
            "t = {t}"
        );
        // fraction == 1 means zero elapsed time.
        assert!(close(nuc.time_to_fraction(1.0).unwrap(), 0.0, 1e-15));
    }

    #[test]
    fn constructors_reject_bad_inputs() {
        assert!(Nuclide::from_decay_constant(0.0).is_err());
        assert!(Nuclide::from_decay_constant(-1.0).is_err());
        assert!(Nuclide::from_half_life(f64::NAN).is_err());
        assert!(Nuclide::from_mean_life(f64::INFINITY).is_err());
    }

    #[test]
    fn evaluation_rejects_bad_inputs() {
        let nuc = Nuclide::from_half_life(1.0).unwrap();
        assert!(nuc.remaining(-1.0, 1.0).is_err()); // n0 <= 0
        assert!(nuc.remaining(1.0, -0.5).is_err()); // t < 0
        assert!(nuc.activity(0.0).is_err()); // n <= 0
        assert!(nuc.time_to_fraction(0.0).is_err()); // fraction not in (0, 1]
        assert!(nuc.time_to_fraction(1.5).is_err());
    }
}
